#include <layout/sim_anneal.h>
#include "random/random.h"

#include <stdlib.h>
#include <math.h>
#include <assert.h>

int lay_sim_anneal(void* context, 
                   const lay_sa_ops ops,
                   const lay_point_ops point_ops, 
                   void* start_stop_point,
                   double* temp, const double temp_min, 
                   const int num_iters, const int max_total_iters, long* rand_seed) {
    
    double old_temp, cur_value, new_value;
    void* point = NULL, *new_point = NULL;
    int i, iters = 0;

    assert(ops.neighbor && ops.eval && ops.reduce_temp);
    assert(start_stop_point && temp && rand_seed);
    assert(point_ops.create && point_ops.copy && point_ops.destroy);

    point = point_ops.create();
    new_point = point_ops.create();
    assert(point && new_point);
    point_ops.copy(start_stop_point, point);
    cur_value = ops.eval(context, point);
    
    while (*temp >= temp_min && iters < max_total_iters) {
        for (i = 0; i < num_iters && iters < max_total_iters; ++i) {
            ++iters;
            
            /* Select a new point in the neighborhood of the current point */
            ops.neighbor(context, point, *temp, new_point);
            new_value = ops.eval(context, new_point);
            
            /* If better, take it.  If not, maybe take it. */
            if (new_value < cur_value || 
                uniform_dev(rand_seed) < exp((cur_value - new_value) / *temp)) {
                point_ops.copy(new_point, point);
                cur_value = new_value;
            }
        }
        
        /* Reduce the temperature. */
        old_temp = *temp;
        *temp = ops.reduce_temp(context, *temp);
        assert(*temp <= old_temp);
    }
    
    point_ops.copy(point, start_stop_point);
    point_ops.destroy(point);
    point_ops.destroy(new_point);
    
    return iters;
}
