#ifndef LAY_SIM_ANNEAL_H
#define LAY_SIM_ANNEAL_H

/** \file sim_anneal.h
    Functional minimization by simulated annealing.
*/

#include <layout/types.h>

#ifdef __cplusplus
extern "C" {
#endif
    
    /** A function that generates a neighbor of \c input and stores it in \c output. 
        The \c context parameter allows the user to store arbitrary information.
    */
    typedef void (*lay_sa_neighbor_func)(void* context, const void const* input, const double cur_temp, void* output);
    
    /** A function to evaluate a particular point. 
        The \c context parameter allows the user to store arbitrary information. 
    */
    typedef double (*lay_sa_eval_func)(void* context, const void const* input);
    
    /** A function to reduce the temperature on some schedule.
        The \c context parameter allows the user to store arbitrary information.
    */
    typedef double (*lay_sa_temp_func)(void* context, const double cur_temp);
    
    /** Operations required for simulated annealing. */
    typedef struct {
        lay_sa_neighbor_func neighbor;      /**< Neighbor-generating function. */
        lay_sa_eval_func eval;              /**< Evaluation function. */
        lay_sa_temp_func reduce_temp;       /**< Temperature-reduction schedule function. */
    } lay_sa_ops;
    
    /** Run the simulated annealing algorithm to search for a minimum point. */
    int lay_sim_anneal(void* context,
                       const lay_sa_ops ops,
                       const lay_point_ops point_ops, 
                       void* start_stop_point,
                       double* temp, const double temp_min, 
                       const int num_iters, const int max_total_iters, long* rand_seed);        
    
#ifdef __cplusplus
}
#endif

#endif