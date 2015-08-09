/* 
    liblayout, an experimental 2D layout library.
    Copyright (C) 2006 Adrian Secord.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

    Contact information for the author is available at http://mrl.nyu.edu/~ajsecord/
    or send an email to ajsecord *at* cs *dot* nyu *dot* edu.
*/

#include <layout/layout.h>
#include <layout/overlap.h>
#include <layout/macopt.h>

#include <float.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

/** Increment a pointer \c count times, skipping \c skip bytes each time. */
#define LAY_INCR_POINTER(pointer, count, skip) (((char*) (pointer)) + (count) * (skip))

/** Get a pointer to the \c count'th position. */
#define LAY_POS_POINTER(state, count) ((lay_coord_t*) LAY_INCR_POINTER(state->pos, count, state->pos_skip))

/** Get a pointer to the next position. */
#define LAY_NEXT_POS(state, pointer) ((lay_coord_t*) (((char*) (pointer)) + state->pos_skip))

/** Get a pointer to the \c count'th size. */
#define LAY_SIZE_POINTER(state, count) ((lay_extent_t*) LAY_INCR_POINTER(state->size, count, state->size_skip))

/** Get a pointer to the next size. */
#define LAY_NEXT_SIZE(state, pointer) ((lay_extent_t*) (((char*) (pointer)) + state->size_skip))

/** Layout state */
struct lay_state {
    /* Rectangle list */
    int num_rects;                  /**< The number of rectangles. */

    lay_coord_t* pos;               /**< Pointer to position data. */
    ptrdiff_t pos_skip;             /**< Number of bytes to skip to get to the next position. */
    
    lay_extent_t* size;             /**< Pointer to size data. */
    ptrdiff_t size_skip;            /**< Number of bytes to skip to get to the next size. */
    
    /* Optimization settings */
    lay_real_t overlap_weight;      /**< The overlap penalty weight. */
    lay_real_t edge_weight;         /**< The edge penalty weight. */
    lay_real_t center_weight;       /**< The center penalty weight. */
    lay_real_t orig_pos_weight;     /**< The original position penalty weight. */
    
    /* Temporary storage */
    float* dof;                     /**< The degrees of freedom, modified by the optimizer. (macopt uses float) */

    /* Optimizer arguments */
    macopt_args opt_args;           /**< Optimizer arguments. */
};

/** Check whether num_rect-based temps have been allocated. */
static int has_num_rect_temps(const lay_statep state) {
    return state->dof != NULL;
}

/** Create any temporary space that is allocated based on the number of rectangles. */
static void create_num_rect_temps(lay_statep state) {
    assert(state && state->num_rects >= 0);
    
    assert(!state->dof);
    if (state->num_rects > 0) {
        state->dof = malloc(state->num_rects * 2 * sizeof(float));
    }
    
    /* Sanity check */
    assert(has_num_rect_temps(state));
}

/** Create any temporary space if it is not already allocated.  Convenience function. */
static void ensure_num_rect_temps(lay_statep state) {
    if (!has_num_rect_temps(state))
        create_num_rect_temps(state);
}

/** Free any temporary space that is allocated based on the number of rectangles. */
static void destroy_num_rect_temps(lay_statep state) {
    assert(state);
    
    if (state->dof) {
        free(state->dof);
        state->dof = NULL;
    }
}

/** Copy user-land positions into a dense array of optimizer floating-point values. */
static void copy_user_pos_to_array(const lay_statep state, float* array) {
    lay_coord_t* p;
    int i;
    
    assert(array && state && state->pos);
    p = state->pos;
    for (i = 0; i < state->num_rects; ++i) {
        array[2*i]   = (float)p[0];
        array[2*i+1] = (float)p[1];
        p = LAY_NEXT_POS(state, p);
    }
}

/** Copy a dense array of optimizer floating-point values into user-land positions. */
static void copy_array_to_user_pos(const float* array, lay_statep state) {
    lay_coord_t* p;
    int i;
    
    assert(array && state && state->pos);
    p = state->pos;
    for (i = 0; i < state->num_rects; ++i) {
        p[0] = (lay_coord_t)array[2*i];
        p[1] = (lay_coord_t)array[2*i+1];
        p = LAY_NEXT_POS(state, p);
    }
}

/** Copy user-land extents into a dense array of lay_extent_t. */
static void copy_user_size_to_array(const lay_statep state, lay_extent_t* array) {
    lay_extent_t* p;
    int i;
    
    assert(array && state && state->pos);
    p = state->pos;
    for (i = 0; i < state->num_rects; ++i) {
        array[2*i]   = p[0];
        array[2*i+1] = p[1];
        p = LAY_NEXT_SIZE(state, p);
    }
}

/** Copy a dense array of lay_extent_t into user-land extents. */
static void copy_array_to_user_size(const lay_extent_t* array, lay_statep state) {
    lay_extent_t* p;
    int i;
    
    assert(array && state && state->pos);
    p = state->pos;
    for (i = 0; i < state->num_rects; ++i) {
        p[0] = array[2*i];
        p[1] = array[2*i+1];
        p = LAY_NEXT_SIZE(state, p);
    }
}


lay_statep lay_create_state() {
    lay_statep state = malloc(sizeof(struct lay_state));
    assert(state);
    
    state->overlap_weight = 1;
    state->edge_weight = 0;
    state->center_weight = 0;
    state->orig_pos_weight = 0;
    
    state->dof = NULL;

    lay_register_rects(state, NULL, 0, NULL, 0, 0);
    
    /* Setup optimizer arguments */
    macopt_defaults(&state->opt_args); 
    state->opt_args.itmax = 400;               /* Maximum interations */
    state->opt_args.verbose = 0;               /* Reporting level */
    state->opt_args.tol = 1e-3;                /* Finishing tolerance */
    state->opt_args.end_if_small_step = 1 ;    /* Finish if step gets small, otherwise grad mag gets small */
    
    /* Sanity check */
    assert(lay_verify_state(state));
    
    return state;
}

int lay_verify_state(const lay_statep state) {
    if (!state)
        return 0;
    
    if (state->num_rects < 0 || state->pos_skip == 0 || state->size_skip == 0)
        return 0;
    
    return 1;
}

void lay_destroy_state(lay_statep state) {
    assert(state);
    
    destroy_num_rect_temps(state);
    
    free(state);
}

void lay_register_rects(lay_statep state, 
                        lay_coord_t* rect_pos, const ptrdiff_t pos_skip,
                        lay_extent_t* rect_size, const ptrdiff_t size_skip,
                        const int count
                        ) {
    assert(state);
    state->pos = rect_pos;
    state->pos_skip = (pos_skip != 0 ? pos_skip : 2 * sizeof(lay_coord_t));
    state->size = rect_size;
    state->size_skip = (size_skip != 0 ? size_skip : 2 * sizeof(lay_extent_t));
    state->num_rects = count;
    
    /* Force reallocation of num_rect-based temps next time they are needed. */
    destroy_num_rect_temps(state);
}

lay_real_t lay_get_overlap_weight(const lay_statep state) {
    assert(state);
    return state->overlap_weight;
}

void lay_set_overlap_weight(lay_statep state, const lay_real_t weight) {
    assert(state);
    state->overlap_weight = weight;
}

lay_real_t lay_get_center_weight(const lay_statep state) {
    assert(state);
    return state->center_weight;
}

void lay_set_center_weight(lay_statep state, const lay_real_t weight) {
    assert(state);
    state->center_weight = weight;
}

lay_real_t lay_get_edge_weight(const lay_statep state) {
    assert(state);
    return state->edge_weight;
}

void lay_set_edge_weight(lay_statep state, const lay_real_t weight) {
    assert(state);
    state->edge_weight = weight;
}

lay_real_t lay_get_orig_pos_weight(const lay_statep state) {
    assert(state);
    return state->orig_pos_weight;
}

void lay_set_orig_pos_weight(lay_statep state, const lay_real_t weight) {
    assert(state);
    state->orig_pos_weight = weight;
}

/** Evaluate the energy and optionally the gradient of the rectangle configuration 
    \c input.  If \c global_grad is not NULL, then it must contain enough space 
    for the number of degrees of freedom per rectangle for *every* rectangle, 
    including fixed rectangles.
    Note that the positions contained in state are *not* used, rather those 
    stored densely in \c cur_pos.
*/
static lay_real_t eval(const lay_statep state, 
                       const lay_coord_t* cur_pos, 
                       lay_real_t* global_grad) {
    const lay_coord_t *p;
    lay_coord_t *q;
    lay_extent_t *size1, *size2;
    lay_real_t layout_energy, dist[4];
    lay_real_t local_grad[4];
    int i, j, grad_num_dof;
    
    assert(lay_verify_state(state));
   
    /* The number of degrees of freedom in the passed-in gradient. */
    grad_num_dof = (global_grad != NULL ? 2 * state->num_rects : 0);
    
    layout_energy = 0;
    for (i = 0; i < grad_num_dof; ++i)
        global_grad[i] = 0;

    for (i = 0; i < state->num_rects; ++i) {
        size1 = LAY_SIZE_POINTER(state, i);
        
        for (j = i+1; j < state->num_rects; ++j) {
            size2 = LAY_SIZE_POINTER(state, j);
            
            layout_energy += lay_overlap_area(cur_pos + 2 * i, size1,
                                              cur_pos + 2 * j, size2, 
                                              local_grad);
            if (global_grad) {
                global_grad[2*i  ] += local_grad[0];
                global_grad[2*i+1] += local_grad[1];
                global_grad[2*j  ] += local_grad[2];
                global_grad[2*j+1] += local_grad[3];
            }
        }
    }
    
    layout_energy *= state->overlap_weight;
    for (i = 0; i < grad_num_dof; ++i)
        global_grad[i] *= state->overlap_weight;
        
#if 0
    /* Add penalty terms to keep rectangles on-screen. */
    if (edge_weight != 0) {
        for (i = 0; i < state->num_rects; ++i) {
            dist[0] = -input->items[i].x;
            dist[1] = -input->items[i].y;
            dist[2] = input->items[i].x + input->items[i].width - screen_width;
            dist[3] = input->items[i].y + input->items[i].height - screen_height;
            
            grad[0][0] = -1; grad[0][1] =  0;
            grad[1][0] =  0; grad[1][1] = -1;
            grad[2][0] =  1; grad[2][1] =  0;
            grad[3][0] =  0; grad[3][1] =  1;
            
            for (j = 0; j < 4; ++j) {
                if (dist[j] > 0) {
                    layout_energy += edge_weight * dist[j] * dist[j];
                    if (global_grad) {
                        global_grad[2*i  ] += edge_weight * 2 * dist[j] * grad[j][0];
                        global_grad[2*i+1] += edge_weight * 2 * dist[j] * grad[j][1];
                    }
                }
            }
        }
    }
#endif
    
#if 0
    /* Penalty terms to keep rectangles near the center of the screen. */
    if (center_weight != 0) {
        for (i = 0; i < state->num_rects; ++i) {
            dist[0] = input->items[i].x - screen_width/2;
            dist[1] = input->items[i].y - screen_height/2;
            dist[2] = dist[0] * dist[0] + dist[1] * dist[1];
            
            layout_energy += center_weight * dist[2];
            
            if (global_grad && dist[2] != 0) {
                global_grad[2*i  ] += center_weight * 2 * dist[0];
                global_grad[2*i+1] += center_weight * 2 * dist[1];
            }
        }
    }
#endif
    
    /* Add terms to keep rectangles near their original positions. */
    if (state->orig_pos_weight != 0) {
        for (i = 0; i < state->num_rects; ++i) {
            p = cur_pos + 2 * i;
            q = LAY_POS_POINTER(state, i);
            
            dist[0] = p[0] - q[0];
            dist[1] = p[1] - q[1];
            dist[2] = dist[0] * dist[0] + dist[1] * dist[1];
            
            layout_energy += state->orig_pos_weight * dist[2];
            
            if (global_grad) {
                global_grad[2*i  ] += state->orig_pos_weight * 2 * dist[0];
                global_grad[2*i+1] += state->orig_pos_weight * 2 * dist[1];
            }
        }
    }
    
    return layout_energy;
}

static float energy(lay_coord_t* x, void* args) {
    return eval((lay_statep) args, x + 1, NULL);    /* Convert one-based array */
}

static void vgrad_energy(lay_coord_t* x, lay_coord_t* grad, void* args) {
    eval((lay_statep) args, x + 1, grad + 1);       /* Convert one-based array */
}

void lay_optimize(lay_statep state) {
    assert(lay_verify_state(state));

    ensure_num_rect_temps(state);
    
    /* Copy the original positions into the minimizer's current state vector. */
    copy_user_pos_to_array(state, state->dof);
    
    /* Check that the gradient_function is the gradient of the function. 
       Note the adjustment for the moronic one-based arrays Numerical Recipes requires.
    */
#if 0
    maccheckgrad(state->dof - 1, 2 * state->num_rects, 1e-3, energy, state, vgrad_energy, state, 0);
#endif
    
    macoptII(state->dof - 1, 2 * state->num_rects, vgrad_energy, state, &state->opt_args);
    
    copy_array_to_user_pos(state->dof, state);
}


