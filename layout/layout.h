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

#ifndef LAY_LAYOUT_H
#define LAY_LAYOUT_H

/** \mainpage The liblayout layout library
* An experimental layout library. 
*/

/** \file layout/layout.h
* Layout function declarations and types.
*/

#include <stddef.h>
#include <layout/types.h>

#ifdef __cplusplus
extern "C" {
#endif

    /** Pointer to the internal liblayout state. */
    typedef struct lay_state* lay_statep;
    
    /** \name Initialization and setup functions */
    /*@{*/
    
    /** Create an internal state structure for use with subsequent liblayout functions.
        When finished, must be destroyed using lay_destroy_state(). 
    */
    lay_statep lay_create_state();
    
    /** Destroy an internal state structure. */
    void lay_destroy_state(lay_statep state);
    
    /** Assert that the internal state is consistent. */
    int lay_verify_state(lay_statep state);
    
    /** Register the positions and sizes of a set of rectangles with liblayout. 
        \param state The internal state structure.
        \param rect_pos A pointer to the x-coordinate of the first size.  The
        y-coordinate of the first size is assumed to exist at <tt>*(rect_pos + 1)</tt>.  
        The pointer can be NULL, but no operations requiring sizes may be used.
        \param pos_skip The number of bytes to add to \c rect_pos to get to the 
        x-coordinate of the next size.  If zero, then the data is assumed to be 
        tightly packed, and pos_skip will be set to <tt>2 * sizeof(lay_coord_t)</tt>.  
        Can be negative.
        \param rect_size A pointer to the width of the first size.  The
        height of the first size is assumed to exist at <tt>*(p + 1)</tt>.  The 
        pointer can be NULL, but no operations requiring sizes may be used.
        \param size_skip The number of bytes to add to \c p to get to the width 
        of the next size.  If zero, then the data is assumed to be tightly packed, 
        and skip will be set to <tt>2 * sizeof(lay_extent_t)</tt>.  Can be negative.
        \param count The number of rectangles.
    */
    void lay_register_rects(lay_statep state, 
                            lay_coord_t* rect_pos, const ptrdiff_t pos_skip,
                            lay_extent_t* rect_size, const ptrdiff_t size_skip,
                            const int count
                            );

    /*@}*/
    
    /** \name Optimization settings */
    /*@{*/
    
    /** Get the overlap penalty weight. */
    lay_real_t lay_get_overlap_weight(const lay_statep state);

    /** Set the overlap penalty weight. */
    void lay_set_overlap_weight(lay_statep state, const lay_real_t weight);
    
    /** Get the edge penalty weight. */
    lay_real_t lay_get_edge_weight(const lay_statep state);

    /** Set the edge penalty weight. */
    void lay_set_edge_weight(lay_statep state, const lay_real_t weight);
    
    /** Get the center penalty weight. */
    lay_real_t lay_get_center_weight(const lay_statep state);

    /** Set the center penalty weight. */
    void lay_set_center_weight(lay_statep state, const lay_real_t weight);
    
    /** Get the original position penalty weight. */
    lay_real_t lay_get_orig_pos_weight(const lay_statep state);
    
    /** Set the original position penalty weight. */
    void lay_set_orig_pos_weight(lay_statep state, const lay_real_t weight);
    
    /*@}*/
    
    /** Optimize the position of the input rectangles. 
        Overwrites the current positions with optimized positions.
    */
    void lay_optimize(lay_statep state);
    
#ifdef __cplusplus
}
#endif

#endif
