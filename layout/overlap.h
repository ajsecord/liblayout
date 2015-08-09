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

#ifndef LAY_OVERLAP_H
#define LAY_OVERLAP_H

/** \file layout/overlap.h
* Overlap computations for rectangles.
*/

#include <layout/types.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Compute the overlap area between two rectangles.
    If \c grad is non-null, then the gradient is saved as an array of four 
    real values, calculated with respect to rect1.x, to rect1.y, to rect2.x, 
    and to rect2.y, in that order.
*/
lay_real_t lay_overlap_area(const lay_coord_t* pos1, const lay_extent_t* size1, 
                            const lay_coord_t* pos2, const lay_extent_t* size2,
                            lay_real_t* grad);

#if 0
/** Compute the total overlap between all rectangles in a list.
    If \c grad is not NULL, then it must have enough space for <tt>2 * num_rects</tt>
    elements and upon return it will contain the gradient of the overlap for 
    rectangle \c i  with respect to \c x and \c y in <tt>grad[2 * i]</tt> and 
    <tt>grad[2 * i + 1]</tt>.
*/
lay_real_t lay_all_overlap_area(const int num_rects, 
                                const lay_coord_t* pos, const lay_extent_t* size, 
                                lay_real_t* grad);
    
/** Return non-zero if rectangle \c index overlaps with any other rectangle in the list. */
int lay_any_overlap(const int num_rects, 
                    const lay_coord_t* pos, const lay_extent_t* size, 
                    const int index);

/** Return non-zero if any rectangle in the list overlaps any other. */
int lay_any_overlap_any(const int num_rects, 
                        const lay_coord_t* pos, const lay_extent_t* size);

#endif

#ifdef __cplusplus
}
#endif

#endif
