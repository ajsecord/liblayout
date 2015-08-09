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

#include <layout/overlap.h>

#include <stdlib.h>
#include <assert.h>
#include <math.h>

lay_real_t lay_overlap_area(const lay_coord_t* pos1, const lay_extent_t* size1, 
                            const lay_coord_t* pos2, const lay_extent_t* size2,
                            lay_real_t* grad) {
    int i;
    lay_coord_t x_len, y_len, x_overlap, y_overlap, x_overlap_grad_x, y_overlap_grad_y;

    assert(pos1 && size1 && pos2 && size2);

    x_len = 2 * (pos2[0] - pos1[0]) + (size2[0] - size1[0]);
    y_len = 2 * (pos2[1] - pos1[1]) + (size2[1] - size1[1]);
    x_overlap = (size1[0] + size2[0]) - LAY_COORD_ABS(x_len);
    y_overlap = (size1[1] + size2[1]) - LAY_COORD_ABS(y_len);
    
    if (x_overlap <= 0 || y_overlap <= 0) {
        if (grad) 
            for (i = 0; i < 4; ++i)
                grad[i] = 0;
        return 0;
    }
    
    if (grad) {
        x_overlap_grad_x = (x_len >= 0 ? 2 : -2);
        y_overlap_grad_y = (y_len >= 0 ? 2 : -2);
        
        grad[0] = x_overlap_grad_x * y_overlap;     /* Grad w.r.t. r1->x */
        grad[1] = y_overlap_grad_y * x_overlap;     /* Grad w.r.t. r1->y */
        grad[2] = -grad[0];                         /* Grad w.r.t. r2->x */
        grad[3] = -grad[1];                         /* Grad w.r.t. r2->y */
    }
    
    return x_overlap * y_overlap;
}

lay_real_t lay_all_overlap_area(const int num_rects, 
                                const lay_coord_t* pos, const lay_extent_t* size, 
                                lay_real_t* grad) {
    int i, j;
    lay_real_t local_grad[4];
    lay_real_t sum = 0;
    for (i = 0; i < num_rects; ++i) {
        for (j = i+1; j < num_rects; ++j) {
            sum += lay_overlap_area(pos + 2 * i, size + 2 * i,
                                    pos + 2 * j, size + 2 * j, 
                                    local_grad);
            if (grad) {
                grad[2*i  ] += local_grad[0];
                grad[2*i+1] += local_grad[1];
                grad[2*j  ] += local_grad[2];
                grad[2*j+1] += local_grad[3];
            }
        }
    }
    
    return sum;
}

int lay_any_overlap(const int num_rects, 
                    const lay_coord_t* pos, const lay_extent_t* size, 
                    const int index) {
    int i;
    const lay_coord_t* this_pos;
    const lay_extent_t* this_size;
    
    assert(num_rects > 0 && pos && size && index >= 0 && index < num_rects);
    this_pos = pos + 2 * index;
    this_size = size + 2 * index;
    
    for (i = 0; i < num_rects; ++i) {
        if (i == index)
            continue;
        
        if (lay_overlap_area(this_pos, this_size, pos + 2 * i, size + 2 * i, NULL) > 0)
            return 1;
    }
    
    return 0;
}

/* Ah, sweet O(N^2). */
int lay_any_overlap_any(const int num_rects, 
                        const lay_coord_t* pos, const lay_extent_t* size) {
    int i;
    
    for (i = 0; i < num_rects; ++i) {
        if (lay_any_overlap(num_rects, pos, size, i))
            return 1;
    }
    
    return 0;
}
