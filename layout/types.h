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

#ifndef LAY_TYPES_H
#define LAY_TYPES_H

/** \file layout/types.h
* Basic type definitions.
*/

#ifdef __cplusplus
extern "C" {
#endif
    
    /** \name Basic data types */
    /*@{*/
    
    /** Define to use floats as the floating-point type. */
    #define LAY_REAL_IS_FLOAT
    
    /** Define to use doubles as the floating-point type. */
    /* #define LAY_REAL_IS_DOUBLE */
    
    /** Define to generate integer coordinates. */
    /*#define LAY_USE_INTEGER_COORDS*/ 
    
    /** Define to generate floating-point coordinates. */
    #define LAY_USE_REAL_COORDS 
    
    
#if defined(LAY_REAL_IS_FLOAT)
    
    /** The floating-point type, when needed. */
    typedef float lay_real_t;
    
    /** The absolute value of a real quantity. */
    #define LAY_REAL_ABS(x) fabsf(x) 
    
#elif defined(LAY_REAL_IS_DOUBLE)
    
    /** The floating-point type, when needed. */
    typedef double lay_real_t;
    
    /** The absolute value of a real quantity. */
    #define LAY_REAL_ABS(x) fabs(x) 
    
#endif
    
#if defined(LAY_USE_INTEGER_COORDS)
    
    /** The basic data type for a position coordinate.  This should be a signed type. */
    typedef int lay_coord_t;
        
    /** The basic data type for a width/height value.   This should be a signed type. */
    typedef int lay_extent_t;
    
    /** The absolute value of a coordinate. */
    #define LAY_COORD_ABS(x) abs(x)
    
    /** The absolute value of an extent. */
    #define LAY_EXTENT_ABS(x) abs(x)
    
#elif defined(LAY_USE_REAL_COORDS)
    
    /** The basic data type for a position coordinate.  This should be a signed type. */
    typedef lay_real_t lay_coord_t;
    
    /** The basic data type for a width/height value.   This should be a signed type. */
    typedef lay_real_t lay_extent_t;
    
    /** The absolute value of a coordinate. */
    #define LAY_COORD_ABS(x) fabs(x)
    
    /** The absolute value of an extent. */
    #define LAY_EXTENT_ABS(x) fabs(x)
#endif
    
    /*@}*/
    
#ifdef __cplusplus
}
#endif

#endif
