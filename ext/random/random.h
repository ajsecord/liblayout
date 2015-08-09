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

#ifndef AJS_RANDOM_H
#define AJS_RANDOM_H

#ifdef __cplusplus
extern "C" {
#endif
    
/** Generate a random seed.  Uses good seed values from /dev/urandom if that
    exists, otherwise the current time.
*/
long rng_gen_random_seed();
    
/** Generate a random shuffle of the integers from 0 to N-1 in the passed array.
    From the comp.lang.c FAQ section 13.19, which is based on 
    Knuth Sec. 3.4.2 pp. 137-8.
*/
void rng_random_shuffle(int* array, const int N, long* randSeed);
    
/** Generate a uniform random deviate in [0,1] */
float rng_uniform_dev(long *seed);
    
/** Generate an exponentially-distributed random deviate. */
float rng_exp_dev(long* seed, float a);

/** Generate a random deviate taken from the uniform Gaussian distribution 
    (mean zero, standard deviation one). 
*/
float rng_gauss_dev(long *seed);

/** Sobol sequence in n dimensions from Numerical Recipes in C, 2nd ed., pp. 312.
    
    Call once with n < 0 and then with 1 <= n <= 6 for each n-dimensional
    point.  Do not modify n without re-initialising. 
*/
void rng_sobol(int *n, float x[]);

/** Generate the Halton sequence in num_dim dimensions. 
    
    Call once with negative num_dim to initialise, and then as many times as
    desired to get additional terms in the sequence.  Do not change 
    num_dim without re-initialising.

    For this implementation, num_dim <= 10, but trivially extended. 
*/
void rng_halton(int *n, float x[]);

/** Generate the Hammersley sequence of N points in num_dimensions.  Note 
    that this sequence is not incremental: changing N generates a completely
    different set of points.  See the comments for halton(). 
    The array x must have num_dim * N floats.

    The one-dimensional Hammersley sequence is not very useful. 
*/
void rng_hammersley(const int num_dim, const int N, float x[]);

void rng_jittered_halton(const int num_dim, const int N, const float jitterMag, 
                         long* seed, float x[]);

void rng_jittered_hammersley(const int num_dim, const int N, const float jitterMag, 
                             long* seed, float x[]);

#ifdef __cplusplus
}
#endif

#endif
