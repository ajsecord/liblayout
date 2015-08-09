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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

/**
"Minimal" random number generator of Park and Miller.  
Returns a uniform random deviate between 0.0 and 1.0.  
Set or reset seed to any integer value (except the unlikely 
value MASK) to initialize the sequence; seed must not be altered 
between calls for successive deviates in a sequence.

The integer constants are chosen by Park and Miller very 
carefully, so don't change them.  MASK is only used to avoid
problems with a seed of zeero.

From Numerical Recipes in C, 2nd ed., pp. 279
*/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

static float rand0(long *seed) {
    long k;
    float result;
    
    *seed ^= MASK;                /* Avoid simple bit patterns */
    
    k = (*seed) / IQ;
    *seed = IA * (*seed - k * IQ) - IR * k;
    
    if (*seed < 0) *seed += IM;
    result = AM * (*seed);
    
    *seed ^= MASK;                /* Unmask */
    return result;
}



long rng_gen_random_seed() {
    long randSeed;
    size_t count;
    FILE* in;
    
    count = 0;
    in = fopen("/dev/urandom", "rb");
    if (in != NULL) {
        count = fread(&randSeed, sizeof(randSeed), 1, in);
        fclose(in);
    }
    
    if (count != 1) {
        randSeed = time(NULL);
    }
    
    return randSeed;
}

void rng_random_shuffle(int* array, const int N, long* randSeed) {
    int i, tmp, c;
    
    assert(array && randSeed);
    
    for(i = 0; i < N; i++)
        array[i] = i;
    
    for(i = 0; i < N-1; i++) {
        c = (int)(rand0(randSeed) * (N - i));
        assert(c != N - i);
        
        tmp = array[i]; 
        array[i] = array[i+c]; 
        array[i+c] = tmp;
    }

}


/*    A uniformly distributed deviate between 0.0 and 1.0.
*/
float rng_uniform_dev(long* seed) {
    return rand0(seed);
}


/*    An exponentially distributed deviate.
*/
float rng_exp_dev(long* seed, float a) {
    float result;
    
    do
        result = rand0(seed);
    while (result == 0.0);
    return -log(result)/a;
}


/*    A normally distributed deviate with zero mean and unit varience,
    using rand0 as the source of uniform deviates.

    From Numerical Recipes in C, 2nd ed., pp. 289
*/
float rng_gauss_dev(long *seed) {
    static int iset = 0;
    static float gset;
    float factor, r_squared, v1, v2;

    /* We need an extra random deviate, but we don't have one, 
       so find a suitable one.
    */
    if (iset == 0) {
        do {
            v1 = 2.0 * rand0(seed) - 1.0;
            v2 = 2.0 * rand0(seed) - 1.0;
            r_squared = v1 * v1 + v2 * v2;
        } while (r_squared >= 1.0 || r_squared == 0.0);

        factor = sqrt(-2.0 * log(r_squared) / r_squared);
    
        /* Now make the Box-Muller transformation to get two normal deviates.
           Return one and save one for next time.
        */
        gset = v1 * factor;
        iset = 1;
        return v2 * factor;
    
    } else {
        iset = 0;
        return gset;
    }
}

/*     Sobol sequence stolen from Numerical Recipes in C, 2nd ed., pp. 312
    Call once with n < 0 and then with 1 <= n <= 6 for each n-dimensional
    point.  Do not modify n without re-initialising. 
*/
#define MAXBIT 30
#define MAXDIM 6
static void sobseq(int *n, float x[]) {
    int j, k, l;
    unsigned long i, im, ipp;
    static float fac;
    static unsigned long in , ix[MAXDIM+1], *iu[MAXBIT+1];
    static unsigned long mdeg[MAXDIM+1]={0,1,2,3,3,4,4};
    static unsigned long ip[MAXDIM+1]={0,0,1,1,2,1,4};
    static unsigned long iv[MAXDIM*MAXBIT+1] = {
        0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};
    
    if (*n < 0) {
        for (k = 1; k <= MAXDIM; ++k) ix[k] = 0;    

        in = 0;
        if (iv[1] != 1) return;
        fac = 1.0/(1L << MAXBIT);
        for (j = 1, k = 0; j <= MAXBIT; j++, k += MAXDIM) iu[j] = &iv[k];
        
        for (k = 1; k <= MAXDIM; ++k) {
            for (j = 1; j <= mdeg[k]; j++) iu[j][k] <<= (MAXBIT-j);
            
            for (j = mdeg[k]+1; j <= MAXBIT; j++) {
                ipp = ip[k];
                i = iu[j-mdeg[k]][k];
                i ^= (i >> mdeg[k]);
                for (l = mdeg[k]-1; l >= 1; l--) {    
                    if (ipp & 1) i ^= iu[j-l][k];
                    ipp >>= 1;
                }
                iu[j][k] = i;
            }
        }

    } else {
        im = in++;
        for (j = 1; j <= MAXBIT; j++) {
            if (!(im & 1)) break;
            im >>= 1;
        }
        if (j > MAXBIT) { 
            fprintf(stderr, "MAXBIT too small in sobseq");        
            exit(-1);
        }
        im = (j-1) * MAXDIM;
        /*for (k = 1; k <= IMIN(*n, MAXDIM); k++) {*/
        for (k = 1; k <= (*n < MAXDIM ? *n : MAXDIM); k++) {
            ix[k] ^= iv[im+k];
            x[k] = ix[k] * fac;
        }
    }
}

/* Convert between a normal 0-based C array and a 1-based array, 
 * which sobseq() requires.  Ugh. 
 */
void rng_sobol(int *n, float x[]) {
    sobseq(n, x - 1);
}


/* Calculate the radical inverse of the number i the given base. 
 * From Alexander Keller's thesis "Quasi-Monte Carlo Methods for 
 * Photorealistic Image Synthesis," from the university 
 * "Vom Fachereich Informatik der Universitat Kaiserslautern"
 * Should be available from any good source concerning low-discrepancy
 * sequences.
 *
 * Note that this modifies the passed-in argument i, which is kind of 
 * bad form, but still safe.
 */
static double radical_inverse(const int base, int i) {
    double x = 0.0;
    double f = 1.0/base;

    while (i) {
        x += f * (double) (i % base);
        i /= base;
        f *= 1.0/base;
    }

    return x;
}

/* Calculate the radical inverse of the number i the given base. 
 * Fast incremental version that takes the radical inverse of i-1 as input.
 * From Alexander Keller's thesis "Quasi-Monte Carlo Methods for 
 * Photorealistic Image Synthesis," from the university 
 * "Vom Fachereich Informatik der Universitat Kaiserslautern"
 * Should be available from any good source concerning low-discrepancy
 * sequences.
 *
 * Note that this modifies the passed-in argument prev_inverse, which is 
 * kind of bad form, but still safe.
 */
static double radical_inverse_inc(const int base, double prev_inverse) {
    double h, hh, r;
    r = 1.0 - prev_inverse - 1e-10;

    if (1.0/base < r) {
        prev_inverse += 1.0/base;
    } else {
        h = 1.0/base;

        do {
            hh = h;
            h *= 1.0/base;
        } while (h >= r);

        prev_inverse += hh + h - 1.0;
    }

    return prev_inverse;        /* It got modified */
}    

/* Generate the Halton sequence in num_dim dimensions. 
 * Call once with negative num_dim to initialise, and then as many times as
 * desired to get additional terms in the sequence.  Do not change 
 * num_dim without re-initialising.
 *
 * For this implementation, num_dim <= 10, but trivially extended.
 */

/*
We follow the suggestion that the incremental version of the radical
inverse function should not be used for long runs of numbers, and so
we call the full version every FULL_VERSION_INTERVAL;  This interval is
a complete guess (although setting to 1000 gives bad results). 
*/
#define HALTON_MAX_DIM 10
#define HALTON_FULL_VERSION_INTERVAL 100
void rng_halton(int* num_dim, float x[]) {
    const int primes[HALTON_MAX_DIM] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 };
    static float last_vector[HALTON_MAX_DIM];
    static int i;
    int j;
    if (*num_dim < 0) {    
        i = 0;
    } else {
        assert(*num_dim <= HALTON_MAX_DIM);

        /* Use the full version of the radical inverse function */ 
        if ((i % HALTON_FULL_VERSION_INTERVAL) == 0) {
            for (j = 0; j < *num_dim; j++) {
                x[j] = radical_inverse(primes[j], i); 
                last_vector[j] = x[j];
            }

        /* Use the fast incremental version */
        } else {
            for (j = 0; j < *num_dim; j++) {
                x[j] = radical_inverse_inc(primes[j], last_vector[j]); 
                last_vector[j] = x[j];
            }
        }

        i++;
    }
}

/* Generate the Hammersley sequence of N points in num_dimensions.  Note 
 * that this sequence is not incremental: changing N generates a completely
 * different set of points.  See the comments for halton(). 
 * The array x must have num_dim * N floats.
 *
 * The one-dimensional Hammersley sequence is not very useful.
 *
 */
#define HAMMER_MAX_DIM 11
#define HAMMER_FULL_VERSION_INTERVAL 100
void rng_hammersley(const int num_dim, const int N, float x[]) {
    const int primes[HAMMER_MAX_DIM-1] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 };
    float last_vector[HAMMER_MAX_DIM];
    int i, j, index;

    assert(num_dim <= HAMMER_MAX_DIM);

    for (i = 0; i < N; i++) {
        index = i * num_dim;

        /* Use the full version of the radical inverse function */ 
        if ((i % HAMMER_FULL_VERSION_INTERVAL) == 0) {
            x[index] = i / (float)N;

            for (j = 1; j < num_dim; j++) {
                x[index + j] = radical_inverse(primes[j-1], i); 
                last_vector[j] = x[index + j];
            }

        /* Use the fast incremental version */
        } else {
            x[index] = i / (float)N;

            for (j = 1; j < num_dim; j++) {
                x[index + j] = radical_inverse_inc(primes[j-1], last_vector[j]);
                last_vector[j] = x[index + j];
            }
        }
    }
}

/* Generate the Halton sequence in num_dim dimensions. 
 *
 * For this implementation, num_dim <= 10, but trivially extended.
 *
 * Since for the jittering we need to determine the resolution of the 
 * grid that the points are selected on, jittered_halton cannot be used
 * incrementally.
 *
 * Each point of the normal Halton sequence is at one corner of a box of 
 * dimensions (1/s_1, 1/s_2, ...) where 
 * s_i = primes[i]^m_i and primes[i]^(m_i - 1) < N < primes[i]^m_i.  We 
 * jitter each Halton point so that it doesn't leave the box.  The 
 * magnitude of the jittering within that box seems to be arbitrary, 
 * however, and Keller suggests limiting it to 1/N.  Jittering with magnitude
 * one within the box would destroy any minimum-distance properties of the 
 * sequence, because two samples could be jittered next to each other on 
 * either side of a box wall.  
 *
 * We follow the suggestion that the incremental version of the radical
 * inverse function should not be used for long runs of numbers, and so
 * we call the full version every FULL_VERSION_INTERVAL;  This interval is
 * a complete guess.
 */
void rng_jittered_halton(const int num_dim, const int N, const float jitterMag, 
                    long* seed, float x[]) {
    const int primes[HALTON_MAX_DIM] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 };
    float last_vector[HALTON_MAX_DIM];
    float jitter[HALTON_MAX_DIM];
    int i, j, index;
    assert(num_dim <= HALTON_MAX_DIM);
    assert(jitterMag >= 0.0 && jitterMag <= 1.0);

    /* Determine the resolution of the underlying grid for each dimension.
     * The resolution is 1 / base^m, where base^(m-1) < N <= base^m.
     */
    for (i = 0; i < num_dim; i++) {
        j = (int)(log(N)/log(primes[i]));
        if (pow(primes[i], j) != (double)N) j++;    /* Hack if N exact power */
        jitter[i] = jitterMag / pow(primes[i], j);
    }

    for (i = 0; i < N; i++) {
        index = i * num_dim;

        /* Use the full version of the radical inverse function */ 
        if ((i % HALTON_FULL_VERSION_INTERVAL) == 0) {
            for (j = 0; j < num_dim; j++) {
                last_vector[j] = radical_inverse(primes[j], i);
                x[index + j] = last_vector[j] + rng_uniform_dev(seed) * jitter[j];
            }

        /* Use the fast incremental version */
        } else {
            for (j = 0; j < num_dim; j++) {
                x[index + j] = radical_inverse_inc(primes[j], last_vector[j]);
                last_vector[j] = x[index + j];
                x[index + j] += rng_uniform_dev(seed) * jitter[j];
            }
        }
    }
}

/* Generate the Hammersley sequence of N points in num_dimensions.  Note 
 * that this sequence is not incremental: changing N generates a completely
 * different set of points.  See the comments for halton(). 
 * The array x must have num_dim * N floats.
 *
 * Jitter the resulting sequence identically to jittered_halton().
 *
 * The one-dimensional Hammersley sequence is not very useful.
 *
 */
void rng_jittered_hammersley(const int num_dim, const int N, const float jitterMag, 
                            long* seed, float x[]) {
    const int primes[HALTON_MAX_DIM] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 };
    float last_vector[HALTON_MAX_DIM];
    float jitter[HALTON_MAX_DIM];
    int i, j, index;
    assert(num_dim <= HALTON_MAX_DIM);
    assert(jitterMag >= 0.0 && jitterMag <= 1.0);

    /* Determine the resolution of the underlying grid for each dimension.
     * The resolution is 1 / base^m, where base^(m-1) < N <= base^m.
     */
    jitter[0] = jitterMag / N;
    for (i = 1; i < num_dim; i++) {
        j = (int)(log(N)/log(primes[i-1]));
        if (pow(primes[i-1], j) != (double)N) j++;    /* Hack if N exact power */
        jitter[i] = jitterMag / pow(primes[i-1], j);
    }

    for (i = 0; i < N; i++) {
        index = i * num_dim;

        x[index] = i / (float)N + rng_uniform_dev(seed) * jitter[0];

        /* Use the full version of the radical inverse function */ 
        if ((i % HAMMER_FULL_VERSION_INTERVAL) == 0) {
            for (j = 1; j < num_dim; j++) {
                last_vector[j] = radical_inverse(primes[j-1], i);  
                x[index + j] = last_vector[j] + rng_uniform_dev(seed) * jitter[j];
            }

        /* Use the fast incremental version */
        } else {
            for (j = 1; j < num_dim; j++) {
                x[index + j] = radical_inverse_inc(primes[j-1], last_vector[j]);
                last_vector[j] = x[index + j];
                x[index + j] += rng_uniform_dev(seed) * jitter[j];
            }
        }
    }
}

