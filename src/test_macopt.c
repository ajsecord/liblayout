#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "macopt.h"
#include "nrutil.h"

double energy(double* x, void* context) {
    assert(x);
    return (x[1] - 3.1415) * (x[1] - 3.1415) + (x[2] + 1) * (x[2] + 1);
}

void vgrad_energy(double* x, double* grad, void* context) {
    assert(x && grad);
    grad[1] = 2 * (x[1] - 3.1415);
    grad[2] = 2 * (x[2] + 1);
}

void print(double* w, int its, void* context) {
    assert(w);
    printf("Iter %i: %g %g\n", its, w[1], w[2]);
}

int main() {
    int dim = 2;
    macopt_args args;
    double* x;
    const double epsilon = 0.001;
    
    x = dvector(1, 2);
    x[1] = -1;
    x[2] = 2;
    
    /* Check that the gradient_function is the gradient of the function  */
    /* You don't have to do this, but it is a good idea when debugging ! */
    maccheckgrad(x, dim, epsilon, energy, NULL, vgrad_energy, NULL, 0) ;
    
    /* initialize the arguments of the optimizer */
    macopt_defaults ( &args ) ; 
    
    /* modify macopt parameters from their default values */
    args.do_newitfunc = 1 ; /* this means that I want to have an auxiliary
        subroutine executed each iteration of the 
        optimization */
    args.newitfunc = &print ; /* setting the function to be performed to
        T_return */
    args.newitfuncarg = NULL;
    args.verbose = 2 ; 
    
    /* Do an optimization */
    macoptII(x, dim, vgrad_energy, NULL, &args) ;
    
    printf("Solution: %g %g: %g\n", x[1], x[2], energy(x, NULL));
    
    free_dvector(x, 1, 2);
    
    return 0;
}