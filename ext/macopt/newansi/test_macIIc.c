/* test_macII.c           example executable that minimizes a quadratic function using macopt

   (c) 2002 David J.C. MacKay
  
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    GNU licenses are here :
    http://www.gnu.org/licenses/licenses.html

    Author contact details are here :
    http://www.inference.phy.cam.ac.uk/mackay/c/macopt.html       mackay@mrao.cam.ac.uk
*/
#include "../newansi/r.h" 
/* #include "../newansi/mynr.h" */
#include "../newansi/macopt.h"

/* 
   test program for macopt solution of equation A x = b. 

   this equation is solved by minimizing the function 1/2 xAx - bx
*/
   
#include "test.h"
void state_printer ( double *w, int its, void *arg ) ;

void main(int argc, char *argv[])
{
  gq_args param;
  macopt_args a ;
  double *x ;
  int n , status ;
  double epsilon=0.001 ;

  /* Load up the parameters of the function that you want to optimize */
  printf("============================================================\n");
  printf("= Demonstration program for macoptIIc                      =\n");
  printf("= Solves A x = b by minimizing the function 1/2 xAx - bx   =\n");
  printf("= A must be positive definite (e.g. 2 1 1 2)               =\n");
  printf("\n  Dimension of A (eg 2)?\n");
  inputi(&(param.n));
  n=param.n; 
  param.A=dmatrix(1,n,1,n);
  param.b=dvector(1,n);
  /* the metric! : */
  a.m=dvector(1,n);
  x=dvector(1,n);
  typeindmatrix(param.A,1,n,1,n);
  printf("  b vector?\n");
  typeindvector(param.b,1,n);
  printf("  Initial condition x?\n");
  typeindvector(x,1,n);
  printf("  Metric m?\n");
  typeindvector(a.m,1,n);

  /* Check that the gradient_function is the gradient of the function  */
  /* You don't have to do this, but it is a good idea when debugging ! */
  maccheckgrad (  x , param.n , epsilon , 
		quadratic , (void *)(&param) , 
		vgrad_quadratic , (void *)(&param) , 
		0
		) ;

  /* initialize the arguments of the optimizer */
  macopt_defaults ( &a ) ; 

  /* modify macopt parameters from their default values */
  a.do_newitfunc = 1 ; /* this means that I want to have an auxiliary
			  subroutine executed each iteration of the 
			  optimization */
  a.newitfunc = &state_printer ; /* setting the function to be performed to
				  T_return */
  a.newitfuncarg = (void *)(&param) ;
  a.metric = 1 ; /* we have put in the metric and want to use it */
  a.verbose = 2 ; /* verbosity */
  a.rich = 0 ; /* verbosity */
  
  /* Do an optimization */
  status = macoptIIc ( x , param.n , 
	    vgrad_quadratic , (void *)(&param) , &a
	    ) ;

  printf("(%d) Solution:\n", status);
  quadratic(x,&param);
}

void state_printer ( double *w, int its, void *arg )
                     /* this routine is called by the
			macopt optimizer ; it reports parameter
			values during the optimization */
{
  int i ; 
  gq_args *param = ( gq_args * ) arg ;

  printf ( "\n----------------------\n");
  printf ( "      state_printer:\n");
  printf ( "Line minimization %d is starting, and the latest params are\n" , its );
  for ( i = 1 ; i <= param->n ; i++ ) printf ( "%g " , w[i] ) ;
  printf ( "\n" ) ;
  printf ( "----------------------\n");

}
