/*   macopt library header file         release 1.1          gradient-based optimizer

     Copyright   (c) 2002   David J.C. MacKay

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
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    GNU licenses are here :
    http://www.gnu.org/licenses/licenses.html

    Author contact details are here :
    http://www.inference.phy.cam.ac.uk/mackay/c/macopt.html       mackay@mrao.cam.ac.uk

    If you find macopt useful, please feel free to make a donation to
    support David MacKay's research group.
*/
/* structure for macopt */
typedef struct {
  double tol ;    /* convergence declared when the gradient vector is smaller
		     in magnitude than this, or when the mean absolute 
		     step is less than this (see above) */
  double grad_tol_tiny ; /* if gradient is less than this, we definitely 
			    stop, even if we are not using a gradient 
			    tolerance */
  double step_tol_tiny ; /* if step is less than this, we stop, even if 
			    we are not using a step tolerance */
  int end_if_small_step ; /* defines the role of tol -- alternative is
			     end_on_small_grad */
  int its ;               /* number of its */
  int itmax ;             /* max */
  int rich ; /* whether to do the extra gradient evaluation at the beginning 
	      of each new line min */
  int verbose ; 
  double stepmax ;        /* largest step permitted (not used in macopt) */

  int linmin_maxits ;     /* in maclinmin */
  double linmin_g1 ;      /* factors for growing and shrinking the interval */
  double linmin_g2 ;
  double linmin_g3 ;
  double lastx     ;      /* keeps track of typical step length */
  double lastx_default ;  /* if maclinmin is reset, lastx is set to this */


  int  do_newitfunc ;  /* whether to run newitfunc each new iteration */
  void (*newitfunc)(double *,int, void *) ; /* this function might for example
					print the current state of the
					simulation */
  void *newitfuncarg ; 

/* These should not be touched by the user. They are handy pointers for macopt
   to use 
*/
  double gtyp ; /* stores the rms gradient for linmin */
  double *pt , *gx , *gy , *gunused ;
  double *xi , *g , *h ;
  double *m ;  /* the metric, used in macoptIIc */
  /* the user is responsible for allocating this and setting it */
  double *mg ;  /* the natural gradient */
  /* g, xi is covariant,  and h and mg are covariant, as is x */
  int metric ; /* whether we are using the metric */
  int n ;                 /* dimension of parameter space */
  int restart ;           /* whether to restart macopt - fresh cg directions */
  /* this is only set to 1 by maclinmin or macopt */
} macopt_args ; 


/* lastx :--- 1.0 might make general sense, (cf N.R.)
				  but the best setting of all is to have 
				  a prior idea of the eigenvalues. If 
				  the objective function is equal to sum of N
				  terms then set this to 1/N, for example 
				  Err on the small side to be conservative. */
int macoptIIc
  (double *,            /* starting vector                                */
   int    ,             /* number of dimensions                           */
   void   (*dfunc)(double *,double *, void *), /* evaluates the gradient   */
   void   *, 
   macopt_args *
   )  ; /* returns a status variable: 1/2 = normal end states 
	   0 = over ran itmax */

void macoptII
  (double *,            /* starting vector                                */
   int    ,             /* number of dimensions                           */
   void   (*dfunc)(double *,double *, void *), /* evaluates the gradient   */
   void   *, 
   macopt_args *
   )  ; 

void maccheckgradc(double *, int, double, 
	       double (*f)(double *, void *), void *,
	       void (*g)(double *,double *, void *), void * , int );	


void maccheckgrad(double *, int, double, 
	       double (*f)(double *, void *), void *,
	       void (*g)(double *,double *, void *), void * , int );	

void macopt_defaults ( macopt_args * ) ;

void macopt_allocate_metric ( macopt_args * , int ) ;
void macopt_allocate ( macopt_args * , int ) ;
void macopt_free ( macopt_args * ) ;

/* the following functions could be declared static within macopt.c */

double maclinminII  /* used by both macoptII and macoptIIc */
(
 double * ,
 void   (*dfunc)(double *,double *, void *), /* evaluates the gradient */
 void   * ,
 macopt_args * ) ;

double macprodII 
( 
 double * , double * , double  ,
 void   (*dfunc)(double *,double *, void *), 
 void   * , 
 macopt_args *
) ;

void macopt_restart ( macopt_args * , int ) ;

void    evaluate_hessian 
( double ** ,
  double * , 
  int  ,     
  double  ,
  void   (*dfunc)(double *,double *, void *), 
  void   *dfunc_arg, 
  int  ) ;

void    evaluate_hessianc 
( double ** ,
  double * , 
  int  ,     
  double  ,
  void   (*dfunc)(double *,double *, void *), 
  void   *dfunc_arg, 
  int  ) ;

