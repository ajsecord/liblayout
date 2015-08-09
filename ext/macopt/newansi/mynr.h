/*   mynr.h         header for various subroutines           release 1.1    

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
    http://www.inference.phy.cam.ac.uk/mackay/
*/
/* for cg.c */
typedef struct {
  int n;
  double *p,*xi,*xt;
  double (*nfunc)(double *,void *);
  void *func_arg;
} f1dim_arg;

typedef struct {
  int n;
  double *p,*xi,*xt;
  void (*ndfunc)(double * , double * , void *);
  void *func_arg;
} df1dim_arg;

typedef struct {
  int n;
  double **m;
  double **in;
  double **lu;
  int *indx;
  double det;
  int in_allocated ; 
  int in_computed ; 
  int lu_computed ; 
} dmatrix_family;


typedef struct {
  int n;
  double *b;
  void (*Afunc)(  double * ,  double * ,double * , void * ) ; 
  void   *Afunc_arg;
} cggq_args;       

void linmin(double *, double *, int n, double *, double (*f)(double *,void *), 
	    void *, double);	
/* void dlinmin();	*/
double brent(double, double, double, double (*f)(double, void *), void *, 
	     double, double *);	
double dbrent(double, double, double, double (*f)(double), 
	      double (*g)(double),
	      double, double *);	
void frprmn(double *, int, double, double, int *, int,
	    double *, double (*f)(double *, void *), void *,
	    void (*g)(double *,double *, void *), void *);	
void checkgrad(double *, int, double, 
	       double (*f)(double *, void *), void *,
	       void (*g)(double *,double *, void *), void *);	
void mnbrak(double *,double *,double *,double *,double *,double *,
	    double (*f)(double,void *),void *);
double f1dim(double, void *);	
/* double df1dim(); */	
void dfpmin(double *, int , double , int *, double *, 
	    double (*f)(double *,void *), 
	    void *, void (*g)(double *,double *,void *), 
	    void *, double **, double);

void macopt
  (double *,
   int    ,
   int ,
   double ,
   int   * ,
   int    ,
   void   (*dfunc)(double *,double *, void *), 
   void   * );
double maclinmin 
(
 double * , double * , int  ,    
 void   (*dfunc)(double *,double *, void *), 
 void   * , double * ,  double * ,  double * );
double macprod 
( 
 double * , double * , double * , double * , double  , int  ,
 void   (*dfunc)(double *,double *, void *), 
 void   * ) ;


void light_speed_cg
  (double *,   int,   int,  double,   int    *,   int ,
   double *,   double (*dfunc)(double *,double *, void *),
   void   *,   void   (*Afunc)(double *,double *,double *, void *),void   *);
void cg_solve  
( double *,  
 int  , 
 double * , 
 void (*f)(double *,double *,double *, void *) , 
 void * , 
 double  ) ;
double cg_grad_quadratic
( double * , double * , void * ) ;

/* for lu */
int allocate_dmatrix_family ( dmatrix_family *s , int n , int control ) ;
int free_dmatrix_family ( dmatrix_family *s ) ;
int lu_invert ( dmatrix_family *s , double *v ) ; 
void show_dmatrix_family ( dmatrix_family *s , FILE * fp ) ;
int write_lu ( dmatrix_family *s , char *file ) ;
int readinlumatrix ( dmatrix_family *s , char *file ) ;
double lu_quadratic_form ( double *v1, dmatrix_family *s ,double *v2 ) ;
int	invert_dmatrix_family( dmatrix_family *mp , int control ) ;
double	lumatrixproduct(double *v1,double **m,double *v2,int n,int *indx) ;
void	symmetrize_dmatrix(double **, int );

int 	ludcmp( double ** , int , int * , double * );
void 	lubksb( double ** , int , int * , double * );

void	find_eigs( double **, int , double *);
void tqli ( double * , double * , int , double ** ) ;
void tred2 ( double ** , int , double * , double * ) ;

/* for matrix.c */
void dmatrixfromdmatrix (double **b,int l1,int h1,int l2,int h2,double **f) ;
void dvectorfromdvector(double *w,int lo,int hi,double *w2) ;

/*
double  evaluate_determinant();
double	quadratic_form();
double	lumatrixproduct();
double	matrixproduct();
double  trace();
void 	find_eigs();
int 	clip_eigs();
void 	det_and_tr_from_eigs();
void 	tred2();
void 	tqli();
double	invert_dmatrix();
double	luinvert_dmatrix();
*/
/* for eig.c */
/*
void	report_eigs();
void	report_eigs2();
void	scale_eigs_by_pow_sqrt_lambda();
void	pos_eigs();
 */

