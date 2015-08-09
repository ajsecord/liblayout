/* test_function.c       example quadratic function for use with macopt

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

#include <stdio.h> 
#include <math.h> 
#include "test.h"

double quadratic(double *x,   void *arg )
{
  int i,j;
  double f=0.0,g;
  gq_args *param = ( gq_args * ) arg ;
 
  printf ( "quadr. at " ) ;
  for ( i = 1 ; i <= param->n ; i++ ) printf( "%g " , x[i] ) ;
  for ( i = 1 ; i <= param->n ; i++ ) {
    g=0.0;
    for(j=1;j<=param->n;j++){
      g += param->A[i][j] * x[j];
    }
    f+=x[i]*(g/2.0 - param->b[i]);
  }
  printf(" : %g\n",f);
  return(f); 
}

double quartic(double *x,   void *arg )
{
  int i,j;
  double f=0.0,g;
  gq_args *param = ( gq_args * ) arg ;
 
  printf ( "quartic at " ) ;
  for ( i = 1 ; i <= param->n ; i++ ) printf( "%g " , x[i] ) ;
  for ( i = 1 ; i <= param->n ; i++ ) {
    g=0.0;
    for(j=1;j<=param->n;j++){
      g += param->A[i][j] * x[j];
    }
    f+=x[i]*(g/2.0 - param->b[i]);
    f += 0.25 * param->q4 * x[i]*x[i]*x[i]*x[i] ; 
    f +=  param->q4 * exp(x[i]) ; 
  }
  printf(" : %g\n",f);
  return(f); 
}

double grad_quadratic(double *x, double *g, void *arg )
{
  int i,j;
  double f=0.0;
  gq_args *param = ( gq_args * ) arg ;
  printf ( "grad_q at " ) ;
  for ( i = 1 ; i <= param->n ; i++ ) printf ( "%g " , x[i] ) ;
  for ( i = 1 ; i <= param->n ; i++ ) {
    g[i] = 0.0 ;
    for ( j = 1 ; j <= param->n ; j++ ) {
      g[i] += param->A[i][j] * x[j] ;
    }
    f += x[i] * ( g[i] / 2.0 - param->b[i] ) ;
    g[i] -= param->b[i] ;
  }
  printf ( " :: %g :: " , f ) ;
  for ( i = 1 ; i <= param->n ; i++ ) printf ( "%g " , g[i] ) ;
  printf ( "\n" ) ;
  return ( f ) ;
}


void vgrad_quadratic(double *x, double *g, void *arg )
{
  int i,j;
  double f=0.0;
  gq_args *param = ( gq_args * ) arg ;
  printf ( "grad_q at " ) ;
  for ( i = 1 ; i <= param->n ; i++ ) printf ( "%g " , x[i] ) ;
  for ( i = 1 ; i <= param->n ; i++ ) {
    g[i] = 0.0 ;
    for ( j = 1 ; j <= param->n ; j++ ) {
      g[i] += param->A[i][j] * x[j] ;
    }
    f += x[i] * ( g[i] / 2.0 - param->b[i] ) ;
    g[i] -= param->b[i] ;
  }
  printf ( " :: %g :: " , f ) ;
  for ( i = 1 ; i <= param->n ; i++ ) printf ( "%g " , g[i] ) ;
  printf ( "\n" ) ;
}

void vgrad_quartic(double *x, double *g, void *arg )
{
  int i,j;
  double f=0.0;
  gq_args *param = ( gq_args * ) arg ;
  printf ( "grad_quartic at " ) ;
  for ( i = 1 ; i <= param->n ; i++ ) printf ( "%g " , x[i] ) ;
  for ( i = 1 ; i <= param->n ; i++ ) {
    g[i] = 0.0 ;
    for ( j = 1 ; j <= param->n ; j++ ) {
      g[i] += param->A[i][j] * x[j] ;
    }
    f += x[i] * ( g[i] / 2.0 - param->b[i] ) ;
    f += 0.25 * param->q4 * x[i]*x[i]*x[i]*x[i] ; 
    f +=  param->q4 * exp(x[i]) ; 
    g[i] -= param->b[i] ;
    g[i] +=  param->q4 * x[i]*x[i]*x[i] ; 
    g[i] +=  param->q4 * exp(x[i]) ; 
  }
  printf ( " :: %g :: " , f ) ;
  for ( i = 1 ; i <= param->n ; i++ ) printf ( "%g " , g[i] ) ;
  printf ( "\n" ) ;
}

void Atimesh( double *x , double *h , double *v , void *arg)
{ /* This ignores the first argument */
  int i,j;
  gq_args *param = ( gq_args * ) arg ;

  for ( i = 1 ; i <= param->n ; i++ ) {
    v[i] = 0 ;
    for ( j = 1 ; j <= param->n ; j++ ) {
      v[i] += param->A[i][j] * h[j] ;
    }
  }
}
