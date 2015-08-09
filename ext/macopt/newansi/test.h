/* test.h           example header file for a quadratic function for use with macopt

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
double quadratic ( double * , void * ) ;
double grad_quadratic ( double *,double *, void *);
void vgrad_quadratic ( double *,double *, void *);
double quartic ( double * , void * ) ;
double grad_quartic ( double *,double *, void *);
void vgrad_quartic ( double *,double *, void *);
void Atimesh ( double *,double *,double *, void *); 
typedef struct {
  int n;
  double **A, *b; /* Used to define f(x) = 1/2 xAx - bx */
  double q4 ; /* coefficient of quartic term if function is quartic */
} gq_args;        /* and g = Ax - b */
