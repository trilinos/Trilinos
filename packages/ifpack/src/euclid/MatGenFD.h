/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef MATGENFD_DH_DH
#define MATGENFD_DH_DH

/*=====================================================================
option summary:
---------------
processor topology
     -px <int> -py <int> -pz <int>
     defaults:  -px 1 -py 1 -pz 0

grid topology
  -m <int>
  if pz=0, each processor has a square grid of dimension m*m,
  hence there are m*m*px*py unknowns.
  if pz > 0, each local grid is of dimension m*m*m, hence
  there are m*m*m*px*py*pz unknowns.


diffusion coefficients (default is 1.0):
    -dx <double> -dy <double> -dz <double>

convection coefficients (default is 0.0)
    -cx <double> -cy <double> -cz <double>

grid dimension; if more than one mpi process, this is
the local size for each processor:
     -m <int>

boundary conditions:
  This is very primitive; boundary conditions can only be generated for
  2D grids; the condition along each side is either dirichlet (constant),
  if bcXX >= 0, or neuman, if bcXX < 0.

   -bcx1 <double>
   -bcx2 <double>
   -bcy1 <double>
   -bcy2 <double>

Misc.
     -debug_matgen
     -striped (may not work?)
=====================================================================*/


#include "euclid_common.h"

#ifdef __cplusplus
extern "C" {
#endif

struct _matgenfd {
  bool allocateMem; 
        /* If true, memory is allocated when run() is called, in which case
         * the caller is responsible for calling FREE_DH for the rp, cval,
         * aval, and rhs arrays.  If false, caller is assumed to have
         * allocated memory when run is called.  
         * Default is "true"
         */
  int px, py, pz;  /* Processor graph dimensions */
  bool threeD;  
  int m;           /* number of matrix rows in local matrix */
  int cc;          /* Dimension of each processor's subgrid */
  double hh;       /* Grid spacing; this is constant,  equal to 1.0/(px*cc-1) */
  int id;          /* the processor whose submatrix is to be generated */
  int np;          /* number of subdomains (processors, mpi tasks) */
  double stencil[8];


  /* derivative coefficients; a,b,c are 2nd derivatives, 
   * c,d,e are 1st derivatives; f,g,h not currently used.
   */
  double a, b, c, d, e, f, g, h;

  int first; /* global number of first locally owned row */
  bool debug;

  /* boundary conditions; if value is < 0, neumen; else, dirichelet */
  double bcX1, bcX2;
  double bcY1, bcY2;
  double bcZ1, bcZ2;
                
  /* The following return coefficients; default is konstant() */
  double (*A)(double coeff, double x, double y, double z);
  double (*B)(double coeff, double x, double y, double z);
  double (*C)(double coeff, double x, double y, double z);
  double (*D)(double coeff, double x, double y, double z);
  double (*E)(double coeff, double x, double y, double z);
  double (*F)(double coeff, double x, double y, double z);
  double (*G)(double coeff, double x, double y, double z);
  double (*H)(double coeff, double x, double y, double z);
};

extern void MatGenFD_Create(MatGenFD *mg);
extern void MatGenFD_Destroy(MatGenFD mg);
extern void MatGenFD_Run(MatGenFD mg, int id, int np, Mat_dh *A, Vec_dh *rhs);

 /* =========== coefficient functions ============== */
extern double konstant(double coeff, double x, double y, double z);
extern double e2_xy(double coeff, double x, double y, double z);



/* 3 boxes nested inside the unit square domain.
   diffusivity constants are: -dd1, -dd2, -dd3.
*/
/* box placement */
#define BOX1_X1 0.1
#define BOX1_X2 0.4
#define BOX1_Y1 0.1
#define BOX1_Y2 0.4

#define BOX2_X1 0.6
#define BOX2_X2 0.9
#define BOX2_Y1 0.1
#define BOX2_Y2 0.4

#define BOX3_X1 0.2
#define BOX3_X2 0.8
#define BOX3_Y1 0.6
#define BOX3_Y2 0.8

/* default diffusivity */
#define BOX1_DD  10
#define BOX2_DD  100
#define BOX3_DD  50

extern double box_1(double coeff, double x, double y, double z);
  /* -bd2 is diffusion coeff outside box;
     -bd1 is diffusion coeff inside box.
  */
     


extern double box_2(double coeff, double x, double y, double z);

#ifdef __cplusplus
}
#endif
#endif
