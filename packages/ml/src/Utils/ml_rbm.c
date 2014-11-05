/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* ******************************************************************** */
/* compute null space from rigid body modes for smoothed aggregation    */
/* -------------------------------------------------------------------- */
/* Original Author : Ray Tuminaro (SNL)                                 */
/* Modified by: Irina Kalashnikova (SNL) to handle case when elasticity */
/* equations are coupled to NscalarDof scalar equations                 */
/*                                                                      */
/* ******************************************************************** */


/* ***************************************************************************
  Input
    Nnodes     - number of nodes in the amalgamated system (i.e., #nodes in physical mesh)
    x,y,z      - coordinates
    rbm        - allocated vector to hold rigid body modes
    Ndof       - total # PDEs (dofs) per node
    NscalarDof - # scalar equations at a node coupled to  elasticity equations (so that there are (Ndof - NscalarDof) elasticity equations)

  Output
    rbm    - vector populated with rigid body modes (array of size Nnodes*(nRBM + NscalarDof)*(Ndof + NscalarDof)
             where nRBM is the number of rigid body modes (6 for 3D elasticity, 3 for 2D elasticity, 1 for 1D elasticity)

    Example (3D shell elements for elasticity + 1 scalar equation per node): local 7 x 7 matrix:

                              translations (elasticity)     scalar   | rotations around (elasticity)
                                      x       y       z      eqn#1   |        x        y         z
                              -----------------------------------------------------------------------
        x-direction (elasticity)      1       0       0        0     |     0      z-zhat    yhat-y
        y-direction (elasticity)      0       1       0        0     |  zhat-z      0       x-xhat
        z-direction (elasticity)      0       0       1        0     |  y-yhat    xhat-x        0
        scalar equation #1            0       0       0        1     |    0        0            0
        x-rot (elasticity)            0       0       0        0     |    1        0            0
        y-rot (elasticity)            0       0       0        0     |    0        1            0
        z-rot (elasticity)            0       0       0        0     |    0        0            1

    3D elasticity with bricks would be the same with the last 3 rows removed. 2D elasticity would
    also remove the 3rd row and columns 3, 5 and 6. 1D elasticty would remove all but rows 1 and
    4 and columns 1 and 4.

******************************************************************************* */
#include <stdio.h>
#include <stdlib.h>
#include "ml_rbm.h"

int ML_Coord2RBM(int Nnodes, double x[], double y[], double z[], double rbm[], int Ndof, int NscalarDof)
{
   int vec_leng, ii, jj, offset, node, dof;
   vec_leng = Nnodes*Ndof;

   for( node = 0 ; node < Nnodes; node++ )
   {
      dof = node*Ndof;
      switch( Ndof - NscalarDof )
      {
         case 6:
            for(ii=3;ii<6+NscalarDof;ii++){ /* lower half = [ 0 I ] */
              for(jj=0;jj<6+NscalarDof;jj++){
                offset = dof+ii+jj*vec_leng;
                rbm[offset] = (ii==jj) ? 1.0 : 0.0;
              }
            }
	    /* There is no break here and that is on purpose */
         case 3:
            for(ii=0;ii<3+NscalarDof;ii++){ /* upper left = [ I ] */
              for(jj=0;jj<3+NscalarDof;jj++){
                offset = dof+ii+jj*vec_leng;
                rbm[offset] = (ii==jj) ? 1.0 : 0.0;
              }
            }
            for(ii=0;ii<3;ii++){ /* upper right = [ Q ] */
              for(jj=3+NscalarDof;jj<6+NscalarDof;jj++){
                offset = dof+ii+jj*vec_leng;
                if(ii == jj-3-NscalarDof) rbm[offset] = 0.0;
                else {
                     if (ii+jj == 4+NscalarDof) rbm[offset] = z[node];
                     else if ( ii+jj == 5+NscalarDof ) rbm[offset] = y[node];
                     else if ( ii+jj == 6+NscalarDof ) rbm[offset] = x[node];
                    else rbm[offset] = 0.0;
                }
              }
            }
            ii = 0; jj = 5+NscalarDof; offset = dof+ii+jj*vec_leng; rbm[offset] *= -1.0;
            ii = 1; jj = 3+NscalarDof; offset = dof+ii+jj*vec_leng; rbm[offset] *= -1.0;
            ii = 2; jj = 4+NscalarDof; offset = dof+ii+jj*vec_leng; rbm[offset] *= -1.0;
            break;
         case 2:
            for(ii=0;ii<2+NscalarDof;ii++){ /* upper left = [ I ] */
              for(jj=0;jj<2+NscalarDof;jj++){
                offset = dof+ii+jj*vec_leng;
                rbm[offset] = (ii==jj) ? 1.0 : 0.0;
              }
            }
            for(ii=0;ii<2+NscalarDof;ii++){ /* upper right = [ Q ] */
              for(jj=2+NscalarDof;jj<3+NscalarDof;jj++){
                offset = dof+ii+jj*vec_leng;
                if (ii == 0) rbm[offset] = -y[node];
                else {
                  if (ii == 1){  rbm[offset] =  x[node];}
                  else rbm[offset] = 0.0;
                }
              }
            }
            break;
         case 1:
             for (ii = 0; ii<1+NscalarDof; ii++) {
               for (jj=0; jj<1+NscalarDof; jj++) {
                  offset = dof+ii+jj*vec_leng;
                  rbm[offset] = (ii == jj) ? 1.0 : 0.0;
                }
             }
            break;

         default:
            printf("ML_Coord2RBM: Ndof = %d not implemented\n",Ndof);
            exit(1);
      } /*switch*/

  } /*for( node = 0 ; node < Nnodes; node++ )*/

  return 1;

} /*ML_Coord2RBM*/
