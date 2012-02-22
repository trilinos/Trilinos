/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* compute null space from rigid body modes for smoothed aggregation    */
/* -------------------------------------------------------------------- */
/* Author : Ray Tuminaro (SNL)                                          */
/* ******************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include "ml_rbm.h"

/* ************************************************************************** *
  Input
    Nnodes - number of nodes in the amalgamated system (i.e., #nodes in physical mesh)
    x,y,z  - coordinates
    rbm    - allocated vector to hold rigid body modes
    Ndof   - #dofs at a node
  Output
    rbm    - vector populated with rigid body modes

    Most general case corresponds to the local 6 x 6 matrix:

                              translations  |         rotations around
                         x       y       z  |    x        y           z
                         ---------------------------------------------------
        x-direction      1       0       0  |    0      z-zhat      yhat-y 
        y-direction      0       1       0  | zhat-z      0         x-xhat
        z-direction      0       0       1  | y-yhat    xhat-x        0
        x-rot            0       0       0  |    1        0           0
        y-rot            0       0       0  |    0        1           0
        z-rot            0       0       0  |    0        0           1

    The above corresponds to 3D elasticity with shell elements. 3D elasticity with
    bricks would be the same with the last 3 rows removed. 2D elasticity would also
    remove the 3rd row and columns 4 and 5. 1D elasticty would remove all but
    the (1,1) entry.

 * ************************************************************************** */

int ML_Coord2RBM(int Nnodes, double x[], double y[], double z[], double rbm[], int Ndof)
{
   int vec_leng, ii, jj, offset, node, dof;

   vec_leng = Nnodes*Ndof;

   for( node = 0 ; node < Nnodes; node++ )
   {
      dof = node*Ndof;
      switch( Ndof )
      {
         case 6: 
            for(ii=3;ii<6;ii++){ /* lower half = [ 0 I ] */
              for(jj=0;jj<6;jj++){
                offset = dof+ii+jj*vec_leng;
                rbm[offset] = (ii==jj) ? 1.0 : 0.0;
              }
            }

         case 3: 
            for(ii=0;ii<3;ii++){ /* upper left = [ I ] */
              for(jj=0;jj<3;jj++){
                offset = dof+ii+jj*vec_leng;
                rbm[offset] = (ii==jj) ? 1.0 : 0.0;        
              }
            }
            for(ii=0;ii<3;ii++){ /* upper right = [ Q ] */
              for(jj=3;jj<6;jj++){
                offset = dof+ii+jj*vec_leng;
                if( ii == jj-3 ) rbm[offset] = 0.0;
                else {
                  if (ii+jj == 4) rbm[offset] = z[node];
                  else if ( ii+jj == 5 ) rbm[offset] = y[node];
                  else if ( ii+jj == 6 ) rbm[offset] = x[node];
                  else rbm[offset] = 0.0;
                }
              }
            }
            ii = 0; jj = 5; offset = dof+ii+jj*vec_leng; rbm[offset] *= -1.0;
            ii = 1; jj = 3; offset = dof+ii+jj*vec_leng; rbm[offset] *= -1.0;
            ii = 2; jj = 4; offset = dof+ii+jj*vec_leng; rbm[offset] *= -1.0;
            break;
         case 2: 
            for(ii=0;ii<2;ii++){ /* upper left = [ I ] */
              for(jj=0;jj<2;jj++){
                offset = dof+ii+jj*vec_leng;
                rbm[offset] = (ii==jj) ? 1.0 : 0.0;        
              }
            }
            for(ii=0;ii<2;ii++){ /* upper right = [ Q ] */
              for(jj=2;jj<3;jj++){
                offset = dof+ii+jj*vec_leng;
                if (ii == 0) rbm[offset] = -y[node];
                else         rbm[offset] =  x[node];
              }
            }
            break;
         case 1: 
             rbm[dof] = 1;
            break;

         default: 
            printf("ML_Coord2RBM: Ndof = %d not implemented\n",Ndof);
            exit(1);
      } /*switch*/

  } /*for( node = 0 ; node < Nnodes; node++ )*/

  return 1;

} /*ML_Coord2RBM*/
