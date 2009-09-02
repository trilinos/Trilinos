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

int ML_Coord2RBM(int Nnodes, double x[], double y[], double z[], 
	double rbm[], int Ndof)
{
   int vec_leng, ii, jj, offset, node, dof;

   vec_leng = Nnodes*Ndof;

   for( node = 0 ; node < Nnodes; node++ ) {
      dof = node*Ndof;
      switch( Ndof ){
         case 1: 
            rbm[node] = 1.0; 
            break;
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
            ii = 0; jj = 5; 
            offset = dof+ii+jj*vec_leng; rbm[offset] *= -1.0;
            ii = 1; jj = 3; 
            offset = dof+ii+jj*vec_leng; rbm[offset] *= -1.0;
            ii = 2; jj = 4; 
            offset = dof+ii+jj*vec_leng; rbm[offset] *= -1.0;
            break;
         default: 
            printf("ML_Coord2RBM: Ndof = %d not implemented\n",Ndof);
            exit(1);
      }
/*

    if( D1_2inv ) { 
      for(ii=kk=0,tt=Ndof*xx;ii<Ndof;ii++,tt++){ 
	temp=1./scale[tt];
	curr->B[kk] = curr->B[kk]*temp;
	curr->B[kk+1] = curr->B[kk+1]*temp;
	curr->B[kk+2] = curr->B[kk+2]*temp;
	curr->B[kk+3] = curr->B[kk+3]*temp;	
	curr->B[kk+4] = curr->B[kk+4]*temp;
	curr->B[kk+5] = curr->B[kk+5]*temp;
	kk += 6;
      }
    }
*/
   }
   return 1;
}

