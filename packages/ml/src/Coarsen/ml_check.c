/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#include <stdlib.h>
#include "ml_struct.h"
#include "ml_check.h"
#include "ml_rap.h"
#include "ml_comm.h"
#include "ml_defs.h"
#include <math.h>

/* ******************************************************************** */
/* check the correctness of the interpolation operator                  */
/* -------------------------------------------------------------------- */

void ML_interp_check(ML *ml, int coarse_level, int fine_level) 
{
   int    ii, jj, ncoarse, nfine;
   double *c_data, *f_data, coords[3], dtemp, d2, dlargest;
   ML_GridFunc *coarse_funs, *fine_funs;
   void   *coarse_data, *fine_data;
   int    nfine_eqn, ncoarse_eqn, stride = 1;

   /* check an interpolated linear function */             

   coarse_data = ml->SingleLevel[coarse_level].Grid->Grid;
   fine_data   = ml->SingleLevel[  fine_level].Grid->Grid;
   coarse_funs = ml->SingleLevel[coarse_level].Grid->gridfcn;
   fine_funs   = ml->SingleLevel[  fine_level].Grid->gridfcn;

   if ( (coarse_data==NULL)||(fine_data==NULL)) {
      printf("ML_interp_check: grid data not found?\n");
      exit(1);
   }
   if ( (coarse_funs==NULL)||(fine_funs==NULL)) {
      printf("ML_interp_check: grid functions not found?\n");
      exit(1);
   }
   if ( (coarse_funs->USR_grid_get_nvertices == 0) ||
        (  fine_funs->USR_grid_get_nvertices == 0)) {
      printf("ML_interp_check: USR_grid_get_nvertices not found?\n");
      exit(1);
   }

   ncoarse     = coarse_funs->USR_grid_get_nvertices(coarse_data);
   nfine       =   fine_funs->USR_grid_get_nvertices(  fine_data);
   nfine_eqn   = ml->SingleLevel[coarse_level].Pmat->outvec_leng;
   ncoarse_eqn = ml->SingleLevel[coarse_level].Pmat->invec_leng;
                                                                    
   c_data  = (double *) ML_allocate(ncoarse_eqn*sizeof(double));    
   f_data  = (double *) ML_allocate(nfine_eqn*sizeof(double));    
   for (ii = 0; ii < ncoarse_eqn; ii++) c_data[ii] = 0.;
   for (ii = 0; ii < nfine_eqn; ii++) f_data[ii] = 0.;

   /* ASSUMING that for each grid point on this processor there are a   */
   /* set of equations and that all points at a grid point are numbered */
   /* consecutively !!!!!!!!!!!!!!                                      */

   stride = nfine_eqn/nfine;
   for (ii = 0 ; ii < ncoarse ; ii++)  {                   
      coarse_funs->USR_grid_get_vertex_coordinate(coarse_data,ii,coords); 
      for (jj = 0; jj < stride; jj++) {
         c_data[ii*stride + jj] = coords[0] + 3.*coords[1] + .5;                
      }
   }                                                       

   ML_Operator_Apply(ml->SingleLevel[coarse_level].Pmat,ncoarse_eqn,c_data,
		     nfine_eqn,f_data);

   dlargest = 0.0;                                         
   for (ii = 0 ; ii < nfine; ii++)  {                      
      fine_funs->USR_grid_get_vertex_coordinate(fine_data , ii, coords);
      dtemp = coords[0] + 3.*coords[1] + .5;
      d2 = ML_dabs(dtemp - f_data[ii*stride])/(ML_dabs(dtemp)+1.e-9);    
      /* Ray debugging
      if ( d2 > 1.e-8)                                      
         printf("%d: f_data[%d] = %e  %e | %e %e\n",ml->comm->ML_mypid,
                ii,f_data[ii*stride],dtemp,coords[0],coords[1]);
      */
      if ( d2 > dlargest) {                                 
            dlargest = d2;   
      }                                                     
   }                                                       
   ML_free(f_data);
   ML_free(c_data);
}

/* For systems derived from time-dependent Maxwell's equations, check the
 * relationship S*T on each level. */

int ML_Reitzinger_Check_Hierarchy(ML *ml, ML_Operator **Tmat_array, int incr_or_decr)
{
  int i,j;
  int finest_level, coarsest_level;
  ML_Operator *Amat, *Tmat;
  double *randvec, *result, *result1;
  double dnorm;

  finest_level = ml->ML_finest_level;
  coarsest_level = ml->ML_coarsest_level;

  if (incr_or_decr == ML_INCREASING) {
    if (ml->comm->ML_mypid == 0) {
      printf("ML_Reitzinger_Check_Hierarchy: ML_INCREASING is not supported ");
      printf(" at this time.  Not checking hierarchy.\n");
    }
    return 1;
  }

  if ( ML_Get_PrintLevel() > 5 ) {
    printf("ML_Reitzinger_Check_Hierarchy: Checking null space\n");
  }

  for (i=finest_level; i>coarsest_level; i--) {

     Amat = ml->Amat+i;
     Tmat = Tmat_array[i];

     /* normalized random vector */
     randvec = (double *) ML_allocate(Tmat->invec_leng * sizeof(double) );
     ML_random_vec(randvec,Tmat->invec_leng, ml->comm);
     dnorm = sqrt( ML_gdot(Tmat->invec_leng, randvec, randvec, ml->comm) );
     for (j=0; j<Tmat->invec_leng; j++) randvec[j] /=  dnorm;

     result = (double *) ML_allocate(Amat->invec_leng * sizeof(double) );
     result1 = (double *) ML_allocate(Amat->outvec_leng * sizeof(double) );

     ML_Operator_Apply(Tmat, Tmat->invec_leng, randvec,
                       Tmat->outvec_leng, result);
     ML_Operator_Apply(Amat, Amat->invec_leng, result,
                       Amat->outvec_leng, result1);

     dnorm = sqrt( ML_gdot(Amat->outvec_leng, result1, result1, ml->comm) );
     if ( (ML_Get_PrintLevel() > 5) && (ml->comm->ML_mypid == 0) ) {
       printf("Level %d: for random v,  ||S*T*v|| = %15.10e\n",i,dnorm);
     }

     ML_free(randvec);
     ML_free(result);
     ML_free(result1);
  }
  if ( (ML_Get_PrintLevel() > 5) && (ml->comm->ML_mypid == 0) ) printf("\n");

  return 0;

}

/*
void ML_check_it(double sol[], double rhs[], ML *ml)
{
   int i,j;
   double *coarse_sol, *coarse_rhs;

   coarse_sol = (double *) ML_allocate((ml->SingleLevel[0].Amat->invec_leng)*
                                       sizeof(double));
   coarse_rhs = (double *) ML_allocate((ml->SingleLevel[0].Amat->invec_leng)*
                                        sizeof(double));
   ML_random_vec(coarse_sol,  ml->SingleLevel[0].Amat->invec_leng, ml->comm);

   if ( ml->SingleLevel[0].BCs->Dirichlet_eqn_list != NULL) {
      for (i = 0; i < ml->SingleLevel[0].BCs->Dirichlet_eqn_length; i++ )
         coarse_sol[ml->SingleLevel[0].BCs->Dirichlet_eqn_list[i]] = 0.0;
   }
   ml->SingleLevel[0].Amat->matvec->func_ptr(ml->SingleLevel[0].Amat->data, coarse_sol, coarse_rhs);

   if ( ml->sl_ptr[0]->Dirichlet_eqn_list != NULL) {
      for (i = 0; i < ml->sl_ptr[0]->Dirichlet_eqn_length; i++ ){
         if (ML_dabs(coarse_rhs[ml->sl_ptr[0]->Dirichlet_eqn_list[i]]) > 1.0e-8)
            printf("disc mat not zero on boundary %e\n",
                          coarse_rhs[ml->sl_ptr[0]->Dirichlet_eqn_list[i]]);
      }
   }

   i = ML_Comm_GsumInt(ml->comm,ml->sl_ptr[1]->Dirichlet_eqn_length);
   j = ML_Comm_GsumInt(ml->comm,ml->sl_ptr[0]->Dirichlet_eqn_length);
   if (ml->comm->ML_mypid == 0)
      printf("boundary length: fine = %d & coarse = %d\n", i,j);

   ML_free(coarse_sol);
   ML_free(coarse_rhs);
   ML_interp_check(ml, 0, 1);
}
*/

/*
ML_Check( ML * ml )
{
   int    i, j, nlevels, myrank, scheme, nprocs, inttmp, status, sum;
   int    flag;
   double dtmp;
   SingleLevel  *sl;

   nlevels = ml->ML_num_levels;
   if ( ml->comm == NULL ) flag = 1; else flag = 0;
   sum = ML_Comm_GsumInt(ml->comm,flag);
   if ( sum > 0 ) printf("Communicator absent. %d\n",sum);
   if ( ml->comm != NULL )
   {
      myrank  = ml->comm->ML_mypid;
      if ( myrank == 0 ) printf("*** Communicator defined. \n");
      nprocs  = ml->comm->ML_nprocs;
      if ( myrank == 0 ) printf("Number of processors = %d.\n",nprocs);

      if ( myrank == 0 ) printf("Checking communicator ... \n");
      if ( ml->comm->USR_sendbytes  == NULL ) flag = 1; else flag = 0;
      if ( ml->comm->USR_irecvbytes == NULL ) flag = 1;
      if ( ml->comm->USR_waitbytes  == NULL ) flag = 1;
      sum = ML_Comm_GsumInt(ml->comm,flag);
      if ( myrank == 0 )
      {
         if ( sum > 0 ) printf("\tCommunicator functions not complete.\n");
         else           printf("\tCommunicator functions complete.\n");
      }
      if ( sum > 0 && myrank == 0 )
      {
         if ( ml->comm->USR_sendbytes == NULL ) 
              printf("\t\tSend function is absent.\n");
         if ( ml->comm->USR_irecvbytes == NULL ) 
              printf("\t\tIrecv function is absent.\n");
         if ( ml->comm->USR_waitbytes == NULL ) 
              printf("\t\tIrecv function is absent.\n");
      }
   } 
   scheme  = ml->ML_scheme;
   if ( myrank == 0 )
   {
      switch ( scheme )
      {
         case ML_NONE :  printf("ML scheme is NONE.\n");
                         break;
         case ML_MGV :   printf("ML scheme is MGV.\n");
                         break;
         case ML_MGFULLV: printf("ML scheme is MGFullV.\n");
                         break;
         case ML_MG2CGC: printf("ML scheme is 2-level MG.\n");
                         break;
         case ML_MGW :   printf("ML scheme is MGW.\n");
                         break;
         default:        printf("ML scheme not defined.\n");
                         break;
      }
   }
   for ( i = 0; i < nlevels; i++ )
   {
      if ( myrank == 0 ) printf("Checking level %d \n", i );
      sl = ml->SingleLevel[i];
      if ( sl->grid == NULL ) flag = 1 ; else flag = 0; 
      sum = ML_Comm_GsumInt(ml->comm,flag);
      if ( myrank == 0 )
      {
         if (sum > 0) printf("\tSL user-defined grid is absent.\n");
         else         printf("\tSL user-defined grid is present.\n");
      }
      if ( sl->gridfcn == NULL ) flag = 1; else flag = 0;
      sum = ML_Comm_GsumInt(ml->comm,flag);
      if ( myrank == 0 )
      {
         if (sum > 0) printf("\tSL grid function object is absent.\n");
         else         printf("\tSL grid function object is present.\n");
      }
      if ( sum == 0 ) 
      {
         status = ML_GridFcns_Check( sl->gridfcn );
         sum = ML_Comm_GsumInt(ml->comm,flag);
         if ( status < 0 && myrank == 0 )
         {
            printf("\tSL grid functions incompletely defined.\n");
         }
      }
      if (myrank==0) 
         printf("\tExamining the Dirichlet list ... \n");
      if ( sl->Dirichlet_list != NULL )
      {
         for ( j = 0; j < sl->Dirichlet_length; j++ )
            inttmp = sl->Dirichlet_list[j];
      }

      if (myrank==0) printf("\tExamining Amat scales ... \n");
      if ( sl->Amat_scales != NULL )
      {
         for ( j = 0; j < sl->Amat_length; j++ )
            dtmp = sl->Amat_scales[j];
      }

      if ( myrank == 0 )
      {
         if ( sl->scaling_lock == 1 )
            printf("\tScaling lock is on. \n");
         else
            printf("\tScaling lock is off. \n");
      }
      if (myrank==0) printf("\tExamining Amat diagonal ... \n");
      if ( sl->Amat_diagonal != NULL )
      {
         for ( j = 0; j < sl->Amat_length; j++ )
         {
            dtmp = sl->Amat_diagonal[j];
            if ( dtmp < 1.0E-15 && dtmp > -1.0E-15 )
               printf("\t%d : diagonal element %d = 0.\n", myrank, j);
         }
      }

      if ( sl->Amat == NULL ) flag = 1; else flag = 0;
      sum = ML_Comm_GsumInt(ml->comm,flag);
      if ( i != 0 && myrank == 0 )
      {
         if (sum == 0) printf("\tAmat is present. \n");
         else          printf("\tAmat is absent. %d\n", sum);
      }

      if ( sl->matvec_struct == NULL ) flag = 1; else flag = 0;
      if ( sl->matvec        == NULL ) flag = 1;
      sum = ML_Comm_GsumInt(ml->comm,flag);
      if ( i != 0 && myrank == 0 )
      {
         if (sum == 0) printf("\tMatvec is complete. \n");
         else          printf("\tMatvec is incomplete. \n");
         if ( sum > 0 )
         {
            if ( sl->matvec == NULL ) 
               printf("\t\tMatvec function is absent. \n");
            if ( sl->matvec_struct == NULL ) 
               printf("\t\tMatvec obsent is absent. \n");
         }
      }

      if ( sl->Rmat       == NULL ) flag = 1; else flag = 0;
      if ( sl->restrictor == NULL ) flag = 1;
      sum = ML_Comm_GsumInt(ml->comm,flag);
      if ( myrank == 0 && i > 0 )
      {
         if (sum == 0) printf("\tRestriction operator is complete.\n");
         else          printf("\tRestriction operator is incomplete.\n");
         if ( sl->restrictor == NULL ) 
            printf("\t\tRestriction function is absent.\n");
         if ( sl->Rmat == NULL ) 
            printf("\t\tRestriction matrix is absent.\n");
      }

      if ( sl->Pmat        == NULL ) flag = 1; else flag = 0;
      if ( sl->prolongator == NULL ) flag = 1;
      sum = ML_Comm_GsumInt(ml->comm,flag);
      if ( myrank == 0 && i < nlevels-1)
      {
         if (sum == 0) printf("\tprolongation operator is complete.\n");
         else          printf("\tprolongation operator is incomplete.\n");
         if ( sl->prolongator == NULL ) 
            printf("\t\tProlongation function is absent.\n");
         if ( sl->Pmat == NULL ) 
            printf("\t\tProlongation matrix is absent.\n");
      }

      if ( sl->relax1_struct == NULL ) flag = 1; else flag = 0;
      if ( sl->relax1        == NULL ) flag = 1;
      sum = ML_Comm_GsumInt(ml->comm,flag);
      if ( myrank == 0 )
      {
         if (sum == 0) printf("\tPre-smoother is complete.\n");
         else          printf("\tPre-smoother is incomplete.\n");
         if ( sl->relax1 == NULL ) 
            printf("\t\tPre-smoothing function is absent.\n");
         if ( sl->relax1_struct == NULL ) 
            printf("\t\tPre-smoothing context is absent.\n");
      }

      if ( sl->relax2_struct == NULL ) flag = 1; else flag = 0;
      if ( sl->relax2        == NULL ) flag = 1;
      sum = ML_Comm_GsumInt(ml->comm,flag);
      if ( myrank == 0 )
      {
         if (sum == 0) printf("\tPost-smoother is complete.\n");
         else          printf("\tPost-smoother is incomplete.\n");
         if ( sl->relax2 == NULL ) 
            printf("\t\tPost-smoothing function is absent.\n");
         if ( sl->relax2_struct == NULL ) 
            printf("\t\tPost-smoothing context is absent.\n");
      }
         
      if ( i == 0 )
      {
         if ( sl->csolve_struct == NULL ) flag = 1; else flag = 0;
         if ( sl->csolve        == NULL ) flag = 1;
         sum = ML_Comm_GsumInt(ml->comm,flag);
         if ( myrank == 0 )
         {
            if (sum == 0) printf("\tCoarse solver is complete.\n");
            else          printf("\tCoarse solver is incomplete.\n");
            if ( sl->csolve == NULL ) 
               printf("\t\tCoarse solver function is absent.\n");
            if ( sl->csolve_struct == NULL ) 
               printf("\t\tCoarse solver context is absent.\n");
         }
      }

      if ( sl->eqn2grid_map == NULL ) flag = 1; else flag = 0;
      sum = ML_Comm_GsumInt(ml->comm,flag);
      if ( myrank == 0 )
      {
         if (sum==0) printf("\tEquation to grid map is present.\n");
         else        printf("\tEquation to grid map is absent.\n");
      }
      if (ML_Mapper_Check(sl->grid2eqn_map) == 0) flag = 1; else flag = 0;
      sum = ML_Comm_GsumInt(ml->comm,flag);
      if ( myrank == 0 )
      {
         if (sum==0) printf("\tGrid to equation map is present.\n");
         else        printf("\tGrid to equation map is absent.\n");
      }
   }
}
*/

