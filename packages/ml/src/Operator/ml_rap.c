/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* ******************************************************************** */
/* Functions for generating Galerkin coarse grid operator               */
/* ******************************************************************** */

#include <stdio.h>
#include <math.h>
#include "ml_struct.h"
#include "ml_rap.h"

/* ******************************************************************** */
/* Given R, A, and P, generate a Galerkin operator at a given level     */
/* i.e. Result = Rmat * Amat * Pmat, where                              */
/*                                                                      */
/*    Rmat, Amat, Pmat - On input, sparse matrices.                     */
/*    Result           - On input, an empty sparse matrix.              */
/*                     - On output, Result = Rmat * Amat * Pmat.        */
/*    matrix_type      - On input, indicates whether the resulting      */
/*                       matrix should be CSR or MSR.                   */
/* -------------------------------------------------------------------- */

void ML_rap(ML_Operator *Rmat, ML_Operator *Amat, 
            ML_Operator *Pmat, ML_Operator *Result, int matrix_type)
{
  int         max_per_proc, i, j, N_input_vector;
   ML_Operator *APmat, *RAPmat, *Pcomm, *RAPcomm, *APcomm, *AP2comm, *tptr;
   ML_CommInfoOP *getrow_comm; 
   double      *scales = NULL;

   /* Check that N_input_vector is reasonable */

   N_input_vector = Pmat->invec_leng;
   getrow_comm = Pmat->getrow->pre_comm;
   if ( getrow_comm != NULL) {
      for (i = 0; i < getrow_comm->N_neighbors; i++) {
         for (j = 0; j < getrow_comm->neighbors[i].N_send; j++) {
            if (getrow_comm->neighbors[i].send_list[j] >= N_input_vector) {
              printf("(%d) Error: N_input_vector (%d) argument to rap() is not \n", Amat->comm->ML_mypid,N_input_vector);
              printf("(%d) Error: larger than %dth element (%d) sent to node %d\n", Amat->comm->ML_mypid,j+1,
                     getrow_comm->neighbors[i].send_list[j],
                     getrow_comm->neighbors[i].ML_id);
printf("(%d) Error: Amat(%d,%d)  Rmat(%d,%d)  Pmat(%d,%d)\n",
                Amat->comm->ML_mypid,
                Amat->outvec_leng,Amat->invec_leng,
                Rmat->outvec_leng,Rmat->invec_leng,
                Pmat->outvec_leng,Pmat->invec_leng);
fflush(stdout);
              exit(1);
            }
         }
      }
   }


   ML_create_unique_col_id(N_input_vector, &(Pmat->getrow->loc_glob_map),
                           getrow_comm, &max_per_proc, Pmat->comm);
   Pmat->getrow->use_loc_glob_map = ML_YES;


 
   if (Amat->getrow->pre_comm != NULL) 
      ML_exchange_rows( Pmat, &Pcomm, Amat->getrow->pre_comm);
   else Pcomm = Pmat;

#ifdef DEBUG
   if ( Pmat->comm->ML_mypid == 0 )
      printf("ML_rap : A * P begins...\n");
#endif
   ML_matmat_mult(Amat, Pcomm , &APmat);
#ifdef DEBUG
   if ( Pmat->comm->ML_mypid == 0 )
      printf("ML_rap : A * P ends.\n");
#endif

   ML_free(Pmat->getrow->loc_glob_map); Pmat->getrow->loc_glob_map = NULL;
   Pmat->getrow->use_loc_glob_map = ML_NO;
   if (Amat->getrow->pre_comm != NULL) {
      tptr = Pcomm;
      while ( (tptr!= NULL) && (tptr->sub_matrix != Pmat)) 
         tptr = tptr->sub_matrix;
      if (tptr != NULL) tptr->sub_matrix = NULL;
      ML_RECUR_CSR_MSRdata_Destroy(Pcomm);
      ML_Operator_Destroy(&Pcomm);
   }

   if (Amat->getrow->post_comm != NULL) {
      ML_exchange_rows(APmat, &APcomm, Amat->getrow->post_comm);
   }
   else APcomm = APmat;

   /* Take into account any scaling in Amat */

   if (Rmat->from != NULL)
      ML_DVector_GetDataPtr(Rmat->from->Amat_Normalization,&scales);
   if (scales != NULL) 
      ML_Scale_CSR(APcomm, scales, 0);

   if (Rmat->getrow->pre_comm != NULL) 
      ML_exchange_rows( APcomm, &AP2comm, Rmat->getrow->pre_comm);
   else AP2comm = APcomm;

#ifdef DEBUG
   if ( Pmat->comm->ML_mypid == 0 )
      printf("ML_rap : R * AP begins...\n");
#endif

   ML_matmat_mult(Rmat,AP2comm, &RAPmat);

#ifdef DEBUG
   if ( Pmat->comm->ML_mypid == 0 )
      printf("ML_rap : R * AP ends.\n");
#endif

   ML_RECUR_CSR_MSRdata_Destroy(AP2comm);
   ML_Operator_Destroy(&AP2comm);

   if (Rmat->getrow->post_comm != NULL) 
      ML_exchange_rows( RAPmat, &RAPcomm, Rmat->getrow->post_comm);
   else RAPcomm = RAPmat;

   scales = NULL;
   if (Rmat->to != NULL)
      ML_DVector_GetDataPtr(Rmat->to->Amat_Normalization,&scales);
   if (scales != NULL) 
      ML_Scale_CSR(RAPcomm, scales, 1);

   RAPcomm->num_PDEs = Amat->num_PDEs;
   RAPcomm->num_rigid = Amat->num_rigid;
   if (matrix_type == ML_MSR_MATRIX)
     ML_back_to_local(RAPcomm, Result, max_per_proc);
   else if (matrix_type == ML_CSR_MATRIX)
     ML_back_to_csrlocal(RAPcomm, Result, max_per_proc);
   else pr_error("ML_RAP: unknown matrix type\n");

   ML_RECUR_CSR_MSRdata_Destroy(RAPcomm);
   ML_Operator_Destroy(&RAPcomm);
#ifdef DEBUG
   if ( Pmat->comm->ML_mypid == 0 )
      printf("ML_rap ends.\n");
#endif
}
