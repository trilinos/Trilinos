#include <stdio.h>
#include <float.h>
#include "az_aztec.h"
void AZ_matvec_mult(double *val, int *indx, int *bindx, int *rpntr, int *cpntr,
                    int *bpntr, double *b, register double *c,
                    int exchange_flag, int *data_org)

/******************************************************************************

  c = Ab:
  Sparse (square) overlapped matrix-vector multiply, using the distributed
  variable block row (DVBR) data structure (A = val).

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  m:               Number of (block) rows in A.

  val:             Array containing the entries of the matrix. The matrix is
                   stored block-row-by-block-row. Each block entry is dense and
                   stored by columns (VBR).

  indx:            The ith element of indx points to the location in val of the
                   (0,0) entry of the ith block entry. The last element is the
                   number of nonzero entries of matrix A plus one.

  bindx:           Contains the block column indices of the non-zero block
                   entries.

  rpntr:           The ith element of rpntr indicates the first point row in
                   the i'th block row. The last element is the number of block
                   rows plus one.

  cpntr:           The jth element of cpntr indicates the first point column in
                   the jth block column. The last element is the number of
  bpntr:           The ith element of bpntr points to the first block entry of
                   the ith row in bindx. The last element is the number of
                   nonzero blocks of matrix A plus one.

  b:               Contains the vector b.

  c:               Contains the result vector c.

  exchange_flag:   Flag which controls call to exchange_bdry.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see Aztec User's Guide).

******************************************************************************/

{

  /* local variables */

   AZ_MATRIX Amat;
   int proc_config[AZ_PROC_SIZE];
   static int first_time = 1;

   if (exchange_flag != 1) {
      printf("Warning: exchange_flag is no longer used in AZ_matvec_mult().\n");
      printf("         Set to '1' to avoid this message.\n");
   }
   Amat.rpntr      = rpntr;   Amat.cpntr    = cpntr;
   Amat.bpntr      = bpntr;   Amat.bindx    = bindx;
   Amat.indx       = indx;    Amat.val      = val;
   Amat.data_org   = data_org;
   Amat.aux_ival   = NULL;
   Amat.aux_dval   = NULL;
   Amat.aux_matrix = NULL;
   Amat.matrix_type = data_org[AZ_matrix_type];
#ifdef AZ_MPI
   AZ_set_comm(proc_config, MPI_COMM_WORLD);
   if (first_time == 1) {
      AZ_set_proc_config(proc_config, MPI_COMM_WORLD);
#else
   if (first_time == 1) {
      AZ_set_proc_config(proc_config, AZ_NOT_MPI);
#endif
      if (proc_config[AZ_node] == 0) {
          printf("Warning: AZ_matvec_mult() should be replaced with either\n");
          printf("          AZ_MSR_matvec_mult or AZ_VBR_matvec_mult()\n");
      }
   }
   first_time = 0;

   if      (Amat.matrix_type == AZ_MSR_MATRIX) Amat.matvec = AZ_MSR_matvec_mult;
   else if (Amat.matrix_type == AZ_VBR_MATRIX) Amat.matvec = AZ_VBR_matvec_mult;

   Amat.matvec(b, c, &Amat, proc_config);
}




