/*******************************************************************************
 * Copyright 1999, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "az_aztec.h"
#include "az_ifpack.h"

void AZ_ifpack_solve(double x[], double b[], int options[], double params[],
              int indx[], int bindx[], int rpntr[], int cpntr[], int bpntr[],
              double val[], int data_org[], double status[], int proc_config[])
{
  AZ_MATRIX *Amat;

   Amat    = AZ_matrix_create(data_org[AZ_N_internal]+data_org[AZ_N_border]);

   options[AZ_output] = 1;
   if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) 
      AZ_set_MSR(Amat, bindx, val, data_org, 0, NULL, AZ_LOCAL);
   else if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX)
      AZ_set_VBR(Amat, rpntr, cpntr, bpntr, indx, bindx, val,
                 data_org, 0, NULL, AZ_LOCAL);
   else {
      fprintf(stderr,"Unknown matrix type (%d)\n",data_org[AZ_matrix_type]);
      fprintf(stderr,"Matrix-free is now available via AZ_iterate()\n");
      exit(1);
   }

  AZ_ifpack_iterate(x, b, options, params, status, proc_config, Amat);

  AZ_matrix_destroy(&Amat);

}
/* AZ_ifpack_solve*/
