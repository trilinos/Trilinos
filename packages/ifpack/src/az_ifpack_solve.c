/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <std::string.h>
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
