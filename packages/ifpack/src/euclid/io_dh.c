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

#include "io_dh.h"
#include "Mat_dh.h"
#include "Vec_dh.h"
#include "Mem_dh.h"
#include "Timer_dh.h"
#include "Parser_dh.h"
#include "mat_dh_private.h"

#undef __FUNC__
#define __FUNC__ "openFile_dh"
FILE *
openFile_dh (const char *filenameIN, const char *modeIN)
{
  START_FUNC_DH FILE * fp = NULL;

  if ((fp = fopen (filenameIN, modeIN)) == NULL)
    {
      sprintf (msgBuf_dh, "can't open file: %s for mode %s\n", filenameIN,
	       modeIN);
      SET_ERROR (NULL, msgBuf_dh);
    }
END_FUNC_VAL (fp)}

#undef __FUNC__
#define __FUNC__ "closeFile_dh"
void
closeFile_dh (FILE * fpIN)
{
  if (fclose (fpIN))
    {
      SET_V_ERROR ("attempt to close file failed");
    }
}

/*----------------------------------------------------------------*/
void
io_dh_print_ebin_mat_private (int m, int beg_row,
			      int *rp, int *cval, double *aval,
			      int *n2o, int *o2n, Hash_i_dh hash,
			      char *filename)
{
}

extern void
io_dh_read_ebin_mat_private (int *m, int **rp, int **cval,
			     double **aval, char *filename)
{
}

void
io_dh_print_ebin_vec_private (int n, int beg_row, double *vals,
			      int *n2o, int *o2n, Hash_i_dh hash,
			      char *filename)
{
}

void
io_dh_read_ebin_vec_private (int *n, double **vals, char *filename)
{
}
