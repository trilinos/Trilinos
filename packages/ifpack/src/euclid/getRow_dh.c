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

#include "getRow_dh.h"
#include "Mat_dh.h"
#include "Euclid_dh.h"
#include "Mem_dh.h"


/*-------------------------------------------------------------------
 *  EPETRA
 *-------------------------------------------------------------------*/

#undef __FUNC__
#define __FUNC__ "EuclidGetRow"
void
EuclidGetRow (void *A, int row, int *len, int **ind, double **val)
{
  START_FUNC_DH int ierr = 0;
  if (ind != NULL)
    ierr += ExtractIndicesView (A, row, len, ind);
  if (ierr != 0)
    {
      sprintf (msgBuf_dh, "ExtractIndicesView(row= %i) returned %i", row + 1,
	       ierr);
      SET_V_ERROR (msgBuf_dh);
    }
  if (val != NULL)
    ierr += ExtractValuesView (A, row, len, val);
  if (ierr != 0)
    {
      sprintf (msgBuf_dh, " ExtractValuesView(row= %i) returned %i", row + 1,
	       ierr);
      SET_V_ERROR (msgBuf_dh);
    }
END_FUNC_DH}

#undef __FUNC__
#define __FUNC__ "EuclidRestoreRow"
void
EuclidRestoreRow (void *A, int row, int *len, int **ind, double **val)
{
  START_FUNC_DH
    /* This is needed in Hypre and Petsc, but not Trilinos. */
END_FUNC_DH}

#undef __FUNC__
#define __FUNC__ "EuclidGetDimensions"
void
EuclidGetDimensions (void *A, int *beg_row, int *rowsLocal, int *rowsGlobal)
{
  START_FUNC_DH int m, n;
  int row_start, row_end, col_start, col_end;

  row_start = MinMaxMyGID (A, true, true);
  row_end = MinMaxMyGID (A, true, false);
  col_start = MinMaxMyGID (A, false, true);
  col_end = MinMaxMyGID (A, false, false);

  m = NumGlobalRowCol (A, false);
  n = NumGlobalRowCol (A, true);
  *beg_row = row_start;
  *rowsLocal = (row_end - row_start + 1);
  *rowsGlobal = n;
END_FUNC_DH}

/*
#undef __FUNC__
#define __FUNC__ "EuclidReadLocalNz"
int EuclidReadLocalNz(void *A)
{
  START_FUNC_DH
  if (ignoreMe) SET_V_ERROR("not implemented");
  return(0);
  END_FUNC_DH
}
*/

#undef __FUNC__
#define __FUNC__ "PrintMatUsingGetRow"
void
PrintMatUsingGetRow (void *A, int beg_row, int m,
		     int *n2o_row, int *n2o_col, char *filename)
{
  START_FUNC_DH FILE *fp;
  int *o2n_col = NULL, pe, i, j, *cval, len;
  int newCol, newRow;
  double *aval;

  /* form inverse column permutation */
  if (n2o_col != NULL)
    {
      o2n_col = (int *) MALLOC_DH (m * sizeof (int));
      CHECK_V_ERROR;
      for (i = 0; i < m; ++i)
	o2n_col[n2o_col[i]] = i;
    }

  for (pe = 0; pe < np_dh; ++pe)
    {

      MPI_Barrier (comm_dh);

      if (myid_dh == pe)
	{
	  if (pe == 0)
	    {
	      fp = fopen (filename, "w");
	    }
	  else
	    {
	      fp = fopen (filename, "a");
	    }
	  if (fp == NULL)
	    {
	      sprintf (msgBuf_dh, "can't open %s for writing\n", filename);
	      SET_V_ERROR (msgBuf_dh);
	    }

	  for (i = 0; i < m; ++i)
	    {

	      if (n2o_row == NULL)
		{
		  EuclidGetRow (A, i + beg_row, &len, &cval, &aval);
		  CHECK_V_ERROR;
		  for (j = 0; j < len; ++j)
		    {
		      fprintf (fp, "%i %i %g\n", i + 1, cval[j], aval[j]);
		    }
		  EuclidRestoreRow (A, i, &len, &cval, &aval);
		  CHECK_V_ERROR;
		}
	      else
		{
		  newRow = n2o_row[i] + beg_row;
		  EuclidGetRow (A, newRow, &len, &cval, &aval);
		  CHECK_V_ERROR;
		  for (j = 0; j < len; ++j)
		    {
		      newCol = o2n_col[cval[j] - beg_row] + beg_row;
		      fprintf (fp, "%i %i %g\n", i + 1, newCol, aval[j]);
		    }
		  EuclidRestoreRow (A, i, &len, &cval, &aval);
		  CHECK_V_ERROR;
		}
	    }
	  fclose (fp);
	}
    }

  if (n2o_col != NULL)
    {
      FREE_DH (o2n_col);
      CHECK_V_ERROR;
    }
END_FUNC_DH}
