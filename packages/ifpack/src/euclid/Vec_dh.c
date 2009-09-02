/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include <stdlib.h>
#include "Vec_dh.h"
#include "Mem_dh.h"
#include "SubdomainGraph_dh.h"
#include "io_dh.h"

#undef __FUNC__
#define __FUNC__ "Vec_dhCreate"
void
Vec_dhCreate (Vec_dh * v)
{
  START_FUNC_DH
    struct _vec_dh *tmp =
    (struct _vec_dh *) MALLOC_DH (sizeof (struct _vec_dh));
  CHECK_V_ERROR;
  *v = tmp;
  tmp->n = 0;
  tmp->vals = NULL;
END_FUNC_DH}

#undef __FUNC__
#define __FUNC__ "Vec_dhDestroy"
void
Vec_dhDestroy (Vec_dh v)
{
  START_FUNC_DH if (v->vals != NULL)
    FREE_DH (v->vals);
  CHECK_V_ERROR;
  FREE_DH (v);
  CHECK_V_ERROR;
END_FUNC_DH}


#undef __FUNC__
#define __FUNC__ "Vec_dhInit"
void
Vec_dhInit (Vec_dh v, int size)
{
  START_FUNC_DH v->n = size;
  v->vals = (double *) MALLOC_DH (size * sizeof (double));
  CHECK_V_ERROR;
END_FUNC_DH}

#undef __FUNC__
#define __FUNC__ "Vec_dhCopy"
void
Vec_dhCopy (Vec_dh x, Vec_dh y)
{
  START_FUNC_DH if (x->vals == NULL)
    SET_V_ERROR ("x->vals is NULL");
  if (y->vals == NULL)
    SET_V_ERROR ("y->vals is NULL");
  if (x->n != y->n)
    SET_V_ERROR ("x and y are different lengths");
  memcpy (y->vals, x->vals, x->n * sizeof (double));
END_FUNC_DH}


#undef __FUNC__
#define __FUNC__ "Vec_dhDuplicate"
void
Vec_dhDuplicate (Vec_dh v, Vec_dh * out)
{
  START_FUNC_DH Vec_dh tmp;
  int size = v->n;
  if (v->vals == NULL)
    SET_V_ERROR ("v->vals is NULL");
  Vec_dhCreate (out);
  CHECK_V_ERROR;
  tmp = *out;
  tmp->n = size;
  tmp->vals = (double *) MALLOC_DH (size * sizeof (double));
  CHECK_V_ERROR;
END_FUNC_DH}

#undef __FUNC__
#define __FUNC__ "Vec_dhSet"
void
Vec_dhSet (Vec_dh v, double value)
{
  START_FUNC_DH int i, m = v->n;
  double *vals = v->vals;
  if (v->vals == NULL)
    SET_V_ERROR ("v->vals is NULL");
  for (i = 0; i < m; ++i)
    vals[i] = value;
END_FUNC_DH}

#undef __FUNC__
#define __FUNC__ "Vec_dhSetRand"
void
Vec_dhSetRand (Vec_dh v)
{
  START_FUNC_DH int i, m = v->n;
  double max = 0.0;
  double *vals = v->vals;

  if (v->vals == NULL)
    SET_V_ERROR ("v->vals is NULL");

#ifdef WIN32
  for (i = 0; i < m; ++i)
    vals[i] = rand ();
#else
  for (i = 0; i < m; ++i)
    vals[i] = rand ();
#endif

  /* find largest value in vector, and scale vector,
   * so all values are in [0.0,1.0]
   */
  for (i = 0; i < m; ++i)
    max = MAX (max, vals[i]);
  for (i = 0; i < m; ++i)
    vals[i] = vals[i] / max;
END_FUNC_DH}


#undef __FUNC__
#define __FUNC__ "Vec_dhPrint"
void
Vec_dhPrint (Vec_dh v, SubdomainGraph_dh sg, char *filename)
{
  START_FUNC_DH double *vals = v->vals;
  int pe, i, m = v->n;
  FILE *fp;

  if (v->vals == NULL)
    SET_V_ERROR ("v->vals is NULL");

  /*--------------------------------------------------------
   * case 1: no permutation information
   *--------------------------------------------------------*/
  if (sg == NULL)
    {
      for (pe = 0; pe < np_dh; ++pe)
	{
	  MPI_Barrier (comm_dh);
	  if (pe == myid_dh)
	    {
	      if (pe == 0)
		{
		  fp = openFile_dh (filename, "w");
		  CHECK_V_ERROR;
		}
	      else
		{
		  fp = openFile_dh (filename, "a");
		  CHECK_V_ERROR;
		}

	      for (i = 0; i < m; ++i)
		fprintf (fp, "%g\n", vals[i]);

	      closeFile_dh (fp);
	      CHECK_V_ERROR;
	    }
	}
    }

  /*--------------------------------------------------------
   * case 2: single mpi task, multiple subdomains
   *--------------------------------------------------------*/
  else if (np_dh == 1)
    {
      int i, j;

      fp = openFile_dh (filename, "w");
      CHECK_V_ERROR;

      for (i = 0; i < sg->blocks; ++i)
	{
	  int oldBlock = sg->n2o_sub[i];
	  int beg_row = sg->beg_rowP[oldBlock];
	  int end_row = beg_row + sg->row_count[oldBlock];

	  printf ("seq: block= %i  beg= %i  end= %i\n", oldBlock, beg_row,
		  end_row);


	  for (j = beg_row; j < end_row; ++j)
	    {
	      fprintf (fp, "%g\n", vals[j]);
	    }
	}
    }

  /*--------------------------------------------------------
   * case 3: multiple mpi tasks, one subdomain per task
   *--------------------------------------------------------*/
  else
    {
      int id = sg->o2n_sub[myid_dh];
      for (pe = 0; pe < np_dh; ++pe)
	{
	  MPI_Barrier (comm_dh);
	  if (id == pe)
	    {
	      if (pe == 0)
		{
		  fp = openFile_dh (filename, "w");
		  CHECK_V_ERROR;
		}
	      else
		{
		  fp = openFile_dh (filename, "a");
		  CHECK_V_ERROR;
		}

	      fprintf (stderr, "par: block= %i\n", id);

	      for (i = 0; i < m; ++i)
		{
		  fprintf (fp, "%g\n", vals[i]);
		}

	      closeFile_dh (fp);
	      CHECK_V_ERROR;
	    }
	}
    }
END_FUNC_DH}


#undef __FUNC__
#define __FUNC__ "Vec_dhPrintBIN"
void
Vec_dhPrintBIN (Vec_dh v, SubdomainGraph_dh sg, char *filename)
{
  START_FUNC_DH if (np_dh > 1)
    {
      SET_V_ERROR ("only implemented for a single MPI task");
    }
  if (sg != NULL)
    {
      SET_V_ERROR ("not implemented for reordered vector; ensure sg=NULL");
    }

  io_dh_print_ebin_vec_private (v->n, 0, v->vals, NULL, NULL, NULL, filename);
  CHECK_V_ERROR;
END_FUNC_DH}

#define MAX_JUNK 200

#undef __FUNC__
#define __FUNC__ "Vec_dhRead"
void
Vec_dhRead (Vec_dh * vout, int ignore, char *filename)
{
  START_FUNC_DH Vec_dh tmp;
  FILE *fp;
  int items, n, i;
  double *v, w;
  char junk[MAX_JUNK];

  Vec_dhCreate (&tmp);
  CHECK_V_ERROR;
  *vout = tmp;

  if (np_dh > 1)
    {
      SET_V_ERROR ("only implemented for a single MPI task");
    }

  fp = openFile_dh (filename, "w");
  CHECK_V_ERROR;

  /* skip over file lines */
  if (ignore)
    {
      printf ("Vec_dhRead:: ignoring following header lines:\n");
      printf
	("--------------------------------------------------------------\n");
      for (i = 0; i < ignore; ++i)
	{
	  fgets (junk, MAX_JUNK, fp);
	  printf ("%s", junk);
	}
      printf
	("--------------------------------------------------------------\n");
    }

  /* count floating point entries in file */
  n = 0;
  while (!feof (fp))
    {
      items = fscanf (fp, "%lg", &w);
      if (items != 1)
	{
	  break;
	}
      ++n;
    }

  printf ("Vec_dhRead:: n= %i\n", n);

  /* allocate storage */
  tmp->n = n;
  v = tmp->vals = (double *) MALLOC_DH (n * sizeof (double));
  CHECK_V_ERROR;

  /* reset file, and skip over header again */
  rewind (fp);
  rewind (fp);
  for (i = 0; i < ignore; ++i)
    {
      fgets (junk, MAX_JUNK, fp);
    }

  /* read values */
  for (i = 0; i < n; ++i)
    {
      items = fscanf (fp, "%lg", v + i);
      if (items != 1)
	{
	  sprintf (msgBuf_dh, "failed to read value %i of %i", i + 1, n);
	}
    }

  closeFile_dh (fp);
  CHECK_V_ERROR;
END_FUNC_DH}

#undef __FUNC__
#define __FUNC__ "Vec_dhReadBIN"
extern void
Vec_dhReadBIN (Vec_dh * vout, char *filename)
{
  START_FUNC_DH Vec_dh tmp;

  Vec_dhCreate (&tmp);
  CHECK_V_ERROR;
  *vout = tmp;
  io_dh_read_ebin_vec_private (&tmp->n, &tmp->vals, filename);
  CHECK_V_ERROR;
END_FUNC_DH}
