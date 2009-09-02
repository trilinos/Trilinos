// @HEADER
// ***********************************************************************
// 
//                 TriUtils: Trilinos Utilities Package
//                 Copyright (2001) Sandia Corporation
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
// @HEADER

#include "Trilinos_Util.h"
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void Trilinos_Util_msr2vbr(double val[], int indx[], int rnptr[], 
                int cnptr[], int bnptr[],
                int bindx[], int msr_bindx[], double msr_val[],
                int total_blk_rows, int total_blk_cols, int blk_space,
                int nz_space, int blk_type)

/*******************************************************************************

  Convert the MSR matrix defined in [msr_val,msr_bindx] to a VBR matrix.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val,
  rnptr,
  bindx,
  indx,
  bnptr,
  cnptr:           Sparse matrix arrays. See User's Guide.
                   On input, the matrix corresponds to the initial ordering
                   (e.g. row i corresponds to global row update[i]).
                   On output, the matrix rows and columns are renumbered to
                   correspond to the ordering given by 'update_index' and
                   'extern_index'. (e.g. global row update[i] now appears
                   as row update_index[i] in the matrix).

  msr_val,
  msr_bindx:       On input, MSR matrix to be converted to VBR.
                   See User's Guide.

  total_blk_rows:  Number of block rows in resulting local VBR matrix.

  total_blk_cols:  Number of block columns in resulting local VBR matrix.

  blk_space:       Length of storage allocated for bindx[] and indx[]. An error
                   message will be printed if we try to write past these arrays.

  nz_space:        Length of storage allocated for val[]. An error message will
                   be printed if we try to write past this array.

  blk_type:        If blk_type > 0, blk_type indicates that all block rows (and
                   colunns) have the same size given by 'blk_type'. If
                   blk_type < 0, the block rows have different sizes.

*******************************************************************************/

{

  /* local variables */

  int therow, thecol;
  int i, j;

  /**************************** execution begins ******************************/

  for (i = 0; i < total_blk_rows; i++) rnptr[i] = cnptr[i];

  Trilinos_Util_convert_values_to_ptrs(rnptr, total_blk_rows, 0);
  Trilinos_Util_convert_values_to_ptrs(cnptr, total_blk_cols, 0);

  indx[0] = bnptr[0] = 0;

  /* go through each block row */

  for (i = 0; i < total_blk_rows; i++) {
    bnptr[i + 1] = bnptr[i];

    for (therow = rnptr[i]; therow < rnptr[i + 1]; therow++) {

      /* add the diagonal entry */

      thecol = therow;
      Trilinos_Util_add_new_ele(cnptr, therow, i, bindx, bnptr, indx, val, therow,
                     msr_val[therow], total_blk_cols, blk_space, nz_space,
                     blk_type);

      /* add off diagonal entries */

      for (j = msr_bindx[therow]; j < msr_bindx[therow + 1]; j++) {
        thecol = msr_bindx[j];
        Trilinos_Util_add_new_ele(cnptr, thecol, i, bindx, bnptr, indx, val, therow,
                       msr_val[j], total_blk_cols, blk_space, nz_space,
                       blk_type);
      }
    }
  }

} /* Trilinos_Util_msr2vbr */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int Trilinos_Util_find_block_col(int cnptr[], int column, int max_blocks, int blk_size)

/*******************************************************************************

  Return the local index of the block column witch contains the point column
  given by local index 'column'.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     int, local index of the block column.
  ============

  Parameter list:
  ===============

  cnptr:           cnptr[0] = 0 and cnptr[i+1] - cnptr[i] gives the column
                   dimension of global block column
                     'update[i]' if  i <  N_update
                     'external[k]' if  i >= N_update
                   where k = i - N_update.

  column:          Local column index of the column for which we are trying to
                   find the block column containing it.

  blk_size:        blk_size > 0 ==> all the blocks are the same size so we can
                   use a shortcut in computing the block column index.
                   blk_size = 0 ==> short cut not used.

*******************************************************************************/

{
  int blk_col;

  if (blk_size > 0)
    blk_col = column / blk_size;
  else
    blk_col = Trilinos_Util_find_closest_not_larger(column, cnptr, max_blocks);

  return blk_col;

} /* find_block_col */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int Trilinos_Util_find_block_in_row(int bindx[], int bnptr[], int blk_row, int blk_col,
                         int indx[], int no_elements, double val[],
                         int blk_space, int nz_space)

/*******************************************************************************

  Search the block row 'blk_row' looking for the block column 'blk_col'. If it
  is not found, create it (and initialize it to all zeros). Return the value
  'index' where indx[index] points to the start of the block (blk_row,blk_col).

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     int, index (see explanation above).
  ============

  Parameter list:
  ===============

  val,
  bindx,
  indx,
  bnptr:           Sparse matrix arrays. See User's Guide.

  blk_row,
  blk_col:         Block indices of the block for which we are looking.

  no_elements:     Number of elements in current block.

  blk_space:       Length of storage allocated for bindx[] and indx[]. An error
                   message will be printed if we try to write past these arrays.

  nz_space:        Length of storage allocated for val[]. An error message will
                   be printed if we try to write past this array.

*******************************************************************************/

{

  /* local variables */

  int   ii, k;

  const char *yo = "find_block_in_row: ";

  /**************************** execution begins ******************************/

  /* look in row 'blk_row' for 'blk_col' */

  for (k = bnptr[blk_row]; k < bnptr[blk_row + 1]; k++) {
    if (bindx[k] == blk_col) return k;
  }

  /* block was not found so let us create a new block */

  if (bnptr[blk_row + 1] + 2 >= blk_space) {
    (void) printf( "%sERROR: not enough space for block ptrs (indx)\n",
                   yo);
    exit(-1);
  }

  if (indx[bnptr[blk_row + 1]] + no_elements >= nz_space) {
    (void) printf( "%sERROR: not enough space for nonzeros (val)\n",
                   yo);
    exit(-1);
  }

  /* create the block (blk_row, blk_col) */

  bindx[bnptr[blk_row + 1]]    = blk_col;
  indx[bnptr[blk_row + 1] + 1] = no_elements + indx[bnptr[blk_row + 1]];

  for (ii = 0; ii < no_elements; ii++) val[ii+indx[bnptr[blk_row + 1]]] = 0.0;
  bnptr[blk_row + 1]++;
  return (bnptr[blk_row + 1] - 1);

} /* find_block_in_row */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void Trilinos_Util_add_new_ele(int cnptr[], int col, int blk_row, int bindx[], int bnptr[],
                    int indx[], double val[], int row, double new_ele,
                    int maxcols, int blk_space, int nz_space, int blk_type)

/*******************************************************************************

  Given a new element 'new_ele' (whose real row and column indices are given by
  'row' and 'col') store it in the VBR matrix given by cnptr[], bindx[],
  bnptr[],indx[], and val[].

  If the new element is in a block that already exists in the data structure,
  then we just add the new entry in 'val[]'. However, if the new element is in a
  block which does not already exist, we must create the new block and return a
  pointer to it (find_block_in_row(...)) before we can add the new entry to
  'val[]'.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val,
  bindx,
  indx,
  bnptr:           Sparse matrix arrays. See User's Guide.

  blk_row:         Block indices of the block for which we are looking.

  row, col:        Point row and column of new entry.

  new_ele:         New value to be placed in matrix.

  maxcols:         Total number of block columns in the local submatrix stored
                   on this processor.

  blk_space:       Length of storage allocated for bindx[] and indx[]. An error
                   message will be printed if we try to write past these arrays.

  nz_space:        Length of storage allocated for val[]. An error message will
                   be printed if we try to write past this array.

  blk_type:        blk_type > 0 ==> all the blocks are the same size so we can
                   use a shortcut in computing the block column index.
                   blk_type = 0 ==> short cut not used.

*******************************************************************************/

{

  /* local variables */

  int  blk_col, no_elements, k, start_location, little_col, little_row, offset;

  /*---------------------- execution begins -----------------------------*/

  /* find block column containing 'col' */

  blk_col = Trilinos_Util_find_block_col(cnptr, col, maxcols, blk_type);

  /* compute number of elements in block containing new point */

  no_elements = (cnptr[blk_col + 1] - cnptr[blk_col]) *
    (cnptr[blk_row + 1] - cnptr[blk_row]);

  /*
   * Search the block row looking for 'blk_col'. If it does not exist, create it
   * (and initialize it to all zeros). Return a ptr (actually an index into
   * indx[]) to the block corresponding to (blk_row,blk_col)
   */

  k = Trilinos_Util_find_block_in_row(bindx, bnptr, blk_row, blk_col, indx, no_elements,
                           val, blk_space, nz_space);

  /* compute the location of the new element in val[] */

  start_location = indx[k];
  little_col     = col - cnptr[blk_col];
  little_row     = row - cnptr[blk_row];

  offset = little_col * (cnptr[blk_row + 1] - cnptr[blk_row]) + little_row;
  val[start_location + offset] = new_ele;

} /* add_new_ele */
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int Trilinos_Util_find_closest_not_larger(int key, int list[], int length)

/*******************************************************************************

  Find the closest number to 'key' in 'list' which is not greater than key and
  return the index number.

  On exit, Trilinos_Util_find_index() returns: i => list[i] = key or list[i] is closest
  number smaller than key.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     int, see explanation above.
  ============

  Parameter list:
  ===============

  key:             Element to be search for in list.

  list:            List (assumed to be in ascending order) to be searched.

  length:          Length of list.

*******************************************************************************/

{

  /* local variables */

  int mid, start, end;

  /**************************** execution begins ******************************/

  if (length == 0) return -1;

  start = 0;
  end   = length - 1;

  while (end - start > 1) {
    mid = (start + end) / 2;
    if (list[mid] > key) end = mid;
    else start = mid;
  }

  if (list[end] > key) return start;
  else return end;

} /* Trilinos_Util_find_closest_not_larger */
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void Trilinos_Util_convert_values_to_ptrs(int array[], int length, int start)

/*******************************************************************************
  Change 'array[]' so that on exit
     1) array[0] = start
     2) array[i+1] - array[i] = value of array[i] on entry to this routine.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  array:           On entry to this routine and array[0] = start.
                   On output, array[i+1] - array[i] = value of array[i].

  length:          Length of array[].

*******************************************************************************/
{

  /* local variables */

  int i;

  /**************************** execution begins ******************************/
  for (i = 1; i < length; i++) array[i] += array[i - 1];
  for (i = length; i > 0; i--)  array[i]  = array[i - 1] + start;

  array[0] = start;

} /* AZ_convert_values_to_ptrs */
