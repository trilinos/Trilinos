/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * This software is distributed under the GNU Lesser General Public License. *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "iohb.h"

#define MAXA 100

/* Program to display a HB-formatted matrix.  
 * Useful mainly for debugging.
 */

main (int argc, char *argv[])
{
int nRow, nCol, nz, nEdge, nPin;
int *rowindex = NULL, *colstart = NULL;
double *val = NULL;
int i, j;
int column;
int A[MAXA][MAXA];

  readHB_newmat_double(argv[1], &nRow, &nCol, &nz, &colstart, &rowindex, &val);

  if (nCol > MAXA || nRow > MAXA) {
    printf("Matrix too big; allocate A bigger.\n");
    exit;
  }

  /* write hypergraph */
  
  for (i = 0; i < MAXA; i++)
    for (j = 0; j < MAXA; j++)
      A[i][j] = 0;
  
  for (column = 0; column < nCol; column++) 
    for (i = colstart[column]; i < colstart[column+1]; i++) 
      A[rowindex[i]][column] = 1;
        
  for (i = 0; i < nRow; i++) {
    for (j = 0; j < nCol; j++)
      if (A[i][j] == 1) 
        printf("1 ");
      else
        printf("  ");
    printf("\n");
  }
  
  free(colstart);
  free(rowindex);
  if (val) free(val);
}
