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

/* Program to create a gnuplot file of an HB-formatted matrix.
 * Useful mainly for pretty viewgraphs.
 */

main (int argc, char *argv[])
{
int nRow, nCol, nz;
int *rowindex = NULL, *colstart = NULL;
double *val = NULL;
int i;
int column;
FILE *fp;

  readHB_newmat_double(argv[1], &nRow, &nCol, &nz, &colstart, &rowindex, &val);

  /* write matrix */
  
  fp = fopen("plotfile", "w");
  for (column = 0; column < nCol; column++) 
    for (i = colstart[column]; i < colstart[column+1]; i++) 
      fprintf(fp, "%d %d\n", column, -rowindex[i]);
        
  fclose(fp);
  
  free(colstart);
  free(rowindex);
  if (val) free(val);
}
