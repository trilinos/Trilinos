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
#include <values.h>
#include "iohb.h"


main (int argc, char *argv[])
    {
    int nRow, nCol, nz, nEdge, nPin;
    int *rowindex = NULL, *colstart = NULL;
    double *val = NULL;
    int i, ii;
    int cnt, maxEdgeCnt = 0, minEdgeCnt = MAXINT;
    int *storage = NULL, *rowstart = NULL, *r = NULL, 
        row, column, *colindex = NULL;
    int *temp = NULL;
    int previous;
    char filename[100], filename2[100], command[200];
    FILE *hg, *g;

    readHB_newmat_double(argv[1], &nRow, &nCol, &nz, &colstart, &rowindex, &val);
    if (val) free (val);

    /* build dual description */
    storage  = calloc (nz,     sizeof (int));
    rowstart = calloc (nRow+3, sizeof (int));
    colindex = calloc (nz,     sizeof (int));

    if (storage == NULL || rowstart == NULL || colindex == NULL)
       {
       printf ("Unable to allocate storage, exiting\n");
       return;
       }

    /* count how many elements are in each row */
    r = rowstart+2;
    for (i = 0; i < nz; i++)
       ++r[rowindex[i]];

    /* partial sum computes storage index */
    for (i = 1; i <= nRow; i++)
       r[i] = r[i] + r[i-1];

    /* now save element ids by row */
    r = rowstart+1;
    for (column = 0; column < nCol; column++)
       for (i = colstart[column]; i < colstart[column+1]; i++)
          {
          colindex [r [rowindex[i]]   ] = column;
          storage  [r [rowindex[i]]++ ] = i;
          }

    /* write graph */
    sprintf (filename, "./%s.graph", argv[1]);
    g  = fopen (filename, "w");
    fprintf (g, "%60s\n", " ");

    if (nRow != nCol)
       {
       printf ("nRow must equal nCol\n");
       return;
       }

    nPin = 0;
    nEdge = 0;
    temp  = calloc (nCol, sizeof (int));
    for (row = 0; row < nRow; row++)
       {
       for (i = 0; i < nCol; i++)
          temp[i] = 0;

       for (i = rowstart[row]; i < rowstart[row+1]; i++)
          temp[colindex[i]] = 1;

       for (i = colstart[row]; i < colstart[row+1]; i++)
          temp[rowindex[i]] = 1;
          
       cnt = 0;
       for (i = 0; i < nRow; i++)
          if (temp[i] != 0) 
             {
             cnt++;
             fprintf (g, "%d ", i+1);
             }
       fprintf (g, "\n");
  
       nPin += cnt;
       if (cnt > maxEdgeCnt) maxEdgeCnt = cnt;
       if (cnt < minEdgeCnt) minEdgeCnt = cnt;
       }
    rewind (g);
    fprintf (g, "%d %d 00", nRow, nPin);    /* header line */
    fclose (g);

    free (storage);
    free (rowstart);
    free (colindex);
    free (colstart);
    free (rowindex);
    free (temp);
    printf("maxEdgeCnt = %d  minEdgeCnt = %d\n", maxEdgeCnt, minEdgeCnt);
    }
