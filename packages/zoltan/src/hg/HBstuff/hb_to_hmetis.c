#include <stdio.h>
#include <stdlib.h>
#include "iohb.h"


main (int argc, char *argv[])
    {
    int nRow, nCol, nz, nEdge;
    int *rowindex, *colstart;
    double *val;
    int i, column;
    char filename[100];
    FILE *hg;

    readHB_newmat_double(argv[1], &nRow, &nCol, &nz, &colstart, &rowindex, &val);

    /* write hypergraph */
    sprintf (filename, "./%s.hgr", argv[1]);
    hg = fopen (filename, "w");
    fprintf (hg, "%60s \n", " ");              /* temporary header line */

    nEdge = 0;
    for (column = 0; column < nCol; column++)
       {
       if  ((colstart[column+1] - colstart[column]) < 2)
          continue;  /* surpress self edges */
       nEdge++;

       for (i = colstart[column]; i < colstart[column+1]; i++)
          fprintf (hg, "%d ", rowindex[i]+1);  /* write vertex */
       fprintf (hg, "\n");
       }
    rewind (hg);
    fprintf (hg, "%d %d", nEdge, nRow);        /* header line for hmetis */
    fclose (hg);
    }
