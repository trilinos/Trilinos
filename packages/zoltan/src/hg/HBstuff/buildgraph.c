#include <stdio.h>
#include <stdlib.h>
#include "iohb.h"


main (int argc, char *argv[])
    {
    int nRow, nCol, nz, count, nEdge, nPin;
    int *rowindex, *colstart;
    double *val;
    int i, ii;
    int *storage, *rowstart, *r, row, column, *colindex;
    int *temp;
    int previous;
    char filename[100], filename2[100], command[200];
    FILE *hg, *g;

     readHB_newmat_double(argv[1], &nRow, &nCol, &nz, &colstart, &rowindex, &val);
/*  free (&val); */


    /* build dual description */
    storage  = calloc (nz,     sizeof (int));
    rowstart = calloc (nRow+3, sizeof (int));
    colindex = calloc (nz,     sizeof (int));

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

    /* write hypergraph */
    sprintf (filename, "./%s.hg.data", argv[1]);
    hg = fopen (filename, "w");

    nEdge = 0;
    nPin = 0;
    for (column = 0; column < nCol; column++)
       {
       count = 0;
       for (i = colstart[column]; i < colstart[column+1]; i++)  /* i is vertex */
          count++;
       if (count < 2)
          continue;
       nEdge++;

       for (i = colstart[column]; i < colstart[column+1]; i++)  /* i is vertex */
          {
          fprintf (hg, "%d ", rowindex[i]+1);
          nPin++;
          }
       fprintf (hg, "\n");
       }
    fclose (hg);

    sprintf (filename2, "./%s.hg.header", argv[1]);
    hg = fopen (filename2, "w");
    fprintf (hg, "%d %d %d 00\n", nRow, nEdge, nPin);        /* header line */
    fclose (hg); fflush (NULL);

    sprintf (command, "cat ./%s ./%s > ./%s.hg", filename2, filename, argv[1]);
    system (command);

    /* write graph */
    sprintf (filename, "./%s.data", argv[1]);
    if (nRow != nCol)
       {
       printf ("nRow must equal nCol\n");
       return;
       }
    g  = fopen (filename, "w");

    nPin = 0;
    nEdge = 0;
    temp  = calloc (nCol, sizeof (int));
    for (row = 0; row < nRow; row++)
       {
       count = 0;
       for (i = 0; i < nCol; i++)
          temp[i] = 0;

       for (i = rowstart[row]; i < rowstart[row+1]; i++)
          temp[colindex[i]] = 1;

       for (i = colstart[row]; i < colstart[row+1]; i++)
          temp[rowindex[i]] = 1;

       for (i = 0; i < nRow; i++)
          if (temp[i] != 0)
             {
             fprintf (g, "%d ", i+1);
             nPin++;
             count++;
             }
       fprintf (g, "\n");
       }
    fclose (g);

    sprintf (filename2, "./%s.header", argv[1]);
    g = fopen (filename2, "w");
    fprintf (g, "%d %d 00\n", nRow, nPin);    /* header line */
    fclose (g);

    sprintf (command, "cat ./%s ./%s > ./%s.graph", filename2, filename, argv[1]);
    system (command);

    free (storage);
    free (rowstart);
    }
