#include <stdio.h>
#include <stdlib.h>
#include "iohb.h"

/* Karen Devine  May 21, 2003 */
/* Program to compute total communication volume for 
 * matrix-vector multiplication given an HB-formatted matrix and 
 * partition information for the rows of that matrix.
 * Matrix is in straight HB-format.  Partition information is in
 * a file that is developed from the partition output files of
 * Zoltan's zdrive file.  Details for generating the partition
 * information file are listed below (near the file read). */

/* to compile:  gcc metric.c iohb.c */

/* usage:  a.out matrix.hb partition_file */

main (int argc, char *argv[]) {
  int nRow, nCol, nz, nEdge, nPin;
  int *counts = NULL, max_count, min_count; 
  float avg_count;
  int *rowindex = NULL, *colstart = NULL;
  double *val = NULL;
  int i, j;
  FILE *fp;
  int nobj;
  int id, junk1, junk2;
  int *parts = NULL;
  int *usedpart = NULL;          /* Array for cut info for a single column. */
  int **partinfo = NULL;         /* 2D Array for cut info for each partition */
  int maxpart;
  int nborpart;
  int colpart;
  int mincut, maxcut, totalcut;  /* min, max, and total cut hyperedges */
  int numcut;
  int minnbor, maxnbor, totalnbor;  /* min, max, and total partition nbors*/
  int numnbor;
  int column;
    

  /* Read in partition information */
  /* Assume partition information is post-processed as follows:
     -  Cat all partition output files together.
     -  remove text (parameters, etc.), leaving only object and part info
        (order info values -1 may be kept in file; they will be read and 
        ignored here).
     -  "sort -n" the files by object IDs.
     -  put number of objects on first line.
   */

  printf("Reading %s\n", argv[2]);
  fp = fopen(argv[2], "r");
  fscanf(fp, "%d", &nobj);
  parts = (int *) malloc(sizeof(int) * nobj);
  maxpart = -1;
  for (i = 0; i < nobj; i++) {
    fscanf(fp, "%d %d %d %d", &id, &(parts[i]), &junk1, &junk2);
    if (parts[i] > maxpart) maxpart = parts[i];
  }
  maxpart++;  /* includes partition 0 */
  fclose(fp);

  /* Compute load imbalance */

  printf("Compute balance\n");
  counts = (int *) calloc(maxpart, sizeof(int));
  for (i = 0; i < nobj; i++)
    counts[parts[i]]++;

  min_count = nobj+1;
  max_count = -1;
  for (i = 0; i < maxpart; i++) {
    if (counts[i] < min_count) min_count = counts[i];
    if (counts[i] > max_count) max_count = counts[i];
  }

  avg_count = nobj / maxpart;
  printf("    Rows per partition:\n"
         "         Min %d  Max %d  Sum %d  Avg %.3f  Imbalance %.3f\n",
         min_count, max_count, nobj, avg_count, 
         (avg_count > 0. ? max_count / avg_count : -1.));

  free(counts);
  
  /* Read matrix information */

  printf("Reading HB matrix\n");
  readHB_newmat_double(argv[1], &nRow, &nCol, &nz, &colstart, &rowindex, &val);
  if (val) free (val);


  printf("Computing communication volume\n");

  partinfo = (int **) malloc(sizeof(int*) * maxpart);
  for (i = 0; i < maxpart; i++) 
    partinfo[i] = (int *) calloc(maxpart, sizeof(int));

  usedpart = (int *) malloc(sizeof(int) * maxpart);
  mincut = maxcut = totalcut = 0;
  for (column = 0; column < nCol; column++) {
    colpart = parts[column];
    for (i = 0; i < maxpart; i++) usedpart[i] = 0;

    for (i = colstart[column]; i < colstart[column+1]; i++) {/* i is vertex */
      nborpart = parts[rowindex[i]];  
      /* KDD  Printing column information for hand-checking small problems.
      printf("col %d on %d nbor %d on %d\n", 
              column, colpart, rowindex[i], nborpart);
      */
      if (nborpart != colpart) {
        if (usedpart[nborpart] == 0) {
          usedpart[nborpart]++;
          partinfo[colpart][nborpart]++;
        }
      }
    }
    
    numcut = 0;
    for (i = 0; i < maxpart; i++) 
      if (usedpart[i] != 0) numcut++;
    /* printf("Column %d  numcut %d\n", column, numcut); */
    
    totalcut += numcut;
  }

  maxnbor = 0;
  minnbor = maxpart+1;
  totalnbor = 0;
  maxcut = 0;
  mincut = totalcut+1;
  for (i = 0; i < maxpart; i++) {
    numnbor = 0;
    numcut = 0;
    for (j = 0; j < maxpart; j++)
      if (partinfo[i][j] != 0) {
        numcut += partinfo[i][j];
        numnbor++;
      }
    totalnbor += numnbor;
    if (numcut > maxcut) maxcut = numcut;
    if (numcut < mincut) mincut = numcut;
    if (numnbor > maxnbor) maxnbor = numnbor;
    if (numnbor < minnbor) minnbor = numnbor;
  }

  printf("    Number of Neighboring partitions:\n"
         "         Min %d  Max %d  Sum %d  Avg %.3f\n",
         minnbor, maxnbor, totalnbor, (double) totalnbor / (double) maxpart);
  printf("    Communication volume:\n"
         "         Min %d  Max %d  Sum %d  Avg %.3f\n",
         mincut, maxcut, totalcut, (double) totalcut / (double) maxpart);

  free(parts);
  free(usedpart);
  for (i = 0; i < maxpart; i++) free(partinfo[i]);
  free(partinfo);

  free(colstart);
  free(rowindex);
}
