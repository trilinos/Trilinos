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

/* Karen Devine  May 21, 2003 */
/* Program to compute total communication volume for 
 * matrix-vector multiplication given an HB-formatted matrix and 
 * partition information for the rows of that matrix.
 * Matrix is in straight HB-format.  Partition information is in
 * a file that is developed from the partition output files of
 * Zoltan's zdrive file.  Details for generating the partition
 * information file are listed below (near the file read). */

/* to compile:  gcc metric.c iohb.c */

/* usage:  a.out matrix.hb zdrive_partition_file 
      or   a.out matrix.hb hmetis_partition_file 1 */

/* Jan 22, 2004:  Added plotting capability with sorting by partition # */

int Plot = 1;
static void Zoltan_quicksort_pointer_inc_int_int(int *, int *, int *, int, int);


main (int argc, char *argv[]) {
  int nRow, nCol, nz, nEdge, nPin, njunk;
  int *counts = NULL, max_count, min_count; 
  float avg_count;
  int *rowindex = NULL, *colstart = NULL;
  double *val = NULL;
  int i, j, k;
  char ch;
  FILE *fp;
  char filename[80];
  int nobj;
  int *ids = NULL;
  int junk1, junk2;
  int *parts = NULL;
  int *usedpart = NULL;          /* Array for cut info for a single column. */
  int **partinfo = NULL;         /* 2D Array for cut info for each partition */
  int numpart;
  int nborpart;
  int colpart;
  int mincut, maxcut, totalcut;  /* min, max, and total cut hyperedges */
  int numcut;
  int minnbor, maxnbor, totalnbor;  /* min, max, and total partition nbors*/
  int numnbor;
  int column;
  int *idx = NULL, *invidx = NULL;
  int hmetis_partition_file = 0;   /* Assume zdrive-formatted partition file. */

  /* Read in partition information */
  /* Assume zdrive partition information is pre-processed as follows:
     -  If zdrive was run on one processor, no pre-processing is needed.
     -  If zdrive was run on multiple processors,
        +  Cat all partition output files together.
        +  remove text (parameters, etc.), leaving only object and part info
           (order info values -1 may be kept in file; they will be read and 
           ignored here).
        +  "sort -n" the files by object IDs.
     No pre-processing of hmetis partition information is needed.
   */

  if (argc == 4) 
    hmetis_partition_file = atoi(argv[3]);
  
  printf("Reading %s\n", argv[2]);

  /* First count vertices and junk lines */

  nobj = 0;  njunk = 0;
  fp = fopen(argv[2], "r");
  while ((ch = fgetc(fp)) != EOF) {
    if (ch >= '0' && ch <= '9')
      nobj++;
    else
      njunk++;
    while ((ch = fgetc(fp)) != '\n');
  }
  fclose(fp);


  /* Next, re-read file to get data */

  fp = fopen(argv[2], "r");
  parts = (int *) malloc(2 * sizeof(int) * nobj);
  ids = parts + nobj;
  numpart = -1;
  for (i = 0; i < njunk; i++)  /* Skip leading text lines */
    while ((ch = fgetc(fp)) != '\n');
  for (i = 0; i < nobj; i++) {
    if (hmetis_partition_file) {
      fscanf(fp, "%d", &(parts[i]));
      ids[i] = i;
    }
    else
      fscanf(fp, "%d %d %d %d", &(ids[i]), &(parts[i]), &junk1, &junk2);
    if (parts[i] > numpart) numpart = parts[i];
  }
  numpart++;  /* includes partition 0 */
  fclose(fp);

  /* Compute load imbalance */

  printf("Compute balance\n");
  counts = (int *) calloc(numpart, sizeof(int));
  for (i = 0; i < nobj; i++)
    counts[parts[i]]++;

  min_count = nobj+1;
  max_count = -1;
  for (i = 0; i < numpart; i++) {
    if (counts[i] < min_count) min_count = counts[i];
    if (counts[i] > max_count) max_count = counts[i];
  }

  avg_count = nobj / numpart;
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

  partinfo = (int **) malloc(sizeof(int*) * numpart);
  for (i = 0; i < numpart; i++) 
    partinfo[i] = (int *) calloc(numpart, sizeof(int));

  usedpart = (int *) malloc(sizeof(int) * numpart);
  mincut = maxcut = totalcut = 0;
  for (column = 0; column < nCol; column++) {
    colpart = parts[column];
    for (i = 0; i < numpart; i++) usedpart[i] = 0;

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
    for (i = 0; i < numpart; i++) 
      if (usedpart[i] != 0) numcut++;
    /* printf("Column %d  numcut %d\n", column, numcut); */
    
    totalcut += numcut;
  }

  maxnbor = 0;
  minnbor = numpart+1;
  totalnbor = 0;
  maxcut = 0;
  mincut = totalcut+1;
  for (i = 0; i < numpart; i++) {
    numnbor = 0;
    numcut = 0;
    for (j = 0; j < numpart; j++)
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
         minnbor, maxnbor, totalnbor, (double) totalnbor / (double) numpart);
  printf("    Communication volume:\n"
         "         Min %d  Max %d  Sum %d  Avg %.3f\n",
         mincut, maxcut, totalcut, (double) totalcut / (double) numpart);

  /* Generate gnuplot files */
  /* Generates one file per partition */
  if (Plot) {

    printf("\nPlotting now...\n");

    idx = (int *) malloc(2 * nRow * sizeof(int));
    invidx = idx + nRow;
    for (i = 0; i < nRow; i++) idx[i] = i;

    Zoltan_quicksort_pointer_inc_int_int(idx, parts, ids, 0, nRow-1);
    for (i = 0; i < nRow; i++) invidx[idx[i]] = i;

    for (i = 0; i < numpart; i++) {
      sprintf(filename, "mplot.%02d", i);
      fp = fopen(filename, "w");

      for (j = 0; j < nCol; j++) {
        for (k = colstart[j]; k < colstart[j+1]; k++) {
          if (parts[rowindex[k]] == i) {
/*
            fprintf(fp, "%d  %d  %d\n", j, -invidx[rowindex[k]], -rowindex[k]);
*/
            fprintf(fp, "%d  %d\n", invidx[j], -invidx[rowindex[k]]);
          }
        }
      }
      fclose(fp);
    }
    free(idx);

    sprintf(filename, "mplot.gnuload");
    fp = fopen(filename, "w");
    fprintf(fp, "set data style dots\n");
    fprintf(fp, "set nokey\n");
    fprintf(fp, "set noxtics\n");
    fprintf(fp, "set noytics\n");
    fprintf(fp, "plot ");
    for (i = 0; i < numpart; i++) {
      fprintf(fp, "\"mplot.%02d\"", i);
      if (i != numpart-1)
        fprintf(fp, ", ");
      else
        fprintf(fp, "\n");
    }
    fclose(fp);

  }

  free(parts);
  free(usedpart);
  for (i = 0; i < numpart; i++) free(partinfo[i]);
  free(partinfo);

  free(colstart);
  free(rowindex);
}

/*****************************************************************************/
/* Sorting pointers in increasing order. Sort key is int. Sub key is int. */

static void quickpart_pointer_inc_int_int (
int *sorted,
int *val1,
int *val2,
int start,
int end,
int *equal,
int *larger
)
{
int i, next, key1, key1_next, key2, key2_next;

  i = (end + start) / 2;
  key1 = val1 ? val1[sorted[i]] : 1;
  key2 = val2 ? val2[sorted[i]] : 1;

  *equal = *larger = start;
  for (i = start; i <= end; i++) {
    next = sorted[i];
    key1_next = val1 ? val1[next] : 1;
    key2_next = val2 ? val2[next] : 1;
    if (key1_next < key1 || (key1_next == key1 && key2_next < key2)) {
      sorted[i]           = sorted[*larger];
      sorted[(*larger)++] = sorted[*equal];
      sorted[(*equal)++]  = next;
    }
    else if (key1_next == key1  &&  key2_next == key2) {
      sorted[i]           = sorted[*larger];
      sorted[(*larger)++] = next;
    }
  }
}

/* 
 * Sorts in increasing order with primary key val1 and secondary key val2.
 * The arrays val1 and val2 are not rearranged; rather the index array
 * sorted is rearranged based on values in val1 and val2. 
 */

/*****************************************************************************/
static void Zoltan_quicksort_pointer_inc_int_int (
int *sorted,   /* index array that is rearranged; should be initialized
                  so that sorted[i] == i. */
int* val1,     /* array of primary key values. */
int *val2,     /* array of secondary key values. */
int start,     /* first array position to be considered for sorting. */
int end        /* last array position to be considered for sorting. */
)
{
int  equal, larger;

  if (start < end) {
    quickpart_pointer_inc_int_int (sorted,val1,val2,start,end,&equal,&larger);
    Zoltan_quicksort_pointer_inc_int_int (sorted, val1, val2, start, equal-1);
    Zoltan_quicksort_pointer_inc_int_int (sorted, val1, val2, larger, end);
  }
}

