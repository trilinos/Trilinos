/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "exodusII.h"

/*
//  Read and EXODUSII database and return a TECPLOT file
*/
void tec(int exoid, const char *filename)
{
  int     i, j, k, idum;
  int     ndim, nnode, nelem, nblk, nnset, neset, nvar, ntime, itime;
  char    title[MAX_LINE_LENGTH + 1];
  char *  nameco[3], **varnames = NULL;
  double *x[3], **q             = NULL, *time;
  int *   elem_id = NULL, *node_per_elem = NULL, *elem_per_blk = NULL, *attr_per_blk = NULL;
  int **  icon = NULL, *ic = NULL, izone;
  char ** elem_type = NULL;
  int     name_size = 0;
  FILE *  tecfile   = NULL;

  void teczone(int, int, int, char *, int, int, int *, int, double **, int, double **, FILE *);

  /*
   * FIRST, READ THE EXODUS DATA BASE
   */

  /*
   *  Open the output file, if we can
   */
  tecfile = fopen(filename, "w");
  if (tecfile == NULL) {
    printf("\nCannot open file %s for writing\n\n", filename);
    exit(1);
  }

  /*
   * Determine max name size used in databsae...
   */
  name_size = ex_inquire_int(exoid, EX_INQ_DB_MAX_USED_NAME_LENGTH);
  ex_set_max_name_length(exoid, name_size);

  /*
   *  Read database size, get coordinates and connectivity
   */
  memset(title, 0, MAX_LINE_LENGTH + 1);
  ex_get_init(exoid, title, &ndim, &nnode, &nelem, &nblk, &nnset, &neset);
  x[0] = x[1] = x[2] = NULL;
  for (i = 0; i < ndim; i++) {
    nameco[i] = (char *)malloc((name_size + 1) * sizeof(char));
    x[i]      = (double *)malloc(nnode * sizeof(double));
  }
  ex_get_coord_names(exoid, nameco);
  if (strlen(nameco[0]) == 0)
    strcpy(nameco[0], "X");
  if (strlen(nameco[1]) == 0)
    strcpy(nameco[1], "Y");
  if (ndim > 2)
    if (strlen(nameco[2]) == 0)
      strcpy(nameco[2], "Z");
  ex_get_coord(exoid, x[0], x[1], x[2]);

  elem_id       = (int *)malloc(nblk * sizeof(int));
  node_per_elem = (int *)malloc(nblk * sizeof(int));
  elem_per_blk  = (int *)malloc(nblk * sizeof(int));
  attr_per_blk  = (int *)malloc(nblk * sizeof(int));
  elem_type     = (char **)malloc(nblk * sizeof(char *));
  icon          = (int **)malloc(nblk * sizeof(int *));
  for (i = 0; i < nblk; i++)
    elem_type[i] = (char *)malloc((name_size + 1) * sizeof(char));
  ex_get_ids(exoid, EX_ELEM_BLOCK, elem_id);
  for (i = 0; i < nblk; i++) {
    ex_get_block(exoid, EX_ELEM_BLOCK, elem_id[i], elem_type[i], &elem_per_blk[i],
                 &node_per_elem[i], NULL, NULL, &attr_per_blk[i]);

    icon[i] = (int *)malloc(elem_per_blk[i] * node_per_elem[i] * sizeof(int));
    ex_get_conn(exoid, EX_ELEM_BLOCK, elem_id[i], icon[i], NULL, NULL);
  }

  /*
   *  Read time step information
   */
  ntime = ex_inquire_int(exoid, EX_INQ_TIME);
  if (ntime > 0) {
    time = (double *)malloc(ntime * sizeof(double));
    ex_get_all_times(exoid, time);
  }

  /*
   *  Read number of nodal variables and save space
   */
  nvar = 0;
  ex_get_variable_param(exoid, EX_NODAL, &nvar);
  if (nvar > 0) {
    varnames = (char **)malloc(nvar * sizeof(char *));
    q        = (double **)malloc(nvar * sizeof(double *));
    for (i = 0; i < nvar; i++) {
      varnames[i] = (char *)malloc((name_size + 1) * sizeof(char));
      q[i]        = (double *)malloc(nnode * sizeof(double));
    }
    ex_get_variable_names(exoid, EX_NODAL, nvar, varnames);
  }

  /* /////////////////////////////////////////////////////////////////////
  //  PROMPT USER FOR INFO AND WRITE TECPLOT FILE
  /////////////////////////////////////////////////////////////////////
  */

  /*
   *  Write the TECPLOT header information
   */

  assert(strlen(title) < (MAX_LINE_LENGTH + 1));
  fprintf(tecfile, "TITLE = \"%s\"\n", title);
  fprintf(tecfile, "VARIABLES = ");
  for (i = 0; i < ndim; i++) {
    fprintf(tecfile, "\"%s\"", nameco[i]);
    if (i < (ndim - 1))
      fprintf(tecfile, ", ");
  }
  if (nvar == 0)
    fprintf(tecfile, "\n");
  else
    fprintf(tecfile, ",\n            ");

  idum = 0;
  for (i = 0; i < nvar; i++) {
    idum += strlen(varnames[i]);
    assert(idum < 1022);

    fprintf(tecfile, "\"%s\"", varnames[i]);
    if (i < (nvar - 1)) {
      if ((i + 1) % 4 == 0) {
        idum = 0;
        fprintf(tecfile, ",\n            ");
      }
      else
        fprintf(tecfile, ", ");
    }
  }
  fprintf(tecfile, "\n");

  /*
   *  Select a time step
   */
  izone = 0;
  if (ntime == 0) {
    printf("\nNo solution variables available, saving mesh only\n\n");
    izone = 1;
  }
  else {
    printf("\nTime step information:\n\n");
    for (i = 0; i < ntime; i++)
      printf("   Time step %5d, time = %e\n", i + 1, time[i]);
    do {
      printf("\nSelect time step number to save,\n");
      printf("  or 0 for zone animation of all time steps: ");
      scanf("%d", &itime);
      printf("\n");
    } while (itime < 0 || itime > ntime);
    printf("\n");
    if (itime == 0)
      izone = 0;
    else
      izone = 1;
  }

  /*
   *  Write time steps
   */

  if (izone == 0) {

    /*
     *  Collapse the zones into one
     */

    /*
     *  Make sure we are using all the same element types
     *  Create one master connectivity array
     */
    for (i = 1; i < nblk; i++)
      if (strcmp(elem_type[0], elem_type[i]) != 0) {
        printf("\nCannot create zone animation because\n");
        ;
        printf("\n  there are multiple element types.");
        exit(1);
      }
    ic = (int *)malloc(nelem * node_per_elem[0] * sizeof(int));
    k  = 0;
    for (j = 0; j < nblk; j++)
      for (i = 0; i < node_per_elem[j] * elem_per_blk[j]; i++)
        ic[k++] = icon[j][i];
    assert(k == nelem * node_per_elem[0]);

    if (itime == 0) {
      for (j = 0; j < ntime; j++) {
        for (i = 0; i < nvar; i++)
          ex_get_var(exoid, j + 1, EX_NODAL, i + 1, 1, nnode, q[i]);

        i = 0;
        teczone(1, nnode, j + 1, elem_type[i], node_per_elem[i], nelem, ic, ndim, x, nvar, q,
                tecfile);
      }
      printf("\n");
    }

    free(ic);
  }
  else if (izone == 1) {

    /*
      ||  Write out each zone individually
    */
    for (i = 0; i < nvar; i++)
      ex_get_var(exoid, itime, EX_NODAL, i + 1, 1, nnode, q[i]);

    for (i = 0; i < nblk; i++)
      teczone(nblk, nnode, elem_id[i], elem_type[i], node_per_elem[i], elem_per_blk[i], icon[i],
              ndim, x, nvar, q, tecfile);
    printf("\n");
  }

  /* /////////////////////////////////////////////////////////////////////
  //  CLEAN UP
  /////////////////////////////////////////////////////////////////////
  */

  fclose(tecfile);

  /*
   *  Free up allocated memory
   */
  for (i = 0; i < ndim; i++) {
    free(nameco[i]);
    free(x[i]);
  }
  free(elem_id);
  free(node_per_elem);
  free(elem_per_blk);
  free(attr_per_blk);
  if (elem_type != NULL) {
    for (i = 0; i < nblk; i++) {
      free(elem_type[i]);
    }
    free(elem_type);
  }
  if (icon != NULL) {
    for (i = 0; i < nblk; i++) {
      free(icon[i]);
    }
    free(icon);
  }
  if (nvar > 0) {
    if (varnames != NULL) {
      for (i = 0; i < nvar; i++) {
        free(varnames[i]);
      }
      free(varnames);
    }
    if (q != NULL) {
      for (i = 0; i < nvar; i++) {
        free(q[i]);
      }
      free(q);
    }
  }
}

/*
 *   Write TECPLOT zonal information
 */
void teczone(int nblk, int nnode, int elem_id, char *elem_type, int node_per_elem, int elem_per_blk,
             int *icon, int ndim, double **x, int nvar, double **q, FILE *tecfile)
{

  int      i, j, k, n, nn, *ic = NULL;
  double * xx[3];
  double **qq = NULL;
  char     zname[80], eltype[16];
  int      ifac;
  int *    isort = NULL, inode;

  void internal_heapsort(int *, int);
  int  hunt(int, int, int *, int *);

  ifac = 1;
  if (strcmp(elem_type, "SHELL") == 0) {
    ifac = 2;
  }
  if (strcmp(elem_type, "BAR") == 0) {
    ifac = 2;
  }
  if (strcmp(elem_type, "TRUSS") == 0) {
    ifac = 4;
  }

#ifdef DEBUG
  printf("...in teczone\n");
  printf("......elem_type %s\n", elem_type);
  printf("......elem_per_blk %d\n", elem_per_blk);
#endif
  if (nblk > 1) {
    /*
     *  Create local block connectivity
     */
    ic    = (int *)malloc(ifac * node_per_elem * elem_per_blk * sizeof(int));
    isort = (int *)malloc(node_per_elem * elem_per_blk * sizeof(int));

    for (j = 0; j < node_per_elem * elem_per_blk; j++)
      isort[j] = icon[j];

#ifdef DEBUG
    printf("...sorting global connectivity\n");
    printf("......total nodes %d\n", elem_per_blk * node_per_elem);
#endif
    /*
     *  Sort nodes in this element block
     */

    internal_heapsort(isort, node_per_elem * elem_per_blk);

#ifdef DEBUG
    printf("...compressing global connectivity\n");
#endif
    /*
     *  Compress node list
     */

    inode = 0;
    for (j = 0; j < node_per_elem * elem_per_blk - 1; j++)
      if (isort[j] != isort[j + 1])
        isort[++inode] = isort[j + 1];
    inode++;

#ifdef DEBUG
    printf("......found %d unique nodes\n", inode);
    printf("...sorting local connectivity\n");
#endif
    /*
     *  Sort local connectivity
     */
    i = (inode + 1) / 2;
    for (k = 0; k < elem_per_blk; k++) {
      for (j = 0; j < node_per_elem; j++) {
        n  = k * node_per_elem + j;
        nn = k * ifac * node_per_elem + j;

        ic[nn] = hunt(icon[n], inode, isort, &i);
      }
      if (ifac == 2)
        for (j = 0; j < node_per_elem; j++) {
          nn                     = k * ifac * node_per_elem + j;
          ic[nn + node_per_elem] = ic[nn];
        }
      else if (ifac == 4)
        for (j = 0; j < node_per_elem; j++) {
          nn                         = k * ifac * node_per_elem + j;
          ic[nn + node_per_elem]     = ic[nn + node_per_elem - 1 - j];
          ic[nn + 2 * node_per_elem] = ic[nn];
          ic[nn + 3 * node_per_elem] = ic[nn + node_per_elem - 1 - j];
        }
    }

#ifdef DEBUG
    printf("...copy local data\n");
#endif
    /*
     *  Copy local data
     */
    for (i = 0; i < ndim; i++)
      xx[i] = (double *)malloc(inode * sizeof(double));

    for (j = 0; j < ndim; j++)
      for (i = 0; i < inode; i++)
        xx[j][i] = x[j][isort[i] - 1];

    if (nvar > 0)
      qq = (double **)malloc(nvar * sizeof(double *));
    for (i = 0; i < nvar; i++)
      qq[i] = (double *)malloc(inode * sizeof(double));

    for (j = 0; j < nvar; j++)
      for (i = 0; i < inode; i++)
        qq[j][i] = q[j][isort[i] - 1];
  }
  else {
    /*
     *  Copy local block connectivity
     */

    if (ifac > 1) {
      printf("\nCannot perform zone animation for degenerate elements\n\n");
      exit(1);
    }

    inode = nnode;
    ic    = icon;

    for (i = 0; i < ndim; i++)
      xx[i] = x[i];

    if (nvar > 0)
      qq = (double **)malloc(nvar * sizeof(double *));
    for (i = 0; i < nvar; i++)
      qq[i] = q[i];
  }

#ifdef DEBUG
  printf("...writing zone data\n");
#endif
  /*
   *  Write zone data
   */
  printf("    Writing Zone %4d, %8s, %6d elements, %6d nodes\n", elem_id, elem_type, elem_per_blk,
         inode);
  sprintf(zname, "Zone %s_%d", elem_type, elem_id);

  if (ndim == 3 && (node_per_elem == 8 || (node_per_elem == 4 && ifac == 2) ||
                    (node_per_elem == 2 && ifac == 4)))
    sprintf(eltype, "BRICK");
  else if (ndim == 3 && node_per_elem == 4)
    sprintf(eltype, "TETRAHEDRON");
  else if (ndim == 2 && (node_per_elem == 4 || (node_per_elem == 2 && ifac == 2)))
    sprintf(eltype, "QUADRILATERAL");
  else if (ndim == 2 && node_per_elem == 3)
    sprintf(eltype, "TRIANGLE");
  else {
    printf("\nBad element type found in teczone\n");
    printf("   Dimensions = %d\n", ndim);
    printf("   Nodes      = %d\n", node_per_elem);
    exit(1);
  }

  fprintf(tecfile, "ZONE T=\"%s\", F=FEBLOCK, ET=%s,\n", zname, eltype);
  fprintf(tecfile, "     N=%d, E=%d\n", inode, elem_per_blk);

#define CHUNK 8

  for (j = 0; j < ndim; j++) {
    for (i = 0; i < inode / CHUNK; i++) {
      for (k = 0; k < CHUNK; k++)
        fprintf(tecfile, "%.16f ", xx[j][i * CHUNK + k]);
      fprintf(tecfile, "\n");
    }
    if (inode % CHUNK != 0) {
      for (i = CHUNK * (inode / CHUNK); i < inode; i++)
        fprintf(tecfile, "%.16f ", xx[j][i]);
      fprintf(tecfile, "\n");
    }
  }

  for (j = 0; j < nvar; j++) {
    for (i = 0; i < inode / CHUNK; i++) {
      for (k = 0; k < CHUNK; k++)
        fprintf(tecfile, "%.16f ", qq[j][i * CHUNK + k]);
      fprintf(tecfile, "\n");
    }
    if (inode % CHUNK != 0) {
      for (i = CHUNK * (inode / CHUNK); i < inode; i++)
        fprintf(tecfile, "%.16f ", qq[j][i]);
      fprintf(tecfile, "\n");
    }
  }

  for (i = 0; i < elem_per_blk; i++) {
    for (j = 0; j < ifac * node_per_elem; j++) {
      k = ic[i * ifac * node_per_elem + j];
      fprintf(tecfile, "%d ", k);
    }
    fprintf(tecfile, "\n");
  }

#ifdef DEBUG
  printf("...freeing allocated space\n");
#endif
  if (nblk > 1) {
    free(ic);
    free(isort);
    for (i = 0; i < ndim; i++)
      free(xx[i]);
    if (qq != NULL) {
      for (i = 0; i < nvar; i++)
        free(qq[i]);
      if (nvar > 0)
        free(qq);
    }
  }
  else {
    if (nvar > 0)
      free(qq);
  }
}

void siftDown(int *a, int start, int count);

#define SWAP(r, s)                                                                                 \
  do {                                                                                             \
    int t = r;                                                                                     \
    r     = s;                                                                                     \
    s     = t;                                                                                     \
  } while (0)

void internal_heapsort(int *a, int count)
{
  int start, end;

  /* heapify */
  for (start = (count - 2) / 2; start >= 0; start--) {
    siftDown(a, start, count);
  }

  for (end = count - 1; end > 0; end--) {
    SWAP(a[end], a[0]);
    siftDown(a, 0, end);
  }
}

void siftDown(int *a, int start, int end)
{
  int root = start;

  while (root * 2 + 1 < end) {
    int child = 2 * root + 1;
    if ((child + 1 < end) && (a[child] < a[child + 1])) {
      child += 1;
    }
    if (a[root] < a[child]) {
      SWAP(a[child], a[root]);
      root = child;
    }
    else
      return;
  }
}

/*
 *   Bisection search with binary hunt,
 *      from Numerical Recipes
 */
int hunt(int node, int n, int *key, int *iguess)
{

  int        i = 0, j = 0, k = 0, inc = 0;
  static int hunt = 1;

  if (node == key[*iguess - 1])
    return *iguess;
  i = *iguess;

  if (hunt) {
    inc = 1;

    if (node > key[i - 1]) { /*  hunt up  */

      if (i == n) {
        printf("\nError in HUNT: node is outside range of table\n");
        printf("      requested node = %d\n", node);
        printf("      last in table  = %d\n", key[n - 1]);
        exit(1);
      }

      j = i + 1;
      while (node > key[j - 1]) {
        i = j;
        inc += inc;
        j = i + inc;
        if (j > n) {
          j = n;
          break;
        }
      }
      if (node == key[j - 1]) {
        *iguess = j;
        return j;
      }
    }
    else if (node < key[i - 1]) { /*  hunt down  */

      if (i == 1) {
        printf("\nError in HUNT: node is outside range of table\n");
        printf("      requested node  = %d\n", node);
        printf("      first in table  = %d\n", key[0]);
        exit(1);
      }

      j = i;
      i -= 1;
      while (node < key[i - 1]) {
        j = i;
        inc += inc;
        i = j - inc;
        if (i < 1) {
          i = 1;
          break;
        }
      }
      if (node == key[i - 1]) {
        *iguess = i;
        return i;
      }
    }
  }
  else {
    i = 1;
    j = n;
  }

  k = (i + j) / 2;
  while (node != key[k - 1]) {
    if (key[k - 1] > node) {
      j = k;
      k = (i + j) / 2;
    }
    else {
      i = k;
      k = (i + j + 1) / 2;
    }
  }

  *iguess = k;
  return k;
}
