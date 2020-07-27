/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/**************************************************************************/
/* Program to perform orthogonal range searches and nearest neighbor      */
/* queries in a more sophisticated k-d tree.  In this implementation the,  */
/* nodes on any given level of the tree do not have the same              */
/* discriminating dimension as the discrimiator is chosen based on the    */
/* dimension with   the "maxspead."                                       */
/*                                                                        */
/* References:  J.H. Friedman, J.L. Bentley, R.A. Finkel  "An Algorithm   */
/* for Finding Best Matches in Logarithmic Expected Time."                */
/* ACM Transactions on Mathematical Software, Vol 3 No. 3 Sept. 1977      */
/* pp. 209-226.                                                           */
/**************************************************************************/

/* Turn off assertions here since SEACAS compiles don't normally do it. */
#ifndef NDEBUG
#define NDEBUG
#endif

#if defined ADDC_
#define KDRECTQUERY kdrectquery_
#define KDKILLTREE kdkilltree_
#define KDBUILDTREE kdbuildtree_
#else
#define KDRECTQUERY kdrectquery
#define KDKILLTREE kdkilltree
#define KDBUILDTREE kdbuildtree
#endif

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <assert.h>

#if defined(Build64)
#define real double
#else
#define real float
#endif

#include "optkd.h"

/* Used to create a new tree in the k-d tree */
#define TESTTREE(PP) ((PP) = (optkdNode *)malloc(sizeof(optkdNode)))
#define NEWTREE(PP)                                                                                \
  if (TESTTREE(PP) == NULL) {                                                                      \
    printf("memory error\n");                                                                      \
    return 0;                                                                                      \
  }

static optkdNode *Root = NULL;
static int *      perm = NULL; /* permutation array */

extern double fabs();

/***************************************************************************/
/* Makes the perm partition the array Values along the element k.          */
/* Adapted from Sedgewick's Algorithms in C (p. 128)                       */
/***************************************************************************/
#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#define sign(x) ((x) >= 0 ? 1 : -1)

/*
 * Floyd and Rivest SELECT Algorithm (CACM Algorithm 489)
 */
void rf_select(real *a, int count, int l, int r, int k, int discrim)
{
  real z, t;
  int  n, i, j, s, sd, ll, rr, tmp;

  while (r > l) {
    if (r - l > 600) {
      /*
       * Use SELECT recursively on a sample of size S to get an estimate
       * for the (K-L+1)-th smallest element into A[K], biased slightly so
       * that the (K-L+1)-th element is expected to lie in the smaller set
       * after partitioning
       */
      n  = r - l + 1;
      i  = k - l + 1;
      z  = log(n);
      s  = 0.5 * exp(2.0 * z / 3.0);
      sd = 0.1 * sqrt(z * s * (n - s) / n) * sign(i - n / 2.0);

      tmp = k - (i * s / n) + sd;
      ll  = max(l, tmp);
      tmp = k + (n - i) * s / n + sd;
      rr  = min(r, tmp);
      rf_select(a, count, ll, rr, k, discrim);
    }
    t = a[discrim * count + perm[k]];

    /* The following code partitions a[L:R] about T. It is similar to
       PARTITION, but will run faster on most machines since subscript
       range checking on I and J has been eliminated.
    */
    i = l;
    j = r;
    /* exchange(x[l], x[k]); */
    tmp     = perm[l];
    perm[l] = perm[k];
    perm[k] = tmp;
    if (a[discrim * count + perm[r]] > t) {
      /* exchange(x[r], x[l]); */
      tmp     = perm[r];
      perm[r] = perm[l];
      perm[l] = tmp;
    }
    while (i < j) {
      /* exchange(x[i], x[j]); */
      tmp     = perm[i];
      perm[i] = perm[j];
      perm[j] = tmp;
      i       = i + 1;
      j       = j - 1;
      while (a[discrim * count + perm[i]] < t)
        i++;
      while (a[discrim * count + perm[j]] > t)
        j--;
    }
    if (a[discrim * count + perm[l]] == t) {
      /* exchange(x[l], x[j]); */
      tmp     = perm[l];
      perm[l] = perm[j];
      perm[j] = tmp;
    }
    else {
      j++;
      /* exchange(x[j], x[r]); */
      tmp     = perm[j];
      perm[j] = perm[r];
      perm[r] = tmp;
    }
    /* Now adjust L, R so they surround the subset containing the (k-l+1)-th element */
    if (j <= k)
      l = j + 1;
    if (k <= j)
      r = j - 1;
  }
}

/****************************************************************************/
int findmaxspread(int l, int u, int dimension, real *points, int N)
{
  int  maxdim    = 0;
  real max       = -FLT_MAX;
  real min       = FLT_MAX;
  real maxspread = -FLT_MAX;
  int  i, j;

  for (i = 0; i < dimension; i++) {
    max = -FLT_MAX;
    min = FLT_MAX;
    for (j = l; j <= u; j++) {
      real val = points[N * i + perm[j]];
      if (max < val) {
        max = val;
      }
      if (min > val) {
        min = val;
      }
    }
    if (maxspread < fabs(max - min)) {
      maxspread = fabs(max - min);
      maxdim    = i;
    }
  }
  return (maxdim);
}

int findmaxvariance(int l, int u, int dimension, real *points, int N)
{
  int  maxdim  = 0;
  real max_var = 0.0;
  int  i, j;

  for (i = 0; i < dimension; i++) {
    real prev_mean = 0.0;
    real mean      = 0.0;
    real variance  = 0.0;
    real count     = 0.0;

    for (j = l; j <= u; j++) {
      real val  = points[N * i + perm[j]];
      prev_mean = mean;
      count += 1.0;
      mean += (val - prev_mean) / count;
      variance += (val - prev_mean) * (val - mean);
    }
    /* True variance would be divided by (count-1.0) */

    if (max_var < variance) {
      max_var = variance;
      maxdim  = i;
    }
  }
  return (maxdim);
}

/*******************************************************************************/
int check(real *points, int N, int l, int u, int m, int discrim)
{
  int i;
  for (i = l; i < m; i++)
    assert(points[discrim * N + perm[i]] <= points[discrim * N + perm[m]]);
  for (i = m; i < u; i++)
    assert(points[discrim * N + perm[i]] >= points[discrim * N + perm[m]]);
  return 1;
}

optkdNode *BuildkdTree(real *points, int N, int l, int u, int dimension)
{
  optkdNode *p;
  int        m;

  NEWTREE(p);
  if (u - l + 1 <= BUCKETSIZE) {
    p->bucket = 1;
    p->lopt   = l;
    p->hipt   = u;
    p->loson  = NULL;
    p->hison  = NULL;
  }
  else {
    p->bucket = 0;

#if 1
    p->discrim = findmaxspread(l, u, dimension, points, N);
#else
    p->discrim = findmaxvariance(l, u, dimension, points, N);
#endif

    m = (l + u) / 2;

    rf_select(points, N, l, u, m, p->discrim);
    assert(check(points, N, l, u, m, p->discrim));

    p->cutval = points[p->discrim * N + perm[m]];
    p->loson  = BuildkdTree(points, N, l, m, dimension);
    p->hison  = BuildkdTree(points, N, m + 1, u, dimension);
  }
  return (p);
}

/*******************************************************************************/
void KDBUILDTREE(real *points, int *numPoints, int *dimension)
{

  int j;

  /* initialize perm array */
  assert(perm == NULL);
  perm = (int *)malloc(*numPoints * sizeof(int));
  assert(perm != NULL);
  for (j = 0; j < *numPoints; j++) {
    perm[j] = j;
  }
  assert(Root == NULL);
  Root = BuildkdTree(points, *numPoints, 0, *numPoints - 1, *dimension);
  assert(Root != NULL);
}

/***************************************************************************/
void KillOptTree(optkdNode *P)
{
  /*  Kills a kd-tree to avoid memory holes.   */
  if (P == NULL) {
    return;
  } /* just to be sure */

  if (P->loson != NULL) {
    KillOptTree(P->loson);
  }

  if (P->hison != NULL) {
    KillOptTree(P->hison);
  }

  free(P);
}

void KDKILLTREE()
{
  if (perm != NULL) {
    free(perm);
    perm = NULL;
  } /* free permutation array */

  KillOptTree(Root);
  Root = NULL;
}

/***************************************************************************/
/* Determines if the treenode P falls inside the rectangular query         */
/* xmin,xmax.  If so, adds the array index of the point to the found       */
/* array.                                                                  */
/***************************************************************************/
void optInRegion(optkdNode *P, int Dimension, real *Points, int N, real *xmin, real *xmax,
                 int *found, int *count)
{
  int index, dc, InsideRange;

  for (index = P->lopt; index <= P->hipt; index++) {
    InsideRange = 1;

    if (Dimension == 3) {
      int inval = perm[index];
      if (Points[0 * N + inval] < xmin[0] || Points[0 * N + inval] > xmax[0] ||
          Points[1 * N + inval] < xmin[1] || Points[1 * N + inval] > xmax[1] ||
          Points[2 * N + inval] < xmin[2] ||
          Points[2 * N + inval] > xmax[2]) { /* P is not in the region */
        InsideRange = 0;
      }
    }
    else {
      for (dc = 0; dc < Dimension; dc++) {
        if (Points[dc * N + perm[index]] < xmin[dc] ||
            Points[dc * N + perm[index]] > xmax[dc]) { /* P is not in the region */
          InsideRange = 0;
          break;
        }
      }
    }
    if (InsideRange) {
      found[(*count)++] = perm[index] + 1;
    }
  }
}

/***************************************************************************/
/* Adds the array index of each point in the bucket the point to the found */
/* array.  There is no need to check if the points are in it because we    */
/* have proven so already.                                                 */
/***************************************************************************/
void optAddRegion(optkdNode *P, int Dimension, int *found, int *count)
{
  int index;

  for (index = P->lopt; index <= P->hipt; index++) {
    found[(*count)++] = perm[index] + 1;
  }
}

/***************************************************************************/
/* Returns true iff the hyper-rectangle defined by bounds array B          */
/* intersects the rectangular query xmin,xmax.                             */
/***************************************************************************/
int optBoundsIntersectRegion(real *B, real *xmin, real *xmax, int Dimension)
{
  if (Dimension == 3) {
    if (B[0] > xmax[0] || B[1] < xmin[0] || B[2] > xmax[1] || B[3] < xmin[1] || B[4] > xmax[2] ||
        B[5] < xmin[2]) {
      return 0;
    }
  }
  else {
    int dc;

    for (dc = 0; dc < Dimension; dc++) {
      if (B[2 * dc] > xmax[dc] || B[2 * dc + 1] < xmin[dc]) {
        return (0);
      }
    }
  }
  return (1);
}

/***************************************************************************/
/* Returns true iff the hyper-rectangle defined by bounds array B          */
/* is completely contained inside the rectangular query xmin,xmax.         */
/***************************************************************************/
int optBoundsContainsRegion(real *B, real *xmin, real *xmax, int Dimension)
{
  if (Dimension == 3) {
    if (xmin[0] <= B[0] && xmax[0] >= B[1] && xmin[1] <= B[2] && xmax[1] >= B[3] &&
        xmin[2] <= B[4] && xmax[2] >= B[5])
      return 1;
    else
      return 0;
  }
  else {
    int dc;

    for (dc = 0; dc < Dimension; dc++) {
      if (!(B[2 * dc] >= xmin[dc] && B[2 * dc + 1] <= xmax[dc])) {
        return (0);
      }
    }
  }
  return (1);
}

#define MAX_DEPTH 60
/***************************************************************************/
void optRangeSearch(optkdNode *P, real *Points, int N, int Dimension, real *xmin, real *xmax,
                    real *B, int *found, int *count, int depth)
{
  static real BLow[MAX_DEPTH][6];
  static real BHigh[MAX_DEPTH][6];

  int disc;

  if (depth >= MAX_DEPTH) {
    fprintf(stderr, "Internal Error in optRangeSearch -- recursion depth too large.\n");
    abort();
  }

  if (P == NULL) {
    fprintf(stderr, "Internal Error in optRangeSearch -- P is NULL.\n");
    abort();
  }

  if (P->bucket) {
    if (optBoundsContainsRegion(B, xmin, xmax, Dimension)) {
      optAddRegion(P, Dimension, found, count);
    }
    else {
      optInRegion(P, Dimension, Points, N, xmin, xmax, found, count);
    }
    return;
  }

  /* Claim: P is not a bucket node */
  disc = P->discrim;
  /* copy the region B into BLow, BHigh */

  if (Dimension == 3) {
    memcpy(&BLow[depth][0], &B[0], 6 * sizeof(real));
    memcpy(&BHigh[depth][0], &B[0], 6 * sizeof(real));
  }
  else {
    memcpy(&BLow[depth][0], &B[0], 2 * Dimension * sizeof(real));
    memcpy(&BHigh[depth][0], &B[0], 2 * Dimension * sizeof(real));
  }

  /* Improve the Bounds for the subtrees */
  BLow[depth][2 * disc + 1] = P->cutval;
  BHigh[depth][2 * disc]    = P->cutval;

  if (optBoundsIntersectRegion(&BLow[depth][0], xmin, xmax, Dimension)) {
    optRangeSearch(P->loson, Points, N, Dimension, xmin, xmax, &BLow[depth][0], found, count,
                   depth + 1);
  }
  if (optBoundsIntersectRegion(&BHigh[depth][0], xmin, xmax, Dimension)) {
    optRangeSearch(P->hison, Points, N, Dimension, xmin, xmax, &BHigh[depth][0], found, count,
                   depth + 1);
  }
}

/***************************************************************************/
void KDRECTQUERY(real *Points, int *N, int *Dimension, real *xmin, real *xmax, int *found,
                 int *count)
{
  static real B[6];

  switch (*Dimension) {
  case 3: B[5] = FLT_MAX; B[4] = -FLT_MAX; /* fall through */
  case 2: B[3] = FLT_MAX; B[2] = -FLT_MAX; /* fall through */
  case 1: B[1] = FLT_MAX; B[0] = -FLT_MAX;
  }

  *count = 0;
  optRangeSearch(Root, Points, *N, *Dimension, xmin, xmax, B, found, count, 0);
}
