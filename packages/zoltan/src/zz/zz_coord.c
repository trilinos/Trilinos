/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <math.h>

#include "zz_const.h"
#include "params_const.h"
#include "rcb.h"
#include "rib.h"
#include "hsfc.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* 
 * This function gets a list of coordinates one way or the other,
 * i.e., by calling either Get_Geom_Multi or Get_Geom for each object.
 * Degenerate geometry detection and transformations may be added here later.
 */

static PARAM_VARS Skip_Dim_Params[] = {
                  { "KEEP_CUTS", NULL, "INT", 0 },
                  { "SKIP_DIMENSIONS", NULL, "INT", 0 },
                  { NULL, NULL, NULL, 0 } };

static void inertial_matrix(ZZ *zz, double *X, int num_obj, 
  double *cm, double (*im)[3]) ;
static int eigenvectors(double (*m)[3], double (*evecs)[3]);
static int tqli(double *d, double *e, double(*z)[3]);
static void tred2(double (*a)[3], double *d, double *e);
static void projected_distances(ZZ *zz, double *coords, int num_obj,
        double *cm, double (*evecs)[3], double *d);
static void order_decreasing(double *d, int *order);

int Zoltan_Get_Coordinates(
  ZZ *zz, 
  int num_obj,               /* Input:  number of objects */
  ZOLTAN_ID_PTR global_ids,  /* Input:  global IDs of objects */
  ZOLTAN_ID_PTR local_ids,   /* Input:  local IDs of objects; may be NULL. */
  int *num_dim,              /* Output: dimension of coordinates */
  double **coords            /* Output: array of coordinates; malloc'ed by
                                        fn if NULL upon input. */
)
{
  char *yo = "Zoltan_Get_Coordinates";
  int i,j,rc;
  int num_gid_entries = zz->Num_GID;
  int num_lid_entries = zz->Num_LID;
  int alloced_coords = 0;
  ZOLTAN_ID_PTR lid;   /* Temporary pointers to local IDs; used to pass 
                          NULL to query functions when NUM_LID_ENTRIES == 0. */
  double cm[3], dist[3];
  double evecs[3][3];
  double M[3][3], im[3][3];
  double (*sav)[3];
  double max_ratio = 10.0;
  double x, y, *cold;
  int order[3];
  int keep_cuts, skip_dimensions;
  int num_skip_dimensions = 0;
  RCB_STRUCT *rcb;
  RIB_STRUCT *rib;
  HSFC_Data *hsfc;
  int ierr = ZOLTAN_OK;
  char msg[256];

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Error check -- make sure needed callback functions are registered. */

  if (zz->Get_Num_Geom == NULL || 
     (zz->Get_Geom == NULL && zz->Get_Geom_Multi == NULL)) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Must register ZOLTAN_NUM_GEOM_FN and "
                       "either ZOLTAN_GEOM_MULTI_FN or ZOLTAN_GEOM_FN");
    goto End;
  }

  /* Get problem dimension. */

  *num_dim = zz->Get_Num_Geom(zz->Get_Num_Geom_Data, &ierr);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                       "Error returned from ZOLTAN_GET_NUM_GEOM_FN");
    goto End;
  }
  if (*num_dim < 0 || *num_dim > 3) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                       "Invalid dimension returned from ZOLTAN_NUM_GEOM_FN");
    goto End;
  }

  /* Get coordinates for object; allocate memory if not already provided. */

  if (*num_dim > 0 && num_obj > 0) {
    if (*coords == NULL) {
      alloced_coords = 1;
      *coords = (double *) ZOLTAN_MALLOC(num_obj * (*num_dim) * sizeof(double));
      if (*coords == NULL) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error");
        goto End;
      }
    }

    if (zz->Get_Geom_Multi != NULL) {
      zz->Get_Geom_Multi(zz->Get_Geom_Multi_Data, zz->Num_GID, zz->Num_LID,
                         num_obj, global_ids, local_ids, *num_dim, *coords,
                         &ierr);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                           "Error returned from ZOLTAN_GET_GEOM_MULTI_FN");
        goto End;
      }
    }
    else {
      for (i = 0; i < num_obj; i++) {
        lid = (num_lid_entries ? &(local_ids[i*num_lid_entries]) : NULL);
        zz->Get_Geom(zz->Get_Geom_Data, num_gid_entries, num_lid_entries,
                     global_ids + i*num_gid_entries, lid, 
                     (*coords) + i*(*num_dim), &ierr);
        if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                             "Error returned from ZOLTAN_GET_GEOM_FN");
          goto End;
        }
      }
    }
  }

  /*
   * For 3D RCB, RIB, and HSFC: if SKIP_DIMENSIONS was selected, compute the
   * center of mass and inertial matrix of the coordinates.  If the
   * geometry is "flat", transform the points so the primary directions
   * lie along the X and Y coordinate axes and project to the Z=0 plane.
   * If in addition the geometry is "skinny", project to the X axis.
   * Return these points to the partitioning algorithm, in effect partitioning
   * in only the 2 (or 1) significant dimensions.  Save the transformation 
   * matrix if KEEP_CUTS is true.
   */

  if ((*num_dim == 3) && 
      ((zz->LB.Method==RCB) || (zz->LB.Method==RIB) || (zz->LB.Method==HSFC))){

    Zoltan_Bind_Param(Skip_Dim_Params, "KEEP_CUTS", (void *)&keep_cuts);
    Zoltan_Bind_Param(Skip_Dim_Params, "SKIP_DIMENSIONS", (void *)&skip_dimensions);

    keep_cuts = 0;
    skip_dimensions = 0;
    Zoltan_Assign_Param_Vals(zz->Params, Skip_Dim_Params, zz->Debug_Level,
      zz->Proc, zz->Debug_Proc);

    if (skip_dimensions){

      /*
       * Get the center of mass and inertial matrix of coordinates.  Ignore
       * vertex weights, we are only interested in geometry.  Global operation.
       */
      inertial_matrix(zz, *coords, num_obj, cm, im); 

      /*
       * The inertial matrix is a 3x3 real symmetric matrix.  Get its
       * three orthonormal eigenvectors.  These usually indicate the 
       * orientation of the geometry.
       */

      rc = eigenvectors(im, evecs);

      if (rc){
        /* eigenvector calculation failed */
        skip_dimensions = 0;
        goto End;
      }

      /*
       * Calculate the extent of the geometry along the three lines defined
       * by the direction of the eigenvectors through the center of mass.
       * (Are the eigenvector magnitudes enough here?)
       */

      projected_distances(zz, *coords, num_obj, cm, evecs, dist); 

      /*
       * Decide whether these distances indicate the geometry is
       * very flat in one or two directions.
       */
  
      order_decreasing(dist, order);
  
      if ((dist[order[0]] / dist[order[2]]) > max_ratio){
        /*
         * We'll rotate geometry so it's aligned with the X/Y plane
         * and project to Z=0.
         */
        num_skip_dimensions = 1;

        if ((dist[order[0]] / dist[order[1]]) > max_ratio){
          /*
           * We'll rotate geometry so it's aligned with the X-axis
           * and project to the X-axis.
           */
          num_skip_dimensions = 2;
        }
      }

      if (num_skip_dimensions > 0){
 
        if ((zz->Debug_Level > 0) && (zz->Proc == 0)){
          sprintf(msg,
           "Geometry approx. %lf x %lf x %lf, we'll treat it as %d dimensional",
           dist[0], dist[1], dist[2],
           ((num_skip_dimensions == 1) ? 2 : 1));

          ZOLTAN_PRINT_INFO(zz->Proc, yo, msg);
        }

        /*
         * Reorder the eigenvectors (they're the columns of evecs) from 
         * longest projected distance to shorted projected distance.  Compute
         * the transpose (the inverse) of the matrix.  This will transform
         * the geometry to align along the X-Y plane.  
         */
        for (i=0; i<2; i++){
          for (j=0; j<3; j++){
            M[i][j] = evecs[j][order[i]];
          }
        }

        for (i=0, cold = *coords; i<num_obj; i++, cold += 3){
          x = M[0][0]*cold[0] + M[0][1]*cold[1] + M[0][2]*cold[2];
          if (num_skip_dimensions == 1){
            /*
             * Orient points to lie along XY plane, major direction along X,
             * secondary along Y.  So it's mostly flat along Z.
             */
            y = M[1][0]*cold[0] + M[1][1]*cold[1] + M[1][2]*cold[2];
          } 
          else{
            /*
             * Orient points to lie along X axis.  Mostly flat along Y and Z.
             */
            y = 0.0;
          }

          cold[0] = x;
          cold[1] = y;
          cold[2] = 0;
        }
         
        if (keep_cuts){
          if (zz->LB.Method == RCB){
            rcb = (RCB_STRUCT *)zz->LB.Data_Structure;
            rcb->Skip_Dimensions = num_skip_dimensions;
            sav = rcb->Transformation;
          }
          else if (zz->LB.Method == RIB){
            rib = (RIB_STRUCT *)zz->LB.Data_Structure;
            rib->Skip_Dimensions = num_skip_dimensions;
  
            sav = rib->Transformation;
          }
          else if (zz->LB.Method == HSFC){
            hsfc = (HSFC_Data *)zz->LB.Data_Structure;
            hsfc->Skip_Dimensions = num_skip_dimensions;
            sav = hsfc->Transformation;
          }

          for (i=0; i<3; i++){
            for (j=0; j<3; j++){
              sav[i][j] = M[i][j];
            }
          }
        }

      } /* If geometry is very flat */
    }  /* If SKIP_DIMENSIONS is true */
  } /* If 3-D rcb, rib or hsfc */

End:
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error found; no coordinates returned.");
    if (alloced_coords) ZOLTAN_FREE(coords);
  }
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

static void inertial_matrix(ZZ *zstruct, double *X, int num_obj, double *cm, double (*im)[3]) 
{
  double    tmp1[6], tmp2[6];
  double    xx, yy, zz, xy, xz, yz;
  double    xdif, ydif, zdif;
  int       j, rank;
  double    cmt[3];
  double    xxt, yyt, zzt, xyt, xzt, yzt;
  double *c, num_coords, total_coords;

  int comm = zstruct->Communicator;
  int proc = zstruct->Proc;
  int nproc = zstruct->Num_Proc;
  int proclower = 0;
  num_coords = (double)num_obj;

  /* Much of this was taken from Zoltan_RIB_inertial3d() */

  cm[0] = cm[1] = cm[2] = 0.0; 
  for (j = 0, c = X; j < num_obj; j++, c += 3) {
    cm[0] += c[0];
    cm[1] += c[1];
    cm[2] += c[2];
  }

  if (zstruct->Tflops_Special) {
     rank = proc - proclower;
     Zoltan_RIB_reduce_double(cm, cmt, 3, comm, nproc, rank, proc, 1);
     Zoltan_RIB_reduce_double(&num_coords, &total_coords, 1, comm, nproc, rank, proc, 1);
  }
  else {
     MPI_Allreduce(cm,cmt,3,MPI_DOUBLE,MPI_SUM,comm);
     MPI_Allreduce(&num_coords,&total_coords,1,MPI_DOUBLE,MPI_SUM,comm);
  }

  /* Global center of mass */
  cm[0] = cmt[0]/total_coords;
  cm[1] = cmt[1]/total_coords;
  cm[2] = cmt[2]/total_coords;

  xx = yy = zz = xy = xz = yz = 0.0;
  for (j = 0, c = X; j < num_obj; j++, c += 3) {
     xdif = c[0] - cm[0];
     ydif = c[1] - cm[1];
     zdif = c[2] - cm[2];
     xx += xdif*xdif;
     yy += ydif*ydif;
     zz += zdif*zdif;
     xy += xdif*ydif;
     xz += xdif*zdif;
     yz += ydif*zdif;
  }

  if (zstruct->Tflops_Special) {
     tmp1[0] = xx; tmp1[1] = yy; tmp1[2] = zz;
     tmp1[3] = xy; tmp1[4] = xz; tmp1[5] = yz;
     Zoltan_RIB_reduce_double(tmp1, tmp2, 6, comm, nproc, rank, proc, 1);
     xxt = tmp2[0]; yyt = tmp2[1]; zzt = tmp2[2];
     xyt = tmp2[3]; xzt = tmp2[4]; yzt = tmp2[5];
  }
  else {
     tmp1[0] = xx; tmp1[1] = yy; tmp1[2] = zz;
     tmp1[3] = xy; tmp1[4] = xz; tmp1[5] = yz;
     MPI_Allreduce(tmp1, tmp2, 6, MPI_DOUBLE, MPI_SUM, comm);
     xxt = tmp2[0]; yyt = tmp2[1]; zzt = tmp2[2];
     xyt = tmp2[3]; xzt = tmp2[4]; yzt = tmp2[5];
  }

  /* Global inertial tensor matrix */

  im[0][0] = xxt;
  im[1][1] = yyt;
  im[2][2] = zzt;
  im[0][1] = im[1][0] = xyt;
  im[0][2] = im[2][0] = xzt;
  im[1][2] = im[2][1] = yzt;
}
static void projected_distances(ZZ *zz, double *coords, int num_obj,
        double *cm, double (*evecs)[3], double *d)
{
int i, j;
double val, min[3], max[3], *c, tmp;
int Tflops_Special, proc, nprocs, proclower;
MPI_Comm local_comm;

  Tflops_Special = zz->Tflops_Special;
  proc = zz->Proc;
  nprocs = zz->Num_Proc;
  proclower = 0;
  local_comm = zz->Communicator;

  for (i=0; i<3; i++){
    for (j=0, c=coords; j<num_obj; j++, c+=3){
      val = (((c[0] - cm[0]) * evecs[0][i]) +
             ((c[1] - cm[1]) * evecs[1][i]) +
             ((c[2] - cm[2]) * evecs[2][i]));

      if (j){
        if (val < min[i]) min[i] = val;
        else if (val > max[i]) max[i] = val;
      }
      else{
        min[i] = max[i] = val;
      }
    }
  }

  if (Tflops_Special){
    for (i=0; i<3; i++){
      Zoltan_RIB_min_max(min+i, max+i, proclower, proc, nprocs,
                       local_comm);
      d[i] = max[i] - min[i];
    }
  }
  else {
    for (i=0; i<3; i++){
      tmp = max[i];
      MPI_Allreduce(&tmp, max + i, 1, MPI_DOUBLE, MPI_MAX, local_comm);
      tmp = min[i];
      MPI_Allreduce(&tmp, min + i, 1, MPI_DOUBLE, MPI_MIN, local_comm);

      d[i] = max[i] - min[i];
    }
  }

  return;
}

static void order_decreasing(double *d, int *order)
{
  if (d[0] > d[1]){     /* Order from longest to shortest direction */
    if (d[0] > d[2]){
      order[0] = 0;
      if (d[1] > d[2]){
        order[1] = 1; order[2] = 2;
      }
      else{
        order[1] = 2; order[2] = 1;
      }
    }
    else{
      order[0] = 2; order[1] = 0; order[2] = 1;
    }
  }
  else{
    if (d[1] > d[2]){
      order[0] = 1;
      if (d[0] > d[2]){
        order[1] = 0; order[2] = 2; }
      else{
        order[1] = 2; order[2] = 0;
      }
    }
    else{
      order[0] = 2; order[1] = 1; order[2] = 0;
    }
  }
}
#define SIGN(a,b) ((b) < 0 ? -fabs(a) : fabs(a))
static int eigenvectors(double (*m)[3], double (*evecs)[3])
{
  /* 
   * Given a real symmetric 3x3 matrix "m", find its eigenvectors.  Put
   * the 3 orthonormal eigenvectors in the columns of the matrix "evecs".
   */

  double d[3], e[3];
  int i, j, rc;

  for (i=0; i<3; i++){
    for (j=0; j<3; j++){
      evecs[i][j] = m[i][j];
    }
  }

  tred2(evecs, d, e); 

  rc = tqli(d, e, evecs);

  return rc;
}

static void tred2(double (*a)[3],  /* Q on output */
                  double *d,       /* on output, the diagonal of T */
                  double *e) /* on output, the off-diagonal (ignore e[0])*/
{
  int l, k, j, i;
  double scale, hh, h, g, f;

  /*
   * Householder reduction from "Numerical Recipes in C". 
   * Take a real symmetric 3x3 matrix "a" and decompose it into
   * an orthogonal Q and a tridiagonal matrix T.
   */

  for (i=2; i>=1; i--){
    l = i - 1;
    h = scale = 0.0;
    if (l > 0){
      for (k=0; k<=l; k++){
        scale += fabs(a[i][k]);
      }
      if (scale == 0.0){
        e[i] = a[i][l];
      }
      else{
        for (k=0; k<=l; k++){
          a[i][k] /= scale;
          h += a[i][k] * a[i][k];
        }
        f = a[i][l];
        g = (f > 0.0 ? -sqrt(h) : sqrt(h));
        e[i] = scale * g;
        h -= f*g;
        a[i][l] = f - g;
        f = 0.0;
        for (j=0; j<=l; j++){
          a[j][i] = a[i][j] / h;
          g = 0.0;
          for (k=0; k<=j; k++){
            g += a[j][k]*a[i][k];
          }
          for (k=j+1; k<=l; k++){
            g += a[k][j]*a[i][k];
          }
          e[j] = g / h;
          f += e[j] * a[i][j];
        }
        hh = f / (h+h);
        for (j=0; j<=l; j++){
          f = a[i][j];
          e[j] = g = e[j] - (hh * f);
          for (k=0; k <= j; k++){
            a[j][k] -= ((f * e[k]) + (g * a[i][k]));
          }
        }
      }
    } else {
      e[i] = a[i][l];
    }

    d[i] = h;
  }

  d[0] = 0.0;
  e[0] = 0.0;

  for (i=0; i<3; i++){
    l = i-1;
    if (d[i]){
      for (j=0; j<=l; j++){  /* skip when i == 0 */
        g = 0.0;
        for (k=0; k <= l; k++){
          g += a[i][k] * a[k][j];
        }
        for (k=0; k <= l; k++){
          a[k][j] -= g * a[k][i];
        }
      }
    }
    d[i] = a[i][i];
    a[i][i] = 1.0;
    for (j=0; j <= l; j++){
      a[j][i] = a[i][j] = 0.0;
    }
  }
}

static int tqli(double *d,     /* input from tred2, output is eigenvalues */
                 double *e,     /* input from tred2, output is garbage */
                 double(*z)[3]) /* input from tred2, output columns are e-vectors */
{
  int m, l, iter, i, k;
  double s, r, p, g, f, dd, c, b;

  /*
  ** QL algorithm with implicit shifts.  Straight out of "Numerical Recipes in C".
  */ 
  e[0] = e[1]; 
  e[1] = e[2];
  e[2] = 0.0;

  for (l=0; l<=2; l++){
    iter = 0;
    do {
      for (m=l; m <= 1; m++){
        dd = fabs(d[m]) + fabs(d[m+1]);
        if ((double)(fabs(e[m]) + dd) == dd) break;
      }
      if (m != l){
        if (iter++ == 300){
          printf("error - too many tqli iterations\n");
          return 1;
        }
        g = (d[l+1] - d[l]) / (2.0 * e[l]);
        r = sqrt((g*g) + 1.0);
        g = d[m] - d[l] + (e[l] / (g + SIGN(r,g)));
        s = c = 1.0;
        p = 0.0;
        for (i=m-1; i >= l; i--){
          f = s * e[i];
          b = c * e[i];
          if (fabs(f) >= fabs(g)){
            c = g / f;
            r = sqrt((c*c) + 1.0);
            e[i+1] = f * r;
            c *= (s = (1.0/r)); 
          }
          else{
            s = f/g;
            r = sqrt((s*s) + 1.0);
            e[i+1] = g * r;
            s *= (c = (1.0/r));
          }
          g = d[i+1] - p;
          r = ((d[i] - g) * s) + (2.0 * c * b);
          p = s * r;
          d[i+1] = g + p;
          g = c * r - b;

          for (k=0; k<=2; k++){
            f = z[k][i+1];
            z[k][i+1] = (s * z[k][i]) + (c * f);
            z[k][i] = (c * z[k][i]) - (s * f);
          }
       }
       d[l] = d[l] - p;
       e[l] = g;
       e[m] = 0.0;
      }
    } while (m != l);
  }
  return 0;
}
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
