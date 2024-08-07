// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


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
#include "inertial.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* 
 * This function gets a list of coordinates one way or the other,
 * i.e., by calling either Get_Geom_Multi or Get_Geom for each object.
 *
 * Note that for 2D or 3D RCB, RIB and HSFC with the REDUCE_DIMENSIONS
 * option on, Zoltan_Get_Coordinates is a global operation.  (A
 * global decision to transform coordinates may be made here.)  So
 * it must be called by all processes in the application or it will
 * hang.
 */

static PARAM_VARS Reduce_Dim_Params[] = {
                  { "KEEP_CUTS", NULL, "INT", 0 },
                  { "REDUCE_DIMENSIONS", NULL, "INT", 0 },
                  { "DEGENERATE_RATIO", NULL, "DOUBLE", 0 },
                  { NULL, NULL, NULL, 0 } };

static void inertial_matrix2D(ZZ *zz, double *X, int num_obj, 
  double *cm, double (*im)[3]) ;
static void inertial_matrix3D(ZZ *zz, double *X, int num_obj, 
  double *cm, double (*im)[3]) ;
static int eigenvectors(double (*m)[3], double (*evecs)[3], int dim);
static void projected_distances(ZZ *zz, double *coords, int num_obj,
        double *cm, double (*evecs)[3], double *d, int dim, int aa, int *order);
static void order_decreasing(double *d, int *order);
static void transform_coordinates(double *coords, int num_obj, int d,
                                  ZZ_Transform *tr);
static int get_target_dimension(double *dist, int *order, 
                                double deg_ratio, int d);

#define NEAR_ONE(x) ((x >= .9999) && (x <= 1.0001))

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
  double dist[3];
  double im[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
  double deg_ratio;
  double x;
  int order[3];
  int reduce_dimensions, d, axis_aligned;
  int target_dim;
  int ierr = ZOLTAN_OK;
  char msg[256];
  ZZ_Transform *tran;

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
   * For RCB, RIB, and HSFC: if REDUCE_DIMENSIONS was selected, compute the
   * center of mass and inertial matrix of the coordinates.  
   *
   * For 3D problems: If the geometry is "flat", transform the points so the
   * two primary directions lie along the X and Y coordinate axes and project 
   * to the Z=0 plane.  If in addition the geometry is "skinny", project to 
   * the X axis.  (This creates a 2D or 1D problem respectively.)
   *
   * For 2D problems: If the geometry is essentially a line, transform it's
   * primary direction to the X axis and project to the X axis, yielding a
   * 1D problem.
   *
   * Return these points to the partitioning algorithm, in effect partitioning
   * in only the 2 (or 1) significant dimensions.  
   */

  if (((*num_dim == 3) || (*num_dim == 2)) && 
      ((zz->LB.Method==RCB) || (zz->LB.Method==RIB) || (zz->LB.Method==HSFC))){

    Zoltan_Bind_Param(Reduce_Dim_Params, "KEEP_CUTS", (void *)&i);
    Zoltan_Bind_Param(Reduce_Dim_Params, "REDUCE_DIMENSIONS", 
                     (void *)&reduce_dimensions);
    Zoltan_Bind_Param(Reduce_Dim_Params, "DEGENERATE_RATIO", (void *)&deg_ratio);

    i = 0;
    reduce_dimensions = 0;
    deg_ratio = 10.0;

    Zoltan_Assign_Param_Vals(zz->Params, Reduce_Dim_Params, zz->Debug_Level,
                             zz->Proc, zz->Debug_Proc);

    if (reduce_dimensions == 0){
      goto End;
    }

    if (deg_ratio <= 1){
      if (zz->Proc == 0){
        ZOLTAN_PRINT_WARN(0, yo, "DEGENERATE_RATIO <= 1, setting it to 10.0");
      }
      deg_ratio = 10.0;
    }

    if (zz->LB.Method == RCB){
      tran = &(((RCB_STRUCT *)(zz->LB.Data_Structure))->Tran);
    } 
    else if (zz->LB.Method == RIB){
      tran = &(((RIB_STRUCT *)(zz->LB.Data_Structure))->Tran);
    }
    else{
      tran = &(((HSFC_Data*)(zz->LB.Data_Structure))->tran);
    }

    d = *num_dim;

    if (tran->Target_Dim >= 0){
      /*
       * On a previous load balancing call, we determined whether
       * or not the geometry was degenerate.  If the geometry was 
       * determined to be not degenerate, then we assume it is still 
       * not degenerate, and we skip the degeneracy calculation.  
       */
      if (tran->Target_Dim > 0){
        /*
         * The geometry *was* degenerate.  We test the extent
         * of the geometry along the principal directions determined
         * last time to determine if it is still degenerate with that
         * orientation.  If so, we transform the coordinates using the
         * same transformation we used last time.  If not, we do the 
         * entire degeneracy calculation again.
         */
 
        if ((tran->Axis_Order[0] >= 0) && 
            (tran->Axis_Order[1] >= 0) && (tran->Axis_Order[2] >= 0)){
          axis_aligned = 1;
        }
        else{
          axis_aligned = 0;
        }

        projected_distances(zz, *coords, num_obj, tran->CM, 
             tran->Evecs, dist, d, axis_aligned, tran->Axis_Order); 

        target_dim = get_target_dimension(dist, order, deg_ratio, d);

        if (target_dim > 0){
          transform_coordinates(*coords, num_obj, d, tran);
        }
        else{
          /* Set's Target_Dim to -1, flag to recompute degeneracy */
          Zoltan_Initialize_Transformation(tran);
        }
      }
    }

    if (tran->Target_Dim < 0){

      tran->Target_Dim = 0;

      /*
       * Get the center of mass and inertial matrix of coordinates.  Ignore
       * vertex weights, we are only interested in geometry.  Global operation.
       */
      if (d == 2){
        inertial_matrix2D(zz, *coords, num_obj, tran->CM, im);
      }
      else{
        inertial_matrix3D(zz, *coords, num_obj, tran->CM, im);
      }

      /*
       * The inertial matrix is a 3x3 or 2x2 real symmetric matrix.  Get its
       * three or two orthonormal eigenvectors.  These usually indicate the 
       * orientation of the geometry.
       */

      rc = eigenvectors(im, tran->Evecs, d);

      if (rc){
        if (zz->Proc == 0){
          ZOLTAN_PRINT_WARN(0, yo, "REDUCE_DIMENSIONS calculation failed");
        }
        goto End; 
      }

      /*
       * Here we check to see if the eigenvectors are very close
       * to the coordinate axes.  If so, we can more quickly
       * determine whether the geometry is degenerate, and also more
       * quickly transform the geometry to the lower dimensional
       * space.
       */

      axis_aligned = 0;

      for (i=0; i<d; i++){
        tran->Axis_Order[i] = -1;
      }

      for (j=0; j<d; j++){
        for (i=0; i<d; i++){
          x = fabs(tran->Evecs[i][j]);

          if (NEAR_ONE(x)){
            tran->Axis_Order[j] = i;  /* e'vector j is very close to i axis */
            break;
          }
        }
        if (tran->Axis_Order[j] < 0){
          break;
        }
      }

      if ((tran->Axis_Order[0] >= 0) && 
          (tran->Axis_Order[1] >= 0) && (tran->Axis_Order[2] >= 0)){
        axis_aligned = 1;
      }

      /*
       * Calculate the extent of the geometry along the three lines defined
       * by the direction of the eigenvectors through the center of mass.
       */

      projected_distances(zz, *coords, num_obj, tran->CM, tran->Evecs, dist, 
                          d, axis_aligned, tran->Axis_Order); 

      /*
       * Decide whether these distances indicate the geometry is
       * very flat in one or two directions.
       */

      target_dim = get_target_dimension(dist, order, deg_ratio, d);

      if (target_dim > 0){
        /*
         * Yes, geometry is degenerate
         */
        if ((zz->Debug_Level > 0) && (zz->Proc == 0)){
          if (d == 2){
            sprintf(msg,
             "Geometry (~%f x %f), exceeds %f to 1.0 ratio",
              dist[order[0]], dist[order[1]], deg_ratio);
          }
          else{
            sprintf(msg,
             "Geometry (~%f x %f x %f), exceeds %f to 1.0 ratio",
              dist[order[0]], dist[order[1]], dist[order[2]], deg_ratio);
          }

          ZOLTAN_PRINT_INFO(zz->Proc, yo, msg);
          sprintf(msg, "We'll treat it as %d dimensional",target_dim);
          ZOLTAN_PRINT_INFO(zz->Proc, yo, msg);
        }

        if (axis_aligned){
          /*
          ** Create new geometry, transforming the primary direction
          ** to the X-axis, and the secondary to the Y-axis.
          */

          tran->Permutation[0] = tran->Axis_Order[order[0]];
          if (target_dim == 2){
            tran->Permutation[1] = tran->Axis_Order[order[1]];
          }
        }
        else{
          /*
           * Reorder the eigenvectors (they're the columns of evecs) from 
           * longest projected distance to shorted projected distance.  Compute
           * the transpose (the inverse) of the matrix.  This will transform
           * the geometry to align along the X-Y plane, or along the X axis. 
           */
  
          for (i=0; i< target_dim; i++){
            tran->Transformation[i][2] = 0.0;
            for (j=0; j<d; j++){
              tran->Transformation[i][j] = tran->Evecs[j][order[i]];

            }
          }
          for (i=target_dim; i< 3; i++){
            for (j=0; j<3; j++){
              tran->Transformation[i][j] = 0.0;
            }
          }
        }

        tran->Target_Dim = target_dim;

        transform_coordinates(*coords, num_obj, d, tran);

      } /* If geometry is very flat */
    }  /* If REDUCE_DIMENSIONS is true */
  } /* If 2-D or 3-D rcb, rib or hsfc */

End:
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error found; no coordinates returned.");
    if (alloced_coords) ZOLTAN_FREE(coords);
  }
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

/*
 * Initialize a coordinate transformation structure.
 */
void Zoltan_Initialize_Transformation(ZZ_Transform *tr)
{
  int i, j;
  tr->Target_Dim = -1;

  for (i=0; i<3; i++){
    for (j=0; j<3; j++){
      tr->Transformation[i][j] = 0.0;
      tr->Evecs[i][j] = 0.0;
    }
    tr->Permutation[i] = -1;
    tr->CM[i] = 0.0;
    tr->Axis_Order[i] = 0;
  } 
}

/*
 * Decide whether the relative lengths of the edges of the oriented
 * bounding box indicate the geometry is very flat in one or two directions.
 */
static int get_target_dimension(double *dist, int *order, 
                                double deg_ratio, int d)
{
  int target_dim = 0;
  double flat;
  
  if (d == 2){
    if (dist[0] < dist[1]){
      order[0] = 1; order[1] = 0;
    }
    else{
      order[0] = 0; order[1] = 1;
    }
    if (dist[order[1]] < (dist[order[0]] / deg_ratio)){
      target_dim = 1;
    }
  }
  else{
    order_decreasing(dist, order);

    flat = dist[order[0]] / deg_ratio;
  
    if (dist[order[2]] < flat){
      /*
       * We'll rotate geometry so it's aligned with the X/Y plane
       * and project to Z=0.
       */
      target_dim = 2;
  
      if (dist[order[1]] < flat){
        /*
         * We'll rotate geometry so it's aligned with the X-axis
         * and project to the X-axis.
         */
        target_dim = 1;
      }
    }
  }

  return target_dim;
}
/*
 * Apply the transformation to the coordinates.
 */
static void transform_coordinates(double *coords, int num_obj, int d,
                                  ZZ_Transform *tr)
{
  int i, a, b;
  double x, y = 0.0;

  if (tr->Permutation[0] >= 0){

    a = tr->Permutation[0];
    b = tr->Permutation[1];

    for (i=0; i<num_obj; i++, coords += d){
      x = coords[a];
      if (tr->Target_Dim == 2){
        y = coords[b];
      }

      coords[0] = x;
      coords[1] = y;
      if (d == 3) coords[2] = 0.0;
    }
  }
  else{
    for (i=0; i<num_obj; i++, coords += d){
  
      x = tr->Transformation[0][0]*coords[0] + 
          tr->Transformation[0][1]*coords[1];

      if (d == 3) x +=  tr->Transformation[0][2]*coords[2];
  
      if (tr->Target_Dim == 2){
        y = tr->Transformation[1][0]*coords[0] + 
            tr->Transformation[1][1]*coords[1]; 
        if (d == 3) y +=  tr->Transformation[1][2]*coords[2];
      } 
  
      coords[0] = x;
      coords[1] = y;
      if (d == 3) coords[2] = 0.0;
    }
  }
}
/*
 * Calculate a 2x2 inertial matrix representing the
 * locations of the 2-dimensional coordinates.
 */
static void inertial_matrix2D(ZZ *zstruct, double *X, 
                            int num_obj, double *cm, double (*im)[3])
{
  double    tmp1[3], tmp2[3];
  double    xx, yy, xy;
  double    xdif, ydif;
  int       j, rank=0;
  double    cmt[2];
  double    xxt, yyt, xyt;
  double *c, num_coords, total_coords;

  MPI_Comm comm = zstruct->Communicator;
  int proc = zstruct->Proc;
  int nproc = zstruct->Num_Proc;
  int proclower = 0;
  num_coords = (double)num_obj;

  cm[0] = cm[1] = 0.0; 
  for (j = 0, c = X; j < num_obj; j++, c += 2) {
    cm[0] += c[0];
    cm[1] += c[1];
  }

  if (zstruct->Tflops_Special) {
     rank = proc - proclower;
     Zoltan_RIB_reduce_double(cm, cmt, 2, comm, nproc, rank, proc, 1);
     Zoltan_RIB_reduce_double(&num_coords, &total_coords, 1, comm, nproc, rank, proc, 1);
  }
  else {
     MPI_Allreduce(cm,cmt,2,MPI_DOUBLE,MPI_SUM,comm);
     MPI_Allreduce(&num_coords,&total_coords,1,MPI_DOUBLE,MPI_SUM,comm);
  }

  /* Global center of mass */
  cm[0] = cmt[0]/total_coords;
  cm[1] = cmt[1]/total_coords;

  xx = yy = xy = 0.0;
  for (j = 0, c = X; j < num_obj; j++, c += 2) {
     xdif = c[0] - cm[0];
     ydif = c[1] - cm[1];
     xx += xdif*xdif;
     yy += ydif*ydif;
     xy += xdif*ydif;
  }

  if (zstruct->Tflops_Special) {
     tmp1[0] = xx; tmp1[1] = yy; tmp1[2] = xy; 
     Zoltan_RIB_reduce_double(tmp1, tmp2, 3, comm, nproc, rank, proc, 1);
     xxt = tmp2[0]; yyt = tmp2[1]; xyt = tmp2[2]; 
  }
  else {
     tmp1[0] = xx; tmp1[1] = yy; tmp1[2] = xy; 
     MPI_Allreduce(tmp1, tmp2, 3, MPI_DOUBLE, MPI_SUM, comm);
     xxt = tmp2[0]; yyt = tmp2[1]; xyt = tmp2[2];
  }

  /* Global inertial tensor matrix */

  im[0][0] = xxt;
  im[1][1] = yyt;
  im[0][1] = im[1][0] = xyt;
}
/*
 * Calculate a 3x3 inertial matrix representing the
 * locations of the 3-dimensional coordinates.
 */
static void inertial_matrix3D(ZZ *zstruct, double *X, 
                            int num_obj, double *cm, double (*im)[3])
{
  double    tmp1[6], tmp2[6];
  double    xx, yy, zz, xy, xz, yz;
  double    xdif, ydif, zdif;
  int       j, rank=0;
  double    cmt[3];
  double    xxt, yyt, zzt, xyt, xzt, yzt;
  double *c, num_coords, total_coords;

  MPI_Comm comm = zstruct->Communicator;
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
/*
 * Calculate the extent of the geometry in the directions indicated
 * by the orthonormal eigenvectors of the inertial matrix.  This is
 * essentially the dimensions of an oriented bounding box around
 * the geometry.
 */
static void projected_distances(ZZ *zz, double *coords, int num_obj,
        double *cm, double (*evecs)[3], double *d, int dim, int aa, int *order)
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

  for (j=0; j<3; j++){
    min[j] = DBL_MAX;
    max[j] = DBL_MIN;
  }

  if (aa){
    /* special case - eigenvectors are axis aligned */

    for (i=0, c = coords; i<num_obj; i++, c += dim){
      for (j=0; j<dim; j++){
        if (c[j] < min[j]) min[j] = c[j];
        if (c[j] > max[j]) max[j] = c[j];
      }
    }

    for (i=0; i<dim; i++){
      tmp = max[i];
      MPI_Allreduce(&tmp, max + i, 1, MPI_DOUBLE, MPI_MAX, local_comm);
      tmp = min[i];
      MPI_Allreduce(&tmp, min + i, 1, MPI_DOUBLE, MPI_MIN, local_comm);
    }

    for (j=0; j<dim; j++){
      /* distance along the j'th eigenvector */
      d[j] = max[order[j]] - min[order[j]];
    }
    
    return;
  }

  for (i=0; i<dim; i++){
    for (j=0, c=coords; j<num_obj; j++, c+=dim){

      val = (((c[0] - cm[0]) * evecs[0][i]) + ((c[1] - cm[1]) * evecs[1][i]));

      if (dim == 3){
        val += ((c[2] - cm[2]) * evecs[2][i]);
      }

      if (j){
        if (val < min[i]) min[i] = val;
        if (val > max[i]) max[i] = val;
      }
      else{
        min[i] = max[i] = val;
      }
    }
  }

  if (Tflops_Special){
    for (i=0; i<dim; i++){
      Zoltan_RIB_min_max(min+i, max+i, proclower, proc, nprocs,
                       local_comm);
      d[i] = max[i] - min[i];
    }
  }
  else {
    for (i=0; i<dim; i++){
      tmp = max[i];
      MPI_Allreduce(&tmp, max + i, 1, MPI_DOUBLE, MPI_MAX, local_comm);
      tmp = min[i];
      MPI_Allreduce(&tmp, min + i, 1, MPI_DOUBLE, MPI_MIN, local_comm);

      d[i] = max[i] - min[i];
    }
  }

  return;
}
/*
 * Order the 2 or 3 lengths from longest to shortest.
 */
static void order_decreasing(double *d, int *order)
{
  if (d[0] > d[1]){
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
/* 
 * Given a real symmetric matrix "m", find its eigenvectors.  Put
 * the orthonormal eigenvectors in the columns of the matrix "evecs".
 * Assume dim is 2 or 3.
 */

#define SIGN(a,b) ((b) < 0 ? -fabs(a) : fabs(a))
static int eigenvectors(double (*m)[3], double (*evecs)[3], int dim)
{
  double temp[2][2];
  double eval1, eval2, eval3;    /* eigenvalue and error in eval calculation */
  double res;
  double tmp;
  int i, j, rc = 0;

  for (i=0; i<3; i++){
    for (j=0; j<3; j++){
      evecs[i][j] = m[i][j];
    }
  }

  if (dim == 3) {
    Zoltan_evals3(m, &eval1, &eval2, &eval3);
    Zoltan_eigenvec3(m, eval1, evecs[0], &res);
    Zoltan_eigenvec3(m, eval2, evecs[1], &res);
    Zoltan_eigenvec3(m, eval3, evecs[2], &res);
  }
  else if (dim == 2) {
    for (i=0; i<2; i++){
      for (j=0; j<2; j++){
          temp[i][j] = m[i][j];
      }
    }
    Zoltan_evals2(temp, &eval1, &eval2);
    Zoltan_eigenvec2(temp, eval1, evecs[0], &res);
    Zoltan_eigenvec2(temp, eval2, evecs[1], &res);
  }

  for (i=0; i<3; i++){
    for (j=0; j<3; j++){
      if (i < j){
        tmp = evecs[i][j];
        evecs[i][j] = evecs[j][i];
        evecs[j][i] = tmp;
      }
    }
  }

  return rc;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
