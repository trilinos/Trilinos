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

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */
/* It was modified by Courtenay T. Vaughan for use in Zoltan */

#include <stdio.h>
#include <math.h>
#include "zz_const.h"
#include "rib.h"

/* macros for routines */ 
#define max(a, b) ((a) < (b) ? (b) : (a)) 
#define min(a, b) ((a) > (b) ? (b) : (a)) 
#define sign(x)   ((x) >= 0 ? 1.0 : -1.0)

/* function prototypes */

static void evals3(double[3][3], double *, double *, double *);
static double determinant(double[3][3]);
static void eigenvec3(double[3][3], double, double *, double *);

int Zoltan_RIB_inertial3d(
     ZZ               *Zz,      /* Zoltan structure Tflops_Special */
     struct Dot_Struct *dotpt,  /* graph data structure */
     int              dotnum,   /* number of vtxs in graph */
     int              wgtflag,  /* are vertex weights being used? */
     double           cm[3],    /* center of mass in each direction */
     double           evec[3],  /* eigenvector */
     double           *value,   /* array for value to sort on */
     MPI_Comm         comm,     /* communicator for partition */
     int proc,          /* Global proc number (Tflops_Special) */
     int nproc,         /* Number of procs in partition (Tflops_Special) */
     int proclower      /* Lowest numbered proc in partition (Tflops_Special)*/
)
{
     double    tensor[3][3];    /* inertia tensor */
     double    xx, yy, zz;      /* elements of inertial tensor */
     double    xy, xz, yz;      /* elements of inertial tensor */
     double    xdif, ydif, zdif;/* deviation from center of mass */
     double    eval, res;       /* eigenvalue and error in eval calculation */
     double    wgt_sum;         /* sum of all the vertex weights */
     int       i;               /* loop counter */
     double    cmt[3], wgtt;    /* temp for center of mass */
     double    xxt, yyt, zzt;   /* temp for inertial tensor */
     double    xyt, xzt, yzt;   /* temp for inertial tensor */
     int       rank;            /* rank in partition (Tflops_Special) */

     /* Compute center of mass and total mass. */
     cm[0] = cm[1] = cm[2] = 0.0;
     if (wgtflag) {
        wgt_sum = 0;
        for (i = 0; i < dotnum; i++) {
           wgt_sum += dotpt[i].Weight;
           cm[0] += dotpt[i].Weight*dotpt[i].X[0];
           cm[1] += dotpt[i].Weight*dotpt[i].X[1];
           cm[2] += dotpt[i].Weight*dotpt[i].X[2];
        }
     }
     else {
        wgt_sum = dotnum;
        for (i = 0; i < dotnum; i++) {
           cm[0] += dotpt[i].X[0];
           cm[1] += dotpt[i].X[1];
           cm[2] += dotpt[i].X[2];
        }
     }

     /* Sum weights across processors */

     if (Zz->Tflops_Special) {
        rank = proc - proclower;
        Zoltan_RIB_reduce_double(cm, cmt, 3, comm, nproc, rank, proc, 1);
        Zoltan_RIB_reduce_double(&wgt_sum, &wgtt, 1, comm, nproc, rank, proc, 1);
     }
     else {
        MPI_Allreduce(cm,cmt,3,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&wgt_sum,&wgtt,1,MPI_DOUBLE,MPI_SUM,comm);
     }

     cm[0] = cmt[0]/wgtt;
     cm[1] = cmt[1]/wgtt;
     cm[2] = cmt[2]/wgtt;

     /* Generate 6 elements of Inertial tensor. */
     xx = yy = zz = xy = xz = yz = 0.0;
     if (wgtflag)
        for (i = 0; i < dotnum; i++) {
           xdif = dotpt[i].X[0] - cm[0];
           ydif = dotpt[i].X[1] - cm[1];
           zdif = dotpt[i].X[2] - cm[2];
           xx += dotpt[i].Weight*xdif*xdif;
           yy += dotpt[i].Weight*ydif*ydif;
           zz += dotpt[i].Weight*zdif*zdif;
           xy += dotpt[i].Weight*xdif*ydif;
           xz += dotpt[i].Weight*xdif*zdif;
           yz += dotpt[i].Weight*ydif*zdif;
        }
     else
        for (i = 0; i < dotnum; i++) {
           xdif = dotpt[i].X[0] - cm[0];
           ydif = dotpt[i].X[1] - cm[1];
           zdif = dotpt[i].X[2] - cm[2];
           xx += xdif*xdif;
           yy += ydif*ydif;
           zz += zdif*zdif;
           xy += xdif*ydif;
           xz += xdif*zdif;
           yz += ydif*zdif;
        }

     /* Sum tensor across processors */

     if (Zz->Tflops_Special) {
        Zoltan_RIB_reduce_double(&xx, &xxt, 1, comm, nproc, rank, proc, 1);
        Zoltan_RIB_reduce_double(&yy, &yyt, 1, comm, nproc, rank, proc, 1);
        Zoltan_RIB_reduce_double(&zz, &zzt, 1, comm, nproc, rank, proc, 1);
        Zoltan_RIB_reduce_double(&xy, &xyt, 1, comm, nproc, rank, proc, 1);
        Zoltan_RIB_reduce_double(&xz, &xzt, 1, comm, nproc, rank, proc, 1);
        Zoltan_RIB_reduce_double(&yz, &yzt, 1, comm, nproc, rank, proc, 1);
     }
     else {
        MPI_Allreduce(&xx,&xxt,1,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&yy,&yyt,1,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&zz,&zzt,1,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&xy,&xyt,1,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&xz,&xzt,1,MPI_DOUBLE,MPI_SUM,comm);
        MPI_Allreduce(&yz,&yzt,1,MPI_DOUBLE,MPI_SUM,comm);
     }

     /* Compute eigenvector with maximum eigenvalue. */

     tensor[0][0] = xxt;
     tensor[1][1] = yyt;
     tensor[2][2] = zzt;
     tensor[0][1] = tensor[1][0] = xyt;
     tensor[0][2] = tensor[2][0] = xzt;
     tensor[1][2] = tensor[2][1] = yzt;
     evals3(tensor, &res, &res, &eval);
     eigenvec3(tensor, eval, evec, &res);

     /* Calculate value to sort/split on for each cell. */
     /* This is inner product with eigenvector. */
     for (i = 0; i < dotnum; i++)
        value[i] = (dotpt[i].X[0] - cm[0])*evec[0] +
                   (dotpt[i].X[1] - cm[1])*evec[1] +
                   (dotpt[i].X[2] - cm[2])*evec[2];

     return(ZOLTAN_OK);
}

/* Find eigenvalues of 3x3 symmetric system by solving cubic. */
static void evals3(
     double H[3][3],            /* 3x3 sym matrix for lowest eigenvalue */
     double *eval1,             /* smallest eigenvalue */
     double *eval2,             /* middle eigenvalue */
     double *eval3              /* largest eigenvalue */
)
{
     double mat[3][3];          /* scaled version of H */
     double a1, a2, a3;         /* coefficents of cubic equation */
     double q, r;               /* intermediate terms */
     double q3, r2;             /* powers of q and r */
     double theta;              /* angular parameter */
     double root1, root2, root3;/* three roots of cubic */
     double tol = 1.0e-6;       /* allowed deviation */
     double xmax;               /* largest matrix element for scaling */
     int    i, j;               /* loop indices */
     double HALFPI = 1.570796327;
     double TWOPI = 6.283185307;

     /* This first requires solving a cubic equation. */
     /* Normalize to avoid any numerical problems. */
     xmax = 0.0;
     for (i = 0; i < 3; i++)
     for (j = i; j < 3; j++)
        if (fabs(H[i][j]) > xmax)
           xmax = fabs(H[i][j]);

     if (xmax != 0)
        for (i = 0; i < 3; i++)
           for (j = 0; j < 3; j++)
              mat[i][j] = H[i][j] / xmax;

     a1 = -(mat[0][0] + mat[1][1] + mat[2][2]);
     a2 = (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]) +
          (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) +
          (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]);
     a3 = -determinant(mat);

     if (a3 == 0) {
        root1 = 0;              /* Solve quadratic. */

        q = -.5 * (a1 + sign(a1) * sqrt(max(0.0, a1 * a1 - 4 * a2)));
        root2 = q;
        root3 = a2 / q;
     }
     else {                     /* solve cubic */
        q = (a1 * a1 - 3 * a2) / 9;
        r = (2 * a1 * a1 * a1 - 9 * a1 * a2 + 27 * a3) / 54;
        q3 = q * q * q;
        r2 = r * r;

        /* To avoid missing a root, check for roundoff. */
        if ((q3 < r2) && fabs(q3 - r2) < tol * (fabs(q3) + fabs(r2))) {
           q3 = r2;
        }

        if (q3 >= r2) {         /* Three real roots. */
           if (r == 0)
           theta = HALFPI;
           else {
           q3 = sqrt(q3);
           if (q3 < fabs(r))
              q3 = fabs(r);
           theta = acos(r / q3);
           }
           q = -2 * sqrt(q);

           root1 = q * cos(theta / 3) - a1 / 3;
           root2 = q * cos((theta + TWOPI) / 3) - a1 / 3;
           root3 = q * cos((theta + 2 * TWOPI) / 3) - a1 / 3;
        }
        else {                  /* Only one real root. */
           theta = sqrt(r2 - q3) + fabs(r);
           theta = pow(theta, 1.0 / 3.0);

           root1 = root2 = root3 = -sign(r) * (theta + q / theta) - a1 / 3;
        }
     }
     root1 *= xmax;
     root2 *= xmax;
     root3 *= xmax;
     *eval1 = min(root1, root2);
     *eval1 = min(*eval1, root3);
     *eval3 = max(root1, root2);
     *eval3 = max(*eval3, root3);
     if (root1 != *eval1 && root1 != *eval3)
        *eval2 = root1;
     else if (root2 != *eval1 && root2 != *eval3)
        *eval2 = root2;
     else
        *eval2 = root3;
}


static double determinant(
     double A[3][3]             /* matrix A */
)
{
     double det;

     det = A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] +
           A[0][2]*A[1][0]*A[2][1] - A[0][1]*A[1][0]*A[2][2] -
           A[0][0]*A[1][2]*A[2][1] - A[0][2]*A[1][1]*A[2][0];

     return det;
}


/* Find the eigenvector of symmetric 3x3 matrix w/ given eigenvalue. */
static void eigenvec3(
     double A[3][3],            /* matrix to find eigenvector of */
     double eval,               /* eigenvalue */
     double evec[3],            /* eigenvector returned */
     double *res                /* returned error estimate */
)
{
     double mat[3][3];          /* copy of A to write over */
     int    ind[3];             /* permutation indices */
     double ex, ey, ez;         /* elements of eigenvector returned */
     double xmax;               /* maximum value in matrix */
     double tmp;                /* intermediate values */
     double norm;               /* norm of eigenvector */
     double res1, res2, res3;   /* elements of residual vector */
     double tol = 1.0e-6;       /* smaller value assumed to be zero */
     int    imax, jmax;         /* indices of max value in matrix */
     int    i, j;               /* loop counters */

     for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
           mat[i][j] = A[i][j];
     for (i = 0; i < 3; i++)
        mat[i][i] -= eval;

     ind[0] = 0;
     ind[1] = 1;
     ind[2] = 2;

     /* Find the largest element in the matrix. */
     xmax = 0.0;
     for (i = 0; i < 3; i++)
        for (j = i; j < 3; j++)
           if (fabs(mat[i][j]) > xmax) {
              imax = i;
              jmax = j;
              xmax = fabs(mat[i][j]);
           }

     if (xmax == 0.0) {         /* Handle completely degenerate case first. */
        evec[0] = 1.0;
        evec[1] = evec[2] = 0.0;
     }
     else {
        /* Scale the matrix so largest value is 1.0 */
        for (i = 0; i < 3; i++)
           for (j = 0; j < 3; j++)
              mat[i][j] /= xmax;

        /* Swap rows if necessary to move max element to first row. */
        if (imax != 0)
           for (j = 0; j < 3; j++) {
              tmp = mat[0][j];
              mat[0][j] = mat[imax][j];
              mat[imax][j] = tmp;
           }

        /* Swap columns if necessary to move max element to first position. */
        if (jmax != 0) {
           for (i = 0; i < 3; i++) {
              tmp = mat[i][0];
              mat[i][0] = mat[i][jmax];
              mat[i][jmax] = tmp;
           }
           ind[0] = jmax;
           ind[jmax] = 0;
        }

        /* Reduce matrix to 2x2 by subtracting off first row. */
        for (i = 1; i < 3; i++)
           for (j = 1; j < 3; j++)
              mat[i][j] = mat[0][0] * mat[i][j] - mat[i][0] * mat[0][j];

        /* Find maximum element in reduced 2x2 matrix. */
        xmax = 0.0;
        for (i = 1; i < 3; i++)
           for (j = i; j < 3; j++)
              if (fabs(mat[i][j]) > xmax) {
                 imax = i;
                 jmax = j;
                 xmax = fabs(mat[i][j]);
              }

        if (xmax < tol) {    /* Handle 2-fold degenerate case - skip to end. */
           ey = 1.0;
           ex = ez = 0;
        }
        else {
           /* Swap rows 2 and 3 to move max element to 2nd row. */
           if (imax != 1)
              for (j = 0; j < 3; j++) {
                 tmp = mat[1][j];
                 mat[1][j] = mat[imax][j];
                 mat[imax][j] = tmp;
              }

           /* Swap columns to move max element to (1,1) position. */
           if (jmax != 1) {
              for (i = 0; i < 3; i++) {
                 tmp = mat[i][1];
                 mat[i][1] = mat[i][2];
                 mat[i][2] = tmp;
              }
              i = ind[1];
              ind[1] = ind[2];
              ind[2] = i;
           }

           /* Compute eigenvector from numerically stabilized matrix. */
           ez = mat[0][0] * mat[1][1];
           ey = -mat[1][2] * mat[0][0];
           ex = mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1];
        }
        /* Reorder the e-vector to undo pivoting - end of 3-D case. */
        evec[ind[0]] = ex;
        evec[ind[1]] = ey;
        evec[ind[2]] = ez;
     }

     /* Normalize eigenvector and calculate a normalized eigen-residual. */
     norm = sqrt(evec[0] * evec[0] + evec[1] * evec[1] + evec[2] * evec[2]);
     for (i = 0; i < 3; i++)
        evec[i] /= norm;
     res1 = (A[0][0] - eval) * evec[0] + A[0][1] * evec[1] + A[0][2] * evec[2];
     res2 = A[1][0] * evec[0] + (A[1][1] - eval) * evec[1] + A[1][2] * evec[2];
     res3 = A[2][0] * evec[0] + A[2][1] * evec[1] + (A[2][2] - eval) * evec[2];
     *res = sqrt(res1 * res1 + res2 * res2 + res3 * res3);

     /* Now normalize the residual */
     res1 = fabs(A[0][0]) + fabs(A[0][1]) + fabs(A[0][2]);
     res2 = fabs(A[1][0]) + fabs(A[1][1]) + fabs(A[1][2]);
     res3 = fabs(A[2][0]) + fabs(A[2][1]) + fabs(A[2][2]);
     res2 = max(res2, res3);
     *res /= max(res1, res2);
}

void Zoltan_RIB_reduce_double(double *in, double *out, int len, MPI_Comm comm,
                      int nproc, int rank, int proc, int n)
{
   int i, m, to, tag = 32107;
   MPI_Status status;
 
   /* this is a recursive function for Tflops_Special that takes a double
      vector in and returns the summed vector out for a subset of processors
      of the entire number of processors.  rank is a processors rank within
      its partition of size nproc while proc is the rank of the processor with
      the entire number of processors being used. */

   m = 2*n;
   if (rank%m) {
      to = proc - n;
      MPI_Send(in, len, MPI_DOUBLE, to, tag, comm);
      MPI_Recv(out, len, MPI_DOUBLE, to, tag, comm, &status);
   }
   else
      if (rank + n < nproc) {
         to = proc + n;
         MPI_Recv(out, len, MPI_DOUBLE, to, tag, comm, &status);
         for (i = 0; i < len; i++)
            in[i] += out[i];
         if (m < nproc)
            Zoltan_RIB_reduce_double(in, out, len, comm, nproc, rank, proc, m);
         else
            for (i = 0; i < len; i++) 
               out[i] = in[i];
         MPI_Send(out, len, MPI_DOUBLE, to, tag, comm);
      }  
      else
         Zoltan_RIB_reduce_double(in, out, len, comm, nproc, rank, proc, m);
}
