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


/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */
/* It was modified by Courtenay T. Vaughan for use in Zoltan */

#include <stdio.h>
#include <math.h>
#include "rib.h"
#include "inertial.h"
#include "zz_const.h"

/* macros for routines */
#define max(a, b) ((a) < (b) ? (b) : (a))
#define min(a, b) ((a) > (b) ? (b) : (a))
#define sign(x)   ((x) >= 0 ? 1.0 : -1.0)


int Zoltan_RIB_inertial2d(
     int Tflops_Special,        /* Use Tflops_Special communication;
                                   should be 0 if called from serial_rib */
     struct Dot_Struct *dotpt,  /* graph data structure for weights */
     int             *dindx,    /* index array into dotpt; if NULL, access dotpt
                                   directly */
     int              dotnum,   /* number of vtxs in graph */
     int              wgtflag,  /* are vertex weights being used? */
     double           cm[3],    /* center of mass in each direction */
     double           evec[3],  /* eigenvector */
     double           *value,   /* array for value to sort on */
     MPI_Comm         comm,     /* communicator for partition */
     int proc,          /* global proc number (Tflops_Special) */
     int nproc,         /* Number of procs in partition (Tflops_Special) */
     int proclower      /* Lowest numbered proc in partition (Tflops_Special)*/
)
{
     double tensor[2][2];       /* inertial tensor */
     double xcm, ycm;           /* center of mass in each direction */
     double tmp1[3], tmp2[3];   /* temporary variables for MPI_Allreduces */
     double xx, yy, xy;         /* elements of inertial tensor */
     double xdif, ydif;         /* deviation from center of mass */
     double eval, res;          /* eigenvalue and error in eval calculation */
     double wgt_sum;            /* sum of all the vertex weights */
     int    i, j;               /* loop counter */
     double xcmt, ycmt, wgtt;   /* temp for center of mass */
     double xxt, yyt, xyt;      /* temp for tensor */
     int    rank = 0;           /* rank in partition (Tflops_Special) */
     double wgt;
     double *x, *y;

     /* Compute center of mass and total mass. */

     x = dotpt->X;
     y = dotpt->Y;
     xcm = ycm = 0.0;

     if (dotpt->nWeights) {
        wgt_sum = 0.0;
        for (j = 0; j < dotnum; j++) {
           i = (dindx ? dindx[j] : j);
           wgt = dotpt->Weight[i * dotpt->nWeights];
           wgt_sum += wgt;
           xcm += wgt*x[i];
           ycm += wgt*y[i];
        }
     }
     else {
        wgt_sum = dotnum;
        for (j = 0; j < dotnum; j++) {
           i = (dindx ? dindx[j] : j);
           xcm += x[i];
           ycm += y[i];
        }
     }

     /* Sum weights across processors */

     if (Tflops_Special) {
        rank = proc - proclower;
        tmp1[0] = xcm; tmp1[1] = ycm; tmp1[2] = wgt_sum;
        Zoltan_RIB_reduce_double(tmp1, tmp2, 3, comm, nproc, rank, proc, 1);
        xcmt = tmp2[0]; ycmt = tmp2[1]; wgtt = tmp2[2];
     }   
     else {
        tmp1[0] = xcm; tmp1[1] = ycm; tmp1[2] = wgt_sum;
        MPI_Allreduce(tmp1, tmp2, 3, MPI_DOUBLE, MPI_SUM, comm);
        xcmt = tmp2[0]; ycmt = tmp2[1]; wgtt = tmp2[2];
     }

     xcm = xcmt/wgtt;
     ycm = ycmt/wgtt;

     /* Generate 3 elements of Inertial tensor. */
     xx = yy = xy = 0.0;
     if (dotpt->nWeights)
        for (j = 0; j < dotnum; j++) {
           i = (dindx ? dindx[j] : j);
           wgt = dotpt->Weight[i * dotpt->nWeights];
           xdif = x[i] - xcm;
           ydif = y[i] - ycm;
           xx += wgt*xdif*xdif;
           yy += wgt*ydif*ydif;
           xy += wgt*xdif*ydif;
        }
     else
        for (j = 0; j < dotnum; j++) {
           i = (dindx ? dindx[j] : j);
           xdif = x[i] - xcm;
           ydif = y[i] - ycm;
           xx += xdif*xdif;
           yy += ydif*ydif;
           xy += xdif*ydif;
        }

     /* Sum tensor across processors */

     if (Tflops_Special) {
        tmp1[0] = xx; tmp1[1] = yy; tmp1[2] = xy;
        Zoltan_RIB_reduce_double(tmp1, tmp2, 3, comm, nproc, rank, proc, 1);
        xxt = tmp2[0]; yyt = tmp2[1]; xyt = tmp2[2];
     }
     else {
        tmp1[0] = xx; tmp1[1] = yy; tmp1[2] = xy;
        MPI_Allreduce(tmp1, tmp2, 3, MPI_DOUBLE, MPI_SUM, comm);
        xxt = tmp2[0]; yyt = tmp2[1]; xyt = tmp2[2];
     }

     /* Compute eigenvector with maximum eigenvalue. */

     tensor[0][0] = xxt;
     tensor[1][1] = yyt;
     tensor[1][0] = tensor[0][1] = xyt;
     Zoltan_evals2(tensor, &res, &eval);
     Zoltan_eigenvec2(tensor, eval, evec, &res);

     /* Calculate value to sort/split on for each cell. */
     /* This is inner product with eigenvector. */
     for (j = 0; j < dotnum; j++) {
        i = (dindx ? dindx[j] : j);
        value[j] = (x[i] - xcm)*evec[0] +
                   (y[i] - ycm)*evec[1];
     }

     cm[0] = xcm;
     cm[1] = ycm;

     /* zero unused third dimension */
     cm[2] = evec[2] = 0.0;

     return(ZOLTAN_OK);
}


/* Find eigenvalues of 2x2 symmetric system by solving quadratic. */
void Zoltan_evals2(
     double H[2][2],            /* symmetric matrix for eigenvalues */
     double *eval1,             /* smallest eigenvalue */
     double *eval2              /* middle eigenvalue */
)
{
     double M[2][2] = {{0.,0.},
                       {0.,0.}};/* normalized version of matrix */
     double b, c;               /* coefficents of cubic equation */
     double root1, root2;       /* roots of quadratic */
     double xmax;               /* largest matrix element */
     int    i, j;               /* loop counters */

     xmax = 0.0;
     for (i = 0; i < 2; i++)
        for (j = i; j < 2; j++)
           if (fabs(H[i][j]) > xmax)
              xmax = fabs(H[i][j]);
     if (xmax != 0)
        for (i = 0; i < 2; i++)
           for (j = 0; j < 2; j++)
              M[i][j] = H[i][j] / xmax;

     b = -M[0][0] - M[1][1];
     c = M[0][0] * M[1][1] - M[1][0] * M[1][0];
     root1 = -.5 * (b + sign(b) * sqrt(max(0.0, b * b - 4 * c)));
     root2 = c / root1;

     root1 *= xmax;
     root2 *= xmax;
     *eval1 = min(root1, root2);
     *eval2 = max(root1, root2);
}


/* Solve for eigenvector of SPD 2x2 matrix, with given eigenvalue. */
void Zoltan_eigenvec2(
     double A[2][2],            /* matrix */
     double eval,               /* eigenvalue */
     double evec[2],            /* eigenvector returned */
     double *res                /* normalized residual */
)
{
     double norm;               /* norm of eigenvector */
     double res1, res2;         /* components of residual vector */
     int    i;                  /* loop counter */

     if (fabs(A[0][0] - eval) > fabs(A[1][1] - eval)) {
        evec[0] = -A[1][0];
        evec[1] = A[0][0] - eval;
     }
     else {
        evec[0] = A[1][1] - eval;
        evec[1] = -A[1][0];
     }

     /* Normalize eigenvector and calculate a normalized eigen-residual. */
     norm = sqrt(evec[0] * evec[0] + evec[1] * evec[1]);
     if (norm == 0) {
        evec[0] = 1;
        evec[1] = 0;
        norm = 1;
     }
     for (i = 0; i < 2; i++)
        evec[i] /= norm;
     res1 = (A[0][0] - eval) * evec[0] + A[1][0] * evec[1];
     res2 = A[1][0] * evec[0] + (A[1][1] - eval) * evec[1];
     *res = sqrt(res1 * res1 + res2 * res2);

     res1 = fabs(A[0][0]) + fabs(A[1][0]);
     res2 = fabs(A[1][1]) + fabs(A[1][0]);
     *res /= max(res1, res2);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
