/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */


#ifndef __RIB_H
#define __RIB_H

#include "zz_const.h"
#include "shared.h"
#include "rib_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* Data structures for parallel recursive inertial bisection method */

struct rib_tree {               /* tree of rib method cuts */
    double    cm[3];            /* center of mass */
    double    ev[3];            /* perpendicular direction from cut */
    double    cut;              /* position of cut */
    int       parent;           /* parent of this node in cut tree */
    int       left_leaf;        /* left child of this node in cut tree */
    int       right_leaf;       /* right child of this node in cut tree */
};

typedef struct RIB_Struct {
    ZOLTAN_ID_PTR Global_IDs;       /* This array is NOT used if Zoltan_RB_Use_IDs returns
                                   FALSE.   */
    ZOLTAN_ID_PTR Local_IDs;        /* This array is NOT used if Zoltan_RB_Use_IDs returns
                                   FALSE.   */
    struct Dot_Struct Dots;
    struct rib_tree   *Tree_Ptr;
    int                Num_Geom;
    ZZ_Transform       Tran;     /* transformation for degenerate geometry */
} RIB_STRUCT;

extern int Zoltan_RIB_Build_Structure(ZZ *, int *, int *, int, double, int,int);
extern void Zoltan_RIB_Print_Structure(ZZ *zz, int howMany);


/* function prototypes */

extern int Zoltan_RIB_inertial1d(struct Dot_Struct *, int *, int, int, double *, double *,
                         double *);
extern int Zoltan_RIB_inertial2d(int, struct Dot_Struct *, int *, int, int, double *,
                         double *, double *, MPI_Comm, int, int, int);
extern int Zoltan_RIB_inertial3d(int, struct Dot_Struct *, int *, int, int, double *,
                         double *, double *, MPI_Comm, int, int, int);
extern void Zoltan_RIB_reduce_double(double *, double *, int, MPI_Comm, int, int, int,
                             int);

extern void Zoltan_RIB_min_max(double *, double *, int, int, int, MPI_Comm);

extern void Zoltan_RIB_inertial3d_all( int, double *,
     int, double *cm, double (*evec)[3], 
     MPI_Comm , int , int , int );

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
