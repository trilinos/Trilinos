/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __RIB_CONST_H
#define __RIB_CONST_H

#include "shared_const.h"

/* Data structures for parallel inertial recursive bisection method */

struct rib_tree {               /* tree of rib method cuts */
    double    cm[3];            /* center of mass */
    double    ev[3];            /* perpendicular direction from cut */
    double    cut;              /* position of cut */
    int       parent;           /* parent of this node in cut tree */
    int       left_leaf;        /* left child of this node in cut tree */
    int       right_leaf;       /* right child of this node in cut tree */
};

typedef struct RIB_Struct {
    LB_ID_PTR Global_IDs;
    LB_ID_PTR Local_IDs;
    struct Dot_Struct *Dots;
    struct rib_tree   *Tree_Ptr;
    int                Num_Geom;
} RIB_STRUCT;

extern int LB_RIB_Build_Structure(LB *, int *, int *, int);
extern void LB_RIB_Free_Structure(LB *);
extern int LB_Set_RIB_Param(char *, char *);

/* function prototypes */

extern int LB_inertial1d(struct Dot_Struct *, int, int, double *, double *,
                         double *);
extern int LB_inertial2d(LB *, struct Dot_Struct *, int, int, double *,
                         double *, double *, MPI_Comm, int, int, int);
extern int LB_inertial3d(LB *, struct Dot_Struct *, int, int, double *,
                         double *, double *, MPI_Comm, int, int, int);
extern void LB_reduce_double(double *, double *, int, MPI_Comm, int, int, int,
                             int);

#endif
