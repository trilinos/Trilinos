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

#ifndef __RIB_H
#define __RIB_H

#include "shared.h"
#include "rib_const.h"

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
    struct Dot_Struct *Dots;
    struct rib_tree   *Tree_Ptr;
    int                Num_Geom;
} RIB_STRUCT;

extern int Zoltan_RIB_Build_Structure(ZZ *, int *, int *, int, int);

/* function prototypes */

extern int Zoltan_RIB_inertial1d(struct Dot_Struct *, int, int, double *, double *,
                         double *);
extern int Zoltan_RIB_inertial2d(ZZ *, struct Dot_Struct *, int, int, double *,
                         double *, double *, MPI_Comm, int, int, int);
extern int Zoltan_RIB_inertial3d(ZZ *, struct Dot_Struct *, int, int, double *,
                         double *, double *, MPI_Comm, int, int, int);
extern void Zoltan_RIB_reduce_double(double *, double *, int, MPI_Comm, int, int, int,
                             int);

#endif
