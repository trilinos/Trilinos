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

#ifndef __IRB_CONST_H
#define __IRB_CONST_H

#include "shared_const.h"

/* Data structures for parallel inertial recursive bisection method */

struct irb_tree {               /* tree of irb method cuts */
    double    cm[3];            /* center of mass */
    double    ev[3];            /* perpendicular direction from cut */
    double    cut;              /* position of cut */
    int       parent;           /* parent of this node in cut tree */
    int       left_leaf;        /* left child of this node in cut tree */
    int       right_leaf;       /* right child of this node in cut tree */
};

typedef struct IRB_Struct {
    LB_ID_PTR Global_IDs;
    LB_ID_PTR Local_IDs;
    struct Dot_Struct *Dots;
    struct irb_tree   *Tree_Ptr;
    int                Num_Geom;
} IRB_STRUCT;

extern int LB_IRB_Build_Structure(LB *, int *, int *, int);
extern void LB_IRB_Free_Structure(LB *);
extern int LB_Set_IRB_Param(char *, char *);

/* function prototypes */

extern int LB_inertial1d(struct Dot_Struct *, int, int, double *, double *,
                         double *);
extern int LB_inertial2d(struct Dot_Struct *, int, int, double *, double *,
                         double *, MPI_Comm);
extern int LB_inertial3d(struct Dot_Struct *, int, int, double *, double *,
                         double *, MPI_Comm);

#endif
