/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/
#ifndef __IRB_CONST_H
#define __IRB_CONST_H

/* Data structures for parallel inertial recursive bisection method
*/

                                /* dot to balance on for inertial method */
struct irb_dot {                /* dot = point in 3-space */
    double    X[3];             /* location of dot */
    double    Weight;           /* weight of dot - if used must be > 0 */
    LB_TAG    Tag;              /* Tag containing IDs for the object.  */
};

struct irb_tree {               /* tree of irb method cuts */
    double    cm[3];            /* center of mass */
    double    ev[3];            /* perpendicular direction from cut */
    double    cut;              /* position of cut */
    int       parent;           /* parent of this node in cut tree */
    int       left_leaf;        /* left child of this node in cut tree */
    int       right_leaf;       /* right child of this node in cut tree */
};

typedef struct IRB_Struct {
    struct irb_dot  *Dots;
    struct irb_tree *Tree_Ptr;
    int               Num_Geom;
} IRB_STRUCT;

extern int LB_IRB_Build_Structure(LB *, int *, int *, int);
extern void LB_IRB_Free_Structure(LB *);
extern int LB_Set_IRB_Param(char *, char *);

#endif
