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
#ifndef __RCB_CONST_H
#define __RCB_CONST_H

#ifndef lint
static char *cvs_rcbconsth_id = "$Id$";
#endif

/* Data structures for parallel RCB

      rcb_dot has 4 required fields as shown below
      other fields can be added by user
      are just carried along as dots migrate to new processors

      examples:

  int       global;                global id # of dot
  int       local;                 local id # (memory loc) of dot before RCB
  int       proc;                  owner of this dot before RCB
*/

/* Steve Plimpton, Sandia National Labs, ABQ, NM  87185
   Dept 9221, MS 1111
   (505) 845-7873
   sjplimp@cs.sandia.gov
*/

                             /* dot to balance on for RCB */ 
struct rcb_dot {	        /* dot = point in 3-space */
  double    X[3];		/* location of dot */
  double    Weight;             /* weight of dot - if used must be > 0 */
  LB_TAG    Tag;                /* Tag containing IDs for the object.  */
};

struct rcb_tree {	     /* tree of RCB cuts */
  double    cut;        	/* position of cut */
  int       dim;	        /* dimension (012) of cut */
};

struct rcb_median {          /* RCB cut info */
  double    totallo, totalhi;   /* weight in each half of active partition */
  double    valuelo, valuehi;	/* position of dot(s) nearest to cut */
  double    wtlo, wthi;         /* total weight of dot(s) at that position */
  int       countlo, counthi;   /* # of dots at that position */
  int       proclo, prochi;	/* unique proc who owns a nearest dot */
};

struct rcb_box {       	     /* bounding box */
  double    lo[3], hi[3];	/* xyz lo/hi bounds */
};

typedef struct RCB_Struct {
  struct rcb_dot *Dots;
  struct rcb_tree *Tree_Ptr;
  struct rcb_box *Box;
  int Dot_Top;
} RCB_STRUCT;

extern void LB_rcb_build_data_structure(LB *, int *, int *, int);
extern int LB_RCB_Set_Param(char *, char *);

#endif
