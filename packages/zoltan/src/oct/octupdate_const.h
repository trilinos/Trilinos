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

#ifndef __OCTUPDATE_CONST_H
#define __OCTUPDATE_CONST_H

#ifndef lint
static char *cvs_octantupdateconsth_id = "$Id$";
#endif

#define MINOCTREGIONS 1              /* minimum number of regions per octant */

extern void    oct_gen_tree_from_input_data(LB *lb, int *c1, int *c2,
					    int *c3, float *c0);
#ifdef LGG_MIGOCT
extern void    oct_roots_in_order(pOctant **roots_ret,int *nroots_ret);
extern void    oct_resetIdCount(int start_count);
extern int     oct_nextId(void);
#endif /* LGG_MIGOCT */
extern void    oct_terminal_refine(pOctant oct,int count);
extern pOctant oct_findId(int i);
extern pOctant oct_global_insert(pRegion region);
extern int     oct_subtree_insert(pOctant oct, pRegion region);

#endif
