/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __MIGOCT_CONST_H
#define __MIGOCT_CONST_H

extern void LB_migrate_regions(LB *lb, pOctant *octs, int *newpids, int nocts,
			       int *nsentags, 
			       pRegion *import_tags, int *nrectags, 
			       float *c2, float *c3, int *counter3, 
			       int *counter4);

extern void LB_fix_tags(LB *lb, LB_ID_PTR *import_global_ids, 
                        LB_ID_PTR *import_local_ids,
                        int **import_procs, int nrectags, pRegion import_regs);

#endif
