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
#ifndef __MIGOCT_CONST_H
#define __MIGOCT_CONST_H

#ifndef lint
static char *cvs_migoctconsth_id = "$Id$";
#endif

extern void Migrate_Objects(pOctant *octs, int *newpids, int nocts,
			    pRegion *export_tags, int *nsentags, 
			    pRegion *import_tags, int *nrectags, 
			    float *c2, float *c3, int *counter3, 
			    int *counter4);

extern void fix_tags(LB_GID **import_global_ids, LB_LID **import_local_ids,
              int **import_procs, int nrectags, pRegion import_regs);

#endif
