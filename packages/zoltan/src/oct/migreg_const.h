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
#ifndef __MIGREG_CONST_H
#define __MIGREG_CONST_H

#ifndef lint
static char *cvs_migregconsth_id = "$Id$";
#endif

extern void migreg_migrate_orphans(pRegion RegionList, int nreg, 
				   int level, Map *array, int *c1, int *c2);

#endif
