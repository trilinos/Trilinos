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
#ifndef __MIGREG_H
#define __MIGREG_H

#ifndef lint
static char *cvs_migregh_id = "$Id$";
#endif

#include "octant_const.h"
#include "octupdate_const.h"
#include "util_const.h"
#include "migreg_const.h"

typedef struct
{
  pRegion region;
  int npid;
} Message;

void migreg_migrate_regions(Region *regions, int *npids, 
			    int nregions, int *c2);
void insert_orphan(Region reg);
void copy_info(pRegion src, pRegion *dest);

#endif
