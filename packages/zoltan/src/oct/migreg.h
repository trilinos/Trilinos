/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __MIGREG_H
#define __MIGREG_H

#include "octant_const.h"
#include "octupdate_const.h"
#include "oct_util_const.h"
#include "migreg_const.h"

typedef struct
{
  pRegion region;
  int npid;
} Message;

void LB_migreg_migrate_regions(LB *lb, Region *regions, int *npids, 
			       int nregions, int *c2);
void LB_insert_orphan(LB *, Region reg);
void LB_copy_info(pRegion src, pRegion *dest);

#endif
