#include "octant_const.h"
#include "octupdate_const.h"
#include "util_const.h"
#include "msg_const.h"
#include "migreg_const.h"
#include <mpi.h>

typedef struct
{
  pRegion region;
  int npid;
} Message;

void migreg_migrate_regions(Region *regions, int *npids, 
			    int nregions, int *c2);
void insert_orphan(Region reg);
void copy_info(pRegion src, pRegion *dest);
