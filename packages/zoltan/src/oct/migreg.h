#include "octant_const.h"
#include "migreg_const.h"
#include "msg_const.h"
#include "octupdate_const.h"
#include "util_const.h"

typedef struct
{
  int octid;
  double min[3];
  double max[3];
} ROOTMSG;

typedef struct
{
  int pid;
  int octid;
  double min[3];
  double max[3];
  double size;
} ROOT;

typedef struct
{
  pRegion region;
  int npid;
} Message;

#define MTYPE_ROOT 100

void migreg_migrate_regions(Message *Array, int nregions);
int copy_info(pRegion *dest, pRegion src);
