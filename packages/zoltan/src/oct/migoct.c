/*** ATTN: name of functions are relics of previous implementations of octree
 *** because of this, they need to be renamed at a later date.....
 *** Thank you for your support.
 ***/

#include <unistd.h>
#include <stdio.h>
#include "msg_const.h"
#include "octant_const.h"
#include "migoct.h"

/* 
 * pOctant *octs;         octs[nocts] 
 * int *newpids;          newpids[nocts]
 * int nocts;             gives size of octs[] and newpids[]
 * LB_TAG **export_tags;  return array of Tags to migrate
 * int *nsentags;         number of Tags to migrate (size of export_tags array)
 * LB_TAG **import_tags;  return array of Tags migrated in
 * int *nrectags;         number of Tags migrated
 */
void Migrate_Objects(pOctant *octs, int *newpids, int nocts,
		     LB_TAG **export_tags, int *nsentags, 
		     LB_TAG **import_tags, int *nrectags)
{
  int i;
  int *tag_pids;

  tag_regions(octs, newpids, nocts, export_tags, nsentags, &tag_pids);
  malloc_new_octants(*nsentags, *export_tags, tag_pids, nrectags, import_tags);

}

/*
 * void tag_regions()
 * Iterates through the list of octants on the processor and finds which
 * are to be migrated. It then looks at the region list for those octants 
 * and stores the migrating regions into the export_tags array.
 */
static void tag_regions(pOctant *octs, int *newpids, int nocts, 
			LB_TAG **export_tags, int *nsentags, int **tag_pids)
{
  int i;                /* index counter */
  pRegion regionlist;   /* list of region on this processor */
  pRegion region;       /* pointer to a region */
  int index;            /* index counter */
  int count;            /* count of number of region on this processor */
  int *export_pids;     /* array of pids where regions are being exported to */
  LB_TAG *mtags;        /* object tags of objects to be migrated */

  /* check for migrating, pointer should not be larger than an int */
  if (sizeof(int)<sizeof(pOctant)) {
    fprintf(stderr,"tag_regions: Fatal error, sizeof(int)<sizeof(ptr)\n");
    abort();
  }

  /* if there are no tags to send, then finished here */
  if (!nsentags)
    return;

  /* find how many objects are being sent */
  count=0;
  for (i=0; i<nocts; i++) {
    /* if newpids != msg_mypid, the it is being migrated */
    if (newpids[i]!=msg_mypid) {
      if(POC_isTerminal(octs[i]))
	count+=POC_nRegions(octs[i]);
      else
	fprintf(stderr, "Warning: %s\n",
		"tried to access nRegions from non-terminal.");
    }
  }

  /* fprintf(stderr, "%d count = %d\n", msg_mypid, count); */
  *nsentags=count;

  /* export_tags should be null at this point */
  if (!export_tags)
    return;

  if (count > 0) {
    mtags = (LB_TAG *)malloc((unsigned)count * sizeof(LB_TAG));
    export_pids = (int *)malloc((unsigned)count * sizeof(int));
    if(export_pids == NULL) {
      fprintf(stderr, "ERROR: unable to malloc export_pids in tag_regions\n");
      abort();
    }
    if(mtags == NULL) {
      fprintf(stderr, "ERROR: unable to malloc mtags in tag_regions\n");
      abort();
    }
  }
  else {
    mtags = NULL;
    export_pids = NULL;
  }
  /* set up return pointers */
  *export_tags=mtags;
  *tag_pids = export_pids;

  index=0;
  for (i=0; i<nocts; i++) {
    if (newpids[i]!=msg_mypid && POC_isTerminal(octs[i])) {
      regionlist=POC_regionlist(octs[i]);
      while(regionlist != NULL) {
	mtags[index]=regionlist->Tag;
	export_pids[index] = newpids[i];
	index++;
	regionlist = regionlist->next;
      }
    }
  }

  if (index!=count)
    abort();
}

static void malloc_new_octants(int nsentags, LB_TAG *export_tags, 
			       int *tag_pids, int *nrectags,
			       LB_TAG **import_tags)
{
  int i;
  int num;
  int from;
  int nsends;
  int nreceives;
  LB_TAG *imp;

  msg_send_init();

  nreceives = 0;
  imp = NULL;
  nsends=0;
  fprintf(stderr, "nsent = %d\n", nsentags);
  for (i=0; i<nsentags; i++) {                 /* Send msg requesting malloc */
    msg_bsend(&export_tags[i], sizeof(LB_TAG), tag_pids[i], MTYPE_DMIGRATE);
    fprintf(stderr, "to whom %d\n", tag_pids[i]);
    nsends++;
  }
  
  nreceives=msg_nreceives();
  *nrectags = nreceives;
  if(nreceives != 0)
    imp = (LB_TAG *)malloc((*nrectags) * sizeof(LB_TAG));
  *import_tags = imp;

  for (i=0; i<nreceives; i++) {                /* Reply to malloc requests */
    msg_breceive(&(imp[i]), sizeof(LB_TAG), &from, MTYPE_DMIGRATE);
  }
}



