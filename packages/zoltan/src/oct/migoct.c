/*** ATTN: name of functions are relics of previous implementations of octree
 *** because of this, they need to be renamed at a later date.....
 *** Thank you for your support.
 ***/

#include <unistd.h>
#include <stdio.h>
#include "msg_const.h"
#include "octant_const.h"
#include "migoct.h"
#include "comm.h"
#include "all_allo_const.h"

void fix_tags(LB_TAG **export_tags, int *nsentags, LB_TAG **import_tags,
	      int *nrectags, LB_TAG *kept_tags, int nkeptags);

/* 
 * void Migrate_Objects(pOctant *octants, int *newpids, int number_of_octants,
 *                      LB_TAG **export_tags, int *number_of_sent_tags,
 *                      LB_TAG **import_tags, int *number_of_received_tags)
 *
 * sets up the export_tags, and import_tags for the application that called
 * this load balancing routine 
 */
void Migrate_Objects(pOctant *octs, int *newpids, int nocts,
		     LB_TAG **export_tags, int *nsentags, 
		     LB_TAG **import_tags, int *nrectags)
{
  int i;                    /* index counter */
  int *tag_pids;            /* array of which processors to send information */
  int nkeptags;             /* number of tags kept on processors */
  LB_TAG *kept_tags;        /* tags of regions being kept on processor */

  /* tag all the regions to be exported */
  tag_regions(octs, newpids, nocts, export_tags, 
	      nsentags, &tag_pids, &kept_tags, &nkeptags);
  /* get all the region tags that are being imported */
  malloc_new_octants(*nsentags, *export_tags, tag_pids, nrectags, import_tags,
		     kept_tags, nkeptags);
  /* fix tags that were previously exported, then imported again */
  /* fix_tags(export_tags, nsentags, import_tags, nrectags, 
   *	   kept_tags, nkeptags);
   */
  
  if(nkeptags > 0)
    LB_safe_free((void **) &kept_tags);
  /*
   *  fprintf(stderr,"(%d) nrectags = %d\n", msg_mypid, (*nrectags));
   *  for(i=0; i<(*nrectags); i++)
   *    fprintf(stderr,"%d\n", (*import_tags)[i].Global_ID);
   */
}

/*
 * void tag_regions()
 * Iterates through the list of octants on the processor and finds which
 * are to be migrated. It then looks at the region list for those octants 
 * and stores the migrating regions into the export_tags array.
 */
static void tag_regions(pOctant *octs, int *newpids, int nocts, 
			LB_TAG **export_tags, int *nsentags, int **tag_pids,
			LB_TAG **kept_tags, int *nkeptags)
{
  int i;                /* index counter */
  pRegion regionlist;   /* list of region on this processor */
  pRegion region;       /* pointer to a region */
  int index;            /* index counter */
  int index2;           /* yet another index counter */
  int count;            /* count of objects exported form this processor */
  int count2;           /* count of objects that are kept on processor */
  int *export_pids;     /* array of pids where regions are being exported to */
  LB_TAG *mtags;        /* object tags of objects to be migrated */
  LB_TAG *ktags;        /* tags of objects that are kept on processor */

  /* check for migrating, pointer should not be larger than an int */
  if (sizeof(int)<sizeof(pOctant)) {
    fprintf(stderr,"tag_regions: Fatal error, sizeof(int)<sizeof(ptr)\n");
    abort();
  }

  if (!nsentags) 
    return;

  /* find how many objects are being sent */
  count = 0;
  count2 = 0;
  for (i=0; i<nocts; i++) {
    /* if newpids != msg_mypid, the it is being migrated */
    if(POC_isTerminal(octs[i])) {
      if (newpids[i]!=msg_mypid)
	count+=POC_nRegions(octs[i]);
      else {
	pRegion regions;

	regions = POC_regionlist(octs[i]);
	while(regions != NULL) {
	  if(regions->Tag.Proc != msg_mypid)
	    count2+=POC_nRegions(octs[i]);
	  regions = regions->next;
	}
      }
    }
  }
  
  /* set up the return pointers */
  *nsentags = count;
  *nkeptags = count2;
  
  if (!export_tags)
    return;
  
  if (count > 0) {
    /* allocate some space */
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
  
  if (count2 > 0) {
    /* allocate some space */
    ktags = (LB_TAG *)malloc((unsigned)count2 * sizeof(LB_TAG));
    if(ktags == NULL) {
      fprintf(stderr, "(%d)ERROR: %d unable to malloc ktags in tag_regions\n",
	      msg_mypid, count2);
      abort();
    }
  }
  else
    ktags = NULL;
  
  /* set up return pointers */
  *kept_tags=ktags;
  
  index = index2 = 0;
  for (i=0; i<nocts; i++) {
    if(POC_isTerminal(octs[i])) {
      if(newpids[i] != msg_mypid) {       /* octant being sent off processor */
	/* get regions associated with the octant */
	regionlist = POC_regionlist(octs[i]);
	while(regionlist != NULL) {
	  /* place information in the appropritate array */
	  mtags[index] = regionlist->Tag;
	  export_pids[index] = newpids[i];
	  index++;                                      /* increment counter */
	  regionlist = regionlist->next;                  /* get next region */
	}
      }
      else if(count2 > 0){               /* octant is kept on this processor */
	/* get regions associated with the octant */
	regionlist=POC_regionlist(octs[i]);
	while(regionlist != NULL) {
	  if(regionlist->Tag.Proc != msg_mypid) {
	    ktags[index2]=regionlist->Tag;	   /* get region information */
	    index2++;                                   /* increment counter */
	  }
	  regionlist = regionlist->next;                  /* get next region */
	}
      }
    }
  }
  
  if (index!=count) {                                         /* error check */
    fprintf(stderr, "%s\n",
	    "ERROR in tag_regions(), inconsistent number of regions.");
    abort();
  }
}

/*
 * void malloc_new_octants(int number_of_sent_tags, LB_TAG *export_tags,
 *                         int *tag_pids, int *number_of_received_tags,
 *                         LB_TAG **import_tags)
 *
 * ATTN: function is misnamed
 * gets the tags being imported into this processor, and sets up the
 * import_tags array, and the nrectags array.
 */
static void malloc_new_octants(int nsentags, LB_TAG *export_tags, 
			       int *tag_pids, int *nrectags,
			       LB_TAG **import_tags, LB_TAG *kept_tags,
			       int nkeptags)
{
  int i;                                  /* index counter */
  int from;                               /* from whom the region originated */
  int nreceives;                          /* number of messages received */
  LB_TAG *imp;                            /* array of tags being imported */
  LB_TAG *tmp;
  int j;
  COMM_OBJ *comm_plan;           /* Communication object returned by 
				    Bruce and Steve's communication routines */

  comm_plan = comm_create(nsentags, tag_pids, MPI_COMM_WORLD, &nreceives);
  tmp = (LB_TAG *)malloc(nreceives * sizeof(LB_TAG));
  
  if((nreceives != 0) && (tmp == NULL)) {
    fprintf(stderr,"ERROR in migreg_migrate_regions: %s\n",
	    "cannot allocate memory for import_objs.");
    abort();
  }
  
  comm_do(comm_plan, (char *) export_tags, sizeof(LB_TAG), (char *) tmp);

  comm_destroy(comm_plan);

  /* get each message sent, and store region in import array */
  j=0;
  for (i=0; i<nreceives; i++) {
    if(tmp[i].Proc != msg_mypid)
      j++;
  }
  
  if((j + nkeptags) != 0) {                   /* malloc import array */
    imp = (LB_TAG *)malloc((j + nkeptags) * sizeof(LB_TAG));
    if(imp == NULL) {
      fprintf(stderr, "ERROR in malloc_new_octants [sic], %s\n",
	      "unable to malloc import array.");
      abort();
    }
  }
  else
    imp == NULL;

  /* setup return pointer */
  *import_tags = imp;
 
  j=0;
  for (i=0; i<nreceives; i++) {
    if(tmp[i].Proc != msg_mypid)
      imp[j++] = tmp[i];
  }
  
  if(nkeptags > 0) {
    for(i=0; i<nkeptags; i++)
      imp[j++] = kept_tags[i];
  }
  *nrectags = j;

  free(tmp);
}

/*
 * void fix_tags(LB_TAG **export_tags, int *number_of_sent_tags,
 *               LB_TAG **import_tags, int *number_of_reveived_tags,
 *               LB_TAG *kept_tags, int number_of_kept_tags)
 *
 * fixes the import tags so that region tags that were previously
 * exported aren't counted when imported back.
 */
void fix_tags(LB_TAG **export_tags, int *nsentags, LB_TAG **import_tags,
	      int *nrectags, LB_TAG *kept_tags, int nkeptags) 
{
  int index;                              /* index counter */
  int i;                                  /* index counter */
  LB_TAG *new_import;                     /* modified array of import tags */

  /* allocate memory */
  new_import = (LB_TAG *)malloc(sizeof(LB_TAG) * (nkeptags + (*nrectags)));
  if(new_import == NULL) {
    fprintf(stderr,"ERROR int fix_tags, unable to allocate space\n");
    abort();
  }
  index = 0;

  /* for each region kept, look at it's origniating processor */
  for(i=0; i<nkeptags; i++) {
    if(kept_tags[i].Proc != msg_mypid)
      new_import[index++] = kept_tags[i];
  }
  /* for each region imported, look at it's origniating processor */
  for(i=0; i<(*nrectags); i++) {
    if((*import_tags)[i].Proc != msg_mypid)
      new_import[index++] = (*import_tags)[i];
  }
  
  /* setup return pointers */
  *import_tags = new_import;
  *nrectags = index;
}
