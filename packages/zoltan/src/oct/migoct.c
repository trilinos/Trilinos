/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include <unistd.h>
#include "lb_const.h"
#include "octant_const.h"
#include "migoct_const.h"
#include "comm_const.h"
#include "all_allo_const.h"

/* function prototypes */
/* static void tag_regions(pOctant *octs, int *newpids, int nocts,
                        LB_TAG **export_tags, int *nsentags, int **tag_pids,
                        LB_TAG **kept_tags, int *nkeptags); */
static void tag_regions(LB *, pOctant *octs, int *newpids, int nocts,
                        Region **export_tags, int *nsentags, int **tag_pids,
                        Region **prev_tags, int *npimtags, float *c2,
                        int *max_objs);

/* static void malloc_new_objects(int ntags, LB_TAG *tags, int *tag_pids,
                               int *nrectags, LB_TAG **import_tags,
                               LB_TAG *kept_tags, int nkeptags); */
static void malloc_new_objects(LB *lb, int nsentags, pRegion export_tags,
                               int *tag_pids, int *nrectags,
                               pRegion *import_tags, pRegion prev_tags,
                               int npimtags, float *c3);

/* void LB_fix_tags(LB_TAG **export_tags, int *nsentags, LB_TAG **import_tags,
 *	            int *nrectags, LB_TAG *prev_regs, int npimregs,
 *	            pRegion import_regs, pRegion export_regs);
 */

/* 
 * void LB_Migrate_Objects(LB *lb, pOctant *octants, int *newpids, 
 *                         int number_of_octants,
 *                         LB_TAG **export_tags, int *number_of_sent_tags,
 *                         LB_TAG **import_tags, int *number_of_received_tags)
 *
 * sets up the export_tags, and import_tags for the application that called
 * this load balancing routine 
 */
void LB_Migrate_Objects(LB *lb, pOctant *octs, int *newpids, int nocts,
		        pRegion *export_regions, int *nsenregs, 
		        pRegion *import_regions, int *nrecregs,
		        float *c2, float *c3, int *counter3, int *counter4)
{
  int i;                    /* index counter */
  int *tag_pids;            /* array of which processors to send information */
  int npimregs;             /* number of regions previously imported */
  Region *pimreg;           /* previously imported regions */
  int max_objs;

  *import_regions = *export_regions = NULL;

  /* tag all the regions to be exported */
  tag_regions(lb, octs, newpids, nocts, export_regions, nsenregs, &tag_pids, 
	      &pimreg, &npimregs, c2, &max_objs);

  /* get all the region tags that are being imported */
  malloc_new_objects(lb, *nsenregs, *export_regions, tag_pids, nrecregs, 
		     import_regions, pimreg, npimregs, c3);

  if(npimregs > 0)
    LB_FREE(&pimreg);

  if(max_objs > (*counter3))
   (*counter3) = max_objs;
   i = (max_objs - (*nsenregs) + (*nrecregs) - npimregs);
   if(i > (*counter3))
     (*counter3) = i;
   (*counter4) = (*nrecregs) - npimregs;

  /*  fprintf(stderr,"(%d) nrectags = %d\n", lb->Proc, (*nrecregs));
   *  for(i=0; i<(*nrecregs); i++)
   *    fprintf(stderr,"%d\n", (*import_regions)[i].Tag.Proc);
   */
  LB_FREE(&tag_pids);
}

/*
 * void tag_regions()
 * Iterates through the list of octants on the processor and finds which
 * are to be migrated. It then looks at the region list for those octants 
 * and stores the migrating regions into the export_tags array.
 */
static void tag_regions(LB *lb, pOctant *octs, int *newpids, int nocts, 
			Region **export_tags, int *nsentags, int **tag_pids,
			Region **prev_tags, int *npimtags, float *c2,
			int *max_objs)
{
  char *yo = "tag_regions";
  int i;                /* index counter */
  pRegion regionlist;   /* list of region on this processor */
  int index;            /* index counter */
  int index2;           /* yet another index counter */
  int count;            /* count of objects exported form this processor */
  int count2;           /* count of objects that are kept on processor */
  int *export_pids;     /* array of pids where regions are being exported to */
  pRegion mtags;        /* object tags of objects to be migrated */
  pRegion ptags;        /* tags of objects that were previously migrated */
  float ex_load;

  ex_load = 0;
  (*max_objs) = 0;

  /* check for migrating, pointer should not be larger than an int */
#ifdef KDDKDD
/* KDD  Don't know why this test is needed   6/2000 */
  if (sizeof(int)<sizeof(pOctant)) {
    fprintf(stderr,"OCT %s: Fatal error, sizeof(int)<sizeof(ptr)\n", yo);
    abort();
  }
#endif

  if (!nsentags) 
    return;

  /* find how many objects are being sent */
  count = 0;
  count2 = 0;
  
  for (i=0; i<nocts; i++) {
    /* if newpids != lb->Proc, the it is being migrated */
    if(POC_isTerminal(octs[i])) {
      (*max_objs) += POC_nRegions(octs[i]);
      if (newpids[i]!=lb->Proc) {
	count+=POC_nRegions(octs[i]);
      }
      else {
	pRegion regions;

	regions = POC_regionlist(octs[i]);
	while(regions != NULL) {
	  if(regions->Tag.Proc != lb->Proc)
	    count2++;
	  regions = regions->next;
	}
      }
    }
  }
  
  /* set up the return pointers */
  *nsentags = count;
  *npimtags = count2;
  
  if (!export_tags) {
    fprintf(stderr, "OCT %s: return code reached\n", yo);
    return;
  }

  if (count > 0) {
    /* allocate some space */
    mtags = (pRegion) LB_MALLOC((unsigned)count * sizeof(Region));
    export_pids = (int *) LB_MALLOC((unsigned)count * sizeof(int));
    if(export_pids == NULL) {
      fprintf(stderr, "OCT ERROR: unable to malloc export_pids in %s\n", yo);
      abort();
    }
    if(mtags == NULL) {
      fprintf(stderr, "OCT ERROR: unable to malloc mtags in %s\n", yo);
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
    ptags = (pRegion) LB_MALLOC((unsigned)count2 * sizeof(Region));
    if(ptags == NULL) {
      fprintf(stderr, "OCT (%d)ERROR: unable to malloc %d ptags in %s\n",
	      lb->Proc, count2, yo);
      abort();
    }
  }
  else
    ptags = NULL;
  
  /* set up return pointers */
  *prev_tags=ptags;
  
  index = index2 = 0;
  for (i=0; i<nocts; i++) {
    if(POC_isTerminal(octs[i])) {
      if(newpids[i] != lb->Proc) {       /* octant being sent off processor */
	/* get regions associated with the octant */
	regionlist = POC_regionlist(octs[i]);
	while(regionlist != NULL) {
	  /* place information in the appropritate array */
	  mtags[index] = *regionlist;
	  ex_load += (float)(regionlist->Weight);
	  export_pids[index] = newpids[i];
	  index++;                                      /* increment counter */
	  regionlist = regionlist->next;                  /* get next region */
	}
      }
      else if(count2 > 0){               /* octant is kept on this processor */
	/* get regions associated with the octant */
	regionlist=POC_regionlist(octs[i]);
	while(regionlist != NULL) {
	  if(regionlist->Tag.Proc != lb->Proc) {
	    ptags[index2] = *regionlist;	   /* get region information */
	    index2++;                                   /* increment counter */
	  }
	  regionlist = regionlist->next;                  /* get next region */
	}
      }
    }
  }
  
  if (index!=count) {                                         /* error check */
    fprintf(stderr, "OCT ERROR in %s, inconsistent number of regions.", yo);
    abort();
  }
  *c2 = ex_load;
}

/*
 * void malloc_new_objects(LB *lb, int number_of_sent_tags, LB_TAG *export_tags,
 *                         int *tag_pids, int *number_of_received_tags,
 *                         LB_TAG **import_tags)
 *
 * gets the tags being imported into this processor, and sets up the
 * import_tags array, and the nrectags array.
 */
static void malloc_new_objects(LB *lb, int nsentags, pRegion export_tags, 
			       int *tag_pids, int *nrectags,
			       pRegion *import_tags, pRegion prev_tags,
			       int npimtags, float *c3)
{
  char *yo = "malloc_new_objects";
  int i;                                  /* index counter */
  int nreceives;                          /* number of messages received */
  pRegion imp;                            /* array of tags being imported */
  pRegion tmp = NULL;
  int msgtag, msgtag2;
  int j;
  int ierr;
  float im_load;
  COMM_OBJ *comm_plan;           /* Object returned by communication routines */

  im_load = 0;
  msgtag = 32767;
  ierr = LB_Comm_Create(&comm_plan, nsentags, tag_pids, lb->Communicator, 
         msgtag, &nreceives);
  if (ierr != COMM_OK && ierr != COMM_WARN) {
    fprintf(stderr, "OCT %s Error %s returned from LB_Comm_Create\n", yo, 
            (ierr == COMM_MEMERR ? "COMM_MEMERR" : "COMM_FATAL"));
    abort();
  }

  if (nreceives > 0) {
    tmp = (pRegion) LB_MALLOC(nreceives * sizeof(Region));
    if(tmp == NULL) {
      fprintf(stderr,"OCT %s ERROR cannot allocate memory for import_objs.",yo);
      abort();
    }
  }
  
  msgtag2 = 32766;
  ierr = LB_Comm_Do(comm_plan, msgtag2, (char *) export_tags, sizeof(Region),
         (char *) tmp);
  if (ierr != COMM_OK && ierr != COMM_WARN) {
    fprintf(stderr, "OCT %s Error %s returned from LB_Comm_Create\n", yo, 
            (ierr == COMM_MEMERR ? "COMM_MEMERR" : "COMM_FATAL"));
    LB_FREE(&tmp);
    abort();
  }
  LB_Comm_Destroy(&comm_plan);

  /* get each message sent, and store region in import array */
  j=0;
  for (i=0; i<nreceives; i++) {
    im_load += tmp[i].Weight;
    if(tmp[i].Tag.Proc != lb->Proc) {
      j++;
    }
  }
  
  if((j + npimtags) != 0) {                   /* malloc import array */
    imp = (pRegion) LB_MALLOC((j + npimtags) * sizeof(Region));
    if(imp == NULL) {
      fprintf(stderr, "OCT %s ERROR unable to malloc import array.", yo);
      abort();
    }
  }
  else
    imp = NULL;

  /* setup return pointer */
  (*import_tags) = imp;

  j=0;
  for (i=0; i<nreceives; i++) {
    if(tmp[i].Tag.Proc != lb->Proc) {
      imp[j++] = tmp[i];
    }
  }
  
  if(npimtags > 0) {
    for(i=0; i<npimtags; i++)
      imp[j++] = prev_tags[i];
  }
  *nrectags = j;

  LB_FREE(&tmp);

  /*
   * fprintf(stderr,
   *     "(%d) nrectags = %d, nreceives = %d, nsentags = %d, nkeptags = %d\n", 
   *     lb->Proc, (*nrectags), nreceives, nsentags, nkeptags);
   */

  if((*nrectags == 0) && (*import_tags != NULL)) {
    fprintf(stderr, "OCT %s (%d) ERROR: import tags not empty but no tags "
                    "received\n", yo, lb->Proc);
    exit(1);
  }

  /*  for(i=0; i<(*nrectags); i++) {
   *    fprintf(stderr,"%d -> %d\n", (*import_tags)[i].Proc, lb->Proc);
   *  }
   */  
  *c3 = im_load;
}

/*****************************************************************************/
/*
 * void LB_fix_tags(LB_TAG **export_tags, int *number_of_sent_tags,
 *                  LB_TAG **import_tags, int *number_of_reveived_tags,
 *                  LB_TAG *kept_tags, int number_of_kept_tags)
 *
 * fixes the import tags so that region tags that were previously
 * exported aren't counted when imported back.
 */
void LB_fix_tags(LB *lb, LB_GID **import_global_ids, LB_LID **import_local_ids,
                 int **import_procs, int nrectags, pRegion import_regs)
{
  char *yo = "LB_fix_tags";
  int i;                                  /* index counter */

  if (nrectags == 0) {
    *import_global_ids = NULL;
    *import_local_ids = NULL;
    *import_procs = NULL;
  }
  else {

    /* allocate memory */

/* KDD  -- added calls to Special_Malloc */

    if (!LB_Special_Malloc(lb,(void **)import_global_ids,nrectags,
                           LB_SPECIAL_MALLOC_GID)) {
      LB_TRACE_EXIT(lb, yo);
      fprintf(stderr,"OCT %s ERROR, unable to allocate space\n", yo);
      /* return LB_MEMERR; */ abort();
    }
    if (!LB_Special_Malloc(lb,(void **)import_local_ids,nrectags,
                           LB_SPECIAL_MALLOC_LID)) {
      LB_Special_Free(lb,(void **)import_global_ids,LB_SPECIAL_MALLOC_GID);
      LB_TRACE_EXIT(lb, yo);
      fprintf(stderr,"OCT %s ERROR, unable to allocate space\n", yo);
      /* return LB_MEMERR; */ abort();
    }
    if (!LB_Special_Malloc(lb,(void **)import_procs,nrectags,
                           LB_SPECIAL_MALLOC_INT)) {
      LB_Special_Free(lb,(void **)import_global_ids,LB_SPECIAL_MALLOC_GID);
      LB_Special_Free(lb,(void **)import_local_ids,LB_SPECIAL_MALLOC_LID);
      LB_TRACE_EXIT(lb, yo);
      fprintf(stderr,"OCT %s ERROR, unable to allocate space\n", yo);
      /* return LB_MEMERR; */ abort();
    }

/* KDD */

    /* for each region imported, look at its originating processor */
    for(i=0; i<nrectags; i++) {
      LB_SET_GID((*import_global_ids)[i], import_regs[i].Tag.Global_ID);
      LB_SET_LID((*import_local_ids)[i], import_regs[i].Tag.Local_ID);
      (*import_procs)[i]      = import_regs[i].Tag.Proc;
    }
  }
#if 0

  /* KDD -- LB_Compute_Destinations will perform this operation for us. */
  new_export = (LB_TAG *) LB_MALLOC((*nsentags) * sizeof(LB_TAG));
  if(((*nsentags) > 0) && (new_export == NULL)) {
    fprintf(stderr,"OCT %s ERROR, unable to allocate space\n", yo);
    abort();
  }

  index = 0;
  for(i=0; i<(*nsentags); i++) {
      new_export[index++] = export_regs[i].Tag;
  }
  /* setup return pointers */
  *export_tags = new_export;

#endif
}

