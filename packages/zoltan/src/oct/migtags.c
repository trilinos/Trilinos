

#include <unistd.h>
#include "lb_const.h"
#include "octree_const.h"
#include "migoct_const.h"
#include "migtags_const.h"
#include "comm_const.h"
#include "all_allo_const.h"

/* function prototypes */

static int tag_regions(LB *lb, pOctant *octs, int *newpids, int nocts,
                        Region **export_tags, 
                        ZOLTAN_ID_PTR *export_gids, ZOLTAN_ID_PTR *export_lids,
                        int *nsentags, int **tag_pids,
                        Region **prev_tags, 
                        ZOLTAN_ID_PTR *prev_gids, ZOLTAN_ID_PTR *prev_lids,
                        int *npimtags, float *c2,
                        int *max_objs);

static int malloc_new_objects(LB *lb, int nsentags, pRegion export_tags, 
                               ZOLTAN_ID_PTR export_gids, ZOLTAN_ID_PTR export_lids,
			       int *tag_pids, int *nrectags,
			       pRegion *import_tags, pRegion prev_tags,
                               ZOLTAN_ID_PTR prev_gids, ZOLTAN_ID_PTR prev_lids,
			       int npimtags, float *c3);

/* 
 * void LB_Migrate_Objects(LB *lb, pOctant *octants, int *newpids, 
 *                         int number_of_octants,
 *                         LB_TAG **export_tags, int *number_of_sent_tags,
 *                         LB_TAG **import_tags, int *number_of_received_tags)
 *
 * sets up the export_tags, and import_tags for the application that called
 * this load balancing routine 
 */
int LB_Migrate_Objects(LB *lb, pOctant *octs, int *newpids, int nocts,
		        int *nsenregs, 
		        pRegion *import_regions, int *nrecregs,
		        float *c2, float *c3, int *counter3, int *counter4)
{
  int i;                    /* index counter */
  int *tag_pids;            /* array of which processors to send information */
  int npimregs;             /* number of regions previously imported */
  Region *pimreg;           /* previously imported regions */
  ZOLTAN_ID_PTR pim_gids;       /* global IDs of previously imported regions */
  ZOLTAN_ID_PTR pim_lids;       /* local IDs of previously imported regions */
  ZOLTAN_ID_PTR export_gids, export_lids;
  Region *export_regions;
  int max_objs;
  int ierr = ZOLTAN_OK;
  char *yo = "LB_Migrate_Objects";
  pimreg = *import_regions = export_regions = NULL;
  export_gids = pim_gids = NULL;
  export_lids = pim_lids = NULL;

  /* tag all the regions to be exported */
  ierr = tag_regions(lb, octs, newpids, nocts, &export_regions, 
                     &export_gids, &export_lids, nsenregs, &tag_pids, 
	      &pimreg, &pim_gids, &pim_lids, &npimregs, c2, &max_objs);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_LB_TRACE_EXIT(lb, yo);
    return (ierr);
  }

  /* get all the region tags that are being imported */
  ierr = malloc_new_objects(lb, *nsenregs, export_regions, 
                     export_gids, export_lids, tag_pids, nrecregs, 
		     import_regions, pimreg, pim_gids, pim_lids, npimregs, c3);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_LB_TRACE_EXIT(lb, yo);
    return (ierr);
  }

  if(npimregs > 0){
    ZOLTAN_FREE(&pimreg);
    ZOLTAN_FREE(&pim_gids);
    ZOLTAN_FREE(&pim_lids);
  }
  ZOLTAN_FREE(&export_regions);
  ZOLTAN_FREE(&export_gids);
  ZOLTAN_FREE(&export_lids);

  if(max_objs > (*counter3))
    (*counter3) = max_objs;
  i = (max_objs - (*nsenregs) + (*nrecregs) - npimregs);
  if(i > (*counter3))
    (*counter3) = i;
  (*counter4) = (*nrecregs) - npimregs;

  /*  fprintf(stderr,"(%d) nrectags = %d\n", lb->Proc, (*nrecregs));
   *  for(i=0; i<(*nrecregs); i++)
   *    fprintf(stderr,"%d\n", (*import_regions)[i].Proc);
   */
  ZOLTAN_FREE(&tag_pids);
  return ierr;
}

/*
 * void tag_regions()
 * Iterates through the list of octants on the processor and finds which
 * are to be migrated. It then looks at the region list for those octants 
 * and stores the migrating regions into the export_tags array.
 */
static int tag_regions(LB *lb, pOctant *octs, int *newpids, int nocts, 
			Region **export_tags, 
                        ZOLTAN_ID_PTR *export_gids, ZOLTAN_ID_PTR *export_lids,
                        int *nsentags, int **tag_pids,
			Region **prev_tags, 
                        ZOLTAN_ID_PTR *prev_gids, ZOLTAN_ID_PTR *prev_lids,
                        int *npimtags, float *c2,
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
  int ierr = ZOLTAN_OK;
  int num_gid_entries = lb->Num_GID;
  int num_lid_entries = lb->Num_LID;

  ex_load = 0;
  (*max_objs) = 0;

  /* check for migrating, pointer should not be larger than an int */
#ifdef KDDKDD
  /* KDDKDD -- DON'T KNOW WHY THIS TEST IS NEEDED  6/2000 */
  if (sizeof(int)<sizeof(pOctant)) {
    ZOLTAN_LB_TRACE_DETAIL(lb, yo, "Fatal error, sizeof(int)<sizeof(ptr)\n");
    return ZOLTAN_FATAL;
  }
#endif /* KDDKDD */

  if (!nsentags) 
    return ierr;

  /* find how many objects are being sent */
  count = 0;
  count2 = 0;
  
  for (i=0; i<nocts; i++) {
    /* if newpids != lb->Proc, the it is being migrated */
    if(LB_Oct_isTerminal(octs[i])) {
      (*max_objs) += LB_Oct_nRegions(octs[i]);
      if (newpids[i]!=lb->Proc) {
	count+=LB_Oct_nRegions(octs[i]);
      }
      else {
	pRegion regions;

	regions = LB_Oct_regionlist(octs[i]);
	while(regions != NULL) {
	  if(regions->Proc != lb->Proc)
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
    return ierr;
  }

  if (count > 0) {
    /* allocate some space */
    if((mtags = (pRegion)ZOLTAN_MALLOC((unsigned)count * sizeof(Region))) == NULL){
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      ZOLTAN_LB_TRACE_EXIT(lb, yo);
      return ZOLTAN_MEMERR;
    }
    if((export_pids = (int *)ZOLTAN_MALLOC((unsigned)count * sizeof(int))) == NULL){
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      ZOLTAN_LB_TRACE_EXIT(lb, yo);
      ZOLTAN_FREE(&mtags);
      return ZOLTAN_MEMERR;
    }
    *export_gids = ZOLTAN_ZOLTAN_MALLOC_GID_ARRAY(lb, count);
    *export_lids = ZOLTAN_ZOLTAN_MALLOC_LID_ARRAY(lb, count);
    if(!(*export_gids) || (num_lid_entries && !(*export_lids))) {
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      ZOLTAN_LB_TRACE_EXIT(lb, yo);
      ZOLTAN_FREE(&mtags);
      ZOLTAN_FREE(&export_pids);
      return ZOLTAN_MEMERR;
    }
  }
  else {
    mtags = NULL;
    export_pids = NULL;
    *export_gids = NULL;
    *export_lids = NULL;
  }
  /* set up return pointers */
  *export_tags=mtags;
  *tag_pids = export_pids;
  
  if (count2 > 0) {
    /* allocate some space */
    if((ptags = (pRegion)ZOLTAN_MALLOC((unsigned)count2 * sizeof(Region))) == NULL){
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient Memory.");
      ZOLTAN_LB_TRACE_EXIT(lb, yo);
      ZOLTAN_FREE(&mtags);
      ZOLTAN_FREE(&export_pids);
      ZOLTAN_FREE(export_gids);
      ZOLTAN_FREE(export_lids);
      return ZOLTAN_MEMERR;
    }
    *prev_gids = ZOLTAN_ZOLTAN_MALLOC_GID_ARRAY(lb, count2);
    *prev_lids = ZOLTAN_ZOLTAN_MALLOC_LID_ARRAY(lb, count2);
    if(!(*prev_gids) || (num_lid_entries && !(*prev_lids))) {
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient Memory.");
      ZOLTAN_LB_TRACE_EXIT(lb, yo);
      ZOLTAN_FREE(&mtags);
      ZOLTAN_FREE(&export_pids);
      ZOLTAN_FREE(export_gids);
      ZOLTAN_FREE(export_lids);
      ZOLTAN_FREE(&ptags);
      ZOLTAN_FREE(prev_gids);
      ZOLTAN_FREE(prev_lids);
      return ZOLTAN_MEMERR;
    }
  }
  else {
    ptags = NULL;
    *prev_gids = NULL;
    *prev_lids = NULL;
  }
  
  /* set up return pointers */
  *prev_tags=ptags;
  
  index = index2 = 0;
  for (i=0; i<nocts; i++) {
    if(LB_Oct_isTerminal(octs[i])) {
      if(newpids[i] != lb->Proc) {       /* octant being sent off processor */
	/* get regions associated with the octant */
	regionlist = LB_Oct_regionlist(octs[i]);
	while(regionlist != NULL) {
	  /* place information in the appropritate array */
	  mtags[index] = *regionlist;
          ZOLTAN_LB_SET_GID(lb, &((*export_gids)[index*num_gid_entries]), regionlist->Global_ID);
          ZOLTAN_LB_SET_LID(lb, &((*export_lids)[index*num_lid_entries]), regionlist->Local_ID);

	  ex_load += (float)(regionlist->Weight);
	  export_pids[index] = newpids[i];
	  index++;                                      /* increment counter */
	  regionlist = regionlist->next;                  /* get next region */
	}
      }
      else if(count2 > 0){               /* octant is kept on this processor */
	/* get regions associated with the octant */
	regionlist=LB_Oct_regionlist(octs[i]);
	while(regionlist != NULL) {
	  if(regionlist->Proc != lb->Proc) {
	    ptags[index2] = *regionlist;	   /* get region information */
            ZOLTAN_LB_SET_GID(lb, &((*prev_gids)[index2*num_gid_entries]), regionlist->Global_ID);
            ZOLTAN_LB_SET_LID(lb, &((*prev_lids)[index2*num_lid_entries]), regionlist->Local_ID);

	    index2++;                                   /* increment counter */
	  }
	  regionlist = regionlist->next;                  /* get next region */
	}
      }
    }
  }
  
  if (index!=count) {                                         /* error check */
    ZOLTAN_LB_TRACE_DETAIL(lb, yo, "Fatal error, inconsistent number of regions.\n");
    return ZOLTAN_FATAL;
  }
  *c2 = ex_load;
  return ierr;
}

/*
 * void malloc_new_objects(LB *lb, int number_of_sent_tags, LB_TAG *export_tags,
 *                         int *tag_pids, int *number_of_received_tags,
 *                         LB_TAG **import_tags)
 *
 * gets the tags being imported into this processor, and sets up the
 * import_tags array, and the nrectags array.
 */
static int malloc_new_objects(LB *lb, int nsentags, pRegion export_tags, 
                               ZOLTAN_ID_PTR export_gids, ZOLTAN_ID_PTR export_lids,
			       int *tag_pids, int *nrectags,
			       pRegion *import_tags, pRegion prev_tags,
                               ZOLTAN_ID_PTR prev_gids, ZOLTAN_ID_PTR prev_lids,
			       int npimtags, float *c3)
{
  char *yo = "malloc_new_objects";
  int i;                                  /* index counter */
  int nreceives;                          /* number of messages received */
  pRegion imp;                            /* array of tags being imported */
  pRegion tmp = NULL;
  ZOLTAN_ID_PTR tmp_gids = NULL;
  ZOLTAN_ID_PTR tmp_lids = NULL;
  int msgtag, msgtag2;
  int j;
  int ierr = ZOLTAN_OK;
  int num_gid_entries = lb->Num_GID;
  int num_lid_entries = lb->Num_LID;
  float im_load;
  ZOLTAN_COMM_OBJ *comm_plan;           /* Object returned by communication routines */

  im_load = 0;
  msgtag = 32767;

  ierr = Zoltan_Comm_Create(&comm_plan, nsentags, tag_pids, lb->Communicator,
				msgtag, &nreceives);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Error returned from Zoltan_Comm_Create.");
    ZOLTAN_LB_TRACE_EXIT(lb, yo);
    return(ierr);
  }


  if (nreceives > 0) {
    tmp = (pRegion) ZOLTAN_MALLOC(nreceives * sizeof(Region));
    tmp_gids = ZOLTAN_ZOLTAN_MALLOC_GID_ARRAY(lb, nreceives);
    tmp_lids = ZOLTAN_ZOLTAN_MALLOC_LID_ARRAY(lb, nreceives);
    if(tmp == NULL || !tmp_gids || (num_lid_entries && !tmp_lids)) {
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      ZOLTAN_FREE(&tmp);
      ZOLTAN_FREE(&tmp_gids);
      ZOLTAN_FREE(&tmp_lids);
      ZOLTAN_LB_TRACE_EXIT(lb, yo);
      return ZOLTAN_MEMERR;
    }
  }
  
  msgtag2 = 32766;
  ierr = Zoltan_Comm_Do(comm_plan, msgtag2, (char *) export_tags, sizeof(Region),
         (char *) tmp);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Error returned from Zoltan_Comm_Do.");
    ZOLTAN_LB_TRACE_EXIT(lb, yo);
    ZOLTAN_FREE(&tmp);
    ZOLTAN_FREE(&tmp_gids);
    ZOLTAN_FREE(&tmp_lids);
    return(ierr);
  }

  msgtag2--;
  ierr = Zoltan_Comm_Do(comm_plan, msgtag2, (char *) export_gids,
                    sizeof(ZOLTAN_ID_TYPE)*num_gid_entries, (char *) tmp_gids);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Error returned from Zoltan_Comm_Do.");
    fprintf(stderr, "OCT %s Error %s returned from Zoltan_Comm_Do\n", yo,
            (ierr == ZOLTAN_MEMERR ? "ZOLTAN_MEMERR" : "ZOLTAN_FATAL"));
    ZOLTAN_FREE(&tmp);
    ZOLTAN_FREE(&tmp_gids);
    ZOLTAN_FREE(&tmp_lids);
    return(ierr);
  }

  if (num_lid_entries > 0) {
    msgtag2--;
    ierr = Zoltan_Comm_Do(comm_plan, msgtag2, (char *) export_lids,
                      sizeof(ZOLTAN_ID_TYPE)*num_lid_entries, (char *) tmp_lids);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Error returned from Zoltan_Comm_Do.");
      ZOLTAN_FREE(&tmp);
      ZOLTAN_FREE(&tmp_gids);
      ZOLTAN_FREE(&tmp_lids);
      return(ierr);
    }
  }

  ierr = Zoltan_Comm_Destroy(&comm_plan);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Error returned from Zoltan_Comm_Destroy.");
    ZOLTAN_LB_TRACE_EXIT(lb, yo);
    ZOLTAN_FREE(&tmp);
    return(ierr);
  }

  /* get each message sent, and store region in import array */
  j=0;
  for (i=0; i<nreceives; i++) {
    im_load += tmp[i].Weight;
    if(tmp[i].Proc != lb->Proc) {
      j++;
    }
  }
  
  if((j + npimtags) != 0) {                   /* malloc import array */
    if((imp = (pRegion) ZOLTAN_MALLOC((j + npimtags) * sizeof(Region))) == NULL) {
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      ZOLTAN_LB_TRACE_EXIT(lb, yo);
      ZOLTAN_FREE(&tmp);
      ZOLTAN_FREE(&tmp_gids);
      ZOLTAN_FREE(&tmp_lids);
      return ZOLTAN_MEMERR;
    }
  }
  else
    imp = NULL;

  /* setup return pointer */
  (*import_tags) = imp;

  j=0;
  for (i=0; i<nreceives; i++) {
    if(tmp[i].Proc != lb->Proc) {
      imp[j] = tmp[i];
      imp[j].Global_ID = ZOLTAN_ZOLTAN_MALLOC_GID(lb);
      imp[j].Local_ID  = ZOLTAN_ZOLTAN_MALLOC_LID(lb);
      if (!(imp[j].Global_ID) || (num_lid_entries && !(imp[j].Local_ID))) {
        ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
        ZOLTAN_LB_TRACE_EXIT(lb, yo);
        ZOLTAN_FREE(&tmp);
        ZOLTAN_FREE(&tmp_gids);
        ZOLTAN_FREE(&tmp_lids);
        return ZOLTAN_MEMERR;
      }
      ZOLTAN_LB_SET_GID(lb, imp[j].Global_ID, &(tmp_gids[i*num_gid_entries]));
      ZOLTAN_LB_SET_LID(lb, imp[j].Local_ID,  &(tmp_lids[i*num_lid_entries]));
      j++;
    }
  }
  
  if(npimtags > 0) {
    for(i=0; i<npimtags; i++) {
      imp[j] = prev_tags[i];
      imp[j].Global_ID = ZOLTAN_ZOLTAN_MALLOC_GID(lb);
      imp[j].Local_ID  = ZOLTAN_ZOLTAN_MALLOC_LID(lb);
      if (!(imp[j].Global_ID) || (num_lid_entries && !(imp[j].Local_ID))) {
        ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
        ZOLTAN_LB_TRACE_EXIT(lb, yo);
        ZOLTAN_FREE(&tmp);
        ZOLTAN_FREE(&tmp_gids);
        ZOLTAN_FREE(&tmp_lids);
        return ZOLTAN_MEMERR;
      }
      ZOLTAN_LB_SET_GID(lb, imp[j].Global_ID, &(prev_gids[i*num_gid_entries]));
      ZOLTAN_LB_SET_LID(lb, imp[j].Local_ID,  &(prev_lids[i*num_lid_entries]));
      j++;
    }
  }
  *nrectags = j;

  ZOLTAN_FREE(&tmp);
  ZOLTAN_FREE(&tmp_gids);
  ZOLTAN_FREE(&tmp_lids);

  /*
   * fprintf(stderr,
   *     "(%d) nrectags = %d, nreceives = %d, nsentags = %d, nkeptags = %d\n", 
   *     lb->Proc, (*nrectags), nreceives, nsentags, nkeptags);
   */

  if((*nrectags == 0) && (*import_tags != NULL)) {
    ZOLTAN_LB_TRACE_DETAIL(lb, yo, "Fatal error, import tags not empty but no tags received\n");
    return ZOLTAN_FATAL;
  }

  /*  for(i=0; i<(*nrectags); i++) {
   *    fprintf(stderr,"%d -> %d\n", (*import_tags)[i].Proc, lb->Proc);
   *  }
   */  
  *c3 = im_load;
  return ierr;
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
int LB_fix_tags(LB *lb, ZOLTAN_ID_PTR *import_global_ids, ZOLTAN_ID_PTR *import_local_ids,
                 int **import_procs, int nrectags, pRegion import_regs)
{
  char *yo = "LB_fix_tags";
  int i;                                  /* index counter */
  int ierr = ZOLTAN_OK;
  int num_gid_entries = lb->Num_GID;
  int num_lid_entries = lb->Num_LID;

    /* allocate memory */

    if (!LB_Special_Malloc(lb,(void **)import_global_ids,nrectags,
                           LB_SPECIAL_MALLOC_GID)) {
      ZOLTAN_PRINT_ERROR(lb->Proc,yo, "Insufficient memory.");
      ZOLTAN_LB_TRACE_EXIT(lb, yo);
      return ZOLTAN_MEMERR;
    }
    if (!LB_Special_Malloc(lb,(void **)import_local_ids,nrectags,
                           LB_SPECIAL_MALLOC_LID)) {
      LB_Special_Free(lb,(void **)import_global_ids,LB_SPECIAL_MALLOC_GID); 
      ZOLTAN_PRINT_ERROR(lb->Proc,yo, "Insufficient memory.");
      ZOLTAN_LB_TRACE_EXIT(lb, yo);
      return ZOLTAN_MEMERR;
    }
    if (!LB_Special_Malloc(lb,(void **)import_procs,nrectags,
                           LB_SPECIAL_MALLOC_INT)) {
      LB_Special_Free(lb,(void **)import_global_ids,LB_SPECIAL_MALLOC_GID);
      LB_Special_Free(lb,(void **)import_local_ids,LB_SPECIAL_MALLOC_LID);
      ZOLTAN_PRINT_ERROR(lb->Proc,yo, "Insufficient memory.");
      ZOLTAN_LB_TRACE_EXIT(lb, yo);
      return ZOLTAN_MEMERR;
    }

    /* for each region imported, look at its originating processor */
    for(i=0; i<nrectags; i++) {
      ZOLTAN_LB_SET_GID(lb, &((*import_global_ids)[i*num_gid_entries]),
                 import_regs[i].Global_ID);
      ZOLTAN_LB_SET_LID(lb, &((*import_local_ids)[i*num_lid_entries]),
                 import_regs[i].Local_ID);
      (*import_procs)[i]      = import_regs[i].Proc;
    }

    return ierr;
}

