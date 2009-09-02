/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include "zz_const.h"
#include "octree_const.h"
#include "migoct_const.h"
#include "migtags_const.h"
#include "all_allo_const.h"

/* function prototypes */

static int tag_regions(ZZ *zz, pOctant *octs, int *newpids, int nocts,
		       Region **export_tags, ZOLTAN_ID_PTR *export_gids,
		       ZOLTAN_ID_PTR *export_lids, int *nsentags,
		       int **tag_pids, Region **prev_tags, 
		       ZOLTAN_ID_PTR *prev_gids, ZOLTAN_ID_PTR *prev_lids,
		       int *npimtags, float *c2, int *max_objs);

static int malloc_new_objects(ZZ *zz, int nsentags, pRegion export_tags, 
			      ZOLTAN_ID_PTR export_gids,
			      ZOLTAN_ID_PTR export_lids, int *tag_pids,
			      int *nrectags, pRegion *import_tags,
			      pRegion prev_tags, ZOLTAN_ID_PTR prev_gids,
			      ZOLTAN_ID_PTR prev_lids, int npimtags, 
			      float *c3);

/* 
 * void Zoltan_Oct_migrate_objects()
 *
 * sets up the export_tags, and import_tags for the application that called
 * this load balancing routine 
 */
int Zoltan_Oct_migrate_objects(ZZ *zz, pOctant *octs, int *newpids, int nocts,
			       int *nsenregs, pRegion *import_regions,
			       int *nrecregs, float *c2, float *c3,
			       int *counter3, int *counter4)
{
  int i;                   /* index counter */
  int *tag_pids;           /* array of which processors to send information */
  int np_regs=0;            /* number of regions previously imported */
  Region *p_reg;          /* previously imported regions */
  ZOLTAN_ID_PTR p_gids;      /* global IDs of previously imported regions */
  ZOLTAN_ID_PTR p_lids;      /* local IDs of previously imported regions */
  ZOLTAN_ID_PTR exported_gids, exported_lids;
  Region *exported_regions;
  int max_objs;
  int ierr = ZOLTAN_OK;
  char *yo = "Zoltan_Oct_migrate_objects";
  p_reg = *import_regions = exported_regions = NULL;
  exported_gids = p_gids = NULL;
  exported_lids = p_lids = NULL;

  /* tag all the regions to be exported */
  ierr = tag_regions(zz, octs, newpids, nocts, &exported_regions, 
                     &exported_gids, &exported_lids, nsenregs, &tag_pids, 
		     &p_reg, &p_gids, &p_lids, &np_regs, c2, &max_objs);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ierr);
  }

  /* get all the region tags that are being imported */
  ierr = malloc_new_objects(zz, *nsenregs, exported_regions, 
			    exported_gids, exported_lids, tag_pids, nrecregs, 
			    import_regions, p_reg, p_gids, p_lids,
			    np_regs, c3);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ierr);
  }

  if(np_regs > 0){
    ZOLTAN_FREE(&p_reg);
    ZOLTAN_FREE(&p_gids);
    ZOLTAN_FREE(&p_lids);
  }
  ZOLTAN_FREE(&exported_regions);
  ZOLTAN_FREE(&exported_gids);
  ZOLTAN_FREE(&exported_lids);

  if(max_objs > (*counter3))
    (*counter3) = max_objs;
  i = (max_objs - (*nsenregs) + (*nrecregs) - np_regs);
  if(i > (*counter3))
    (*counter3) = i;
  (*counter4) = (*nrecregs) - np_regs;

  /*  fprintf(stderr,"(%d) nrectags = %d\n", zz->Proc, (*nrecregs));
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
static int tag_regions(ZZ *zz,
		       pOctant *octs,
		       int *newpids,
		       int nocts, 
		       Region **exported_tags,
		       ZOLTAN_ID_PTR *exported_gids,
		       ZOLTAN_ID_PTR *exported_lids,
		       int *nsentags, 
		       int **tag_pids,
		       Region **p_tags, 
		       ZOLTAN_ID_PTR *p_gids,
		       ZOLTAN_ID_PTR *p_lids,
		       int *npimtags,
		       float *c2,
		       int *max_objs)
{
  char *yo = "tag_regions";
  int i;               /* index counter */
  pRegion regionlist;  /* list of region on this processor */
  int index;           /* index counter */
  int index2;          /* yet another index counter */
  int count;           /* count of objects exported form this processor */
  int count2;          /* count of objects that are kept on processor */
  int count3;
  int *exported_pids;  /* array of pids where regions are being exported to */
  pRegion mtags;       /* object tags of objects to be migrated */
  pRegion ptags;       /* tags of objects that were previously migrated */
  float ex_load;
  int ierr = ZOLTAN_OK;
  int num_gid_entries = zz->Num_GID;
  int num_lid_entries = zz->Num_LID;

  ex_load = 0;
  (*max_objs) = 0;

  if (!nsentags) 
    return ierr;

  /* find how many objects have been exported */
  count = 0;
  /* find number of local objs to export */
  count2 = 0;
  count3 = 0;

  for (i=0; i<nocts; i++) {
    if(Zoltan_Oct_isTerminal(octs[i])) {
      (*max_objs) += Zoltan_Oct_nRegions(octs[i]);

      regionlist = Zoltan_Oct_regionlist(octs[i]);
      while(regionlist != NULL) {
	count3++;
	if(regionlist->Proc != zz->Proc) {
	  count++;
	  if(newpids[i] != zz->Proc)
	    regionlist->newProc = newpids[i];
	  else
	    regionlist->newProc = zz->Proc;
	}
	else {
	  if(newpids[i] != zz->Proc) {
	    count2++;
	    regionlist->newProc = newpids[i];
	  }
	  else
	    regionlist->newProc = zz->Proc;
	}
	regionlist = regionlist->next;                 /* get next region */
      }
    }
  }

#if 0
  {
    {
      if (newpids[i]!=zz->Proc) {
	count+=Zoltan_Oct_nRegions(octs[i]);
      }
      else {
	pRegion regions;

	regions = Zoltan_Oct_regionlist(octs[i]);
	while(regions != NULL) {
	  if(regions->Proc != zz->Proc)
	    count2++;
	  regions = regions->next;
	}
      }
    }
  }
#endif

  /* set up the return pointers */
  *nsentags = count;
  *npimtags = count2;

  if (!exported_tags) {
    return ierr;
  }

  if (count > 0) {
    /* allocate some space */
    if((mtags=(pRegion)ZOLTAN_MALLOC((unsigned)count*sizeof(Region)))==NULL){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    if((exported_pids = (int *)ZOLTAN_MALLOC((unsigned)count*sizeof(int))) ==
       NULL){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      ZOLTAN_FREE(&mtags);
      return ZOLTAN_MEMERR;
    }
    *exported_gids = ZOLTAN_MALLOC_GID_ARRAY(zz, count);
    *exported_lids = ZOLTAN_MALLOC_LID_ARRAY(zz, count);
    if(!(*exported_gids) || (num_lid_entries && !(*exported_lids))) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      ZOLTAN_FREE(&mtags);
      ZOLTAN_FREE(&exported_pids);
      return ZOLTAN_MEMERR;
    }
  }
  else {
    mtags = NULL;
    exported_pids = NULL;
    *exported_gids = NULL;
    *exported_lids = NULL;
  }
  /* set up return pointers */
  *exported_tags=mtags;
  *tag_pids = exported_pids;
  
  if (count2 > 0) {
    /* allocate some space */
    if((ptags=(pRegion)ZOLTAN_MALLOC((unsigned)count2*sizeof(Region)))==NULL){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient Memory.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      ZOLTAN_FREE(&mtags);
      ZOLTAN_FREE(&exported_pids);
      ZOLTAN_FREE(exported_gids);
      ZOLTAN_FREE(exported_lids);
      return ZOLTAN_MEMERR;
    }
    *p_gids = ZOLTAN_MALLOC_GID_ARRAY(zz, count2);
    *p_lids = ZOLTAN_MALLOC_LID_ARRAY(zz, count2);
    if(!(*p_gids) || (num_lid_entries && !(*p_lids))) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient Memory.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      ZOLTAN_FREE(&mtags);
      ZOLTAN_FREE(&exported_pids);
      ZOLTAN_FREE(exported_gids);
      ZOLTAN_FREE(exported_lids);
      ZOLTAN_FREE(&ptags);
      ZOLTAN_FREE(p_gids);
      ZOLTAN_FREE(p_lids);
      return ZOLTAN_MEMERR;
    }
  }
  else {
    ptags = NULL;
    *p_gids = NULL;
    *p_lids = NULL;
  }
  
  /* set up return pointers */
  *p_tags=ptags;
  
  index = index2 = 0;
  for (i=0; i<nocts; i++) {
    if(Zoltan_Oct_isTerminal(octs[i])) {
      regionlist = Zoltan_Oct_regionlist(octs[i]);
      while(regionlist != NULL) {
	if(regionlist->Proc != zz->Proc) {
	  /* place information in the appropritate array */
	  mtags[index] = *regionlist;
          ZOLTAN_SET_GID(zz, &((*exported_gids)[index*num_gid_entries]), 
			 regionlist->Global_ID);
          ZOLTAN_SET_LID(zz, &((*exported_lids)[index*num_lid_entries]),
			 regionlist->Local_ID);

	  /*ex_load += (float)(regionlist->Weight);*/
	  exported_pids[index] = regionlist->Proc;
	  index++;                                     /* increment counter */
	}
	else if(newpids[i] != zz->Proc) { 
	  ptags[index2] = *regionlist;	  /* get region information */
	  ZOLTAN_SET_GID(zz, &((*p_gids)[index2*num_gid_entries]),
			 regionlist->Global_ID);
	  ZOLTAN_SET_LID(zz, &((*p_lids)[index2*num_lid_entries]),
			 regionlist->Local_ID);
	  
	  index2++;                                  /* increment counter */
	}
	regionlist = regionlist->next;                 /* get next region */
      }
    }
  }
  
  if (index!=count) {                                        /* error check */
    ZOLTAN_TRACE_DETAIL(zz, yo, 
			"Fatal error, inconsistent number of regions.\n");
    return ZOLTAN_FATAL;
  }
  *c2 = ex_load;
  return ierr;
}

/*
 * void malloc_new_objects();
 *
 * gets the tags being imported into this processor, and sets up the
 * import_tags array, and the nrectags array.
 */
static int malloc_new_objects(ZZ *zz, int nsentags, pRegion exported_tags, 
			      ZOLTAN_ID_PTR exported_gids,
			      ZOLTAN_ID_PTR exported_lids, int *tag_pids,
			      int *nstags, pRegion *ex_tags,
			      pRegion prev_tags, ZOLTAN_ID_PTR prev_gids,
			      ZOLTAN_ID_PTR prev_lids, int npimtags,
			      float *c3)
{
  char *yo = "malloc_new_objects";
  int i;                                  /* index counter */
  int nreceives;                          /* number of messages received */
  pRegion t_b_exp;                        /* array of tags to be exported */
  pRegion tmp = NULL;
  ZOLTAN_ID_PTR tmp_gids = NULL;
  ZOLTAN_ID_PTR tmp_lids = NULL;
  int msgtag, msgtag2;
  int j;
  int ierr = ZOLTAN_OK;
  int num_gid_entries = zz->Num_GID;
  int num_lid_entries = zz->Num_LID;
  float im_load;
  ZOLTAN_COMM_OBJ *comm_plan;  /* Object returned by communication routines */

  im_load = 0;
  msgtag = 32767;

  ierr = Zoltan_Comm_Create(&comm_plan, nsentags, tag_pids, zz->Communicator,
				msgtag, &nreceives);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
		       "Error returned from Zoltan_Comm_Create.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ierr);
  }


  if (nreceives > 0) {
    tmp = (pRegion) ZOLTAN_MALLOC(nreceives * sizeof(Region));
    tmp_gids = ZOLTAN_MALLOC_GID_ARRAY(zz, nreceives);
    tmp_lids = ZOLTAN_MALLOC_LID_ARRAY(zz, nreceives);
    if(tmp == NULL || !tmp_gids || (num_lid_entries && !tmp_lids)) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_FREE(&tmp);
      ZOLTAN_FREE(&tmp_gids);
      ZOLTAN_FREE(&tmp_lids);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
  }
  
  msgtag2 = 32766;
  ierr = Zoltan_Comm_Do(comm_plan, msgtag2, (char *) exported_tags,
			sizeof(Region), (char *) tmp);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_Comm_Do.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    ZOLTAN_FREE(&tmp);
    ZOLTAN_FREE(&tmp_gids);
    ZOLTAN_FREE(&tmp_lids);
    return(ierr);
  }

  msgtag2--;
  ierr = Zoltan_Comm_Do(comm_plan, msgtag2, (char *) exported_gids,
			sizeof(ZOLTAN_ID_TYPE)*num_gid_entries, 
			(char *) tmp_gids);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_Comm_Do.");
    fprintf(stderr, "OCT %s Error %s returned from Zoltan_Comm_Do\n", yo,
            (ierr == ZOLTAN_MEMERR ? "ZOLTAN_MEMERR" : "ZOLTAN_FATAL"));
    ZOLTAN_FREE(&tmp);
    ZOLTAN_FREE(&tmp_gids);
    ZOLTAN_FREE(&tmp_lids);
    return(ierr);
  }

  if (num_lid_entries > 0) {
    msgtag2--;
    ierr = Zoltan_Comm_Do(comm_plan, msgtag2, (char *) exported_lids,
			  sizeof(ZOLTAN_ID_TYPE)*num_lid_entries,
			  (char *) tmp_lids);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_Comm_Do.");
      ZOLTAN_FREE(&tmp);
      ZOLTAN_FREE(&tmp_gids);
      ZOLTAN_FREE(&tmp_lids);
      return(ierr);
    }
  }

  ierr = Zoltan_Comm_Destroy(&comm_plan);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
		       "Error returned from Zoltan_Comm_Destroy.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    ZOLTAN_FREE(&tmp);
    return(ierr);
  }

  /* get each message sent, and store region in import array */
  j=0;
  for (i=0; i<nreceives; i++) {
    im_load += tmp[i].Weight;
    if(tmp[i].newProc != zz->Proc) {
      j++;
    }
  }
  
  if((j + npimtags) != 0) {                   /* malloc import array */
    if((t_b_exp = (pRegion)ZOLTAN_MALLOC((j+npimtags)*sizeof(Region)))==NULL){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      ZOLTAN_FREE(&tmp);
      ZOLTAN_FREE(&tmp_gids);
      ZOLTAN_FREE(&tmp_lids);
      return ZOLTAN_MEMERR;
    }
  }
  else
    t_b_exp = NULL;

  /* setup return pointer */
  (*ex_tags) = t_b_exp;

  j=0;
  for (i=0; i<nreceives; i++) {
    if(tmp[i].newProc != zz->Proc) {
      t_b_exp[j] = tmp[i];
      t_b_exp[j].Global_ID = ZOLTAN_MALLOC_GID(zz);
      t_b_exp[j].Local_ID  = ZOLTAN_MALLOC_LID(zz);
      if (!(t_b_exp[j].Global_ID) ||
	  (num_lid_entries && !(t_b_exp[j].Local_ID))) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        ZOLTAN_TRACE_EXIT(zz, yo);
        ZOLTAN_FREE(&tmp);
        ZOLTAN_FREE(&tmp_gids);
        ZOLTAN_FREE(&tmp_lids);
        return ZOLTAN_MEMERR;
      }
      ZOLTAN_SET_GID(zz, t_b_exp[j].Global_ID,
		     &(tmp_gids[i*num_gid_entries]));
      ZOLTAN_SET_LID(zz, t_b_exp[j].Local_ID, &(tmp_lids[i*num_lid_entries]));
      j++;
    }
  }
  
  if(npimtags > 0) {
    for(i=0; i<npimtags; i++) {
      t_b_exp[j] = prev_tags[i];
      t_b_exp[j].Global_ID = ZOLTAN_MALLOC_GID(zz);
      t_b_exp[j].Local_ID  = ZOLTAN_MALLOC_LID(zz);
      if (!(t_b_exp[j].Global_ID) ||
	  (num_lid_entries && !(t_b_exp[j].Local_ID))) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        ZOLTAN_TRACE_EXIT(zz, yo);
        ZOLTAN_FREE(&tmp);
        ZOLTAN_FREE(&tmp_gids);
        ZOLTAN_FREE(&tmp_lids);
        return ZOLTAN_MEMERR;
      }
      ZOLTAN_SET_GID(zz, t_b_exp[j].Global_ID,
		     &(prev_gids[i*num_gid_entries]));
      ZOLTAN_SET_LID(zz, t_b_exp[j].Local_ID, 
		     &(prev_lids[i*num_lid_entries]));
      j++;
    }
  }
  *nstags = j;

  ZOLTAN_FREE(&tmp);
  ZOLTAN_FREE(&tmp_gids);
  ZOLTAN_FREE(&tmp_lids);

  if((*nstags == 0) && (*ex_tags != NULL)) {
    ZOLTAN_TRACE_DETAIL(zz, yo,
                 "Fatal error, import tags not empty but no tags received\n");
    return ZOLTAN_FATAL;
  }

  *c3 = im_load;
  return ierr;
}

/****************************************************************************/
/*
 * void Zoltan_Oct_fix_tags()
 *
 * fixes the import tags so that region tags that were previously
 * exported aren't counted when imported back.
 */
int Zoltan_Oct_fix_tags(ZZ *zz, ZOLTAN_ID_PTR *import_global_ids, 
			ZOLTAN_ID_PTR *import_local_ids, int **import_procs, 
			int **import_to_part, int nrectags, 
			pRegion import_regs) {
  char *yo = "Zoltan_Oct_fix_tags";
  int i;                                  /* index counter */
  int ierr = ZOLTAN_OK;
  int num_gid_entries = zz->Num_GID;
  int num_lid_entries = zz->Num_LID;

    /* allocate memory */

    if (!Zoltan_Special_Malloc(zz,(void **)import_global_ids,nrectags,
                           ZOLTAN_SPECIAL_MALLOC_GID)) {
      ZOLTAN_PRINT_ERROR(zz->Proc,yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    if (!Zoltan_Special_Malloc(zz,(void **)import_local_ids,nrectags,
                           ZOLTAN_SPECIAL_MALLOC_LID)) {
      Zoltan_Special_Free(zz,(void **)import_global_ids,
			  ZOLTAN_SPECIAL_MALLOC_GID); 
      ZOLTAN_PRINT_ERROR(zz->Proc,yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    if (!Zoltan_Special_Malloc(zz,(void **)import_procs,nrectags,
                           ZOLTAN_SPECIAL_MALLOC_INT)) {
      Zoltan_Special_Free(zz,(void **)import_global_ids,
			  ZOLTAN_SPECIAL_MALLOC_GID);
      Zoltan_Special_Free(zz,(void **)import_local_ids,
			  ZOLTAN_SPECIAL_MALLOC_LID);
      ZOLTAN_PRINT_ERROR(zz->Proc,yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    if (!Zoltan_Special_Malloc(zz,(void **)import_to_part,nrectags,
                           ZOLTAN_SPECIAL_MALLOC_INT)) {
      Zoltan_Special_Free(zz,(void **)import_global_ids,
			  ZOLTAN_SPECIAL_MALLOC_GID);
      Zoltan_Special_Free(zz,(void **)import_local_ids,
			  ZOLTAN_SPECIAL_MALLOC_LID);
      Zoltan_Special_Free(zz,(void **)import_procs,ZOLTAN_SPECIAL_MALLOC_INT);
      ZOLTAN_PRINT_ERROR(zz->Proc,yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }

    /* for each region imported, look at its originating processor */
    for(i=0; i<nrectags; i++) {
      ZOLTAN_SET_GID(zz, &((*import_global_ids)[i*num_gid_entries]),
                 import_regs[i].Global_ID);
      ZOLTAN_SET_LID(zz, &((*import_local_ids)[i*num_lid_entries]),
                 import_regs[i].Local_ID);
      (*import_procs)[i]   = import_regs[i].newProc;
      /*(*import_to_part)[i] = zz->Proc;*/
      (*import_to_part)[i] = import_regs[i].newProc;
    }

    return ierr;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
