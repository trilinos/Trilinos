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

#include "migreg.h"
#include "hilbert_const.h"
#include "comm_const.h"
#include "dfs_const.h"
#include "all_allo_const.h"

/*
 * LB_migreg_migrate_regions(Message *message_Array, int number_of_regions)
 *
 * migrate regions to the processor that have the octant owning their centroid.
 * receive regions from other processors and try to insert them into one of
 * the local subtree.
 */

void LB_migreg_migrate_regions(LB *lb, Region *regions, int *npids, 
			       int nregions, int *c2) {
  int i;                         /* index counter */
  int n_import;
<<<<<<< migreg.c
  int msgtag, msgtag2;
  COMM_OBJ *comm_plan;           /* Object returned by communication routines */
  Region *import_objs;          /* Array of import objects used to request 
=======
  COMM_OBJ *comm_plan;           /* Communication object returned by 
				    Bruce and Steve's communication routines */
  Region *import_objs = NULL;    /* Array of import objects used to request 
>>>>>>> 1.16
				    the objs from other processors. */

  msgtag = 32768;
  LB_Comm_Create(&comm_plan, nregions, npids, lb->Communicator, msgtag, &n_import);
  *c2 = n_import;
  if (n_import > 0) {
    import_objs = (Region *) LB_Array_Alloc(__FILE__, __LINE__, 1, n_import,
                                            sizeof(Region));

    if(import_objs == NULL) {
      fprintf(stderr,"ERROR in LB_migreg_migrate_regions: %s\n",
  	    "cannot allocate memory for import_objs.");
      abort();
    }
  }

  msgtag2 = 32767;
  LB_Comm_Do(comm_plan, msgtag2, (char *) regions, sizeof(Region), 
          (char *) import_objs);

  for (i=0; i<n_import; i++) {
    LB_insert_orphan(lb, import_objs[i]);
  }

  LB_FREE(&import_objs);
  LB_Comm_Destroy(&comm_plan);
}

/*
 * void LB_insert_orphan(pRegion region)
 *
 * Insert orphan regions migrated from off processors, or to insert
 * regions that lie on the boundary.
 */
void LB_insert_orphan(LB *lb, Region reg) {
  pRList rootlist;                 /* list of all local roots */
  int rflag;                       /* flag to indicate region fits in octant */
  int i, j;                        /* index counters */
  double upper,                    /* upper bounds of the octant */
         lower;                    /* lower bounds of the octant */
  OCT_Global_Info *OCT_info = (OCT_Global_Info *)(lb->Data_Structure);

  rootlist = POC_localroots(OCT_info);                 /* get a list all local roots */
  if(rootlist == NULL)                                        /* error check */
    fprintf(stderr,"ERROR in LB_insert_orphans(), rootlist is NULL\n");

  if (OCT_info->OCT_dimension == 2)
    i = 2;                                           /* ignore z coordinates */
  else
    i = 3;

  rflag = 0;
  while(rootlist != NULL) {
    rflag = 1;
    /* if(LB_Proc == 1)
     *   fprintf(stderr, "(%d) %lf %lf %lf in (%lf %lf %lf, %lf %lf %lf)\n",
     *     LB_Proc, reg.Coord[0], reg.Coord[1], reg.Coord[2],
     *     (rootlist->oct)->min[0], (rootlist->oct)->min[1], 
     *     (rootlist->oct)->min[2], (rootlist->oct)->max[0],
     *     (rootlist->oct)->max[1], (rootlist->oct)->max[2]);
     */     
    /* for each of the x,y,z-coordinates check if region fits in the octant */
    for (j=0; j<i; j++) {
      lower = rootlist->oct->min[j];
      upper = rootlist->oct->max[j];
      if (reg.Coord[j]<lower || reg.Coord[j]>upper) {
	/* if region coord lie outside bounds, then cannot fit */
	rflag = 0;
	break;
      }
    }
    if(rflag == 1)                              /* region fits inside octant */
      break;
    rootlist = rootlist->next;
  }
  
  if(rootlist == NULL) {                                      /* error check */
    fprintf(stderr,
	    "(%d) ERROR: LB_migreg_migrate_regions could not insert region\n",
	    lb->Proc);
    abort();
  }
  else                                     /* found a place to insert region */
    LB_oct_subtree_insert(lb, rootlist->oct, &reg);
}


/*
 * void LB_migreg_migrate_orphans(pRegion RegionList, int number_of_regions,
 *                                int level_of_refinement, Map *map_array)
 *
 * Migrate regions that are not attached to an octant to the processor owning
 * the octant where their centroid is.
 */
void LB_migreg_migrate_orphans(LB *lb, pRegion RegionList, int nregions,
                               int level, Map *array, int *c1, int *c2) {
  int     i, j, k;                    /* index counters */
  pRegion ptr;                        /* region in the mesh */
  COORD   origin;                     /* centroid coordinate information */
  pRegion *regions;                   /* an array of regions */
  int     *npids;
  Region  *regions2;                  /* an array of regions */
  int     *npids2;
  int     nreg;                       /* number of regions */
  COORD   min,                        /* minimum bounds of an octant */
          max;                        /* maximum bounds of an octant */
  COORD   cmin,                       /* minimum bounds of a child octant */
          cmax;                       /* maximum bounds of a child octant */
  int     new_num;
  int     n;
  OCT_Global_Info *OCT_info = (OCT_Global_Info *)(lb->Data_Structure);

  /* create the array of messages to be sent to other processors */
  /* Array = (Message *) LB_Array_Alloc(__FILE__, __LINE__, 1, nregions,
                                        sizeof(Message)); */

  regions = (pRegion *) LB_Array_Alloc(__FILE__, __LINE__, 1, nregions,
                                       sizeof(pRegion));
  npids = (int *) LB_Array_Alloc(__FILE__, __LINE__, 1, nregions,
                                 sizeof(int));

  ptr = RegionList;
  n = nreg = 0;
  while((ptr != NULL) && (nregions > 0)) {
    if(ptr->attached == 1) {
      /* if region already attached to an octant, then skip to next region */
      ptr = ptr->next;
      continue;
    }

    /*
     * fprintf(stderr,"(%d) %lf %lf %lf\n", LB_Proc,
     *	    ptr->Coord[0], ptr->Coord[1], ptr->Coord[2]);
     */
    /* region not attached, have to find which processor to send to */
    j=0;
    vector_set(min, OCT_info->OCT_gmin);
    vector_set(max, OCT_info->OCT_gmax);
    /* 
     * for each level of refinement, find which child region belongs to.
     * translate which child to which entry in map array.
     */
    for(i=0; i<level; i++) {
      LB_bounds_to_origin(min, max, origin);
      if(OCT_info->OCT_dimension == 2)
	j = j * 4;
      else
	j = j * 8;
      k = LB_child_which(OCT_info,origin, ptr->Coord);
      if(OCT_info->HILBERT) {
	if(OCT_info->OCT_dimension == 3) 
	  new_num = LB_change_to_hilbert(OCT_info,min, max, origin, k);
	else 
	  new_num = LB_change_to_hilbert2d(OCT_info,min, max, origin, k);
      }
      else if(OCT_info->GRAY)
	new_num = LB_convert_to_gray(k);
      else
	new_num = k;

      j += new_num;
      LB_child_bounds(min, max, origin, k, cmin, cmax);
      vector_set(min, cmin);
      vector_set(max, cmax);
    }
    
    /* inform message which processor to send to */
    if(array[j].npid != -1) {
      npids[n] = array[j].npid;
      LB_copy_info(ptr, &(regions[n++]));
    }
    else {
      LB_insert_orphan(lb, *ptr);
    }    
    /* copy region info to the message array */
    /* LB_copy_info(&(Array[nreg].region), ptr); */
  
    nreg++;                                      /* increment region counter */
    ptr = ptr->next;                                  /* look at next region */
  }

  /*
   * if regions looked at != number of regions in region list, 
   * then there is an error
   */
  if (nreg!=nregions) {
    fprintf(stderr,"%d LB_migreg_migrate_orphans: "
	    "%d regions found != %d expected\n",
	    lb->Proc,nreg,nregions);
    abort();
  }

  regions2 = (Region *) LB_Array_Alloc(__FILE__, __LINE__, 1, n,
                                       sizeof(Region));
  npids2 = (int *) LB_Array_Alloc(__FILE__, __LINE__, 1, n, sizeof(int));
  
  /* fprintf(stderr,"(%d) n = %d\n", LB_Proc, n); */
  for(i=0; i<n; i++) {
    npids2[i] = npids[i];
    vector_set(regions2[i].Coord, regions[i]->Coord);
    regions2[i].Weight = regions[i]->Weight;
    LB_SET_GID(regions2[i].Tag.Global_ID, regions[i]->Tag.Global_ID);
    LB_SET_LID(regions2[i].Tag.Local_ID, regions[i]->Tag.Local_ID);
    regions2[i].Tag.Proc = regions[i]->Tag.Proc;
    regions2[i].attached = 0;
  }

  *c1 = n;
  /* migrate the orphan regions according to the message array */
  LB_migreg_migrate_regions(lb, regions2, npids2, n, c2);
  
  for (i=0; i < n; i++) 
    LB_FREE(&(regions[i]));
  LB_FREE(&regions);
  LB_FREE(&npids);
  LB_FREE(&regions2);
  LB_FREE(&npids2);
}

/*
 * int LB_copy_info(pRegion *destination, pRegion source)
 *
 * Copies region information from the source to the destination
 */
void LB_copy_info(pRegion src, pRegion *dest) {
  pRegion copy;

  /* mallloc space for destination */
  copy = (pRegion) LB_MALLOC(sizeof(Region));
  if(copy == NULL) {
    fprintf(stderr,"ERROR in LB_copy_info, cannot allocate memory\n");
    abort();
  }
  
  /* set up return pointer */
  *dest = copy;

  /* copy all important information */
  vector_set(copy->Coord, src->Coord);
  copy->Weight = src->Weight;
  LB_SET_GID(copy->Tag.Global_ID, src->Tag.Global_ID);
  LB_SET_LID(copy->Tag.Local_ID, src->Tag.Local_ID);
  copy->Tag.Proc = src->Tag.Proc;
  copy->attached = 0;
}
