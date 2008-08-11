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


#include "migreg.h"
#include "zoltan_comm.h"
#include "dfs_const.h"
#include "oct_util_const.h"
#include "octree_const.h"

#define MIGMIGREGCommCreate 32767
#define MIGMIGREGCommDo     32766

static int Zoltan_Oct_insert_orphan(ZZ *zz, Region reg);
static int Zoltan_Oct_copy_info(ZZ *zz, pRegion src, pRegion *dest);
/*
 * Zoltan_Oct_migreg_migrate_regions(Message *message_Array, int number_of_regions)
 *
 * migrate regions to the processor that have the octant owning their centroid.
 * receive regions from other processors and try to insert them into one of
 * the local subtree.
 */

static int Zoltan_Oct_migreg_migrate_regions(ZZ *zz, Region *regions, 
                              ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
                              int *npids, int nregions, int *c2) 
{
  char *yo = "Zoltan_Oct_migreg_migrate_regions";
  int i;                         /* index counter */
  int ierr = ZOLTAN_OK;
  int n_import;
  ZOLTAN_COMM_OBJ *comm_plan;           /* Object returned by communication routines */
  Region *import_objs = NULL;    /* Array of import objects used to request 
				    the objs from other processors. */
  ZOLTAN_ID_PTR import_gids = NULL;  /* Array of global IDs of import_objs. */
  ZOLTAN_ID_PTR import_lids = NULL;  /* Array of local IDs of import_objs. */
  int num_gid_entries = zz->Num_GID;
  int num_lid_entries = zz->Num_LID;

  ierr = Zoltan_Comm_Create(&comm_plan, nregions, npids, zz->Communicator, 
                        MIGMIGREGCommCreate, &n_import);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_Comm_Create");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ierr);
  }
  *c2 = n_import;
  if (n_import > 0) {
    import_objs = (Region *) ZOLTAN_MALLOC(n_import * sizeof(Region));
    import_gids = ZOLTAN_MALLOC_GID_ARRAY(zz, n_import);
    import_lids = ZOLTAN_MALLOC_LID_ARRAY(zz, n_import);

    if (!import_objs || !import_gids || (num_lid_entries && !import_lids)) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
  }
  ierr = Zoltan_Comm_Do(comm_plan, MIGMIGREGCommDo, (char *) regions, sizeof(Region), 
                   (char *) import_objs);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_Comm_Do.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    ZOLTAN_FREE(&import_objs);
    ZOLTAN_FREE(&import_gids);
    ZOLTAN_FREE(&import_lids);
    return (ierr);
  }

  ierr = Zoltan_Comm_Do(comm_plan, MIGMIGREGCommDo-1, (char *) gids, 
                    sizeof(ZOLTAN_ID_TYPE)*num_gid_entries, (char *) import_gids);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_Comm_Do.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    ZOLTAN_FREE(&import_objs);
    ZOLTAN_FREE(&import_gids);
    ZOLTAN_FREE(&import_lids);
    return (ierr);
  }

  if (num_lid_entries > 0) {
    ierr = Zoltan_Comm_Do(comm_plan, MIGMIGREGCommDo-2, (char *) lids, 
                      sizeof(ZOLTAN_ID_TYPE)*num_lid_entries, (char *) import_lids);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_Comm_Do.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      ZOLTAN_FREE(&import_objs);
      ZOLTAN_FREE(&import_gids);
      ZOLTAN_FREE(&import_lids);
      return (ierr);
    }
  }
  for (i=0; i<n_import; i++) {
    import_objs[i].Global_ID = &(import_gids[i*num_gid_entries]);
    import_objs[i].Local_ID = (num_lid_entries 
                                 ? &(import_lids[i*num_lid_entries]) 
                                 : NULL);
    Zoltan_Oct_insert_orphan(zz, import_objs[i]);
  }
  
  ZOLTAN_FREE(&import_objs);
  ZOLTAN_FREE(&import_gids);
  ZOLTAN_FREE(&import_lids);

  ierr = Zoltan_Comm_Destroy(&comm_plan);
  if(ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ierr);
  }
  return ierr;
}

/*
 * void Zoltan_Oct_insert_orphan(pRegion region)
 *
 * Insert orphan regions migrated from off processors, or to insert
 * regions that lie on the boundary.
 */
static int Zoltan_Oct_insert_orphan(ZZ *zz, Region reg) {
  pRList  RootList;                /* list of all local roots */
  pOctant RootOct;
  int rflag;                       /* flag to indicate region fits in octant */
  int i, j;                        /* index counters */
  double upper,                    /* upper bounds of the octant */
         lower;                    /* lower bounds of the octant */
  OCT_Global_Info *OCT_info = (OCT_Global_Info *)(zz->LB.Data_Structure);
  char *yo = "Zoltan_Oct_insert_orphan";
  int ierr = ZOLTAN_OK;

  if (OCT_info->OCT_dimension == 2)
    i = 2;                                           /* ignore z coordinates */
  else
    i = 3;

  rflag = 0;
  RootList = Zoltan_Oct_POct_localroots(OCT_info);               /* get a list all local roots */
  while((RootOct = RL_nextRootOctant(&RootList))) {
    rflag = 1;
    for (j=0; j<i; j++) {
      lower = RootOct->min[j];
      upper = RootOct->max[j];
      if (reg.Coord[j]<lower || reg.Coord[j]>upper) {
	/* if region coord lie outside bounds, then cannot fit */
	rflag = 0;
	break;
      }
    }
    if(rflag == 1) { 
      /* region fits inside octant */
      /* found a place to insert region */
      Zoltan_Oct_subtree_insert(zz, RootOct, &reg);
      return ierr;
    }
  }
  ierr = ZOLTAN_WARN;
  ZOLTAN_TRACE_DETAIL(zz, yo, "could not insert region");
  
  fprintf(stderr,"%s failed to insert %f %f %f on proc %d\n",yo, reg.Coord[0], reg.Coord[1], reg.Coord[2], zz->Proc);

  RootList = Zoltan_Oct_POct_localroots(OCT_info); 
  RL_printRootOctants(RootList);
  return ierr;
}

int Zoltan_Oct_migreg_migrate_orphans(ZZ *zz, pRegion RegionList, int nregions,
                               int level, Map *array, int *c1, int *c2) {
  int     i, j, k;                    /* index counters */
  pRegion ptr;                        /* region in the mesh */
  COORD   origin;                     /* centroid coordinate information */
  pRegion *regions = NULL;            /* an array of regions */
  int     *npids = NULL;
  Region  *regions2 = NULL;           /* an array of regions */
  int     *npids2 = NULL;
  int     nreg;                       /* number of regions */
  COORD   min,                        /* minimum bounds of an octant */
          max;                        /* maximum bounds of an octant */
  COORD   cmin,                       /* minimum bounds of a child octant */
          cmax;                       /* maximum bounds of a child octant */
  COORD   rmin,                       /* minimum bounds of a remote octant */
          rmax;                       /* maximum bounds of a remote octant */
  int     new_num;
  int     n;
  int     dir = 0;
  pRList  RootList;              
  pOctant RootOct;
  OCT_Global_Info *OCT_info = (OCT_Global_Info *)(zz->LB.Data_Structure);
  char *yo = "Zoltan_Oct_migreg_migrate_orphans_static";
  int ierr = ZOLTAN_OK;
  ZOLTAN_ID_PTR gids2, lids2;
  int num_gid_entries = zz->Num_GID;
  int num_lid_entries = zz->Num_LID;

  if(nregions > 0) {
    /* create the array of messages to be sent to other processors */
    /* Array = (Message *) ZOLTAN_MALLOC(nregions * sizeof(Message)); */
    
    if((regions = (pRegion *) ZOLTAN_MALLOC(nregions * sizeof(pRegion))) == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    if((npids = (int *) ZOLTAN_MALLOC(nregions * sizeof(int))) == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      ZOLTAN_FREE(&regions);
      return ZOLTAN_MEMERR;
    }
  }
  ptr = RegionList;
  n = nreg = 0;
  while((ptr != NULL) && (nregions > 0)) {
    if(ptr->attached == 1) {
      /* if region already attached to an octant, then skip to next region */
      ptr = ptr->next;
      continue;
    }

    /* region not attached, have to find which processor to send to */
    j=0;
    dir = 0;
    vector_set(min, OCT_info->OCT_gmin);
    vector_set(max, OCT_info->OCT_gmax);
    /* 
     * for each level of refinement, find which child region belongs to.
     * translate which child to which entry in map array.
     */
    for(i=0; i<level; i++) {
      Zoltan_Oct_bounds_to_origin(min, max, origin);
      if(OCT_info->OCT_dimension == 2)
	j = j * 4;
      else
	j = j * 8;
      k = Zoltan_Oct_child_which(OCT_info,origin, ptr->Coord);
      new_num = Zoltan_Oct_convert_idx_from_map(OCT_info, dir, k);
      dir = Zoltan_Oct_get_child_dir(OCT_info, dir, new_num);
      j += new_num;
      Zoltan_Oct_child_bounds(min, max, origin, k, cmin, cmax);
      vector_set(min, cmin);
      vector_set(max, cmax);
    }
    /* inform message which processor to send to */
    npids[n] = array[j].npid;
    RootList = array[j].list;
    while((RootOct = RL_nextRootOctant(&RootList))) {
      Zoltan_Oct_bounds(RootOct,rmin,rmax);
      if (Zoltan_Oct_in_box_closure(OCT_info, ptr->Coord ,rmin, rmax)) {
	npids[n] = RootOct->npid;
	break;
      }
    }
    if((npids[n] != -1) && (npids[n] != zz->Proc)) {
      Zoltan_Oct_copy_info(zz, ptr, &(regions[n++]));
    }
    else {
      Zoltan_Oct_insert_orphan(zz, *ptr);
    }
    nreg++;                                      /* increment region counter */
    ptr = ptr->next;                                  /* look at next region */
  }

  /*
   * if regions looked at != number of regions in region list, 
   * then there is an error
   */
  if (nreg!=nregions) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "regions found != to expected number of regions");
    return ZOLTAN_FATAL;
  }

  regions2 = (Region *) ZOLTAN_MALLOC(n * sizeof(Region));
  gids2 = ZOLTAN_MALLOC_GID_ARRAY(zz, n);
  lids2 = ZOLTAN_MALLOC_LID_ARRAY(zz, n);
  npids2 = (int *) ZOLTAN_MALLOC(n * sizeof(int));
  
  for(i=0; i<n; i++) {
    npids2[i] = npids[i];
    vector_set(regions2[i].Coord, regions[i]->Coord);
    regions2[i].Weight = regions[i]->Weight;
    regions2[i].Global_ID = &(gids2[i*num_gid_entries]);
    regions2[i].Local_ID = (num_lid_entries 
                              ? &(lids2[i*num_lid_entries]) 
                              : NULL);
    ZOLTAN_SET_GID(zz, &(gids2[i*num_gid_entries]), regions[i]->Global_ID);
    ZOLTAN_SET_LID(zz, &(lids2[i*num_lid_entries]), regions[i]->Local_ID);
    regions2[i].Proc = regions[i]->Proc;
    regions2[i].attached = 0;
  }

  *c1 = n;
  /* migrate the orphan regions according to the message array */
  Zoltan_Oct_migreg_migrate_regions(zz, regions2, gids2, lids2, npids2, n, c2);
  
  for (i=0; i < n; i++) {
    ZOLTAN_FREE(&(regions[i]->Global_ID));
    ZOLTAN_FREE(&(regions[i]->Local_ID));
    ZOLTAN_FREE(&(regions[i]));
  }
  ZOLTAN_FREE(&regions);
  ZOLTAN_FREE(&npids);
  ZOLTAN_FREE(&regions2);
  ZOLTAN_FREE(&gids2);
  ZOLTAN_FREE(&lids2);
  ZOLTAN_FREE(&npids2);

  return ierr;
}

/*
 * int Zoltan_Oct_copy_info(pRegion *destination, pRegion source)
 *
 * Copies region information from the source to the destination
 */
static int Zoltan_Oct_copy_info(ZZ *zz, pRegion src, pRegion *dest) {
  pRegion copy;
  char *yo = "Zoltan_Oct_copy_info";
  int ierr = ZOLTAN_OK;

  /* mallloc space for destination */
  copy = (pRegion) ZOLTAN_MALLOC(sizeof(Region));
  if(copy == NULL) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_MEMERR;
  }
  copy->Global_ID = ZOLTAN_MALLOC_GID(zz);
  copy->Local_ID  = ZOLTAN_MALLOC_LID(zz);
  if (copy->Global_ID == NULL || (zz->Num_LID && copy->Local_ID == NULL)) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_MEMERR;
  }
  
  /* set up return pointer */
  *dest = copy;

  /* copy all important information */
  vector_set(copy->Coord, src->Coord);
  copy->Weight = src->Weight;
  ZOLTAN_SET_GID(zz, copy->Global_ID, src->Global_ID);
  ZOLTAN_SET_LID(zz, copy->Local_ID, src->Local_ID);
  copy->Proc = src->Proc;
  copy->attached = 0;
  return ierr;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
