#include "migreg.h"
#include "hilbert_const.h"
#include "comm.h"

typedef double Coord[3];

/*
 * migreg_migrate_regions(Message *message_Array, int number_of_regions)
 *
 * migrate regions to the processor that have the octant owning their centroid.
 * receive regions from other processors and try to insert them into one of
 * the local subtree.
 */

void migreg_migrate_regions(Region *regions, int *npids, 
			    int nregions, int *c2) {
  int i;                         /* index counter */
  int nreceives,                 /* number of messages received */
      from;                      /* from whom the message originated */
  Region reg;                    /* region to be inserted into the subtree */
  int dest;                      /* where the message is being sent */
  int n_import;
  COMM_OBJ *comm_plan;           /* Communication object returned by 
				    Bruce and Steve's communication routines */
  Region *import_objs;          /* Array of import objects used to request 
				    the objs from other processors. */

  comm_plan = comm_create(nregions, npids, MPI_COMM_WORLD, &n_import);
  *c2 = n_import;
  import_objs = (Region *)malloc(n_import * sizeof(Region));

  if((n_import != 0) && (import_objs == NULL)) {
    fprintf(stderr,"ERROR in migreg_migrate_regions: %s\n",
	    "cannot allocate memory for import_objs.");
    abort();
  }

  comm_do(comm_plan, (char *) regions, sizeof(Region), 
          (char *) import_objs);

  for (i=0; i<n_import; i++) {
    insert_orphan(import_objs[i]);
  }

  free(import_objs);
  comm_destroy(comm_plan);
}

/*
 * void insert_orphan(pRegion region)
 *
 * Insert orphan regions migrated from off processors, or to insert
 * regions that lie on the boundary.
 */
void insert_orphan(Region reg) {
  pRList rootlist;                 /* list of all local roots */
  int rflag;                       /* flag to indicate region fits in octant */
  int i, j;                        /* index counters */
  double upper,                    /* upper bounds of the octant */
         lower;                    /* lower bounds of the octant */

  rootlist = POC_localroots();                 /* get a list all local roots */
  if(rootlist == NULL)                                        /* error check */
    fprintf(stderr,"ERROR in insert_orphans(), rootlist is NULL\n");

  if (dimension == 2)
    i = 2;                                           /* ignore z coordinates */
  else
    i = 3;

  rflag = 0;
  while(rootlist != NULL) {
    rflag = 1;
    /* fprintf(stderr, "(%d) attempting to insert %lf %lf %lf in octant %d\n",
     *       msg_mypid, reg.Coord[0], reg.Coord[1], reg.Coord[2],
     *       POC_id(rootlist->oct));
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
	    "ERROR: migreg_migrate_regions could not insert region\n");
    abort();
  }
  else                                     /* found a place to insert region */
    oct_subtree_insert(rootlist->oct, &reg);
}


/*
 * void migreg_migrate_orphans(pRegion RegionList, int number_of_regions,
 *                             int level_of_refinement, Map *map_array)
 *
 * Migrate regions that are not attached to an octant to the processor owning
 * the octant where their centroid is.
 */
void migreg_migrate_orphans(pRegion RegionList, int nregions, int level,
			    Map *array, int *c1, int *c2) {
  int     i, j, k;                    /* index counters */
  pRegion ptr;                        /* region in the mesh */
  double  origin[3];                  /* centroid coordinate information */
  pRegion *regions;                   /* an array of regions */
  int     *npids;
  Region  *regions2;                   /* an array of regions */
  int     *npids2;
  int     nreg;                       /* number of regions */
  double  min[3],                     /* minimum bounds of an octant */
          max[3];                     /* maximum bounds of an octant */
  double  cmin[3],                    /* minimum bounds of a child octant */
          cmax[3];                    /* maximum bounds of a child octant */
  int     new_num;
  int     n;

  /* create the array of messages to be sent to other processors */
  /* Array = (Message *)malloc(nregions * sizeof(Message)); */

  regions = (pRegion *)malloc(nregions * sizeof(pRegion));
  npids = (int *)malloc(nregions * sizeof(int));

  ptr = RegionList;
  n = nreg = 0;
  while(ptr != NULL && nregions > 0) {
    if(ptr->attached == 1) {
      /* if region already attached to an octant, then skip to next region */
      ptr = ptr->next;
      continue;
    }

    /* region not attached, have to find which processor to send to */
    j=0;
    vector_set(min, gmin);
    vector_set(max, gmax);
    /* 
     * for each level of refinement, find which child region belongs to.
     * translate which child to which entry in map array.
     */
    for(i=0; i<level; i++) {
      bounds_to_origin(min, max, origin);
      if(dimension == 2)
	j = j * 4;
      else
	j = j * 8;
      k = child_which(origin, ptr->Coord);
      if(HILBERT) {
	new_num = change_to_hilbert(min, max, origin, k);
      }
      else if(GRAY)
	new_num = convert_to_gray(k);
      else
	new_num = k;
      j += new_num;
      child_bounds(min, max, origin, k, cmin, cmax);
      vector_set(min, cmin);
      vector_set(max, cmax);
    }
    
    /* inform message which processor to send to */
    if(array[j].npid != -1) {
      npids[n] = array[j].npid;
      copy_info(ptr, &(regions[n++]));
    }
    else {
      insert_orphan(*ptr);
    }    
    /* copy region info to the message array */
    /* copy_info(&(Array[nreg].region), ptr); */
  
    nreg++;                                      /* increment region counter */
    ptr = ptr->next;                                  /* look at next region */
  }

  /*
   * if regions looked at != number of regions in region list, 
   * then there is an error
   */
  if (nreg!=nregions) {
    fprintf(stderr,"%d migreg_migrate_orphans: "
	    "%d regions found != %d expected\n",
	    msg_mypid,nreg,nregions);
    abort();
  }

  regions2 = (Region *)malloc(n * sizeof(Region));
  npids2 = (int *)malloc(n * sizeof(int));
  
  /* fprintf(stderr,"(%d) n = %d\n", msg_mypid, n); */
  for(i=0; i<n; i++) {
    npids2[i] = npids[i];
    vector_set(regions2[i].Coord, regions[i]->Coord);
    regions2[i].Weight = regions[i]->Weight;
    regions2[i].Tag.Global_ID = regions[i]->Tag.Global_ID;
    regions2[i].Tag.Local_ID = regions[i]->Tag.Local_ID;
    regions2[i].Tag.Proc = regions[i]->Tag.Proc;
    regions2[i].attached = 0;
  }

  *c1 = n;
  /* migrate the orphan regions according to the message array */
  migreg_migrate_regions(regions2, npids2, n, c2);
  
  free(regions);
  free(npids);
  free(regions2);
  free(npids2);
}

/*
 * int copy_info(pRegion *destination, pRegion source)
 *
 * Copies region information from the source to the destination
 */
void copy_info(pRegion src, pRegion *dest) {
  pRegion copy;

  /* mallloc space for destination */
  copy = (pRegion)malloc(sizeof(Region));
  if(copy == NULL) {
    fprintf(stderr,"ERROR in copy_info, cannot allocate memory\n");
    abort();
  }
  
  /* set up return pointer */
  *dest = copy;

  /* copy all important information */
  vector_set(copy->Coord, src->Coord);
  copy->Weight = src->Weight;
  copy->Tag.Global_ID = src->Tag.Global_ID;
  copy->Tag.Local_ID = src->Tag.Local_ID;
  copy->Tag.Proc = src->Tag.Proc;
  copy->attached = 0;
}
