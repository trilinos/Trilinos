#include "migreg.h"

/*
 * migreg_migrate_regions()
 *
 */
void migreg_migrate_regions(Message *Array, int nregions) {
  int i,j;
  pOctant poct;
  pRegion region;
  int dest;
  int num_sends;             /* num of procs to send to               */
  int *dest_pids;            /* which pids to send to                 */
  int num_recvs;             /* num of procs we received from         */
  long *num_recv_regions;    /* how many regions received from each   */
  int *src_pids;             /* which pids received from              */
  double centroid[3];
  int nreceives, from;
  Region reg;
  Message recmsg;

  num_sends=0;

  msg_send_init();
  for (i=0; i<nregions; i++) {
    dest = Array[i].npid;
    
    if (dest==msg_mypid) {
      fprintf(stderr,"%s %d %s\n", 
	      "migreg_migrate_regions: tried to send region",
	      "to self.", Array[i].region->Tag.Local_ID);
      abort();
    }

#if 0
    fprintf(stderr, "%d sending coord:(%lf, %lf, %lf) to %d\n",
	    msg_mypid, Array[i].region->Coord[0], 
	    Array[i].region->Coord[1], Array[i].region->Coord[2],
	    Array[i].npid);
#endif

    /* send the regions from here */
    msg_bsend(Array[i].region, sizeof(Region), 
	      Array[i].npid, MTYPE_DMIGRATE);
  }

  nreceives=msg_nreceives();
  for(i=0; i<nreceives; i++) {
    msg_breceive(&(reg), sizeof(Region), &from, MTYPE_DMIGRATE);

#if 0
    fprintf(stderr, "(%d) %d coord:(%lf, %lf, %lf) from %d\n",
	    msg_mypid, reg.Tag.Global_ID,
	    reg.Coord[0], reg.Coord[1], reg.Coord[2], from);
#endif

    if (!oct_global_insert(&reg)) {
      abort();
    }
  }
  fprintf(stderr,"finished migrating orphans...\n");
}

/*
 * void migreg_migrate_orphans(pMeshPB pmeshpb, pMesh mesh, int nregions)
 *
 * Migrate regions without octant pointers to the processor owning
 * the octant where their centroid is.
 *
 */
void migreg_migrate_orphans(pRegion RegionList, int nregions, int level,
			    Map *array) {
  int     i, j;                   /* index counters */
  pRegion ptr;                    /* region in the mesh */
  double  origin[3];              /* centroid coordinate information */
  pRegion *regions = NULL;        /* an array of regions */
  int     nreg;                   /* number of regions */
  int     sender;                 /* who sent the region */
  pOctant root;                   /* root of a subtree */
  ROOTMSG rmsg;                   /* message sent about the root */
  ROOT    *roots = NULL;          /* root information */
  int     parent;                 /* parent info */
  Message *Array;
  double  min[3], max[3];
  double  cmin[3], cmax[3];
  int k;

  Array = (Message *)malloc(nregions * sizeof(Message));
  ptr = RegionList;

  nreg = 0;
  while(ptr != NULL && nregions > 0) {
    if(ptr->attached == 1) {
      ptr = ptr->next;
      continue;
    }

#if 0
    ptr2 = NULL;
    bounds_to_origin(map->min, map->max, origin);
    i = child_which(origin, ptr->Coord);
    ptr2 = POC_child(map, i);
    while(!POC_isTerminal(ptr2)) {
      bounds_to_origin(ptr2->min, ptr2->max, origin);
      i = child_which(origin, ptr->Coord);
      ptr2 = POC_child(ptr2, i);
    }
#endif

    j=0;
    vector_set(min, gmin);
    vector_set(max, gmax);
    for(i=0; i<level; i++) {
      bounds_to_origin(min, max, origin);
      j = j * 8;
      k = child_which(origin, ptr->Coord);
      j += k;
      child_bounds(min, max, origin, k, cmin, cmax);
      vector_set(min, cmin);
      vector_set(max, cmax);
    }
    
    Array[nreg].npid = array[j].npid;
    copy_info(&(Array[nreg].region), ptr);

    nreg++;
    ptr = ptr->next;
  }

  /* free(roots); */
  
  if (nreg!=nregions) {
    fprintf(stderr,"%d migreg_migrate_orphans: "
	    "%d regions found != %d expected\n",
	    msg_mypid,nreg,nregions);
    abort();
  }

#if 0
  PRINT_IN_ORDER()
    printf("Searching... regions: %5d  roots %5d\n",nreg,nroots);
#endif

  migreg_migrate_regions(Array, nreg);
  /* free(regions); */
  
#if 0
  /* Check that it worked right */
  temp=NULL;
  while (region=M_nextRegion(mesh,&temp)) {
    elt_centroid(region,centroid);
    
    if (!oct_global_find(centroid)) {
      fprintf(stderr,"%d migreg_migrate_orphans: region outside of"
	      "local octree\n",msg_mypid);
      abort();
    }
  }
#endif
}


int copy_info(pRegion *dest, pRegion src) {
  (*dest) = (pRegion)malloc(sizeof(Region));
  vector_set((*dest)->Coord, src->Coord);
  (*dest)->Weight = src->Weight;
  (*dest)->Tag.Global_ID = src->Tag.Global_ID;
  (*dest)->Tag.Local_ID = src->Tag.Local_ID;
  (*dest)->Tag.Proc = src->Tag.Proc;
  (*dest)->attached = 0;
}
