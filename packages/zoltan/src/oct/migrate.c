
int migrate_regions() {
  msg_send_init();

  for(i=0; i<nreg; i++) {
    rmsg.octid=array[i].id;
    rmsg.region = region;
    msg_bsend(&rmsg,sizeof(REGIONMSG), rmsg[i].npid, MTYPE_REGION);
  }
  
  nroots=msg_nreceives();
  /* replaces: SAFE_MALLOC(roots,ROOT *,sizeof(ROOT)*nroots); */
  if(nroots)
    regs=(REGION *)malloc(nroots * sizeof(REGION)); 
  
  for (i=0; i<nroots; i++) {
    msg_breceive(&rmsg,sizeof(ROOTMSG),&sender,MTYPE_ROOT);
    regs[i].pid=sender;
    regs[i].octid=rmsg.octid;
    roots[i].region=rmsg.region;
  }
  

  /* replaces : SAFE_MALLOC(regions,pRegion *,sizeof(pRegion) * nregions); */
  if(nregions)
    regions=(pRegion *)malloc(sizeof(pRegion) * nregions);
  nreg=0;

  temp=NULL;
  ptr = RegionList;
  while(ptr != NULL) {
    if(ptr->attached == 1)
      continue;
    
    regions[nreg]=region;

    parent= -1;
    for (i=0; i<nroots; i++)
      if (in_box(region->Coord,roots[i].min,roots[i].max) &&
	  (parent<0 || roots[i].size<roots[parent].size))
	parent=i;
    
    if (parent<0) {
      fprintf(stderr,"migreg_migrate_orphans: "
	      "no neighbor contains region %d \n",(int)region->Tag.Local_ID);
      abort();
    }
    
    array[nreg].npid = roots[parent].pid;
    array[nreg].octid = roots[parent].octid;
    array[nreg].region = region;

    /*
     * dataI_set_new(region,"DEST",roots[parent].pid);
     * dataI_set_new(region,"OCTP",(int)roots[parent].ptr);
     */
    nreg++;
    ptr = ptr->next;
  }
  free(roots);
}
