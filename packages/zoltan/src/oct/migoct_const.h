
extern void Migrate_Objects(pOctant *octs, int *newpids, int nocts,
			    pRegion *export_tags, int *nsentags, 
			    pRegion *import_tags, int *nrectags, 
			    float *c2, float *c3, int *counter3, 
			    int *counter4);

void fix_tags(LB_TAG **export_tags, int *nsentags, LB_TAG **import_tags,
	      int *nrectags, pRegion import_regs, pRegion export_regs);

