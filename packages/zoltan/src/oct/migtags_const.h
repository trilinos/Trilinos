#ifndef __MIGTAGS_CONST_H
#define __MIGTAGS_CONST_H

extern int LB_Migrate_Objects(LB *lb, pOctant *octs, int *newpids, int nocts,
			       int *nsentags, 
			       pRegion *import_tags, int *nrectags, 
			       float *c2, float *c3, int *counter3, 
			       int *counter4);

extern int LB_fix_tags(LB *lb, LB_ID_PTR *import_global_ids, LB_ID_PTR *import_local_ids,
                        int **import_procs, int nrectags, pRegion import_regs);

#endif /* __MIGTAGS_CONST_H */
