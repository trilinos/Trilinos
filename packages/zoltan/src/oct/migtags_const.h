/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __MIGTAGS_CONST_H
#define __MIGTAGS_CONST_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


extern int Zoltan_Oct_migrate_objects(ZZ *zz, pOctant *octs, int *newpids, int nocts,
			       int *nsentags, 
			       pRegion *import_tags, int *nrectags, 
			       float *c2, float *c3, int *counter3, 
			       int *counter4);

extern int Zoltan_Oct_fix_tags(ZZ *zz, ZOLTAN_ID_PTR *import_global_ids, ZOLTAN_ID_PTR *import_local_ids,
                        int **import_procs, int nrectags, pRegion import_regs);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif /* __MIGTAGS_CONST_H */
