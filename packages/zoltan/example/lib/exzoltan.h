/*
** $Id$
**
** Functions to support writing simple Zoltan examples 
*/

#include "zoltan.h"

#ifdef __cplusplus
extern "C" {
#endif

void exSetDivisions(int div);
int exInitializePoints(float **retPts, int **retIds, int rank, int size);
void exShowIds(int *ids, int n);
void exShowPoints(float *pts, int n, int *ids);

void exShowError(int val, char *s, int me);
int exGlobalSuccess(int rc, int nprocs, int me, int verbose);
int exPrintGlobalResult(char *s, int nprocs, int me,
  int begin, int import, int export, int change);

int exGetNumberOfAssignedObjects(void *userDefinedData, int *err);
int exGetObjectSize(void *userDefinedData, int *err);
void exGetObjectList(void *userDefinedData, int numGlobalIds, int numLids,
  ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int wgt_dim, float *obj_wgts,
  int *err);
void exGetObject(void *userDefinedData, int numGlobalIds, int numLids, int numObjs,
  ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int numDim, double *pts, int *err);

#ifdef __cplusplus
}
#endif





