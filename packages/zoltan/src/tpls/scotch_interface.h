

#ifndef __SCOTCH_INTERFACE_H
#define __SCOTCH_INTERFACE_H

#include <limits.h>
#include "zoltan_comm.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

int Zoltan_Scotch(
  ZZ *, float *, int *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *,
  int **, int **, int *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *,
  int **, int **);


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
