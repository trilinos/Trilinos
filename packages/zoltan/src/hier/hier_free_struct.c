

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zz_const.h"
#include "hier.h"

void Zoltan_Hier_Free_Structure(
  ZZ *zz                 /* Zoltan structure */
) {
  /* currently, no persistent data is defined by hierarchical balancing,
     so nothing needs to happen here */
}

int Zoltan_Hier_Copy_Structure(
  ZZ *newzz, ZZ const *oldzz                 /* Zoltan structure */
) {
  /* currently, no persistent data is defined by hierarchical balancing,
     so nothing needs to happen here */
  return ZOLTAN_OK;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
