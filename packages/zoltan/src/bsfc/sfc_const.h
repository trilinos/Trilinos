#ifndef _LB_SFC_CONST_H
#define _LB_SFC_CONST_H

/* should change this later on but keylength is okay for now */
#define SFC_KEYLENGTH 3
#define SFC_NO_CUT 0
#define SFC_CUT 1
#define SFC_NOT_BALANCED 1
#define SFC_BALANCED 0
#define SFC_BOUNDING_BOX_EPSILON 0.0000001
#define SFC_FIRST_LEVEL_FLAG 1 /* first level after coarse level */
#define SFC_COARSE_LEVEL_FLAG 2

extern int LB_Set_SFC_Param(char *, char *);

#endif
