#ifndef _LB_SFC_CONST_H
#define _LB_SFC_CONST_H

/* define some constants that are used in multiple files */
#define SFC_KEYLENGTH 3
#define SFC_NO_CUT 0
#define SFC_CUT 1
#define SFC_NOT_BALANCED 1
#define SFC_BALANCED 0
#define SFC_COARSE_LEVEL_FLAG 2

extern int LB_Set_SFC_Param(char *, char *);

#endif
