/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "phg.h"


char *uMe(PHGComm *hgc)
{
    static char msg[1024];

    sprintf(msg, "<%d/%d>: (%d,%d)/[%d,%d] ->", hgc->myProc, hgc->nProc, hgc->myProc_x, hgc->myProc_y, hgc->nProc_x, hgc->nProc_y);
    return msg;
}

void uprintf(PHGComm *hgc, char *f_str,...)
{
va_list argp;

fflush(stdout);
fflush(stderr);
printf("%s", uMe(hgc)); 
va_start(argp, f_str);
vfprintf(stdout, f_str, argp);
va_end(argp);
fflush(stdout);
}

/*************************************************************************
* -------------------------- Error Exit ----------------------------------
**************************************************************************/
void errexit(char *f_str,...)
{
va_list argp;

fflush(stdout);
fflush(stderr);
fprintf(stderr, "\n****** Error:\n");
va_start(argp, f_str);
vfprintf(stderr, f_str, argp);
va_end(argp);

fprintf(stderr," ******\n");
fflush(stderr);
exit(1);
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
