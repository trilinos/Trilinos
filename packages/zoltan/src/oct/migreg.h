/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __MIGREG_H
#define __MIGREG_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include "octant_const.h"
#include "octupdate_const.h"
#include "oct_util_const.h"
#include "migreg_const.h"

typedef struct
{
  pRegion region;
  int npid;
} Message;


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
