/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/
#ifndef lint
static char *cvs_gr_list_const_h = "$Id$";
#endif

/*****************************************************************************/

#ifndef __LIST_CONST_H
#define __LIST_CONST_H

#include "gr_hash_const.h"

extern void LB_add_to_bucket(LIST_ENTRY **, VERTEX *);
extern void LB_remove_from_bucket(LIST_ENTRY **, VERTEX *);
extern VERTEX *LB_search_bucket(LIST_ENTRY *, ID *);
extern void LB_free_bucket(LIST_ENTRY **);

#endif
