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

extern void add_to_bucket(LIST_ENTRY **, VERTEX *);
extern void remove_from_bucket(LIST_ENTRY **, VERTEX *);
extern VERTEX *search_bucket(LIST_ENTRY *, ID *);
extern void free_bucket(LIST_ENTRY **);

#endif
