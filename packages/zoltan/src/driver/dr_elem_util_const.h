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
 * $Name$
 *====================================================================*/

#ifndef _DR_ELEM_UTIL_CONST_H_
#define _DR_ELEM_UTIL_CONST_H_
#ifndef lint
static char *cvs_elem_util_ch_id = "$Id$";
#endif

#include "dr_const.h"

/* Function prototypes */
extern void initialize_element(ELEM_INFO *elem);
extern void free_element_arrays(ELEM_INFO *elem);

#endif
