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

#ifndef __ID_UTIL_CONST_H
#define __ID_UTIL_CONST_H

#ifndef lint
static char *cvs_id_util_const_h = "$Id$";
#endif

struct ID_Util_Struct {
  ID_INT_FN  *Compare;
  ID_VOID_FN *Assign;
  ID_NEW_FN  *New_ID;
  ID_PRINT_FN *Print_ID;
};

extern ID_UTIL BL_ID_Util;

#endif
