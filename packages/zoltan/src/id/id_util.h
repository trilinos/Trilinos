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

#ifndef __ID_UTIL_H
#define __ID_UTIL_H

#ifndef lint
static char *cvs_id_util_h = "$Id$";
#endif

#include "id_const.h"

ID_UTIL BL_ID_Util = {compare_combo_id, assign_combo_id, 
                      new_combo_id, print_combo_id};

#endif
