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

#include "dr_const.h"

/* Function prototypes */
extern void initialize_element(ELEM_INFO *elem);
extern void free_mesh_arrays(MESH_INFO_PTR mesh);
extern void free_element_arrays(ELEM_INFO *elem, MESH_INFO_PTR mesh);

#endif
