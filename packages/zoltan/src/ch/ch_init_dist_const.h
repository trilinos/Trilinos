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

#ifndef CH_INIT_DIST_CONST_H
#define CH_INIT_DIST_CONST_H

#ifndef lint
static char *cvs_ch_init_dist_const_id = "$Id$";
#endif

#include "dr_const.h"
#include "dr_input_const.h"

/* define the Chaco initial distribution types */
#define INITIAL_FILE   0
#define INITIAL_LINEAR 1
#define INITIAL_CYCLIC 2

extern void ch_dist_init(int, int, PARIO_INFO_PTR);
extern int ch_dist_num_vtx(int);
extern int ch_dist_max_num_vtx();
extern void ch_dist_vtx_list(int*, int*, int);
extern int ch_dist_proc(int);

#endif  /* CH_INIT_DIST_CONST_H */
