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

#ifndef __ID_CONST_H
#define __ID_CONST_H

#ifndef lint
static char *cvs_id_const_h = "$Id$";
#endif

struct ID_Struct {
  int Proc;
  int Number;
};

extern ID_INT_FN LB_compare_combo_id;
extern ID_VOID_FN LB_assign_combo_id;
extern ID_NEW_FN LB_new_combo_id;
extern ID_PRINT_FN LB_print_combo_id;

#endif
