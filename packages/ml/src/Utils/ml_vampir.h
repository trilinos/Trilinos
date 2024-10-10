/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifndef __MLVAMPIREH__
#define __MLVAMPIREH__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "ml_common.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef ML_VAMPIR
/*vampir header file*/
#include <VT.h>

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

/* variables for different classes */
extern int ml_vt_genhierarchy_class;

/* variables for different states */
extern int /*ml_vt_reitzinger_state,*/
           ml_vt_aggregating_state,
           ml_vt_building_coarse_T_state,
           ml_vt_build_Pe_state,
           ml_vt_smooth_Pe_state,
           ml_vt_make_coarse_A_with_rap_state,
           ml_vt_reitzinger_cleanup_state;

extern void ML_Vampir_Setup(void);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif /*ifdef ML_VAMPIR*/
#endif
