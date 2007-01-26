/*
#@HEADER
# ************************************************************************
#
#               ML: A Multilevel Preconditioner Package
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
#
# ************************************************************************
#@HEADER
*/

#ifndef __MLVAMPIREH__
#define __MLVAMPIREH__

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
