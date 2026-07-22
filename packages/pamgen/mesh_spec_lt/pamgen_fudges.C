// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "pamgen_code_types.h"

namespace PAMGEN_NEVADA {

extern const Real ROUND_OFF_FACTOR  = 1.0e-12;
extern const Real TINY = REAL_MIN;
extern const Real a_lot_positive     = 1.0e-6;
extern const Real a_little_positive  = (100.0 * ROUND_OFF_FACTOR);
extern const Real a_little_negative  = -a_little_positive;

} // end namespace PAMGEN_NEVADA {

