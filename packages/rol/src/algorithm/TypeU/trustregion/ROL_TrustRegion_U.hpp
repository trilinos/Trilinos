// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TRUSTREGION_U_H
#define ROL_TRUSTREGION_U_H

/** \class ROL::TrustRegion_U
    \brief Provides interface for and implements trust-region subproblem solvers.
*/

#include "ROL_Vector.hpp"
#include "ROL_TrustRegionModel_U.hpp"

namespace ROL {

template<typename Real>
class TrustRegion_U {
public:
  virtual ~TrustRegion_U() {}

  virtual void initialize(const Vector<Real> &x, const Vector<Real> &g) {}

  virtual void solve(Vector<Real>             &s,          // Step (to be computed)
                     Real                     &snorm,      // Step norm (to be computed)
                     Real                     &pRed,       // Predicted reduction (to be computed)
                     int                      &iflag,      // Exit flag (to be computed)
                     int                      &iter,       // Iteration count (to be computed)
                     const Real                del,        // Trust-region radius
                     TrustRegionModel_U<Real> &model) = 0; // Trust-region model
};

} // namespace ROL

#endif
