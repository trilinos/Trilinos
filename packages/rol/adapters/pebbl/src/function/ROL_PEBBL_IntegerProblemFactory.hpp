// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PEBBL_INTEGERPROBLEMFACTORY_H
#define ROL_PEBBL_INTEGERPROBLEMFACTORY_H

#include "ROL_Ptr.hpp"
#include "ROL_PEBBL_IntegerProblem.hpp"

/** @ingroup func_group
    \class ROL::PEBBL::IntegerProblemFactory
    \brief Defines the pebbl IntegerProblemFactory interface.

    ROL's IntegerProblemFactory constructs a new (identical)
    instance of an optimization problem for use in pebbl.

    ---
*/

namespace ROL {
namespace PEBBL {

template <class Real>
class IntegerProblemFactory {
public:
  virtual ~IntegerProblemFactory(void) {}

  virtual Ptr<IntegerProblem<Real>> build(void) = 0;

#ifdef HAVE_MPI
  virtual void setCommunicator(MPI_Comm comm = MPI_COMM_WORLD) {};
#endif

}; // class IntegerProblemFactory

} // namespace PEBBL
} // namespace ROL

#endif
