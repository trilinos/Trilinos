// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_ABSTRACT_IMPLICIT_WEIGHTING_H
#define NOX_ABSTRACT_IMPLICIT_WEIGHTING_H

#include "NOX_Common.H"  // for NOX_Config.h

namespace NOX {
namespace Abstract {

/*!
  \brief A pure virtual interface for enabling/disabling any user defined implicit weighting of the concrete object.

  This class is used to toggle on/off any implicit weighting that
  might be handled internal to a particular user defined class.  The
  classical use case is to disable an implicitly defined weighting on
  the inner product and norm computations of the abstract vector class
  so that classes like the NormWRMS can apply their own weighting.
  This prevents a "double weight" effect being applied to specific
  status tests.

  The typical use case is the get the current weighting flag, then
  disable the weighting, perform the computations on the unweighted
  object and then set the flag back to its original value.  By storing
  the original value, one can return it to its prior state regardless
  if an implicit weighting is defined.

  \code
  NOX::ImplicitWeighting& v;

  bool saved_status = v.getImplicitWeighting();
  v.setImplicitWeighting(false);
  .
  <perform unweighted computaiton>
  .
  v.setImplicitWeighting(saved_status);
  \endcode

  CAUTION: This is a power user feature and should only be used in
  concert with specialized NOX::Abstract::Vector implementations.

  \author Roger Pawlowski (SNL 1444)
*/
class ImplicitWeighting {

public:

  virtual ~ImplicitWeighting() {}

  /** \brief Return true if implicit weighting is currently enabled

      CAUTION: This is a power user feature and should only be used in
      concert with specialized NOX::Abstract::Vector implementations.
  */
  virtual bool getImplicitWeighting() const = 0;

  /** \brief Set to true to enable implicit weighting, false disables

      CAUTION: This is a power user feature and should only be used in
      concert with specialized NOX::Abstract::Vector implementations.
  */
  virtual void setImplicitWeighting(bool do_implicit_weighting) = 0;

};

} // namespace Abstract
} // namespace NOX

#endif
