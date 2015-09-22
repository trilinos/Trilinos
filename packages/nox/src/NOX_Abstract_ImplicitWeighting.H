// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

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
