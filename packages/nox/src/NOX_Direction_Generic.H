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

#ifndef NOX_DIRECTION_GENERIC_H
#define NOX_DIRECTION_GENERIC_H

#include "Teuchos_RCP.hpp"

// Forward declarations
namespace Teuchos {
  class ParameterList;
}
namespace NOX {
  class GlobalData;
  namespace Abstract {
    class Vector;
    class Group;
  }
  namespace Solver {
    class Generic;
    class LineSearchBased;
  }
}

namespace NOX {

//! Search direction strategies
/*!
  The Generic class spells out the interface. The Manager class can
  instantiate any other concrete implementation. The remaining classes
  are concrete implementations.
*/
namespace Direction {

//! %Generic direction interface
/*!
  %Generic interface for calculating a search direction, \f$d\f$, to
  be used in updating the iterate.
 */
class Generic {

public:

  //! Constructor.
  /*!
    Constructors of derived objects should look like reset().
   */
  Generic() {};

  //! Destructor
  virtual ~Generic() {};

  //! Reset direction based on possibly new parameters.
  virtual bool reset(const Teuchos::RCP<NOX::GlobalData>& gd,
             Teuchos::ParameterList& params) = 0;

  /*!
    \brief Compute the direction vector, \c dir, for a specific method given
    the current group, \c grp.

    The \c grp is not const so that we can compute the F vector, the
    Jacobian matrix, the %Newton vector, and so on.

    Const access to the solver is used for getting additional
    information such as the past solution, the iteration number, and
    so on.
  */
  virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
               const NOX::Solver::Generic& solver) = 0;

  /*!

    \brief Same as compute(NOX::Abstract::Vector&,
    NOX::Abstract::Group&, const NOX::Solver::Generic&)

    Enables direct support for line search based solvers for the
    purpose of efficiency since the LineSearchBased object has a
    getStep() function that some directions require.

    If it is not redefined in the derived class, it will just call the
    compute with the NOX::Solver::Generic argument.
  */
  virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
               const NOX::Solver::LineSearchBased& solver);
};
} // namespace Direction
} // namespace NOX

#endif
