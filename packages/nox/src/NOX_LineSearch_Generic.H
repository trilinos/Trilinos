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

#ifndef NOX_LINESEARCH_GENERIC_H
#define NOX_LINESEARCH_GENERIC_H

#include "Teuchos_RCP.hpp"

// Forward Declarations
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
  }
}

namespace NOX {

//! NOX linesearches
namespace LineSearch {

/*!  \brief Base class line search interface

  Every line search should respect the following Parameter:

  - "Max Iters" - maximum number of iterations (i.e., RHS computations)
*/
class Generic {

public:

  //! Default constructor
  Generic() {};

  //! Destructor
  virtual ~Generic() {};

  //! Perform a line search
  /*!

  On input:
  @param grp The initial solution vector, \f$x_{\rm old}\f$.
  @param dir A vector of directions to be used in the line search, \f$d\f$.
  @param s The nonlinear solver.

  On output:
  @param step The distance the direction was scaled, \f$ \lambda \f$.
  @param grp The \c grp is updated with a new solution, \f$ x_{\rm new} \f$, resulting from the linesearch.  Normally, for a single direction line search, this is computed as:

    \f[
    x_{\rm new} = x_{\rm old} + \lambda d.
    \f]

  Ideally, \f$ \|F(x_{\rm new})\| < \|F(x_{\rm old})\| \f$ (e.g the final direction is a descent direction).

  Note that the \c dir object is a std::vector. For typical line
  searches as described in the above equation, this vector is of size
  one.  We have used a std::vector to allow for special cases of
  multi-directional line searches such as the Bader/Schnabel
  curvillinear line search.

  Return value is true for a successful line search computation.

   */
  virtual bool compute(NOX::Abstract::Group& grp,
               double& step,
               const NOX::Abstract::Vector& dir,
               const NOX::Solver::Generic& s) = 0;

};
} // namespace LineSearch
} // namespace NOX

#endif
