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

#ifndef NOX_LINESEARCH_UTILS_SLOPE_H
#define NOX_LINESEARCH_UTILS_SLOPE_H

#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "Teuchos_RCP.hpp"

// Forward declaration
namespace NOX {
  class GlobalData;
  namespace Abstract {
    class Vector;
    class Group;
  }
}

namespace NOX {

namespace LineSearch {

namespace Utils {

//! %Common line search utilites for computing the slope of a function.
/*!
  This class provides routines for computing the slope of a give function.  There are two methods, one that uses a Jacobian and the other that estimates the action of the Jacobian by directional derivatives.
*/

class Slope {

public:

  //! Default constructor
  Slope(const Teuchos::RCP<NOX::GlobalData>& gd);

  //! Destructor
  virtual ~Slope();

  //! Reset method.
  void reset(const Teuchos::RCP<NOX::GlobalData>& gd);

  //! Compute the inner product of the given direction and the gradient associated with the given group.
  /*!
    Calculates and returns
    \f[
    \zeta = d^T \nabla f(x).
    \f]

    Here \f$d\f$ represents the input parameter \c dir and \f$\nabla
    f(x)\f$ is the gradient associated with the given group.
  */
  double computeSlope(const NOX::Abstract::Vector& dir,
              const NOX::Abstract::Group& grp);

  //! This is a variant of the computeSlope() method above optimized to work with out having to compute an explicit Jacobian.
  /*!
    Calculates and returns
    \f[
    \zeta = d^T \nabla f(x) = d^TJ^TF
    \f]

    Here \f$d\f$ represents the input parameter \c dir \f$\nabla
    f(x)\f$ is the gradient associated with the given group (for nonlinear solves this equates to \f$ J^TF \f$ where \f$ J \f$ is the Jacobian and \f$ F \f$ is the original nonlinear function).

    We can rewrite this equation as:

    \f[ d^TJ^TF = F^TJd \f]

    which allows us to use directional derivatives to estimate \f$ J^TF \f$:

    \f[ F^TJd = F^T \frac{F(x + \eta d) - F(x)}{\eta} \f]

    This may allow for faster computations of the slope if the Jacobian is expensive to evaluate.

where \f$\eta\f$ is a scalar perturbation calculated by:

\f[ \eta = \lambda * (\lambda + \frac{\| x\|}{\| d\|} ) \f]

\f$ \lambda \f$ is a constant fixed at 1.0e-6.

  */
  double computeSlopeWithOutJac(const NOX::Abstract::Vector& dir,
                const NOX::Abstract::Group& grp);

private:
  //! Disallow default constructor.
  Slope() {};

private:

  //! Printing object
  NOX::Utils utils;

  //! A vector that *may* be used in computeSlope()
  Teuchos::RCP<NOX::Abstract::Vector> vecPtr;

  //! An extra group that will only be allocated if computeSlopeWithOutJac() is called.
  Teuchos::RCP<NOX::Abstract::Group> grpPtr;

};
} // namespace Utils
} // namespace LineSearch
} // namespace NOX

#endif
