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

#ifndef NOX_DIRECTION_STEEPESTDESCENT_H
#define NOX_DIRECTION_STEEPESTDESCENT_H

#include "NOX_Direction_Generic.H" // base class

// Forward declaration
namespace NOX {
class Utils;
  class GlobalData;
  namespace MeritFunction {
    class Generic;
  }
}

namespace NOX {
namespace Direction {

//! %Calculates the steepest descent direction.
/*! 

  Calculates the direction
  \f[
  d = - \gamma \nabla f(x) = -\gamma J(x)^T F(x)
  \f]

  This is the (scaled) gradient of the function \f$f(x) = \frac{1}{2}
  \|F(x)\|^2\f$. The valued of \f$\gamma\f$ depends on the choice of
  "Scaling Type" below.

<B>Parameters</B>

"Direction":
<ul>
<li>"Method" = "Steepest Descent"
</ul>

"Direction"/"Steepest Descent":

<ul>
<li>"Scaling Type" can be either of:

    <ul>
    <li>"2-Norm" - \f$ \gamma = \displaystyle\frac{1}{\|d\|_2} \f$ [default]

    <li>"Quadratic Model Min" - \f$ \gamma = \displaystyle\frac{\|d\|_2^2}{d^T J^T J d} \f$

    <li>"F 2-Norm" - \f$ \gamma = \displaystyle\frac{1}{\|F(x)\|_2} \f$

    <li>"None" - No scaling

    </ul>
</ul>

 */
class SteepestDescent : public Generic {

public:

  //! Scaling types
  enum ScalingType {
    //! Scale by the 2-norm of the gradient
    TwoNorm, 
    //! Scale by the step that minimizes the 1-D quadratic
    QuadMin, 
    //! Scale by the 2-norm of F
    FunctionTwoNorm, 
    //! No scaling at all
    None
  };

  //! Constructor
  SteepestDescent(const Teuchos::RCP<NOX::GlobalData>& gd, 
		  Teuchos::ParameterList& params);

  //! Destructor
  ~SteepestDescent();

  // derived
  bool reset(const Teuchos::RCP<NOX::GlobalData>& gd,
	     Teuchos::ParameterList& params);

  // derived
  bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp, 
	       const NOX::Solver::Generic& solver);

  // derived
  virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp, 
		       const NOX::Solver::LineSearchBased& solver);

private:

  //! Print error message and throw error
  void throwError(const std::string& functionName, const std::string& errorMsg);

private:

  //! Global data pointer.  Keep this so any stored parameter list remains valid.
  Teuchos::RCP<NOX::GlobalData> globalDataPtr;

  //! Printing Utils
  Teuchos::RCP<NOX::Utils> utils;

  //! Temporary vector used in Cauchy scaling factor calculation
  Teuchos::RCP<NOX::Abstract::Vector> tmpVecPtr;

  //! Type of scaling of the Cauchy direction
  ScalingType scaleType;

  //! Pointer to a user defined merit function
  Teuchos::RCP<NOX::MeritFunction::Generic> meritFuncPtr;
  
};
} // namespace Direction
} // namespace NOX
#endif
