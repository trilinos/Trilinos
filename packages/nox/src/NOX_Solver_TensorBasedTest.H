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

#ifndef NOX_SOLVER_TENSORBASEDTEST_H
#define NOX_SOLVER_TENSORBASEDTEST_H

#include "NOX_Common.H"

#ifdef WITH_PRERELEASE

#include "NOX_Solver_Generic.H"        // base class
#include "NOX_Direction_Tensor.H"   // class data element
#include "Teuchos_ParameterList.hpp"        // class data element
#include "NOX_LineSearch_Tensor.H"  // class data element

// Forward declarations
namespace NOX {
  class GlobalData;
  class Observer;
}

namespace NOX {
namespace Solver {

/*!
  \brief Linesearch-based solver that uses tensor methods.

  The following parameters are valid for this solver:
<ul>
<li>   "%Direction" - Sublist of the direction parameters, passed to the
       NOX::Direction::Tensor constructor. Defaults to an empty list.

<li> "Linear %Solver" - Sublist of "Direction" and contains the linear
     solver parameters, passed to
     NOX::Abstract::Group::computeNewton() and to
     NOX::Direction::Tensor::compute().  Valid choices used in tensor
     step calculation are listed below.  More valid choices are
     available for the first step when the Newton step is computed by
     computeNewton().

<ul>
<li>    "Max Iterations" - Maximum number of Arnoldi iterations (also
        max Krylov space dimension)
<li>    "Output Frequency" - Print output at every number of iterations
        [default = 20]
<li>    "Preconditioning" - Name of preconditioning type.
<li>    "PreconditioningSide" - Side to apply preconditioning
        (None, Left, Right)
<li>    "Reorthogonalize" - When to reorthogonalize Krylov basis vectors
        (None, As Needed, Always)
<li>    "Tolerance" - Relative tolerance for solving local model
        [default = 1e-4]
</ul>
<br>

<li>   "Line Search" - Sublist of the line search parameters, passed to
       the NOX::LineSearch::Manager constructor. Defaults to an empty list.

    The following parameters can be specified for this line search:

<ul>
<li> "Method" - Name of the line search available to tensor methods
     Valid choices are
   <ul>
   <li> "Standard" - Backtrack along tensor direction unless it is not a
   descent direction, in which case backtrack along Newton direction.
   <li> "Dual" - Backtrack along both the Newton and tensor directions
   and choose the better of the two.
   <li> "Curvilinear" - Backtrack along the "curvilinear" path that
   spans the tensor direction and the Newton direction and that maintains
   monotonicity on the tensor model. [Default]
   </ul>
<li> "Lambda selection" - flag for how to calculate next lambda
     (quadratic, halving)
<li> "Default Step" - starting step parameter (defaults to 1.0)
<li> "Minimum Step" - minimum acceptable step parameter (defaults to 1.0e-12)
<li> "Recovery Step" - step parameter to take when the line search fails
     (defaults to value for "Default Step")
<li> "Max Iters" - maximum number of iterations (i.e., backtracks)
</ul>
<li> "Solver Options" - Sublist of general solver options.
     <ul>
     <li> "User Defined Pre/Post Operator" is supported.  See NOX::Parameter::PrePostOperator for more details.
     </ul>
</ul>

    \author Brett W. Bader (CU-Boulder) bader&nbsp;\@&nbsp;cs.colorado.edu

*/

class TensorBasedTest : public Generic {

public:

  //! %TensorBasedTest-method constructor.
  TensorBasedTest(const Teuchos::RCP<NOX::Abstract::Group>& grp,
          const Teuchos::RCP<NOX::StatusTest::Generic>& tests,
          const Teuchos::RCP<Teuchos::ParameterList>& params);

  //! Destructor
  virtual ~TensorBasedTest();

  virtual void reset(const NOX::Abstract::Vector& initialGuess,
                     const Teuchos::RCP<NOX::StatusTest::Generic>& tests);
  virtual void reset(const NOX::Abstract::Vector& initialGuess);
  virtual void reset();
  virtual NOX::StatusTest::StatusType step();
  virtual NOX::StatusTest::StatusType solve();
  virtual const NOX::Abstract::Group& getSolutionGroup() const;
  virtual const NOX::Abstract::Group& getPreviousSolutionGroup() const;
  virtual NOX::StatusTest::StatusType getStatus() const;
  virtual int getNumIterations() const;
  virtual const Teuchos::ParameterList& getList() const;
  virtual Teuchos::RCP< const NOX::Abstract::Group > getSolutionGroupPtr() const;
  virtual Teuchos::RCP< const NOX::Abstract::Group > getPreviousSolutionGroupPtr() const;
  virtual Teuchos::RCP< const Teuchos::ParameterList > getListPtr() const;
  virtual Teuchos::RCP<const NOX::SolverStats> getSolverStatistics() const;

  //!  Return the direction object
  virtual const NOX::Direction::Tensor& getDirection() const;

protected:

  //! Print out initialization information and calcuation the RHS.
  virtual void init();

  //! Prints the current iteration information.
  virtual void printUpdate();

protected:

  //! Pointer to the global data object.
  Teuchos::RCP<NOX::GlobalData> globalDataPtr;

  //! Utils
  Teuchos::RCP<NOX::Utils> utilsPtr;

  //! Current solution.
  Teuchos::RCP<NOX::Abstract::Group> solnptr;

  //! Previous solution pointer.
  Teuchos::RCP<NOX::Abstract::Group> oldsolnptr;

  //! Current search direction.pointer.
  Teuchos::RCP<NOX::Abstract::Vector> dirptr;

  //! Stopping test.
  Teuchos::RCP<NOX::StatusTest::Generic> testptr;

  //! Type of check to use for status tests.  See NOX::StatusTest for more details.
  NOX::StatusTest::CheckType checkType;

  //! Input parameters.
  Teuchos::RCP<Teuchos::ParameterList> paramsPtr;

  //! Linesearch
  NOX::LineSearch::Tensor lineSearch;

  //! %Search %Direction.
  NOX::Direction::Tensor direction;  // This is Direction::Manager in LSBased

  //! Current step.
  double stepSize;

  //! Number of nonlinear iterations.
  int niter;

  //! %Status of nonlinear solver.
  NOX::StatusTest::StatusType status;

  //! Pointer to a user defined NOX::Observer object.
  Teuchos::RCP<NOX::Observer> observer;

};
} // namespace Solver
} // namespace NOX

#endif

#endif
