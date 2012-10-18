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

#ifndef NOX_SOLVER_PSEUDO_TRANSIENT_HPP
#define NOX_SOLVER_PSEUDO_TRANSIENT_HPP

#include "NOX_Solver_Generic.H"	         // base class
#include "NOX_Solver_PrePostOperator.H"  // class data element
#include "Teuchos_ParameterList.hpp"	         // class data element
#include "NOX_Utils.H"		         // class data element
#include "Teuchos_RCP.hpp"       // class data element

// Forward declarations
namespace NOX {
  class GlobalData;
  namespace LineSearch {
    class Generic;
  }
  namespace Direction {
    class Generic;
  }
  namespace Thyra {
    class Group;
  }
}

namespace Thyra {
  template <typename ScalarT> class VectorBase;
}

namespace NOX {
namespace Solver {

/*!
  \brief Pseudo-transient solver.

  Requires that the Thyra::ModelEvaluator implement transient support
  (uses alpha, beta, Jacobian evaluation and optionally x_dot in
  residual).

  Based on the 1998 Kelley Keyes paper, with minor modifications.
*/
class PseudoTransient : public NOX::Solver::Generic {

public:
  
  PseudoTransient(const Teuchos::RCP<NOX::Abstract::Group>& xGrp, 
		  const Teuchos::RCP<NOX::StatusTest::Generic>& tests, 
		  const Teuchos::RCP<Teuchos::ParameterList>& params);
  
  virtual void reset(const NOX::Abstract::Vector& initialGuess, 
		     const Teuchos::RCP<NOX::StatusTest::Generic>& tests);
  virtual void reset(const NOX::Abstract::Vector& initialGuess);
  virtual NOX::StatusTest::StatusType getStatus();
  virtual NOX::StatusTest::StatusType step();
  virtual NOX::StatusTest::StatusType solve();
  virtual const NOX::Abstract::Group& getSolutionGroup() const;
  virtual const NOX::Abstract::Group& getPreviousSolutionGroup() const;
  virtual int getNumIterations() const;
  virtual const Teuchos::ParameterList& getList() const;
  virtual double getStepSize() const;

  inline virtual Teuchos::RCP< const NOX::Abstract::Group > getSolutionGroupPtr() const {return solnPtr;};
  inline virtual Teuchos::RCP< const NOX::Abstract::Group > getPreviousSolutionGroupPtr() const {return oldSolnPtr;};
  inline virtual Teuchos::RCP< const Teuchos::ParameterList > getListPtr() const {return paramsPtr;};

protected:
  
  virtual void init();

  virtual void printUpdate();

protected:
  
  //! Pointer to the global data object.
  Teuchos::RCP<NOX::GlobalData> globalDataPtr;

  //! Utils
  Teuchos::RCP<NOX::Utils> utilsPtr;
  
  //! Current solution.
  Teuchos::RCP<NOX::Abstract::Group> solnPtr;		

  //! Previous solution pointer. 
  Teuchos::RCP<NOX::Abstract::Group> oldSolnPtr;	

  //! Group used to evaluate a transient residual
  Teuchos::RCP<NOX::Abstract::Group> transientResidualGroup;	

  //! Current search direction pointer.
  Teuchos::RCP<NOX::Abstract::Vector> dirPtr;

  //! Stopping test.
  Teuchos::RCP<NOX::StatusTest::Generic> testPtr;		

  //! Input parameters.
  Teuchos::RCP<Teuchos::ParameterList> paramsPtr;	

  //! Linesearch. 
  Teuchos::RCP<NOX::LineSearch::Generic> lineSearchPtr; 

  //! %Search %Direction. 
  Teuchos::RCP<NOX::Direction::Generic> directionPtr; 

  //! Current step.
  double stepSize;			

  //! Number of nonlinear iterations.
  int nIter;                    

  //! %Status of nonlinear solver.
  NOX::StatusTest::StatusType status;

  //! Type of check to use for status tests.  See NOX::StatusTest for more details.
  NOX::StatusTest::CheckType checkType;

  //! Pointer to a user defined NOX::Abstract::PrePostOperator object.
  NOX::Solver::PrePostOperator prePostOperator;




  //! Pointer to solnPtr casted back to a thyra group 
  Teuchos::RCP<NOX::Thyra::Group> thyraSolnGroup; 
  //! Pointer to oldSolnPtr casted back to a thyra group 
  Teuchos::RCP<NOX::Thyra::Group> thyraOldSolnGroup; 
  //! Group used to evaluate a transient residual
  Teuchos::RCP<NOX::Thyra::Group> thyraTransientResidualGroup; 

  //! Step size for pseudo-transient stepping
  double delta;
  //! Inverse step size for pseudo-transient stepping
  double inv_delta;
  //! Initial step size
  double deltaInit;
  //! Maximum step size
  double deltaMax;
  //! Minimum step size
  double deltaMin;
  //! Step size from previous iteration
  double deltaOld;
  //! Pseudo-transient time
  double time;
  //! solution time derivative used for scaling operator V in pseudo-transient paper
  Teuchos::RCP< ::Thyra::VectorBase<double> > x_dot;

  //! If set to true, the candidate direction will use the transient residual instead of the steady-state residual.  This is a modification of the Kelley-Keyes paper.
  bool use_transient_residual;

  //! Maximum number of iterations before pseudo-transient is disabled and the algorithm switches to a line search-based direct to steady state solve. 
  int max_pseudo_transient_iterations;
};
} // namespace Solver
} // namespace NOX

#endif

