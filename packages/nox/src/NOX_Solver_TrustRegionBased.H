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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef NOX_SOLVER_TRUSTREGIONBASED_H
#define NOX_SOLVER_TRUSTREGIONBASED_H

#include "NOX_Solver_Generic.H"	    // base class
#include "NOX_LineSearch_Manager.H" // class data element
#include "NOX_Direction_Manager.H"  // class data element
#include "NOX_Parameter_List.H"	    // class data element

namespace NOX {
namespace Solver {

/*!
  \brief %Newton-like solver using a trust region.

Our goal is to solve:
\f[
F(x) = 0,
\f]
where \f$ F:\Re^n \rightarrow \Re^n \f$. 
Alternatively, we might say that we wish to solve
\f[
\min f(x) \equiv \frac{1}{2} \|F(x)\|^2_2.
\f]
The trust region subproblem at iteration \f$k\f$ is given by
\f[
\min m_k(s) \equiv f_k + g_k^T s + \frac{1}{2} s^T B_k s,
\mbox{ s.t. } \|s\| \leq \Delta_k
\f]
where 
\f$ f_k = f(x_k) = \frac{1}{2} \|F(x_k)\|^2_2 \f$,
\f$ g_k = \nabla f(x_k) = J(x_k)^T F(x_k) \f$,
\f$ B_k =  J(x_k)^T J(x_k) \approx \nabla^2 f(x_k) \f$,
\f$ J(x_k)\f$ is the Jacobian of \f$F\f$ at \f$x_k\f$, and
\f$ \Delta_k \f$ is the trust region radius.

  
The "improvement ratio" for a given step \f$ s \f$ is defined as
\f[

\frac{ f(x) - f(x + s) } { m(0) - m(s) }

\f]


  The following parameters are valid for this solver:

  - "%Direction" - Sublist of the direction parameters for the
    %Newton point, passed to the NOX::Direction::Manager
    constructor. If "Method" is not specified, it defaults to
    "%Newton".

  - "%Cauchy %Direction" - Sublist of the direction parameters for the
    Cauchy point, passed to the NOX::Direction::Manager
    constructor. If "Method" is not specified, it defaults to
    "Steepest Descent".

- "Minimum Trust Region Radius" - Minimum allowable trust region
  radius. Defaults to 1.0e-6.

- "Maximum Trust Region Radius" - Minimum allowable trust region
  radius. Defaults to 1.0e+10.

- "Minimum Improvement Ratio" - Minimum improvement ratio to accept
  the step. Defaults to 1.0e-4.

- "Contraction Trigger Ratio" - If the improvement ratio is less than
  this value, then the trust region is contracted by the amount
  specified by the "Contraction Factor". Must be larger than "Minimum
  Improvement Ratio". Defaults to 0.1.

- "Contraction Factor" - See above. Defaults to 0.25.

- "Expansion Trigger Ratio" - If the improvement ratio is greater than this
  value, then the trust region is contracted by the amount specified
  by the "Expansion Factor". Defaults to 0.75.

- "Expansion Factor" - See above. Defaults to 4.0.

- "Recovery Step" - Defaults to 1.0.

A sublist for output parameters called "Output" will be created and contain the following parameters:

- "Nonlinear Iterations" - Number of nonlinear iterations

- "2-Norm or Residual" - Two-norm of final residual

  \author Tammy Kolda (SNL 8950), Roger Pawlowski (SNL 9233)
*/

class TrustRegionBased : public Generic {

public:

  //! %TrustRegionBased-method constructor. 
  TrustRegionBased(NOX::Abstract::Group& xgrp, NOX::StatusTest::Generic& t, const NOX::Parameter::List& p);

  //! Destructor
  virtual ~TrustRegionBased();

  virtual bool reset(NOX::Abstract::Group& xgrp, NOX::StatusTest::Generic& t, const NOX::Parameter::List& p);
  virtual NOX::StatusTest::StatusType getStatus();
  virtual NOX::StatusTest::StatusType iterate();
  virtual NOX::StatusTest::StatusType solve();
  virtual const NOX::Abstract::Group& getSolutionGroup() const;
  virtual const NOX::Abstract::Group& getPreviousSolutionGroup() const;
  virtual int getNumIterations() const;
  virtual const NOX::Parameter::List& getParameterList() const;

protected:
  
  //! Print out initialization information and calcuation the RHS.
  virtual void init();

  //! Prints the current iteration information.
  virtual void printUpdate();

protected:
  
  //! Current solution.
  NOX::Abstract::Group* solnPtr;		

  //! Previous solution pointer. 
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
  NOX::Abstract::Group* oldSolnPtr;	
  //! Previous solution reference.
  NOX::Abstract::Group& oldSoln;	

  //! Current search direction.pointer.
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
  NOX::Abstract::Vector* newtonVecPtr;
  //! Current Newton direction.reference.
  NOX::Abstract::Vector& newtonVec;

  //! Current search direction.pointer.
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
  NOX::Abstract::Vector* cauchyVecPtr;
  //! Current Cauchy direction.reference.
  NOX::Abstract::Vector& cauchyVec;

  //! Extra vector used in computations
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
  NOX::Abstract::Vector* aVecPtr;
  //! Extra vector
  NOX::Abstract::Vector& aVec;

  //! Extra vector used in computations
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
  NOX::Abstract::Vector* bVecPtr;
  //! Extra vector
  NOX::Abstract::Vector& bVec;

  //! Stopping test.
  NOX::StatusTest::Generic* testPtr;		

  //! Input parameters.
  NOX::Parameter::List params;	

  //! %Newton %Search %Direction. 
  NOX::Direction::Manager newton; 

  //! Cauchy %Search %Direction. 
  NOX::Direction::Manager cauchy; 

  //! Radius of the trust region
  double radius;

  //! Minimum improvement ratio to accept step
  double minRatio;

  //! Minimum trust region radius
  double minRadius;

  //! Maximum trust region radius
  double maxRadius;

  //! ratio < alpha triggers contraction
  double contractTriggerRatio;

  //! ratio > beta triggers expansion
  double expandTriggerRatio;

  //! Expansion factor
  double expandFactor;

  //! Constraction factor
  double contractFactor;

  //! Take a step of this length in the Newton direction if the
  //! trust-region search fails
  double recoveryStep;

  //! Value of \f$ f \f$ at current solution
  double newF;
  //! Value of \f$ f \f$ at previous solution
  double oldF;

  //! norm(xnew - xold)
  double dx;

  //! Number of nonlinear iterations.
  int nIter;                    

  //! %Status of nonlinear solver.
  NOX::StatusTest::StatusType status;

  //! Enumerated list for each direction that may be required in the Trust region computation.
  enum StepType {
    //! Use the Newton direction
    Newton, 
    //! Use the Cauchy direction
    Cauchy, 
    //! Use the doglog direction
    Dogleg
  };

  //! Type of step to be taken.
  StepType stepType;


};
} // namespace Solver
} // namespace NOX

#endif

