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

#ifndef NOX_STATUS_NORM_WRMS_H
#define NOX_STATUS_NORM_WRMS_H

#include "NOX_StatusTest_Generic.H"	// base class

// Forward declarations
namespace NOX {
namespace Abstract {
class Vector;
}
}

namespace NOX {

namespace StatusTest {

//! Convergence test based on the weighted root mean square norm fo the solution update between iterations.
/*!
  Convergence of the nonlinear solver is achieved if the following three criteria are met:

  <ol>
  <li> Weigthed root mean square norm is less than a specified tolerance:

  \f[
  ||\delta x^k||_{wrms} < \mbox{tolerance}
  \f]
 
  where \f$ ||\delta x^k||_{wrms} \f$ is the weigthed root mean square norm 
  at nonlinear iteration \f$ k \f$.  "tolerance" is a user supplied parameter 
  that defaults to 1.0 (1.0 is used almost universally). The user can change 
  this value in the constructor by specifying the "tolerance" variable.

  The weighted root mean square norm at iteration \f$ k \f$ is defined by:

  \f[
  ||\delta x^k||_{wrms} =   \sqrt{ \frac{1}{N} \sum_{i=1}^N \left( \frac {C*(x^k_i-x^{k-1}_i)}{RTOL |x^{k-1}_i| + ATOL_i} \right) ^2  }
  \f]

  Where \f$ N \f$ is the total number of unknowns in the problem. 
  \f$ x^k_i \f$ denotes the \f$ i \f$ -th component of the solution 
  vector \f$ x \f$ at nonlinear iteration \f$ k \f$. \f$ RTOL \f$ is the 
  relative error tolerance (a scalar value).  \f$ ATOL_i \f$
  is the absolute error tolerance and can be a scalar or
  a vector of the same size as the solution vector. 
  \f$ C \f$ is a scalar value that is typically set to 1.0.  When 
  NormWRMS is used in conjunction with time stepping algorithms such as
  Backward Differentiation Formulas (BDF) a different value may be required 
  based on the specific BDF algorithm utlilized.  The user is free to set 
  this parameter in the consturctor with the "BDFMultiplier" variable.
  Defaults to 1.0.

  NOTE:  If no iterations have been performed, then there is no information for\f$ x_{k-1} \f$ and thus the WRMS norm can not be computed.  If the test is called while the number of iterations is zero (i.e. before an iteration takes place), the code returns the status test as Unconverged and sets the value of the WRMS norm to 1.0e+12.  Once an iteration has taken place the WRMS Norm value is calculated as above.  

  <li> If a line search is used, the line search step size, \f$ \lambda \f$, must be greater than a specified step size value, \f$ \alpha \f$:

  \f[ \lambda > \alpha \f]

  The WRMS norm is essentially a stagnation test: i.e. termination of the convergence test is achieved when the norm of the <B>change</B> in the solution vector between iterations (\f$ \delta x = x^k - x^{k-1} \f$) is below a specified value.  If a line search is used, the line search algorithm could cut the step length so that the change in \f$ x \f$ is small enough to prematurely claim convergence.  This second criteria will prevent this premature convergence by insisting that the line search step size be larger than a specified value \f$ \alpha \f$.  The value \f$ \alpha \f$ can be set in the constructor with the "alpha" variable.  Defaults to 1.0.  Setting this value > 1.0 essentially eliminates this criteria from the status test.  NOTE: This criteria will only be used if a "Line Search Based" solver is being used as the "Nonlinear %Solver" (see NOX::Solver::LineSearchBased).

  <li> The achieved linear solver tolerance, \f$ \eta^k \f$ for nonlinear iteration \f$ k \f$, must be less than a specified tolerance value, \f$ \beta \f$:

  \f[ \eta^k < \beta \f]
The WRMS norm is essentially a stagnation test: i.e. termination of the convergence test is achieved when the norm of the <B>change</B> in the solution vector between iterations (\f$ \delta x = x^k - x^{k-1} \f$) is below a specified value.  If an adjustable forcing term algorithm is used (see NOX::Direction::Newton for details), then far from the solution when the linear solve tolerance, \f$ \eta^k \f$, is large, the linear solver may take so few iterations (or possibly zero iterations) that the change in \f$ x \f$ is small enough to prematurely claim convergence.  This third criteria prevents premature convergence by requiring that the linear solve tolerance be below a specified value.  The value of the linear solve tolerance, \f$ \beta \f$, can be specified in the constructor with the "beta" variable.  Defaults to 0.5 (heuristic).  Setting this value >= 1.0 essentially eliminates this criteria from the status test.  NOTE: This criteria will only be used if the "Achieved Tolerance" parameter (the value of \f$ \eta \f$) is set by the linear solver in the "Output" sublist of the "%Linear %Solver" sublist.

</ol>

  References: 

  (1) "Classics in Applied Mathematics #14: Numerical Solution of 
  Initial-Value Problems in Differential-Algebraic Equations", 
  Brennan, K. E., Campbell, S. L., Petzold, L. R., SIAM 1996. 

  (2) "PVODE, an ODE %Solver for Parallel Computers", Byrne,
  G. D., and Hindmarsh, A. C., Center for Applied Scientific Computing
  (CASC), UCRL-JC-132361, Rev. 1, May 1999.

*/
class NormWRMS : public Generic {

public:

  //! Constructor where ATOL is a scalar
  NormWRMS(double rtol, 
	   double atol, 
	   double BDFMultiplier = 1.0, 
	   double tolerance = 1.0,
	   double alpha = 1.0,
	   double beta = 0.5);

  //! Constructor where ATOL is a vector
  NormWRMS(double rtol,
	   Abstract::Vector& atol,
	   double BDFMultiplier = 1.0, 
	   double tolerance = 1.0,
	   double alpha = 1.0,
	   double beta = 0.5);

  //! Destructor.
  virtual ~NormWRMS();

  // Derived classes
  virtual StatusType checkStatus(const NOX::Solver::Generic& problem);
  virtual StatusType getStatus() const;
  virtual ostream& print(ostream& stream, int indent = 0) const;

private:

  //! Current value of the 
  double value;

  //! Relative tolerance for convergence.
  double rtol;

  //! flag to tell whether atol is a scalar or a vector.
  bool atolIsScalar;

  //! Absolute tolerance for convergence (scalar).
  double atol;

  //! Absolute tolerance for convergence (vector).
  Abstract::Vector* atolVec;

  //! Time integration method multiplier (BDF Multiplier).
  double factor;

  //! Required tolerance for the NormWRMS to be declared converged.
  double tolerance;

  //! Minimum step size allowed during a line search for WRMS norm to be flagged as converged.
  double alpha;

  //! Actual step size used during line search.
  double computedStepSize;
  
  //! Maximum linear solve tolerance allowed for WRMS norm to be flagged as converged.
  double beta;

  //! Actual tolerance achieved by the linear solver during the last linear solve.
  double achievedTol;

  //! %Status
  StatusType status;

  //! Temporary vector used in computation.
  Abstract::Vector* u;

  //! Temporary vector used in computation.
  Abstract::Vector* v;

  //! Flag that tells the print method whether to print the criteria 2 information.
  bool printCriteria2Info;

  //! Flag that tells the print method whether to print the criteria 3 information.
  bool printCriteria3Info;

};

} // namespace Status
} // namespace NOX

#endif
