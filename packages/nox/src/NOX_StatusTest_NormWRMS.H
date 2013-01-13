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

#ifndef NOX_STATUS_NORM_WRMS_H
#define NOX_STATUS_NORM_WRMS_H

#include "NOX_StatusTest_Generic.H"	// base class
#include "Teuchos_RCP.hpp"

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

` If the number of iterations is zero, then the status is set to
  NOX::StatusTest::Unconverged and returned. (Also, #value is set to
  1.0e+12.)

  Otherwise, returns NOX::StatusTest::Converged if the three criteria
  listed below are satisfied. Note that use of Criteria #2 and #3
  depend on the options set in the solver.

  <ol>
  <li> Weigthed root mean square norm is less than a specified tolerance:

  \f[
  ||\delta x^k||_{wrms} < \mbox{tolerance}
  \f]
 
  where

 \f[
  ||\delta x^k||_{wrms} \equiv C \sqrt{ \frac{1}{N} \sum_{i=1}^N \left( \frac
    {(x^k_i-x^{k-1}_i)}{RTOL |x^{k-1}_i| + ATOL_i} \right) ^2 }
  \f]

  Here:

  - \f$x_i^k\f$ denotes component \f$i\f$ of nonlinear iterate \f$k\f$.
  - \f$N\f$ denotes the number of unknowns
  - \f$RTOL\f$ denotes the relative error tolerance, specified via \c
    rtol in the constructor
  - \f$ATOL\f$ denotes the absolution error tolerance, specified via
    \c atol in the constructor. This can be a vector or a scalar.
  - \f$C\f$ denotes a weight, specified via the parameter \c
    BDFMultiplier in the constructor.

  <li> <em>If a line search based solver is used</em>, the line search
  step size, \f$ \lambda \f$, must be greater than a specified step
  size value, \f$ \alpha \f$:

  \f[ \lambda > \alpha \f]

  The motivation for this test is to avoid detecting stagnation
  when in fact the true problem is that the step size is just small.

  The value of \f$\alpha\f$ is set in the constructor via the argument
  \c alpha. Setting \f$\alpha\f$ to zero effectively eliminates this
  part of the test.

  <li> The achieved linear solver tolerance, \f$ \eta^k \f$ for
  nonlinear iteration \f$ k \f$, must be less than a specified
  tolerance value, \f$ \beta \f$; i.e.,

  \f[ \eta^k < \beta \f]

  The motivation for this test is to avoid detecting stagnation when
  in fact the true problem is that the linear solve tolerance was not
  accurate enough. 

  The value of \f$\beta\f$ is set in the constructor via the argument
  \c beta. Setting \f$\beta\f$ to 1.0 effectively eliminates this
  part of the test.

  \note This criteria will only be used if the "Achieved Tolerance"
  parameter (the value of \f$ \eta^k \f$) is set by the linear solver
  in the "Newton"/"Linear Solver"/"Output" sublist.  The
  checkStatus() method will search for this parameter.  

</ol>

  References: 

  <ol>

  <li> K. E. Brennam, S. L. Cambell, L. R. Petzold, <em>Numerical
  Solution of Initial-Value Problems in Differential-Algebraic
  Equations</em>, Classics in Applied Mathematics 14, SIAM 1996.

  <li> G. D. Byrne and A. C. Hindmarch, <em>PVODE, an ODE %Solver for
  Parallel Computers</em>, Technical Report UCRL-JC-132361, Rev. 1,
  Center for Applied Scientific Computing (CASC), Lawrence Livermore
  National Lab, May 1999.

  </ol>

*/
class NormWRMS : public Generic {

public:

  //! Constructor where ATOL is a scalar
  NormWRMS(double rtol, 
	   double atol, 
	   double BDFMultiplier = 1.0, 
	   double tolerance = 1.0,
	   double alpha = 1.0,
	   double beta = 0.5,
	   bool disable_implicit_weighting = true);

  //! Constructor where ATOL is a vector
  NormWRMS(double rtol,
	   const Teuchos::RCP<const NOX::Abstract::Vector>& atol,
	   double BDFMultiplier = 1.0, 
	   double tolerance = 1.0,
	   double alpha = 1.0,
	   double beta = 0.5,
	   bool disable_implicit_weighting = true);

  //! Destructor.
  virtual ~NormWRMS();

  // Derived classes
  virtual StatusType 
  checkStatus(const NOX::Solver::Generic& problem,
	      NOX::StatusTest::CheckType checkType);
  virtual StatusType getStatus() const;
  virtual std::ostream& print(std::ostream& stream, int indent = 0) const;


  /* @name Accessor Functions
     Used to query current values of variables in the status test.
  */
  //@{

  //! Returns the value of WRMS norm.
  virtual double getNormWRMS() const;
   
  //! Returns the requested tolerance set in the constructor.
  virtual double getTolerance() const;

  //! Returns the realative tolerance set in the constructor.
  virtual double getRTOL() const;

  //! Returns the absolute tolerance set in the constructor.  If ATOL is a vector, this will return a value of -1.0.
  virtual double getATOL() const;

  //! Returns the value of the BDFMultiplier set in the constructor.
  virtual double getBDFMultiplier() const;

  //! Returns the value of 'alpha' set in the constructor.
  virtual double getAlpha() const;

  //! Returns the value of 'beta' set in the constructor.
  virtual double getBeta() const;

  //! Returns true if implicit weighting is disabled during norm computation
  virtual bool getDisableImplicitWeighting() const;

  //@}

private:

  //! Current value of the WRMS norm (Must be less than 'tolerance' for convergence).  
  double value;

  //! Relative tolerance for convergence.
  double rtol;

  //! flag to tell whether atol is a scalar or a vector.
  bool atolIsScalar;

  //! Absolute tolerance for convergence (scalar).
  double atol;

  //! Absolute tolerance for convergence (vector).
  Teuchos::RCP<const NOX::Abstract::Vector> atolVec;

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
  Teuchos::RCP<NOX::Abstract::Vector> u;

  //! Temporary vector used in computation.
  Teuchos::RCP<NOX::Abstract::Vector> v;

  //! Flag that tells the print method whether to print the criteria 2 information.
  bool printCriteria2Info;

  //! Flag that tells the print method whether to print the criteria 3 information.
  bool printCriteria3Info;

  //! If true, the implicit weighting will be disabled during norm calculation
  bool m_disable_implicit_weighting;
};

} // namespace Status
} // namespace NOX

#endif
