// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_CONTINUATION_STATUSTEST_PARAMETERRESIDUALNORM_H
#define LOCA_CONTINUATION_STATUSTEST_PARAMETERRESIDUALNORM_H

#include "NOX_StatusTest_Generic.H"	// base class

namespace LOCA {
  
  namespace Continuation {

    //! Convergence tests relating to continuation groups
    namespace StatusTest {

      /*!
       * \brief A convergence test based on the parameter component of the
       * residual for continuation.
       */
      /*!
       * Consider a continuation method with parameter equation \f$g = 0\f$
       * (see LOCA::Continuation::ExtendedGroup).  This convergence test
       * defines convergence of the parameter equation when the following
       * is true
       * \f[ 
            \frac{|g|}{\epsilon_r|\Delta s| + \epsilon_a} < \tau
       * \f]
       * where \f$\Delta s\f$ is the current step size, \f$\epsilon_r\f$
       * is the relative tolerance, \f$\epsilon_a\f$ is the absolute
       * tolerance, and \f$\tau\f$ is an overall scale factor (typically
       * \f$\tau = 1\f$).
       *
       * Note that this status test deals only with the parameter component
       * of the continuation equations.  This status test should be combined
       * with other status tests for the solution component (using 
       * NOX::StatusTest::Combo and LOCA::StatusTest::Wrapper) to build
       * a composite status test for the entire system.
       *
       * Also note that if the group returned by the getSolutionGroup()
       * method of the solver supplied in checkStatus() is not a 
       * continuation group (i.e., not derived from 
       * LOCA::Continuation::ExtendedGroup), checkStatus() returns 
       * NOX::StatusTest::Converged.  This allows the status test to 
       * be used in situations other than continuation, e.g., steady-state
       * solves, without raising error conditions.
       */
      class ParameterResidualNorm : public NOX::StatusTest::Generic {

      public:

	//! Constructor
	/*!
	 * \em rtol is the relative tolerance \f$\epsilon_r\f$, \em atol
	 * is the absolute tolerance \f$\epsilon_a\f$, and \em tol is the
	 * overall scale factor \f$\tau\f$ defined above.
	 */
	ParameterResidualNorm(double rtol, double atol, double tol);

	//! Destructor.
	virtual ~ParameterResidualNorm();

	//! Evaluates convergence criteria specified above
	virtual NOX::StatusTest::StatusType 
	checkStatus(const NOX::Solver::Generic& problem);
	
	//! Returns status as defined above
	virtual NOX::StatusTest::StatusType 
	getStatus() const;

	//! Prints current status
	virtual ostream& print(ostream& stream, int indent = 0) const;


	/* @name Accessor Functions
	 * Used to query current values of variables in the status test.
	 */
	//@{
   
	//! Returns the value of scaled parameter residual norm
	double getResidualNorm() const;

	//! Returns the realative tolerance set in the constructor.
	double getRTOL() const;

	//! Returns the absolute tolerance set in the constructor.  
	double getATOL() const;

	//! Returns the tolerance set in the constructor
	double getTOL() const;

	//@}

      protected:

	//! Relative tolerance for convergence.
	double rtol;

	//! Absolute tolerance for convergence.
	double atol;

	//! Overall scale factor
	double tol;

	//! Current value of the scaled parameter residual norm
	double paramResidualNorm;

	//! Status of entire status test
	NOX::StatusTest::StatusType status;
      };
    } // namespace StatusTest
  } // namespace Continuation
} // namespace LOCA

#endif
