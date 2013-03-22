// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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

#ifndef LOCA_BIFURCATION_TPBORD_STATUSTEST_NULLVECTORNORMWRMS_H
#define LOCA_BIFURCATION_TPBORD_STATUSTEST_NULLVECTORNORMWRMS_H

#include "NOX_StatusTest_Generic.H"	// base class

namespace LOCA {
  
  namespace Bifurcation {
    
    namespace TPBord {

      namespace StatusTest {

	/*!
	 * \brief A convergence test based on the weighted root-mean-square
	 * norm of the update to the null vector component for turning point 
	 * location.
	 */
	/*!
	 * Let \f$n\f$ be the approximation to the null vector for turning 
	 * point tracking (see LOCA::Bifurcation::TPBord::ExtendedGroup).  
	 * This convergence test defines convergence for the null vector when 
	 * the following is true
	 * \f[ 
	 \sqrt{\frac{1}{N}\sum_{i=1}^N\left(\frac{|n_i-(n_0)_i|}{\epsilon_r|n_i| + \epsilon_a}\right)} < \tau
	 * \f]
	 * where \f$n_0\f$ is the previous approximation to the null vector, 
	 * \f$N\f$ is the length of \f$n\f$, \f$\epsilon_r\f$
	 * is the relative tolerance, \f$\epsilon_a\f$ is the absolute
	 * tolerance, and \f$\tau\f$ is an overall scale factor (typically
	 * \f$\tau = 1\f$).
	 *
	 * Note that this status test deals only with the null vector component
	 * of the turning point equations.  This status test should be combined
	 * with other status tests for the solution and parameter components
	 * (using NOX::StatusTest::Combo and LOCA::StatusTest::Wrapper) to 
	 * build a composite status test for the entire system.
	 *
	 * Also note that if the group returned by the getSolutionGroup()
	 * method of the solver supplied in checkStatus() is not a 
	 * turning point group (i.e., not derived from 
	 * LOCA::Bifurcation::TPBord::ExtendedGroup), checkStatus() returns 
	 * NOX::StatusTest::Converged.  This allows the status test to 
	 * be used in situations other than turning point tracking, e.g., 
	 * steady-state solves, without raising error conditions.
	 */
	class NullVectorNormWRMS : public NOX::StatusTest::Generic {

	public:

	  //! Constructor
	  /*!
	   * \em rtol is the relative tolerance \f$\epsilon_r\f$, \em atol
	   * is the absolute tolerance \f$\epsilon_a\f$, and \em tol is the
	   * overall scale factor \f$\tau\f$ defined above.
	   */
	  NullVectorNormWRMS(double rtol, double atol, double tol);

	  //! Destructor.
	  virtual ~NullVectorNormWRMS();

	  //! Evaluates convergence criteria specified above
	  virtual NOX::StatusTest::StatusType 
	  checkStatus(const NOX::Solver::Generic& problem);
	
	  //! Returns status as defined above
	  virtual NOX::StatusTest::StatusType 
	  getStatus() const;

	  //! Prints current status
	  virtual std::ostream& print(ostream& stream, int indent = 0) const;


	  /* @name Accessor Functions
	   * Used to query current values of variables in the status test.
	   */
	  //@{

	  //! Returns the value of weighted parameter update norm.
	  double getNullVectorNormWRMS() const;

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

	  //! Tolerance for convergence
	  double tol;

	  //! Current value of the parameter update norm
	  double normWRMS;

	  //! %Status of entire %status test
	  NOX::StatusTest::StatusType status;
	};
      } // namespace StatusTest
    } // namespace TPBord
  } // namespace Bifurcation
} // namespace LOCA

#endif
