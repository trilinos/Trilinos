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

#ifndef LOCA_CONTINUATION_STATUSTEST_PARAMETERUPDATENORM_H
#define LOCA_CONTINUATION_STATUSTEST_PARAMETERUPDATENORM_H

#include "NOX_StatusTest_Generic.H"	// base class

namespace LOCA {
  
  namespace Continuation {

    namespace StatusTest {

      class ParameterUpdateNorm : public NOX::StatusTest::Generic {

      public:

	//! Constructor
	ParameterUpdateNorm(double rtol, double atol, double tol);

	//! Destructor.
	virtual ~ParameterUpdateNorm();

	//! check status
	virtual NOX::StatusTest::StatusType 
	checkStatus(const NOX::Solver::Generic& problem);
	
	//! get status
	virtual NOX::StatusTest::StatusType 
	getStatus() const;

	//! print
	virtual ostream& print(ostream& stream, int indent = 0) const;


	/* @name Accessor Functions
	 * Used to query current values of variables in the status test.
	 */
	//@{

	//! Returns the value of weighted parameter update norm.
	double getParameterUpdateNorm() const;

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
	double paramUpdateNorm;

	//! %Status of entire %status test
	NOX::StatusTest::StatusType status;
      };
    } // namespace StatusTest
  } // namespace Continuation
} // namespace LOCA

#endif
