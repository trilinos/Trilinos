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

#ifndef LOCA_STATUSTEST_WRAPPER_H
#define LOCA_STATUSTEST_WRAPPER_H

#include "NOX_StatusTest_Generic.H" // base class

namespace LOCA {
  namespace StatusTest {

    //! Adaptor between NOX and LOCA status tests.
    class Wrapper : public NOX::StatusTest::Generic {

    public:

      //! Constructor
      Wrapper(NOX::StatusTest::Generic& s);

      //! Destructor
      virtual ~Wrapper();

      //! Calls checkStatus of underlying status test
      virtual NOX::StatusTest::StatusType 
      checkStatus(const NOX::Solver::Generic& problem);

      //! Calls checkStatus of underlying status test
      virtual NOX::StatusTest::StatusType
      checkStatus(const NOX::Solver::Generic& problem, 
		 NOX::StatusTest::CheckType checkType);

      //! Calls getStatus of underlying status test
      virtual NOX::StatusTest::StatusType getStatus() const;
      
      //! Calls print of underlying status test
      virtual ostream& print(ostream& stream, int indent = 0) const;

      //! Returns underlying status test
      NOX::StatusTest::Generic&
      getUnderlyingStatusTest();

      //! Returns underlying status test
      const NOX::StatusTest::Generic&
      getUnderlyingStatusTest() const;
  
    protected:
      NOX::StatusTest::Generic *statusTestPtr;
      
    }; // class Wrapper

  } // namespace StatusTest
} // namespace LOCA


#endif
