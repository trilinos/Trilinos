// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
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

#ifndef LOCA_STATUSTEST_WRAPPER_H
#define LOCA_STATUSTEST_WRAPPER_H

#include "Teuchos_RCP.hpp"

#include "NOX_StatusTest_Generic.H" // base class

namespace LOCA {

  //! %LOCA %StatusTest classes
  namespace StatusTest {

    /*!
     * \brief A %wrapper class for wrapping a %NOX status test.
     */
    /*!
     * The %LOCA::StatusTest::Wrapper class provides a wrapper for %NOX
     * status tests to change the solver passed to the wrapped status test.
     * The solver passed through the checkStatus() method is wrapped via 
     * the LOCA::Solver::Wrapper class and then forwarded to the checkStatus()
     * method of the wrapped status test.  The purpose of this is to allow
     * status tests that use concrete group data to function correctly when
     * the group is stored in an extended continuation or bifurcation group.
     * (See LOCA::Solver::Wrapper for more details or the %LOCA status tests
     * page for examples on how to effectively use this class.)
     */
    class Wrapper : public NOX::StatusTest::Generic {

    public:

      //! Constructor
      Wrapper(const Teuchos::RCP<NOX::StatusTest::Generic>& s);

      //! Destructor
      virtual ~Wrapper();

      //! Calls checkStatus of underlying status test
      virtual NOX::StatusTest::StatusType
      checkStatus(const NOX::Solver::Generic& problem, 
		  NOX::StatusTest::CheckType checkType);

      //! Calls getStatus of underlying status test
      virtual NOX::StatusTest::StatusType getStatus() const;
      
      //! Calls print of underlying status test
      virtual ostream& print(ostream& stream, int indent = 0) const;

      //! Returns underlying status test
      Teuchos::RCP<NOX::StatusTest::Generic>
      getUnderlyingStatusTest();

      //! Returns underlying status test
      Teuchos::RCP<const NOX::StatusTest::Generic>
      getUnderlyingStatusTest() const;

    private:

      //! Private to prohibit copying
      Wrapper(const Wrapper&);

      //! Private to prohibit copying
      Wrapper& operator=(const Wrapper&);
  
    protected:

      //! Stores pointer to underlying status test
      Teuchos::RCP<NOX::StatusTest::Generic> statusTestPtr;
      
    }; // class Wrapper

  } // namespace StatusTest
} // namespace LOCA


#endif
