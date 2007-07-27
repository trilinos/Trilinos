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

#ifndef NOX_STATUSTEST_FACTORY_H
#define NOX_STATUSTEST_FACTORY_H

#include "Teuchos_RCP.hpp" 

namespace Teuchos {
  class ParameterList;
}

namespace NOX {

  class Utils;

namespace StatusTest {

//! Factory to build a set of status tests from a parameter list
/*! 
    This object takes either an xml file name or a
    Teuchos::ParameterList and generates an entire set (a tree) of
    status tests for use in a NOX::Solver derived object.

    The tagged_tests field in the constructors allows users to store
    tests from the tree in a flat list in case they want to change the
    tolerance values during a run.  The tagged_tests flag is optional.

    \author Roger Pawlowski (SNL 1416)
*/
class Factory {

public:

  //! Constructor.
  Factory();

  //! Destructor.
  virtual ~Factory();

  //! Returns a status test set from a parameter list xml file.
  Teuchos::RCP<NOX::StatusTest::Generic> 
  buildStatusTests(const std::string& file_name, const NOX::Utils& utils,
		   std::map<std::string, Teuchos::RCP<NOX::StatusTest::Generic> >*
		   tagged_tests=0) const;

  //! Returns a status test set from a parameter list.
  Teuchos::RCP<NOX::StatusTest::Generic> 
  buildStatusTests(Teuchos::ParameterList& p, const NOX::Utils& utils,
		   std::map<std::string, Teuchos::RCP<NOX::StatusTest::Generic> >*
		   tagged_tests=0) const;

private:

  Teuchos::RCP<NOX::StatusTest::Generic> 
  buildComboTest(Teuchos::ParameterList& p, 
		 const NOX::Utils& utils,
		 std::map<std::string, Teuchos::RCP<NOX::StatusTest::Generic> >*
		 tagged_tests=0) const;

  Teuchos::RCP<NOX::StatusTest::Generic> 
  buildNormFTest(Teuchos::ParameterList& p, 
		 const NOX::Utils& utils) const;

  Teuchos::RCP<NOX::StatusTest::Generic> 
  buildNormUpdateTest(Teuchos::ParameterList& p, 
		      const NOX::Utils& utils) const;

  Teuchos::RCP<NOX::StatusTest::Generic> 
  buildNormWRMSTest(Teuchos::ParameterList& p, 
		    const NOX::Utils& utils) const;

  Teuchos::RCP<NOX::StatusTest::Generic> 
  buildMaxItersTest(Teuchos::ParameterList& p, 
		    const NOX::Utils& utils) const;

  Teuchos::RCP<NOX::StatusTest::Generic> 
  buildFiniteValueTest(Teuchos::ParameterList& p, 
		       const NOX::Utils& utils) const;

  Teuchos::RCP<NOX::StatusTest::Generic> 
  buildDivergenceTest(Teuchos::ParameterList& p, 
		      const NOX::Utils& utils) const;

  Teuchos::RCP<NOX::StatusTest::Generic> 
  buildStagnationTest(Teuchos::ParameterList& p, 
		      const NOX::Utils& utils) const;

  //! Checks if a tag is present in the param list and adds the test to the tagged_test std::map if true.  Returns true if a tag was present. 
  bool checkAndTagTest(const Teuchos::ParameterList& p,
		       const Teuchos::RCP<NOX::StatusTest::Generic>& test,
         std::map<std::string, Teuchos::RCP<NOX::StatusTest::Generic> >* 
		       tagged_tests) const;
};

} // namespace StatusTest
} // namespace NOX

#endif
