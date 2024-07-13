// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_STATUSTEST_FACTORY_H
#define LOCA_STATUSTEST_FACTORY_H

#include "Teuchos_RCP.hpp"

#include "LOCA_StatusTest_Abstract.H"
#include "LOCA_GlobalData.H"

namespace Teuchos {
  class ParameterList;
}

namespace LOCA {

  class Utils;

namespace StatusTest {

  // TODO Update documentation.

//! Factory to build a set of status tests from a parameter list
/*!
    This object takes either an XML file name or a
    Teuchos::ParameterList and generates an entire set (a tree) of
    status tests for use in a LOCA::Stepper derived object.

    The tagged_tests field in the constructors allows users to store
    tests from the tree in a flat list in case they want to change the
    tolerance values during a run.  The tagged_tests flag is optional.

    Please use the related nonmember functions instead of calling the
    factory directly (See example below).

    Valid parameters are as follows:

    <ul>
      <li> "Test Type" <std::string> Type of test this list contains.  Valid tests include:
      <ul>
         <li> "Combo" - NOX::StatusTest::Combo
     <li> "MaxIters" - LOCA::StatusTest::MaxIters
     <li> "User Defined" - A user constructed test, derived from NOX::StatusTest::Generic.
       </ul>

       <li> "Tag" <std::string> A unique identifier that will place the test in the map for tagged_tests.  This allows users to access individual tests to change tolerances on the fly or query values while still using the factory to build objects.

     </ul>

     Additional parameters valid for a Combo test (LOCA::StatusTest::Combo):
     <ul>
       <li> "Combo Type" <std:string> Type of combo to use.  Valid options are:
       <ul>
         <li> "AND"
         <li> "OR"
       </ul>

       <li> "Number of Tests" <int> Number of sublists that contain tests to be added to this combo test.  The sublists must be named "Test X" where "X" represents the test number starting with 0 and preceeding to "Number of Tests - 1".

       <li> "Test X" <Teuchos::ParameterList> A sublist containing a test to add to the current combo test.  The "X" represents the number of the test.  the numbering starts with 0 and is valid through "Number of Tests - 1" tests.
     </ul>

     Additional parameters valid for a Maximum Iterations test
     (NOX::StatusTest::MaxIters):
     <ul>
       <li> "Maximum Iterations" <int>
     </ul>

     Additional parameters valid for a "User Defined" test:
     <ul>
       <li> "User Status Test" < Teuchos::RCP<LOCA::StatusTest::Abstract> > A status test suppied by the user.  It is very important that when registering this status test, that the user set it as a "Generic" object since there is no implicit casting on the ParameterList's get method.  See the example below.
     </ul>

     Example usage:
     \code
     Teuchos::ParameterList stl;
     stl.set("Test Type", "Combo");
     stl.set("Combo Type", "OR");
     stl.set("Number of Tests", 5);
     Teuchos::ParameterList& conv = stl.sublist("Test 0");
     Teuchos::ParameterList& fv = stl.sublist("Test 1");
     Teuchos::ParameterList& divergence = stl.sublist("Test 2");
     Teuchos::ParameterList& stagnation = stl.sublist("Test 3");
     Teuchos::ParameterList& maxiters = stl.sublist("Test 4");

     conv.set("Test Type", "Combo");
     conv.set("Combo Type", "AND");
     conv.set("Number of Tests", 4);
     Teuchos::ParameterList& normF = conv.sublist("Test 0");
     Teuchos::ParameterList& normWRMS = conv.sublist("Test 1");
     Teuchos::ParameterList& normUpdate = conv.sublist("Test 2");
     Teuchos::ParameterList& userDefined = conv.sublist("Test 3");
     normF.set("Test Type", "NormF");
     normF.set("Tolerance", 1.0e-12);
     normF.set("Norm Type", "Two Norm");
     normF.set("Scale Type", "Unscaled");
     normWRMS.set("Test Type", "NormWRMS");
     normWRMS.set("Absolute Tolerance", 1.0e-8);
     normWRMS.set("Relative Tolerance", 1.0e-5);
     normWRMS.set("Tolerance", 1.0);
     normWRMS.set("BDF Multiplier", 1.0);
     normWRMS.set("Alpha", 1.0);
     normWRMS.set("Beta", 0.5);
     normUpdate.set("Test Type", "NormUpdate");
     normUpdate.set("Norm Type", "One Norm");
     normUpdate.set("Scale Type", "Scaled");
     userDefined.set("Test Type", "User Defined");
     Teuchos::RCP<NOX::StatusTest::Generic> myTest =
       Teuchos::rcp(new MyTest(1.0e-3));
     userDefined.set("User Status Test", myTest);


     fv.set("Test Type", "FiniteValue");
     fv.set("Vector Type", "F Vector");
     fv.set("Norm Type", "Two Norm");

     divergence.set("Test Type", "Divergence");
     divergence.set("Tolerance", 1.0e+20);
     divergence.set("Consecutive Iterations", 3);

     stagnation.set("Test Type", "Stagnation");
     stagnation.set("Tolerance", 1.0);
     stagnation.set("Consecutive Iterations", 5);

     maxiters.set("Test Type", "MaxIters");
     maxiters.set("Maximum Iterations", 20);

     Teuchos::RCP<NOX::StatusTest::Generic> status_tests =
       NOX::StatusTest::buildStatusTests(stl, utils);
     \endcode

    \author Nico Schloemer
*/
class Factory {

public:

  //! Constructor.
  Factory();

  //! Destructor.
  virtual ~Factory();

  //! Returns a status test set from a parameter list xml file.
  Teuchos::RCP<LOCA::StatusTest::Abstract>
  buildStatusTests(const std::string& file_name,
                   const Teuchos::RCP<const LOCA::GlobalData> & globalData,
           std::map<std::string, Teuchos::RCP<LOCA::StatusTest::Abstract> >*
           tagged_tests=0) const;

  //! Returns a status test set from a parameter list.
  Teuchos::RCP<LOCA::StatusTest::Abstract>
  buildStatusTests(Teuchos::ParameterList& p,
                   const Teuchos::RCP<const LOCA::GlobalData> & globalData,
           std::map<std::string, Teuchos::RCP<LOCA::StatusTest::Abstract> >*
           tagged_tests=0) const;

private:

  Teuchos::RCP<LOCA::StatusTest::Abstract>
  buildComboTest(Teuchos::ParameterList& p,
                 const Teuchos::RCP<const LOCA::GlobalData> & globalData,
         std::map<std::string, Teuchos::RCP<LOCA::StatusTest::Abstract> >*
         tagged_tests=0) const;

  Teuchos::RCP<LOCA::StatusTest::Abstract>
  buildMaxItersTest(Teuchos::ParameterList& p,
                    const Teuchos::RCP<const LOCA::GlobalData> & globalData ) const;

  //! Checks if a tag is present in the param list and adds the test to the tagged_test std::map if true.  Returns true if a tag was present.
  bool checkAndTagTest(const Teuchos::ParameterList& p,
               const Teuchos::RCP<LOCA::StatusTest::Abstract>& test,
         std::map<std::string, Teuchos::RCP<LOCA::StatusTest::Abstract> >*
               tagged_tests) const;
};

/*! \brief Nonmember helper function for the NOX::StatusTest::Factory.

\relates NOX::StatusTest::Factory

*/
Teuchos::RCP<LOCA::StatusTest::Abstract>
buildStatusTests(const std::string& file_name,
                 const Teuchos::RCP<const LOCA::GlobalData> & globalData,
                 std::map<std::string, Teuchos::RCP<LOCA::StatusTest::Abstract> >*
         tagged_tests=0);

/*! \brief Nonmember helper function for the NOX::StatusTest::Factory.

\relates NOX::StatusTest::Factory

*/
Teuchos::RCP<LOCA::StatusTest::Abstract>
buildStatusTests(Teuchos::ParameterList& p,
                 const Teuchos::RCP<const LOCA::GlobalData> & globalData,
             std::map<std::string, Teuchos::RCP<LOCA::StatusTest::Abstract> >*
         tagged_tests=0);

} // namespace StatusTest
} // namespace NOX

#endif
