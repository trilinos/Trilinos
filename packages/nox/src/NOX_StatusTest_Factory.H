// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

    Please use the related nonmember functions instead of calling the
    factory directly (See example below).

    Valid parameters are as follows:

    <ul>
      <li> "Test Type" <std::string> Type of test this list contains.  Valid tests include:
      <ul>
         <li> "Combo" - NOX::StatusTest::Combo
     <li> "NormF" - NOX::StatusTest::NormF
     <li> "NormUpdate" - NOX::StatusTest::NormUpdate
     <li> "NomrWRMS" - NOX::StatusTest::NormWRMS
     <li> "MaxIters" - NOX::StatusTest::MaxIters
     <li> "FiniteValue" - NOX::StatusTest::FiniteValue
     <li> "Divergence" - NOX::StatusTest::Divergence
     <li> "Stagnation" - NOX::StatusTest::Stagnation
     <li> "RelativeNormF" - NOX::StatusTest::RelativeNormF
     <li> "User Defined" - A user constructed test, derived from NOX::StatusTest::Generic.
       </ul>

       <li> "Tag" <std::string> A unique identifier that will place the test in the map for tagged_tests.  This allows users to access individual tests to change tolerances on the fly or query values while still using the factory to build objects.

     </ul>

     Additional parameters valid for a Combo test (NOX::StatusTest::Combo):
     <ul>
       <li> "Combo Type" <std:string> Type of combo to use.  Valid options are:
       <ul>
         <li> "AND"
         <li> "OR"
       </ul>

       <li> "Number of Tests" <int> Number of sublists that contain tests to be added to this combo test.  The sublists must be named "Test X" where "X" represents the test number starting with 0 and preceeding to "Number of Tests - 1".

       <li> "Test X" <Teuchos::ParameterList> A sublist containing a test to add to the current combo test.  The "X" represents the number of the test.  the numbering starts with 0 and is valid through "Number of Tests - 1" tests.
     </ul>

     Additional parameters valid for a Norm F test
     (NOX::StatusTest::NormF):
     <ul>
       <li> "Tolerance" <double> required tolerance for the test to return a converged status. (default = 1.0e-8)
       <li> "Norm Type" <std::string> Type of norm to use. Valid options are:
         <ul>
         <li> "Two Norm" (default)
         <li> "One Norm"
         <li> "Max Norm"
         </ul>
       <li> "Scale Type" <std::string> Type of scaling to use.  Valid options are:
         <ul>
         <li> "Unscaled" (default)
         <li> "Scaled"
         </ul>
       <li> "Initial Guess" < Teuchos::RCP<NOX::Abstract::Group> > If present, a relative tolerance will be used where the group passed in will be used to compute \f$F_0\f$.
     </ul>

     Additional parameters valid for a Norm Update test
     (NOX::StatusTest::NormUpdate):
     <ul>
       <li> "Tolerance" <double> required tolerance for the test to return a converged status. (default = 1.0e-3)
       <li> "Norm Type" <std::string> Type of norm to use. Valid options are:
         <ul>
         <li> "Two Norm" (default)
         <li> "One Norm"
         <li> "Max Norm"
         </ul>
       <li> "Scale Type" <std::string> Type of scaling to use.  Valid options are:
         <ul>
         <li> "Unscaled" (default)
         <li> "Scaled"
         </ul>
     </ul>

     Additional parameters valid for a Norm WRMS test
     (NOX::StatusTest::NormWRMS):
     <ul>
       <li> "Tolerance" <double> (default = 1.0)
       <li> "Relative Tolerance" <double> (default = 1.0e-5)
       <li> "Absolute Tolerance" <double>  or < Teuchos::RCP<const NOX::Abstract::Vector> > (default = 1.0e-8)
       <li> "BDF Multiplier" <double> (default = 1.0)
       <li> "Alpha" <double> (default = 1.0)
       <li> "Beta" <double> (default = 0.5)
     </ul>

     Additional parameters valid for a Maximum Iterations test
     (NOX::StatusTest::MaxIters):
     <ul>
       <li> "Maximum Iterations" <int>
     </ul>

     Additional parameters valid for a Finite Value test
     (NOX::StatusTest::FiniteValue):
     <ul>
       <li> "Vector Type" <std::string> Type of vector to check. Valid options are:
         <ul>
         <li> "F Vector" (default)
         <li> "Solution Vector"
         </ul>
       <li> "Norm Type" <std::string> Type of norm to use. Valid options are:
         <ul>
         <li> "Two Norm" (default)
         <li> "One Norm"
         <li> "Max Norm"
         </ul>
     </ul>

     Additional parameters valid for a Divergence test
     (NOX::StatusTest::Divergence):
     <ul>
       <li> "Tolerance" <double> (default = 1.0e+12)
       <li> "Consecutive Iterations" <int> (default = 1)
     </ul>

     Additional parameters valid for a Stagnation test
     (NOX::StatusTest::Stagnation):
     <ul>
       <li> "Tolerance" <double> (default = 1.0e+12)
       <li> "Consecutive Iterations" <int> (default = 1)
     </ul>

     Additional parameters valid for a RelativeNormF test
     (NOX::StatusTest::RelativeNormF):
     <ul>
       <li> "Tolerance" <double> required tolerance for the test to return a converged status. (default = 1.0e-8)
       <li> "Scale Norms by Length" <bool> if set to true, the norms are scaled by the length of the vector. (default = false)
     </ul>

     Additional parameters valid for a NStep test
     (NOX::StatusTest::NStep):
     <ul>
       <li> "Number of Nonlinear Iterations" <int> Number of nonlinear iterations to take until convergence (default = 1)
       <li> "Number of Initial Ramping Steps" <int> Number of nonlinear solves in an initial ramping phase (default = 0)
       <li> "Number of Nonlinear Iterations in Ramping Phase" <int> Number of nonlinear iterations to take until convergence in ramping phase (default = 10)
     </ul>

     Additional parameters valid for a "User Defined" test:
     <ul>
       <li> "User Status Test" < Teuchos::RCP<NOX::StatusTest::Generic> > A status test suppied by the user.  It is very important that when registering this status test, that the user set it as a "Generic" object since there is no implicit casting on the ParameterList's get method.  See the example below.
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
                   std::map<std::string, Teuchos::RCP<NOX::StatusTest::Generic> >* tagged_tests=nullptr) const;

  //! Returns a status test set from a parameter list.
  Teuchos::RCP<NOX::StatusTest::Generic>
  buildStatusTests(Teuchos::ParameterList& p, const NOX::Utils& utils,
                   std::map<std::string, Teuchos::RCP<NOX::StatusTest::Generic> >* tagged_tests=nullptr) const;

private:

  Teuchos::RCP<NOX::StatusTest::Generic>
  buildComboTest(Teuchos::ParameterList& p,
                 const NOX::Utils& utils,
                 std::map<std::string, Teuchos::RCP<NOX::StatusTest::Generic> >* tagged_tests=nullptr) const;

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

  Teuchos::RCP<NOX::StatusTest::Generic>
  buildRelativeNormFTest(Teuchos::ParameterList& p,
                         const NOX::Utils& utils) const;

  Teuchos::RCP<NOX::StatusTest::Generic>
  buildNStepTest(Teuchos::ParameterList& p,
                 const NOX::Utils& utils) const;

  Teuchos::RCP<NOX::StatusTest::Generic>
  buildUserDefinedTest(Teuchos::ParameterList& p,
                       const NOX::Utils& utils) const;

  //! Checks if a tag is present in the param list and adds the test to the tagged_test std::map if true.  Returns true if a tag was present.
  bool checkAndTagTest(const Teuchos::ParameterList& p,
                       const Teuchos::RCP<NOX::StatusTest::Generic>& test,
                       std::map<std::string, Teuchos::RCP<NOX::StatusTest::Generic> >* tagged_tests) const;
};

/*! \brief Nonmember helper function for the NOX::StatusTest::Factory.

\relates NOX::StatusTest::Factory

*/
Teuchos::RCP<NOX::StatusTest::Generic>
buildStatusTests(const std::string& file_name, const NOX::Utils& utils,
                 std::map<std::string, Teuchos::RCP<NOX::StatusTest::Generic> >* tagged_tests=nullptr);

/*! \brief Nonmember helper function for the NOX::StatusTest::Factory.

\relates NOX::StatusTest::Factory

*/
Teuchos::RCP<NOX::StatusTest::Generic>
buildStatusTests(Teuchos::ParameterList& p, const NOX::Utils& utils,
                 std::map<std::string, Teuchos::RCP<NOX::StatusTest::Generic> >* tagged_tests=nullptr);

} // namespace StatusTest
} // namespace NOX

#endif
