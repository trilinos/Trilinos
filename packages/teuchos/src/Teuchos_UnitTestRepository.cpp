// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER


#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_UnitTestBase.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_StandardCatchMacros.hpp"


namespace Teuchos {


struct UnitTestData {

  Teuchos::UnitTestBase *unitTest;
  std::string groupName;
  std::string testName;

  UnitTestData(Teuchos::UnitTestBase *unitTest_in,
    const std::string groupName_in, const std::string testName_in)
    :unitTest(unitTest_in), groupName(groupName_in), testName(testName_in)
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_ASSERT(unitTest_in);
#endif
    }

private:
  UnitTestData(); // Not defined!
};



bool operator<(const UnitTestData &a, const UnitTestData &b)
{
  if (a.groupName < b.groupName) {
    return true;
  }
  else if (a.groupName > b.groupName) {
    return false;
  }
  return a.testName < b.testName;
}



std::string getUnitTestName(const std::string groupName,
  const std::string testName)
{
  std::ostringstream oss;
  oss << groupName<<"_"<<testName<<"_UnitTest";
  return oss.str();
}


enum EShowTestDetails {
  SHOW_TEST_DETAILS_ALL,
  SHOW_TEST_DETAILS_TEST_NAMES,
  SHOW_TEST_DETAILS_FINAL_RESULTS
};


bool strMatch( const std::string &fullMatchStr, const std::string &str )
{

  const std::string::size_type npos = std::string::npos;

  const int strLen = str.length();
  const int fullMatchStrLen = fullMatchStr.length();

  if (fullMatchStrLen == 0) {
    return true;
  }

  const bool beginGlob = fullMatchStr[0] == '*';
  const bool endGlob = fullMatchStr[fullMatchStrLen-1] == '*';

  const int matchStrLen =
	fullMatchStrLen + (beginGlob ? -1 : 0) + (endGlob ? -1 : 0);

  if (matchStrLen == 0) {
    return true;
  }

  if (matchStrLen > strLen) {
    return false;
  }

  if (beginGlob && endGlob) {
    return str.find(fullMatchStr.substr(1, matchStrLen)) != npos;
  }

  if (endGlob) {
    return fullMatchStr.substr(0, matchStrLen) == str.substr(0, matchStrLen);
  }

  if (beginGlob) {
    return fullMatchStr.substr(1, matchStrLen) ==
      str.substr(strLen-matchStrLen, matchStrLen);
  }

  return fullMatchStr == str;

}


} // namespace Teuchos




namespace Teuchos {


// Implementation class


class UnitTestRepository::InstanceData {
public:

  typedef Teuchos::Array<UnitTestData> unitTests_t;

  unitTests_t unitTests;
  CommandLineProcessor clp;
  EShowTestDetails showTestDetails;
  bool showSrcLocation;
  bool noOp;
  std::string groupName;
  std::string testName;
  std::string notUnitTestName;

  InstanceData()
    :clp(false),
     showTestDetails(SHOW_TEST_DETAILS_TEST_NAMES),
     showSrcLocation(false),
     noOp(false)
    {}

};


// public


CommandLineProcessor& UnitTestRepository::getCLP()
{
  return getData().clp;
}


bool UnitTestRepository::runUnitTests(FancyOStream &out)
{

  typedef InstanceData::unitTests_t unitTests_t;

  out << "\n***\n*** Unit test suite ...\n***\n\n";

  InstanceData &data = getData();

  const bool showAll = data.showTestDetails == SHOW_TEST_DETAILS_ALL;
  const bool showTestNames = data.showTestDetails == SHOW_TEST_DETAILS_TEST_NAMES || showAll;

  bool success = true;
  int testCounter = 0;
  int numTestsRun = 0;
  int numTestsFailed = 0;

  try {
    
    out << "\nSorting tests by group name then by test name ...\n";
    std::sort( data.unitTests.begin(), data.unitTests.end() );

    out << "\nRunning unit tests ...\n\n";
    unitTests_t::iterator iter = data.unitTests.begin();
    for ( ; iter != data.unitTests.end(); ++iter, ++testCounter ) {

      const UnitTestData &utd = (*iter);

      const std::string unitTestName = getUnitTestName(utd.groupName, utd.testName);

      if (
        (
          strMatch(data.groupName, utd.groupName)
          &&
          strMatch(data.testName, utd.testName)
          )
        &&
        (
          data.notUnitTestName.length() == 0
          ||
          !strMatch(data.notUnitTestName, unitTestName)
          )
        )
      {

        ++numTestsRun;

        std::ostringstream testHeaderOSS;
        testHeaderOSS <<testCounter<<". "<<unitTestName<<" ... ";
        const std::string testHeader = testHeaderOSS.str();

        if (showAll)
          out <<"\n";

        if (showTestNames)
          out <<testHeader;

        {

          RCP<std::ostringstream> oss;
          RCP<FancyOStream> localOut;
          if (showAll) {
            out << "\n";
            localOut = rcpFromRef(out);
          }
          else {
            oss = rcp(new std::ostringstream);
            localOut = fancyOStream(rcp_implicit_cast<std::ostream>(oss));
          }

          OSTab tab(out);

          if (!data.noOp) {

            const bool result = utd.unitTest->runUnitTest(*localOut);

            if (!result) {
              
              if (!showTestNames)
                out <<testHeader<<"\n";
              else if (!showAll)
                out <<"\n";
              
              if (!is_null(oss))
                out << oss->str();
              
              out
                << "[FAILED]\n"
                << "Location: "<<utd.unitTest->unitTestFile()<<":"
                <<utd.unitTest->unitTestFileLineNumber()<<"\n";
              
              if (!is_null(oss))
                out << "\n";
              
              success = false;
              
              ++numTestsFailed;
              
            }
            else {
              
              if (showTestNames)
                out << "[Passed]\n";
              
              if (showAll && data.showSrcLocation)
                out
                  << "Location: "<<utd.unitTest->unitTestFile()<<":"
                  <<utd.unitTest->unitTestFileLineNumber()<<"\n";
              
            }

          }
          else {

            if (showTestNames)
              out << "[Not Run]\n";
            
          }

        }
   
      }

    }

    TEUCHOS_ASSERT_EQUALITY(testCounter, as<int>(data.unitTests.size()));

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, out, success);

  out
    << "\nSummary: total = " << testCounter
    << ", run = " << numTestsRun;

  if (!data.noOp) {
    out
      << ", passed = " << (numTestsRun-numTestsFailed)
      << ", failed = " << numTestsFailed << "\n";
  }
  else {
    out
      << ", passed = ???"
      << ", failed = ???\n";
  }
    
  return success;

}


int UnitTestRepository::runUnitTestsFromMain( int argc, char* argv[] )
{

  const RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();

  CommandLineProcessor &clp = getData().clp;
  setUpCLP(outArg(clp));
  CommandLineProcessor::EParseCommandLineReturn parse_return =
    clp.parse(argc,argv);
  if ( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) {
    *out << "\nEnd Result: TEST FAILED" << std::endl;
    return parse_return;
  }


  const bool success = runUnitTests(*out);

  if (success)
    *out << "\nEnd Result: TEST PASSED" << std::endl;
  else
    *out << "\nEnd Result: TEST FAILED" << std::endl;

  return (success ? 0 : 1);

}


void UnitTestRepository::addUnitTest( UnitTestBase *unitTest,
  const std::string groupName, const std::string testName )
{
  getData().unitTests.push_back(UnitTestData(unitTest, groupName, testName));
}


// private:


UnitTestRepository::UnitTestRepository()
{}


void UnitTestRepository::setUpCLP(const Ptr<CommandLineProcessor>& clp)
{

  const int numShowTestDetails = 3;
  const EShowTestDetails showTestDetailsValues[numShowTestDetails] =
    { SHOW_TEST_DETAILS_ALL,
      SHOW_TEST_DETAILS_TEST_NAMES,
      SHOW_TEST_DETAILS_FINAL_RESULTS
    };
  const char* showTestDetailsNames[numShowTestDetails] =
    { "ALL",
      "TEST_NAMES",
      "FINAL_RESULTS"
    };
  clp->setOption(
    "show-test-details", &getData().showTestDetails,
    numShowTestDetails, showTestDetailsValues, showTestDetailsNames,
    "Level of detail to show in the tests"
    );

  clp->setOption(
    "show-src-location", "no-show-src-location", &getData().showSrcLocation,
    "If true, then the location of the unit test source code is shown."
    "  Only meaningfull if --show-test-details=ALL."
    );

  clp->setOption(
    "group-name", &getData().groupName,
    "If specified, selects only tests that match the group name glob." );

  clp->setOption(
    "test-name", &getData().testName,
    "If specified, selects only tests that match the test name glob." );

  clp->setOption(
    "not-unit-test", &getData().notUnitTestName,
    "If specified, full unit tests with glob matches will *not* be run." );

  clp->setOption(
    "no-op", "do-op", &getData().noOp,
    "If --no-op, then only the names of the tests that would be run are run."
    );
  
}


UnitTestRepository::InstanceData& UnitTestRepository::getData()
{
  static UnitTestRepository::InstanceData data;
  return data;
}


} // namespace Teuchos
