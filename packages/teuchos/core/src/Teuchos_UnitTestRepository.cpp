// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_UnitTestBase.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_StandardCatchMacros.hpp"


namespace Teuchos {


struct UnitTestData {

  const Teuchos::UnitTestBase * unitTest;
  std::string groupName;
  std::string testName;
  int insertionIndex;

  UnitTestData(
    Teuchos::UnitTestBase *unitTest_in,
    const std::string groupName_in,
    const std::string testName_in
    )
    : unitTest(unitTest_in), groupName(groupName_in), testName(testName_in),
      insertionIndex(insersionIndexCounter_++)
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_ASSERT(unitTest_in);
#endif
    }

private:
  UnitTestData(); // Not defined!
  static int insersionIndexCounter_;
};


int UnitTestData::insersionIndexCounter_ = 0;


bool operator<(const UnitTestData &a, const UnitTestData &b)
{
  if (a.groupName < b.groupName) {
    return true;
  }
  else if (a.groupName > b.groupName) {
    return false;
  }
  return a.insertionIndex < b.insertionIndex;
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

  const size_t strLen = str.length();
  const size_t fullMatchStrLen = fullMatchStr.length();

  if (fullMatchStrLen == 0) {
    return true;
  }

  const bool beginGlob = fullMatchStr[0] == '*';
  const bool endGlob = fullMatchStr[fullMatchStrLen-1] == '*';

  const size_t matchStrLen =
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
  bool globallyReduceUnitTestResult;
  bool showSrcLocation;
  bool showFailSrcLocation;
  bool noOp;
  std::string groupName;
  std::string testName;
  std::string notUnitTestName;
  int testCounter;

  InstanceData()
    :clp(false),
     showTestDetails(SHOW_TEST_DETAILS_TEST_NAMES),
#if defined(HAVE_TEUCHOS_GLOBALLY_REDUCE_UNITTEST_RESULTS)
     globallyReduceUnitTestResult(true),
#else
     globallyReduceUnitTestResult(false),
#endif
     showSrcLocation(false),
     showFailSrcLocation(true),
     noOp(false),
     testCounter(0)
    {}

};


// public


CommandLineProcessor& UnitTestRepository::getCLP()
{
  return getData().clp;
}


void UnitTestRepository::setGloballyReduceTestResult(
  const bool globallyReduceUnitTestResult)
{
  getData().globallyReduceUnitTestResult = globallyReduceUnitTestResult;
}


bool UnitTestRepository::getGloballyReduceTestResult()
{
  return getData().globallyReduceUnitTestResult;
}


bool UnitTestRepository::runUnitTests(FancyOStream &out)
{

  typedef InstanceData::unitTests_t unitTests_t;

  using std::setprecision;

  Time overallTimer("overallTimer", true);
  Time timer("timer");

  const int timerPrec = 3;

  out << "\n***\n*** Unit test suite ...\n***\n\n";

  InstanceData &data = getData();

  const bool showAll = data.showTestDetails == SHOW_TEST_DETAILS_ALL;
  const bool showTestNames = data.showTestDetails == SHOW_TEST_DETAILS_TEST_NAMES || showAll;

  showTestFailureLocation(data.showFailSrcLocation);

  bool success = true;
  int testCounter = 0;
  int numTestsRun = 0;
  int numTestsFailed = 0;

  Array<std::string> failedTests;

  try {

    out << "\nSorting tests by group name then by the order they were added ...";
    timer.start(true);
    std::sort( data.unitTests.begin(), data.unitTests.end() );
    timer.stop();
    out << " (time = "<<setprecision(timerPrec)<<timer.totalElapsedTime()<<")\n";

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
          out <<testHeader<<std::flush;

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

            timer.start(true);
            const bool result = runUnitTestImpl(*utd.unitTest, *localOut);
            timer.stop();

            if (!result) {

              failedTests.push_back(testHeader);

              if (!showTestNames)
                out <<testHeader<<"\n"<<std::flush;
              else if (!showAll)
                out <<"\n";

              if (!is_null(oss))
                out << oss->str();

              out
                <<"[FAILED] "
                <<" "<<setprecision(timerPrec)<<"("<<timer.totalElapsedTime()<< " sec)"
                <<" "<<unitTestName<<"\n"
                <<"Location: "<<utd.unitTest->unitTestFile()<<":"
                <<utd.unitTest->unitTestFileLineNumber()<<"\n";

              if (!is_null(oss))
                out << "\n";

              success = false;

              ++numTestsFailed;

            }
            else {

              if (showTestNames)
                out << "[Passed] "
                    << setprecision(timerPrec)<<"("<<timer.totalElapsedTime()<<" sec)\n";

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

  if (failedTests.size()) {
    out << "\nThe following tests FAILED:\n";
    for (Teuchos_Ordinal i = 0; i < failedTests.size(); ++i)
      out << "    " << failedTests[i] << "\n";
  }

  overallTimer.stop();
  out << "\nTotal Time: " << setprecision(timerPrec)
      << overallTimer.totalElapsedTime() << " sec\n";

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

  clp.printFinalTimerSummary(out.ptr());

  return (success ? 0 : 1);

}


void UnitTestRepository::addUnitTest( UnitTestBase *unitTest,
  const std::string groupName, const std::string testName_in )
{
  InstanceData &data = getData();
  std::string testName = testName_in;
  data.unitTests.push_back(UnitTestData(unitTest, groupName, testName));
}


bool UnitTestRepository::verboseUnitTests()
{
  return (getData().showTestDetails == SHOW_TEST_DETAILS_ALL);
}


// private:


UnitTestRepository::UnitTestRepository()
{}


void UnitTestRepository::setUpCLP(const Ptr<CommandLineProcessor>& clp)
{

  clp->addOutputSetupOptions(true);

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
    "details", &getData().showTestDetails,
    numShowTestDetails, showTestDetailsValues, showTestDetailsNames,
    "Short for --show-test-details"
    );

  clp->setOption(
    "show-src-location", "no-show-src-location", &getData().showSrcLocation,
    "If true, then the location of the unit test source code is shown."
    "  Only meaningful if --show-test-details=ALL."
    );

  clp->setOption(
    "show-fail-src-location", "no-show-fail-src-location", &getData().showFailSrcLocation,
    "If true, then the location of every failed unit test check is printed."
    );

  clp->setOption(
    "globally-reduce-test-result", "no-globally-reduce-test-result",
    &getData().globallyReduceUnitTestResult,
    "If true, individual unit test pass/fail is globally reduced across MPI processes."
    );

  clp->setOption(
    "group-name", &getData().groupName,
    "If specified, selects only tests that match the group name glob." );
  clp->setOption(
    "group", &getData().groupName,
    "Short for --group-name." );

  clp->setOption(
    "test-name", &getData().testName,
    "If specified, selects only tests that match the test name glob." );
  clp->setOption(
    "test", &getData().testName,
    "Short for --test-name." );

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


bool UnitTestRepository::runUnitTestImpl(const UnitTestBase &unitTest,
  FancyOStream &out)
{
  const bool result = unitTest.runUnitTest(out);
  if (getData().globallyReduceUnitTestResult) {
    const int globalSum = GlobalMPISession::sum(result ? 0 : 1);
    if (globalSum == 0) {
      return true;
    }
    else {
      // Only print that there are failures on processes where the local
      // unit test actally passed.  On processes where the local unit test
      // fails, users already know that test failed so there is no need to
      // exlain it.
      if (result) {
        out << "NOTE: Global reduction shows failures on other processes!\n"
            << "(rerun with --output-to-root-rank-only=-1 to see output\n"
            << "from other processes to see what process failed!)\n";
      }
      else {
        // The test failed on the root process so the user already knows it failed!
      }
      // Determine what processes have failing tests
      const int numProcs = GlobalMPISession::getNProc();
      Array<int> passFailFlags(numProcs);
      GlobalMPISession::allGather( result ? 0 : 1, passFailFlags());
      Array<int> procsThatFailed;
      for ( int proc_k = 0; proc_k < numProcs; ++proc_k ) {
        if (passFailFlags[proc_k] != 0) {
          procsThatFailed.push_back(proc_k);
        }
      }
      // Print what processes have the failing tests.  If there is only one
      // processes, don't print anything.
      if (numProcs > 1) {
        if (procsThatFailed.size() == numProcs) {
          out << "NOTE: Unit test failed on all processes!\n";
          // NOTE: when all the processes are failing it is useless to print
          // out a list of all of the processes.
        }
        else {
          out << "NOTE: Unit test failed on processes = " << procsThatFailed << "\n"
              << "(rerun with --output-to-root-rank-only=<procID> to see output\n"
              << "from individual processes where the unit test is failing!)\n";
        }
      }
      return false;
    }
  }
  return result;
}


} // namespace Teuchos
