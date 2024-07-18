// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_UNIT_TEST_REPOSITORY_HPP
#define TEUCHOS_UNIT_TEST_REPOSITORY_HPP


/*! \file Teuchos_UnitTestRepository.hpp
    \brief Unit testing support.
*/


#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_Ptr.hpp"


namespace Teuchos {


class UnitTestBase;


class CommandLineProcessor;


/** \brief Singleton unit testing repository.
 *
 * This class is the universal driver for unit testing.  This class should be
 * almost completely invisible to a user of the test harness.  The main
 * interaction is through command-line arguments set and processed by calling
 * runUnitTestsFromMain() from a main() function.  See --help for details of
 * the options.  For a more general overview, see \ref Teuchos_UnitTest_grp.
 */
class TEUCHOSCORE_LIB_DLL_EXPORT UnitTestRepository {
public:

  /** \brief Return the CLP to add options to. */
  static CommandLineProcessor& getCLP();

  /** \brief Set if the unit tests should reduce pass/fail across processes. */
  static void setGloballyReduceTestResult(const bool globallyReduceUnitTestResult);

  /** \brief Get if the unit tests should reduce across processes or not. */
  static bool getGloballyReduceTestResult();

  /** \brief Run the registered unit tests */
  static bool runUnitTests(FancyOStream &out);

  /** \brief Run the unit tests from main() passing in (argc, argv).
   *
   * \returns Returns the appropriate int for main()
   */
  static int runUnitTestsFromMain(int argc, char* argv[]);

  /** \brief Add an unit test (called indirectly through macros.
   *
   * unittest [in] The unit test.  NOTE: the memory of *unittest must be persistant.
   */
  static void addUnitTest(UnitTestBase *unitTest, const std::string groupName,
    const std::string testName);

  /** \brief Returns if unit tests are verbose or not.
   *
   * This can be used in individual unit tests that need to know if the unit
   * test harness is running in verbose mode.  This is useful when the unit
   * test's std::ostream 'out' can not be printed to (for example, when
   * Fortran is testing is running).
   */
  static bool verboseUnitTests();

private:

  UnitTestRepository();

  static void setUpCLP(const Ptr<CommandLineProcessor>& clp);

  class InstanceData;

  static InstanceData& getData();

  static bool runUnitTestImpl(const UnitTestBase &unitTest, FancyOStream &out);

};


} // namespace Teuchos


#endif  // TEUCHOS_UNIT_TEST_REPOSITORY_HPP
