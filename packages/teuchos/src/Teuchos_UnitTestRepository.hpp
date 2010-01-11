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

#ifndef TEUCHOS_UNIT_TEST_REPOSITORY_HPP
#define TEUCHOS_UNIT_TEST_REPOSITORY_HPP


/*! \file Teuchos_UnitTestHarness.hpp
    \brief Unit testing support.
*/


#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_Ptr.hpp"


namespace Teuchos {


class UnitTestBase;


class CommandLineProcessor;


/** \brief Singleton unit testing repository. */
class UnitTestRepository {
public:

  /** \brief . */
  static CommandLineProcessor& getCLP();

  /** \brief Run the registered unit tests */
  static bool runUnitTests(FancyOStream &out);

  /** \brief Run the unit tests from main() passing in (argc, argv).
   *
   * \returns Returns the appropriate int for main()
   */
  static int runUnitTestsFromMain(int argc, char* argv[]);

  /** \brief . */
  static void addUnitTest(UnitTestBase *unitTest, const std::string groupName,
    const std::string testName);

private:

  UnitTestRepository();

  static void setUpCLP(const Ptr<CommandLineProcessor>& clp);

  class InstanceData;

  static InstanceData& getData();

};


} // namespace Teuchos


#endif  // TEUCHOS_UNIT_TEST_REPOSITORY_HPP
