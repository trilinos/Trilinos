About unit tests in Zoltan2
===========================

Algorithms and Problems are tested with system tests rather 
than unit tests.  System tests are in directories like "partitioning" 
and "ordering".  All other classes should be tested with unit tests.

The names of the unit test subdirectories mirror the names of
the source file subdirectories.  The name of each class unit test
is the name of the class it is testing.  

A certain structure and a minimum of Doxygen markup is needed 
in your test in order to document the unit tests.

Suppose I'm writing the unit test models/SomeModel.cpp which tests the
class SomeModel which is defined in src/models/Zoltan2_SomeModel.hpp.

It would look like this:

////////////////////////////////////////////////////////////
// Standard copyright message here

/*! \file SomeModel.cpp
 *
 *  \brief Tests the SomeModel class.
 *
 *  \todo Here I describe some part of the test that still needs to be done.
 */

#include <Zoltan2_SomeModel.hpp>
#include <Zoltan2_TestHelpers.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

void SomeModelTest(Teuchos::RCP<const Teuchos::Comm<int> > &comm)
{
  some code

  /*! \test feature 1 of SomeModel
   */

  some code // If an error, call TEST_FAIL_AND_EXIT

  /*! \test feature 2 of SomeModel
   */

  some code // If an error, call TEST_FAIL_AND_EXIT

  return;
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

  SomeModelTest(comm);

  if (comm->getRank() > 0)
    std::cout << "PASS" << std::endl;

  return 0;
}

//
////////////////////////////////////////////////////////////

Zoltan2_TestHelpers.hpp has error handling macros and methods.  It also 
includes a class that can generate user input for tests.  And it
typedefs gno_t, lno_t and scalar_t based on whether Tpetra was
explicitly instantiated with LocalOrdinal, GlobalOrdinal and Scalar
data types.

Naming the test after the class that is being tested, and including
a brief description of what is being tested with "\test", results
in helpful test documentation when doxygen is run.
