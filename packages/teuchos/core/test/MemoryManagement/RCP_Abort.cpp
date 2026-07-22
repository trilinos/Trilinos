// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_RCP.hpp"

// This test will abort for TEUCHOS_DEBUG only. RCPNode will detect an
// exception in the destructor of A and then execute the abort because handling
// the exception would not be thread safe.
int main(int argc, char* argv[])
{
  // create a test class which throws on destruct
  class A {
    public:
      ~A() TEUCHOS_NOEXCEPT_FALSE {
       throw std::runtime_error( "Test Class A throws on destructor..." );
      }
  };

  // create RCP of the test class
  Teuchos::RCP<A> a = Teuchos::rcp(new A);

  // release RCP - will trigger the exception from the destructor of A
  a = Teuchos::null;

  // we should not get here because this test should only be running
  // in debug mode and it should have aborted already.
  std::cout << "Unexpected RCP_Abort should have aborted" << std::endl;
  return 0;
}