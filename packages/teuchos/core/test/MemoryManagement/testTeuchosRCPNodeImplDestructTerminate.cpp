// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCPNode.hpp"

int main(int argc, char* argv[])
{

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  {
    int obj = 1;
    using Dealloc_T = Teuchos::DeallocDelete<int>;
    Teuchos::RCPNodeTmpl<int, Teuchos::DeallocDelete<int> >
      rcpNode(&obj, Dealloc_T(), false);
  }
  // When the above destructor executes it should terminate the program with
  // an error message!

  return 1; // Will never be called!

}
