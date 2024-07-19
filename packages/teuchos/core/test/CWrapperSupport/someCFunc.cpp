// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_CWrapperSupport_Cpp.hpp"
#include "Teuchos_Assert.hpp"


extern "C" {


int someCFunc(int input, int *ierr)
{
  int output = -1;
  TEUCHOS_CWRAPPER_TRY(ierr) {
    TEUCHOS_ASSERT_INEQUALITY(input, >=, 0);
    if (input > 10) {
      TEUCHOS_CWRAPPER_SET_ERROR_CODE(ierr, -2);
    }
    else {
      output = input;
    }
  } TEUCHOS_CWRAPPER_CATCH_ERROR_CODE(ierr);
  return output;
}


} // extern "C"
