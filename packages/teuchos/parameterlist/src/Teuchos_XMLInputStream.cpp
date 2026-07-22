// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_XMLInputStream.hpp"
#include "Teuchos_Assert.hpp"

using namespace Teuchos;


unsigned int XMLInputStream::curPos() const
{
  // NOTE (mfh 15 Sep 2014): Most compilers have figured out that the
  // return statement below is unreachable.  Some older compilers
  // might not realize this.  That's why the return statement was put
  // there, so that those compilers don't warn that this function
  // doesn't return a value.  If it's a choice between one warning and
  // another, I would prefer the choice that produces less code and
  // doesn't have unreachable code (which never gets tested).

  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
    "XMLInputStream::curPos() should never be called. It exists only for "
    "compatibility with Xerces.");
  // return 0; // unreachable
}
