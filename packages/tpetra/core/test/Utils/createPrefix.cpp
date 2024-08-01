// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Util.hpp"
#include <sstream>

namespace { // (anonymous)

  TEUCHOS_UNIT_TEST( TpetraUtil, CreatePrefix )
  {
    using Tpetra::TestingUtilities::getDefaultComm;
    using Tpetra::Details::createPrefix;
    using std::endl;

    out << "Test Tpetra::Details::createPrefix" << endl;
    Teuchos::OSTab tab1 (out);

    auto comm = getDefaultComm();
    TEST_ASSERT( comm.getRawPtr() != nullptr );
    const int myRank = comm->getRank();

    {
      // This overload takes any integer and any string.
      std::unique_ptr<std::string> pfx =
        createPrefix(418, "HI FOLKS");
      TEST_ASSERT( pfx.get() != nullptr );
      std::ostringstream os;
      os << "Proc 418: HI FOLKS: ";
      TEST_EQUALITY( *pfx, os.str() );
    }

    {
      std::unique_ptr<std::string> pfx =
        createPrefix(comm.getRawPtr(), "myFunc");
      TEST_ASSERT( pfx.get() != nullptr );
      std::ostringstream os;
      os << "Proc " << myRank << ": Tpetra::myFunc: ";
      TEST_EQUALITY( *pfx, os.str() );
    }

    {
      std::unique_ptr<std::string> pfx =
        createPrefix(nullptr, "myFunc");
      TEST_ASSERT( pfx.get() != nullptr );
      std::ostringstream os;
      os << "Proc -1: Tpetra::myFunc: ";
      TEST_EQUALITY( *pfx, os.str() );
    }

    {
      std::unique_ptr<std::string> pfx =
        createPrefix(comm.getRawPtr(), "MyClass", "myMethod");
      TEST_ASSERT( pfx.get() != nullptr );
      std::ostringstream os;
      os << "Proc " << myRank << ": Tpetra::MyClass::myMethod: ";
      TEST_EQUALITY( *pfx, os.str() );
    }

    {
      std::unique_ptr<std::string> pfx =
        createPrefix(nullptr, "MyClass", "myMethod");
      TEST_ASSERT( pfx.get() != nullptr );
      std::ostringstream os;
      os << "Proc -1: Tpetra::MyClass::myMethod: ";
      TEST_EQUALITY( *pfx, os.str() );
    }
  }

} // namespace (anonymous)
