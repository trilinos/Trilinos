/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************
// @HEADER
*/

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
