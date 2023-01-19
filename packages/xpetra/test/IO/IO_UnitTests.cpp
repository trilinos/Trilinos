// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_as.hpp>
#include <Xpetra_UnitTestHelpers.hpp>
#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DefaultPlatform.hpp"
#include <Xpetra_IO.hpp>


namespace {

  TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( IO, MMMissingRows, M, MA, Scalar, LO, GO, Node )
  {
    using Teuchos::as;

    // get a comm and node
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
    TEUCHOS_ASSERT_EQUALITY(comm->getSize(), 1);

    if (Teuchos::ScalarTraits<Scalar>::isComplex)
      return;

    M testMap(1,0,comm);
    Xpetra::UnderlyingLib lib = testMap.lib();

    auto A = Xpetra::IO<Scalar,LO,GO,Node>::Read("test.mtx", lib, comm, false);
    TEUCHOS_ASSERT_EQUALITY(A->getGlobalNumRows(), 5);
    TEUCHOS_ASSERT_EQUALITY(A->getGlobalNumCols(), 5);
    TEUCHOS_ASSERT_EQUALITY(A->getGlobalNumEntries(), 3);

    auto colmap = A->getColMap();
    auto crsA = Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node> >(A, true)->getCrsMatrix();
    Teuchos::ArrayView< const LO > indices;
    Teuchos::ArrayView< const Scalar > values;
    crsA->getLocalRowView(0, indices, values);
    TEST_EQUALITY(indices.size(), 2);
    TEST_EQUALITY(colmap->getGlobalElement(indices[0]), 0);
    TEST_EQUALITY(colmap->getGlobalElement(indices[1]), 3);
    TEST_EQUALITY(values[0], as<Scalar>(2.));
    TEST_EQUALITY(values[1], as<Scalar>(3.));

    crsA->getLocalRowView(1, indices, values);
    TEST_EQUALITY(indices.size(), 1);
    TEST_EQUALITY(colmap->getGlobalElement(indices[0]), 4);
    TEST_EQUALITY(values[0], as<Scalar>(4.));

  }

  TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( IO, BinaryMissingRows, M, MA, Scalar, LO, GO, Node )
  {
    using Teuchos::as;

    // get a comm and node
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
    TEUCHOS_ASSERT_EQUALITY(comm->getSize(), 1);

    M testMap(1,0,comm);
    Xpetra::UnderlyingLib lib = testMap.lib();

    auto A = Xpetra::IO<Scalar,LO,GO,Node>::Read("test.mtx.bin", lib, comm, true);
    TEUCHOS_ASSERT_EQUALITY(A->getGlobalNumRows(), 5);
    TEUCHOS_ASSERT_EQUALITY(A->getGlobalNumCols(), 5);
    TEUCHOS_ASSERT_EQUALITY(A->getGlobalNumEntries(), 3);

    auto colmap = A->getColMap();
    auto crsA = Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node> >(A, true)->getCrsMatrix();
    Teuchos::ArrayView< const LO > indices;
    Teuchos::ArrayView< const Scalar > values;
    crsA->getLocalRowView(0, indices, values);
    TEST_EQUALITY(indices.size(), 2);
    TEST_EQUALITY(colmap->getGlobalElement(indices[0]), 0);
    TEST_EQUALITY(colmap->getGlobalElement(indices[1]), 3);
    TEST_EQUALITY(values[0], as<Scalar>(2.));
    TEST_EQUALITY(values[1], as<Scalar>(3.));

    crsA->getLocalRowView(1, indices, values);
    TEST_EQUALITY(indices.size(), 1);
    TEST_EQUALITY(colmap->getGlobalElement(indices[0]), 4);
    TEST_EQUALITY(values[0], as<Scalar>(4.));

  }


//
// INSTANTIATIONS
//

  #define XPETRA_TPETRA_TYPES( S, LO, GO, N) \
    typedef typename Xpetra::TpetraMap<LO,GO,N> M##LO##GO##N; \
    typedef typename Xpetra::TpetraCrsMatrix<S,LO,GO,N> MA##S##LO##GO##N;

#ifdef HAVE_XPETRA_EPETRA

  #define XPETRA_EPETRA_TYPES( S, LO, GO, N) \
    typedef typename Xpetra::EpetraMapT<GO,N> M##LO##GO##N; \
    typedef typename Xpetra::EpetraCrsMatrixT<GO,N> MA##S##LO##GO##N;

#endif


  //list of all tests which run both with Epetra and Tpetra
#define XP_IO_INSTANT(S,LO,GO,N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( IO, MMMissingRows,     M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( IO, BinaryMissingRows, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N )


#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XPETRA_TPETRA_TYPES )
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XP_IO_INSTANT )


#if defined(HAVE_XPETRA_EPETRA)

#include "Xpetra_Map.hpp" // defines EpetraNode
typedef Xpetra::EpetraNode EpetraNode;
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
XPETRA_EPETRA_TYPES(double,int,int,EpetraNode)
XP_IO_INSTANT(double,int,int,EpetraNode)
#endif
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
typedef long long LongLong;
XPETRA_EPETRA_TYPES(double,int,LongLong,EpetraNode)
XP_IO_INSTANT(double,int,LongLong,EpetraNode)
#endif

#endif

}
