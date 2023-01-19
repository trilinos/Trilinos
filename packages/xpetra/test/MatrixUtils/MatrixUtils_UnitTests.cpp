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
#include <Xpetra_UnitTestHelpers.hpp>

#include "Xpetra_DefaultPlatform.hpp"

#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_Map.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_MatrixUtils.hpp"

#include "Xpetra_TpetraCrsMatrix.hpp"
#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraCrsMatrix.hpp"
#endif
namespace {
  using Xpetra::viewLabel_t;

  bool testMpi = true;
  double errorTolSlack = 1e+1;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }

  Teuchos::RCP<const Teuchos::Comm<int> > getDefaultComm()
  {
    if (testMpi) {
      return Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
    }
    return Teuchos::rcp(new Teuchos::SerialComm<int>());
  }


  //
  // UNIT TESTS
  //


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( MatrixUtils, SwitchMatrixToStridedMaps, M, MA, Scalar, LO, GO, Node )
  {
    using CrsMatrixWrap = Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node>;
    using Map = Xpetra::Map<LO,GO,Node>;
    using MapFactory = Xpetra::MapFactory<LO,GO,Node>;
    using MatrixUtils = Xpetra::MatrixUtils<Scalar,LO,GO,Node>;
    using StridedMap = Xpetra::StridedMap<LO,GO,Node>;
    using StridedMapFactory = Xpetra::StridedMapFactory<LO,GO,Node>;

    using Teuchos::Array;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;

    Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();
    M testMap(1,0,comm);
    Xpetra::UnderlyingLib lib = testMap.lib();

    const size_t numLocalElements = 10;
    const GO INVALID = Teuchos::OrdinalTraits<GO>::invalid();
    RCP<const Map> map = MapFactory::createContigMapWithNode(lib, INVALID, numLocalElements, comm);

    std::vector<size_t> stridingInfo;
    stridingInfo.push_back(1);
    RCP<const StridedMap> stridedMap = StridedMapFactory::Build(map, stridingInfo);

    RCP<CrsMatrixWrap> matrix = rcp(new CrsMatrixWrap(stridedMap, stridedMap, 1));
    TEST_EQUALITY_CONST(matrix->GetCurrentViewLabel(), matrix->GetDefaultViewLabel());
    TEST_EQUALITY_CONST(matrix->GetCurrentViewLabel(), matrix->SwitchToView(matrix->GetCurrentViewLabel()));

    for (LO lRow = 0; lRow < Teuchos::as<LO>(numLocalElements); ++lRow) {
      Array<LO> lCols (1, lRow);
      Array<Scalar> values (1, Teuchos::ScalarTraits<Scalar>::one());
      matrix->insertLocalValues(lRow, lCols, values);
    }

    matrix->fillComplete();
    TEST_ASSERT(matrix->isFillComplete());

    TEST_EQUALITY_CONST(matrix->IsView("stridedMaps"), false)
    MatrixUtils::convertMatrixToStridedMaps(matrix, stridingInfo, stridingInfo);
    TEST_EQUALITY_CONST(matrix->IsView("stridedMaps"), true)

    TEST_ASSERT(!matrix->getRangeMap().is_null());
    TEST_ASSERT(!matrix->getDomainMap().is_null());
    TEST_ASSERT(!matrix->getRowMap().is_null());
    TEST_ASSERT(!matrix->getColMap().is_null());

    matrix->SwitchToView("stridedMaps");
    RCP<const StridedMap> stridedRowMap = rcp_dynamic_cast<const StridedMap>(matrix->getRowMap("stridedMaps"));
    RCP<const StridedMap> stridedColMap = rcp_dynamic_cast<const StridedMap>(matrix->getColMap("stridedMaps"));

    TEST_ASSERT(!stridedRowMap.is_null());
    TEST_ASSERT(!stridedColMap.is_null());
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
#define XP_MATRIX_INSTANT(S,LO,GO,N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( MatrixUtils, SwitchMatrixToStridedMaps, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N )


#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XPETRA_TPETRA_TYPES )
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XP_MATRIX_INSTANT )


#if defined(HAVE_XPETRA_EPETRA)

#include "Xpetra_Map.hpp" // defines EpetraNode
typedef Xpetra::EpetraNode EpetraNode;
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
XPETRA_EPETRA_TYPES(double,int,int,EpetraNode)
XP_MATRIX_INSTANT(double,int,int,EpetraNode)
#endif
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
typedef long long LongLong;
XPETRA_EPETRA_TYPES(double,int,LongLong,EpetraNode)
XP_MATRIX_INSTANT(double,int,LongLong,EpetraNode)
#endif

#endif


}
