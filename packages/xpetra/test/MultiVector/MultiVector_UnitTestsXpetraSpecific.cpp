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
#include <Teuchos_Comm.hpp>

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DefaultPlatform.hpp"

#include "Xpetra_MultiVector.hpp"
#include "Xpetra_Vector.hpp"

#include "Xpetra_MapFactory.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include "Xpetra_MapExtractor.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraMultiVector.hpp"
#include "Xpetra_TpetraVector.hpp"
#endif

#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraMultiVector.hpp"
#include "Xpetra_EpetraVector.hpp"
#endif

namespace {

  TEUCHOS_STATIC_SETUP()
  {
  }

  Teuchos::RCP<const Teuchos::Comm<int> > getDefaultComm()
  {
    Teuchos::RCP<const Teuchos::Comm<int> > ret;
    ret = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
    return ret;
  }

  // Get an instance of the given Kokkos Node type.
  //
  // \warning This function is NOT reentrant, and therefore NOT thread safe.
  template <class Node>
  Teuchos::RCP<Node> getNode () {
    Teuchos::ParameterList pl;
    return Teuchos::rcp (new Node (pl));
  }

  // Test getVector() / getVectorNonConst()
  // More specifically, this test verifies that the newly created vector will remain valid after the disappearance of the references to the multivector in user code.
  TEUCHOS_UNIT_TEST_TEMPLATE_7_DECL( MultiVector, XpetraSpecific_GetVector, M, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node )
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    EXTRACT_LIB(comm,M) // returns mylib

    const Xpetra::global_size_t INVALID = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();

    const size_t numLocal = 4;

    RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::createContigMap(mylib, INVALID, numLocal, comm);
    RCP< Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > mv = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(map, 3, false);
    for(size_t k=0; k < 3; k++) {
      Teuchos::ArrayRCP<Scalar> mvData = mv->getDataNonConst(k);

      for(size_t i=0; i < numLocal; i++) {
        mvData[i] = Teuchos::as<Scalar>(i*(k+1));
      }
    }

    Teuchos::RCP< const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > v         = mv->getVector(1);         // second vector
    Teuchos::RCP<       Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > vNonConst = mv->getVectorNonConst(2); // third vector

    mv = Teuchos::null;

    {
      Teuchos::ArrayRCP<const Scalar> vData = v->getData(0);
      for(size_t i=0; i< numLocal; i++) {
        TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(2*i));
      }
    }

    {
      Teuchos::ArrayRCP<Scalar> vData = vNonConst->getDataNonConst(0);
      for(size_t i=0; i< numLocal; i++) {
        TEST_EQUALITY(vData[i], Teuchos::as<Scalar>(3*i));
      }
    }
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_7_DECL( MultiVector, XpetraSpecific_GetHostLocalView, M, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node )
  {
#ifdef HAVE_XPETRA_KOKKOS_REFACTOR
    using Teuchos::RCP;
    using Teuchos::rcp;

    RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    EXTRACT_LIB(comm,M) // returns mylib

    typedef typename Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> mv_type;
    typedef typename mv_type::dual_view_type dual_view_type;

    const Xpetra::global_size_t INVALID = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
    const size_t numLocal = 4;

    RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::createContigMap(mylib, INVALID, numLocal, comm);

    // create new vector and fill it with data
    RCP< mv_type > mv = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(map, 3, false);
    for(size_t k=0; k < 3; k++) {
      Teuchos::ArrayRCP<Scalar> mvData = mv->getDataNonConst(k);

      for(size_t i=0; i < numLocal; i++) {
        mvData[i] = i*(k+1) + comm->getRank();
      }
    }

    // get a view of the multivector data on the host memory
    typename dual_view_type::t_host_um hostView = mv->getHostLocalView ();

    TEST_EQUALITY(hostView.extent(0), numLocal);
    TEST_EQUALITY(hostView.extent(1), 3);
    TEST_EQUALITY(hostView.size(), numLocal * 3);
    for(size_t k=0; k < 3; k++) {
      for(size_t i = 0; i < hostView.extent(0); i++) {
        TEST_EQUALITY(Teuchos::as<Scalar>(hostView(i,k)), Teuchos::as<Scalar>(i*(k+1) + comm->getRank()));
      }
    }

    // overwrite data in hostView
    for(size_t r = 0; r < hostView.extent(0); r++) {
      for(size_t c = 0; c < hostView.extent(1); c++) {
        hostView(r,c) = comm->getRank() + c*hostView.extent(1) + r + 42.0;
      }
    }

    // check data in multivector

    for(size_t k=0; k < 3; k++) {
      Teuchos::ArrayRCP<const Scalar> vData = mv->getData(k);
      for(size_t i=0; i< numLocal; i++) {
        TEST_EQUALITY(Teuchos::as<Scalar>(vData[i]), Teuchos::as<Scalar>(comm->getRank() + k*hostView.extent(1) + i + 42.0));
      }
    }

    // change data in multivector
    for(size_t k=0; k < 3; k++) {
      Teuchos::ArrayRCP<Scalar> vData = mv->getDataNonConst(k);
      for(size_t i=0; i< numLocal; i++) {
        vData[i] = k * numLocal + i;
      }
    }

    // check updated data in view
    for(size_t r = 0; r < hostView.extent(0); r++) {
      for(size_t c = 0; c < hostView.extent(1); c++) {
        TEST_EQUALITY(Teuchos::as<Scalar>(hostView(r,c)), Teuchos::as<Scalar>(c * numLocal + r));
      }
    }

    // delete vector
    mv = Teuchos::null;
#endif // HAVE_XPETRA_KOKKOS_REFACTOR
  }



  //
  // INSTANTIATIONS
  //
  #ifdef HAVE_XPETRA_TPETRA

    #define XPETRA_TPETRA_TYPES( S, LO, GO, N) \
      typedef typename Xpetra::TpetraMap<LO,GO,N> M##LO##GO##N; \
      typedef typename Xpetra::TpetraMultiVector<S,LO,GO,N> MV##S##LO##GO##N; \
      typedef typename Xpetra::TpetraVector<S,LO,GO,N> V##S##LO##GO##N;       \

  #endif

  #ifdef HAVE_XPETRA_EPETRA

    #define XPETRA_EPETRA_TYPES( S, LO, GO, N) \
      typedef typename Xpetra::EpetraMapT<GO,N> M##LO##GO##N; \
      typedef typename Xpetra::EpetraMultiVectorT<GO,N> MV##S##LO##GO##N; \
      typedef typename Xpetra::EpetraVectorT<GO,N> V##S##LO##GO##N;       \

  #endif

  // List of tests which run only with Tpetra
  #define XP_MULTIVECTOR_INSTANT(S,LO,GO,N) \
        TEUCHOS_UNIT_TEST_TEMPLATE_7_INSTANT( MultiVector, XpetraSpecific_GetHostLocalView , M##LO##GO##N , MV##S##LO##GO##N , V##S##LO##GO##N , S, LO, GO, N ) \
        TEUCHOS_UNIT_TEST_TEMPLATE_7_INSTANT( MultiVector, XpetraSpecific_GetVector        , M##LO##GO##N , MV##S##LO##GO##N , V##S##LO##GO##N , S, LO, GO, N ) \

#if defined(HAVE_XPETRA_TPETRA)

#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
// no ordinal types as scalar for testing as some tests use ScalarTraits::eps...
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XPETRA_TPETRA_TYPES )
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XP_MULTIVECTOR_INSTANT )

#endif


#if defined(HAVE_XPETRA_EPETRA)

#include "Xpetra_Map.hpp" // defines EpetraNode
typedef Xpetra::EpetraNode EpetraNode;
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
XPETRA_EPETRA_TYPES(double,int,int,EpetraNode)
XP_MULTIVECTOR_INSTANT(double,int,int,EpetraNode)
#endif
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
typedef long long LongLong;
XPETRA_EPETRA_TYPES(double,int,LongLong,EpetraNode)
XP_MULTIVECTOR_INSTANT(double,int,LongLong,EpetraNode)
#endif
#endif

} // namespace
