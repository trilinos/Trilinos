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

#include "Xpetra_Map.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraCrsMatrix.hpp"
#endif
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

  // Get an instance of the given Kokkos Node type.
  //
  // \warning This function is NOT reentrant, and therefore NOT thread safe.
  template <class Node>
  Teuchos::RCP<Node> getNode () {
    //static Teuchos::RCP<Node> node_ = Teuchos::null; // Defaults to null
    //if (node_.is_null ()) {
    //  Teuchos::ParameterList pl;
    //  pl.set<int> ("Num Threads", 0);
    //  pl.set<int> ("Verbose", 1);
    //  node_ = Teuchos::rcp (new Node (pl));
    //}
    //return node_;
    Teuchos::ParameterList pl;
    return Teuchos::rcp (new Node (pl));
  }

  //
  // UNIT TESTS
  //


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( Matrix, ViewSwitching, M, MA, Scalar, LO, GO, Node )
  {
#ifdef HAVE_XPETRA_TPETRA
    typedef Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node> CrsMatrixWrap;
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

    const size_t numLocal = 10;
    const size_t INVALID = Teuchos::OrdinalTraits<size_t>::invalid(); // TODO: global_size_t instead of size_t
    //Teuchos::RCP<const Xpetra::Map<LO,GO,Node> > map = Xpetra::useTpetra::createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    Teuchos::RCP<const Xpetra::Map<LO,GO,Node> > map =
        Xpetra::MapFactory<LO,GO,Node>::createContigMapWithNode (Xpetra::UseTpetra,INVALID,numLocal,comm);
     {
       Xpetra::TpetraCrsMatrix<Scalar, LO, GO, Node> t =  Xpetra::TpetraCrsMatrix<Scalar,LO,GO,Node>(map, numLocal);

       // Test of constructor
       CrsMatrixWrap op(map,1);
       TEST_EQUALITY_CONST(op.GetCurrentViewLabel(), op.GetDefaultViewLabel());
       TEST_EQUALITY_CONST(op.GetCurrentViewLabel(), op.SwitchToView(op.GetCurrentViewLabel()));

       // Test of CreateView
       TEST_THROW(op.CreateView(op.GetDefaultViewLabel(),op.getRowMap(),op.getColMap()), Xpetra::Exceptions::RuntimeError); // a
       op.CreateView("newView",op.getRowMap(),op.getColMap());                                                              // b
       TEST_THROW(op.CreateView("newView",op.getRowMap(),op.getColMap()), Xpetra::Exceptions::RuntimeError);                // c

       // Test of SwitchToView
       // a
       viewLabel_t viewLabel    = op.GetCurrentViewLabel();
       viewLabel_t oldViewLabel = op.SwitchToView("newView");
       TEST_EQUALITY_CONST(viewLabel, oldViewLabel);
       TEST_EQUALITY_CONST(op.GetCurrentViewLabel(), "newView");
       // b
       TEST_THROW(op.SwitchToView("notAView"), Xpetra::Exceptions::RuntimeError);

       // Test of SwitchToDefaultView()
       // a
       viewLabel    = op.GetCurrentViewLabel();
       oldViewLabel = op.SwitchToDefaultView();
       TEST_EQUALITY_CONST(viewLabel, oldViewLabel);
       TEST_EQUALITY_CONST(op.GetCurrentViewLabel(), op.GetDefaultViewLabel());

       // Test of RemoveView()
       TEST_THROW(op.RemoveView(op.GetDefaultViewLabel()), Xpetra::Exceptions::RuntimeError); // a
       TEST_THROW(op.RemoveView("notAView"), Xpetra::Exceptions::RuntimeError);               // b
       op.RemoveView("newView");                                                               // c
       TEST_THROW(op.RemoveView("newView"), Xpetra::Exceptions::RuntimeError);

       op.fillComplete();
     }
#endif
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( Matrix, StridedMaps_Tpetra, M, MA, Scalar, LO, GO, Node )
  {
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    const size_t numLocal = 10;
    const size_t INVALID = Teuchos::OrdinalTraits<size_t>::invalid(); // TODO: global_size_t instead of size_t


#ifdef HAVE_XPETRA_TPETRA
    typedef Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node> CrsMatrixWrap;
    Teuchos::RCP<const Xpetra::Map<LO,GO,Node> > map =
        Xpetra::MapFactory<LO,GO,Node>::createContigMapWithNode (Xpetra::UseTpetra,INVALID,numLocal,comm);
     {
       Xpetra::TpetraCrsMatrix<Scalar, LO, GO, Node> t =  Xpetra::TpetraCrsMatrix<Scalar,LO,GO,Node>(map, numLocal);

       // Test of constructor
       CrsMatrixWrap op(map,1);
       op.fillComplete();

       TEST_EQUALITY_CONST(op.GetCurrentViewLabel(), op.GetDefaultViewLabel());
       TEST_EQUALITY_CONST(op.GetCurrentViewLabel(), op.SwitchToView(op.GetCurrentViewLabel()));

       op.SetFixedBlockSize(2);
       TEST_EQUALITY(op.IsView("stridedMaps"),true);
       TEST_EQUALITY(op.IsView("StridedMaps"),false);
       int blkSize = op.GetFixedBlockSize();
       TEST_EQUALITY_CONST(blkSize, 2);
     }
#endif
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( Matrix, StridedMaps_Epetra, M, MA, Scalar, LO, GO, Node )
  {
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    const size_t numLocal = 10;
    const size_t INVALID = Teuchos::OrdinalTraits<size_t>::invalid();

#ifdef HAVE_XPETRA_EPETRA
    typedef Xpetra::CrsMatrixWrap<double, int, GO, Xpetra::EpetraNode> EpCrsMatrix;

    Teuchos::RCP<const Xpetra::Map<int,GO,Node> > epmap = Xpetra::MapFactory<int,GO,Node>::createContigMap(Xpetra::UseEpetra, INVALID, numLocal, comm);
     {
       Xpetra::EpetraCrsMatrixT<GO,Node> t =  Xpetra::EpetraCrsMatrixT<GO,Node>(epmap, numLocal);

       // Test of constructor
       EpCrsMatrix op(epmap,1);
       op.fillComplete();

       TEST_EQUALITY_CONST(op.GetCurrentViewLabel(), op.GetDefaultViewLabel());
       TEST_EQUALITY_CONST(op.GetCurrentViewLabel(), op.SwitchToView(op.GetCurrentViewLabel()));

       op.SetFixedBlockSize(2);
       TEST_EQUALITY(op.IsView("stridedMaps"),true);
       TEST_EQUALITY(op.IsView("StridedMaps"),false);
       int blkSize = op.GetFixedBlockSize();
       TEST_EQUALITY_CONST(blkSize, 2);
     }
#endif
  }

//
// INSTANTIATIONS
//
#ifdef HAVE_XPETRA_TPETRA

  #define XPETRA_TPETRA_TYPES( S, LO, GO, N) \
    typedef typename Xpetra::TpetraMap<LO,GO,N> M##LO##GO##N; \
    typedef typename Xpetra::TpetraCrsMatrix<S,LO,GO,N> MA##S##LO##GO##N;

#endif

#ifdef HAVE_XPETRA_EPETRA

  #define XPETRA_EPETRA_TYPES( S, LO, GO, N) \
    typedef typename Xpetra::EpetraMapT<GO,N> M##LO##GO##N; \
    typedef typename Xpetra::EpetraCrsMatrixT<GO,N> MA##S##LO##GO##N;

#endif

// List of tests which run only with Tpetra
#define XP_TPETRA_MATRIX_INSTANT(S,LO,GO,N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( Matrix, StridedMaps_Tpetra,  M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( Matrix, ViewSwitching, M##LO##GO##N , MA##S##LO##GO##N , S, LO, GO, N )

// List of tests which run only with Epetra
#define XP_EPETRA_MATRIX_INSTANT(S,LO,GO,N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( Matrix, StridedMaps_Epetra, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N )

// list of all tests which run both with Epetra and Tpetra
//#define XP_MATRIX_INSTANT(S,LO,GO,N)




#if defined(HAVE_XPETRA_TPETRA)

#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XPETRA_TPETRA_TYPES )
//TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XP_MATRIX_INSTANT )
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XP_TPETRA_MATRIX_INSTANT )

#endif


#if defined(HAVE_XPETRA_EPETRA)

#include "Xpetra_Map.hpp" // defines EpetraNode
typedef Xpetra::EpetraNode EpetraNode;
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
XPETRA_EPETRA_TYPES(double,int,int,EpetraNode)
//XP_MATRIX_INSTANT(double,int,int,EpetraNode)
XP_EPETRA_MATRIX_INSTANT(double,int,int,EpetraNode)
#endif
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
typedef long long LongLong;
XPETRA_EPETRA_TYPES(double,int,LongLong,EpetraNode)
//XP_MATRIX_INSTANT(double,int,LongLong,EpetraNode)
XP_EPETRA_MATRIX_INSTANT(double,int,LongLong,EpetraNode)
#endif

#endif


}
