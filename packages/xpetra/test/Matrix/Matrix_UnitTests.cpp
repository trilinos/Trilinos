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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include "Xpetra_ConfigDefs.hpp" //TODO
#include "Xpetra_DefaultPlatform.hpp" //TODO
#include "Teuchos_as.hpp"

//#include "Xpetra_ConfigDefs.hpp"
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
  using Teuchos::Array;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::arcp;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Teuchos::Tuple;
  using Teuchos::tuple;
  using std::sort;
  using std::find;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::Comm;

  using Xpetra::DefaultPlatform;
  using Xpetra::Matrix;
  using Xpetra::CrsMatrixWrap;
#ifdef HAVE_XPETRA_TPETRA
  using Xpetra::TpetraCrsMatrix; //TMP
#endif
#ifdef HAVE_XPETRA_EPETRA
  using Xpetra::EpetraCrsMatrix; //TMP
#endif

  using Xpetra::Map;

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

  RCP<const Comm<int> > getDefaultComm()
  {
    if (testMpi) {
      return DefaultPlatform::getDefaultPlatform().getComm();
    }
    return rcp(new Teuchos::SerialComm<int>());
  }

  //
  // UNIT TESTS
  //


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( Matrix, ViewSwitching, Scalar, LO, GO, Node ) //TODO: add template parameter <Node,...>
  {
#ifdef HAVE_XPETRA_TPETRA
    typedef Teuchos::ScalarTraits<Scalar> ST;
    typedef Matrix<Scalar, LO, GO, Node> Matrix;
    typedef CrsMatrixWrap<Scalar, LO, GO, Node> CrsMatrixWrap;
    RCP<const Comm<int> > comm = getDefaultComm();

    const size_t numLocal = 10;
    const size_t INVALID = OrdinalTraits<size_t>::invalid(); // TODO: global_size_t instead of size_t
    RCP<const Map<LO,GO,Node> > map = Xpetra::useTpetra::createContigMap<LO,GO>(INVALID,numLocal,comm);
     {
       TpetraCrsMatrix<Scalar, LO, GO, Node> t =  TpetraCrsMatrix<Scalar,LO,GO,Node>(map, numLocal);

       // Test of constructor
       CrsMatrixWrap op(map,1);
       TEST_EQUALITY_CONST(op.GetCurrentViewLabel(), op.GetDefaultViewLabel());
       TEST_EQUALITY_CONST(op.GetCurrentViewLabel(), op.SwitchToView(op.GetCurrentViewLabel()));

       // Test of CreateView
       TEST_THROW(op.CreateView(op.GetDefaultViewLabel(),op.getRowMap(),op.getColMap()), Xpetra::Exceptions::RuntimeError); // a
       op.CreateView("newView",op.getRowMap(),op.getColMap());                                                               // b
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
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( Matrix, StridedMaps, Scalar, LO, GO, Node )
  {
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numLocal = 10;
    const size_t INVALID = OrdinalTraits<size_t>::invalid(); // TODO: global_size_t instead of size_t

#ifdef HAVE_XPETRA_TPETRA
    //typedef Teuchos::ScalarTraits<Scalar> ST;
    //typedef Matrix<Scalar, LO, GO, Node> Matrix;
    typedef CrsMatrixWrap<Scalar, LO, GO, Node> CrsMatrixWrap;


    RCP<const Map<LO,GO,Node> > map = Xpetra::useTpetra::createContigMap<LO,GO>(INVALID,numLocal,comm);
     {
       TpetraCrsMatrix<Scalar, LO, GO, Node> t =  TpetraCrsMatrix<Scalar,LO,GO,Node>(map, numLocal);

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

#ifdef HAVE_XPETRA_EPETRA
    typedef Xpetra::CrsMatrixWrap<double, int, int, Node> EpCrsMatrix;

    RCP<const Map<int,int,Node> > epmap = Xpetra::MapFactory<int,int,Node>::createContigMap(Xpetra::UseEpetra, INVALID, numLocal, comm);
     {
       EpetraCrsMatrix t =  EpetraCrsMatrix(epmap, numLocal);

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

#   define UNIT_TEST_GROUP_ORDINAL( SC, LO, GO, Node )                       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Matrix, ViewSwitching, SC, LO, GO, Node ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Matrix, StridedMaps, SC, LO, GO, Node )

  typedef Kokkos::DefaultNode::DefaultNodeType DefaultNodeType;
  UNIT_TEST_GROUP_ORDINAL(double, int, int, DefaultNodeType)

}
