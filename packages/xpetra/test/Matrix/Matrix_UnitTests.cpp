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
#include "Xpetra_MatrixUtils.hpp"
#include "Xpetra_MatrixFactory.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
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
  TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( Matrix, ViewSwitching, M, MA, Scalar, LO, GO, Node )
  {
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
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( Matrix, StridedMaps_Tpetra, M, MA, Scalar, LO, GO, Node )
  {
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    const size_t numLocal = 10;
    const size_t INVALID = Teuchos::OrdinalTraits<size_t>::invalid(); // TODO: global_size_t instead of size_t


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
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( Matrix, StridedMaps_Epetra, M, MA, Scalar, LO, GO, Node )
  {
#ifdef HAVE_XPETRA_EPETRA
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    const size_t numLocal = 10;
    const size_t INVALID = Teuchos::OrdinalTraits<size_t>::invalid();

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


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( Matrix, BuildCopy_StridedMaps_Tpetra, M, MA, Scalar, LO, GO, Node )
  {
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    const size_t numLocal = 10;
    const size_t INVALID = Teuchos::OrdinalTraits<size_t>::invalid(); // TODO: global_size_t instead of size_t


    typedef Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node> CrsMatrixWrap;
    Teuchos::RCP<const Xpetra::Map<LO,GO,Node> > map =
        Xpetra::MapFactory<LO,GO,Node>::createContigMapWithNode (Xpetra::UseTpetra,INVALID,numLocal,comm);
     {
       Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > s0(new Xpetra::TpetraCrsMatrix<Scalar, LO, GO, Node>(map, numLocal));
       s0->fillComplete();

       Teuchos::RCP<Xpetra::Matrix<Scalar, LO, GO, Node> > s(new CrsMatrixWrap(s0));
       s->SetFixedBlockSize(2);

       Teuchos::RCP<Xpetra::Matrix<Scalar, LO, GO, Node> > t = Xpetra::MatrixFactory2<Scalar,LO,GO,Node>::BuildCopy(s);
       
       int blkSize = t->GetFixedBlockSize();
       TEST_EQUALITY_CONST(blkSize, 2);
     }
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( Matrix, BlockDiagonalUtils_Tpetra, M, MA, Scalar, LO, GO, Node )
  {
    Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
    const size_t numLocal = 10;
    const size_t INVALID = Teuchos::OrdinalTraits<size_t>::invalid(); // TODO: global_size_t instead of size_t
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_const_cast;
    Scalar SC_one = Teuchos::ScalarTraits<Scalar>::one();


    using MV  = Xpetra::MultiVector<Scalar,LO,GO,Node>;
    typedef Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node> CrsMatrixWrap;
    RCP<const Xpetra::Map<LO,GO,Node> > map =
      Xpetra::MapFactory<LO,GO,Node>::createContigMapWithNode (Xpetra::UseTpetra,INVALID,numLocal,comm);   
     {
       // create the identity matrix, via three arrays constructor
       ArrayRCP<size_t> rowptr(numLocal+1);
       ArrayRCP<LO>     colind(numLocal); // one unknown per row
       ArrayRCP<Scalar> values(numLocal); // one unknown per row
       
       for(size_t i=0; i<numLocal; i++){
         rowptr[i] = i;
         colind[i] = Teuchos::as<LO>(i);
         values[i] = SC_one + SC_one;
       }
       rowptr[numLocal]=numLocal;
       RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > eye2  = Xpetra::CrsMatrixFactory<Scalar,LO,GO,Node>::Build(map,map,0);
       TEST_NOTHROW( eye2->setAllValues(rowptr,colind,values) );
       TEST_NOTHROW( eye2->expertStaticFillComplete(map,map) );

       RCP<const Xpetra::Matrix<Scalar, LO, GO, Node> > eye2x(new CrsMatrixWrap(eye2));
       
       // Just extract & scale; don't test correctness (Tpetra does this)
       RCP<MV> diag5c = Xpetra::MultiVectorFactory<Scalar,LO,GO,Node>::Build(map,5);  
       RCP<MV> diag5 = rcp_const_cast<MV>(diag5c);
       diag5->putScalar(SC_one);

       Xpetra::MatrixUtils<Scalar, LO, GO, Node>::extractBlockDiagonal(*eye2x,*diag5);

       RCP<MV> toScale5 = Xpetra::MultiVectorFactory<Scalar,LO,GO,Node>::Build(map,2); toScale5->putScalar(SC_one);
       
       Xpetra::MatrixUtils<Scalar, LO, GO, Node>::inverseScaleBlockDiagonal(*diag5c,false,*toScale5);
      

     }
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

// List of tests which run only with Tpetra
#define XP_TPETRA_MATRIX_INSTANT(S,LO,GO,N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( Matrix, StridedMaps_Tpetra,  M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( Matrix, BuildCopy_StridedMaps_Tpetra,  M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( Matrix, ViewSwitching, M##LO##GO##N , MA##S##LO##GO##N , S, LO, GO, N ) 

// Tpetra, but no complex
#define XP_TPETRA_MATRIX_INSTANT_NO_COMPLEX(S,LO,GO,N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( Matrix, BlockDiagonalUtils_Tpetra, M##LO##GO##N , MA##S##LO##GO##N , S, LO, GO, N )

// List of tests which run only with Epetra
#define XP_EPETRA_MATRIX_INSTANT(S,LO,GO,N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT( Matrix, StridedMaps_Epetra, M##LO##GO##N , MA##S##LO##GO##N, S, LO, GO, N )

// list of all tests which run both with Epetra and Tpetra
//#define XP_MATRIX_INSTANT(S,LO,GO,N)




#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XPETRA_TPETRA_TYPES )
//TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XP_MATRIX_INSTANT )
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XP_TPETRA_MATRIX_INSTANT )

#if !defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE) && !defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XP_TPETRA_MATRIX_INSTANT_NO_COMPLEX )
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
