// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
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
#include "MueLu_Utilities_def.hpp"

namespace MueLu {

  RCP<Xpetra::Matrix<double, int, int> > Utils2<double, int, int>::Transpose(RCP<Xpetra::Matrix<double, int, int> > const &Op, bool const & optimizeTranspose)
  {
   typedef double Scalar;
   typedef int LocalOrdinal;
   typedef int GlobalOrdinal;
   typedef Kokkos::DefaultNode::DefaultNodeType Node;
   typedef Kokkos::DefaultKernels<double,int,NO>::SparseOps LocalMatOps;

   Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("ZZ Entire Transpose"));

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
    std::string TorE = "epetra";
#else
    std::string TorE = "tpetra";
#endif

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
    RCP<Epetra_CrsMatrix> epetraOp;
    try {

      epetraOp = Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2NonConstEpetraCrs(Op);
    }
    catch (...) {
      TorE = "tpetra";
    }
#endif

#ifdef HAVE_MUELU_TPETRA
    RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tpetraOp;
    if (TorE=="tpetra") {
      try {
        tpetraOp = Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2TpetraCrs(Op);
      }
      catch (...) {
        throw(Exceptions::RuntimeError("Utils::Transpose: Can only transpose Crs matrices"));
      }
    } //if
#endif

    if (TorE == "tpetra") {
#ifdef HAVE_MUELU_TPETRA
      // Compute the transpose A of the Tpetra matrix tpetraOp.
      RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > A;
      {
        Teuchos::TimeMonitor tmm (*Teuchos::TimeMonitor::getNewCounter ("ZZ Tpetra Transpose Only"));
        Tpetra::RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node,
          LocalMatOps> transposer (tpetraOp);
        A = transposer.createTranspose();
      }
      RCP<Xpetra::TpetraCrsMatrix<SC> > AA = rcp(new Xpetra::TpetraCrsMatrix<SC>(A) );
      RCP<Xpetra::CrsMatrix<SC> > AAA = rcp_implicit_cast<Xpetra::CrsMatrix<SC> >(AA);
      RCP<Xpetra::CrsMatrixWrap<SC> > AAAA = rcp( new Xpetra::CrsMatrixWrap<SC> (AAA) );
      return AAAA;
#else
      throw(Exceptions::RuntimeError("Tpetra"));
#endif // HAVE_MUELU_TPETRA
    } else {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
      //epetra case
      Epetra_RowMatrixTransposer et(&*epetraOp);
      Epetra_CrsMatrix *A;
      int rv = et.CreateTranspose(false,A);
      if (rv != 0) {
        std::ostringstream buf;
        buf << rv;
        std::string msg = "Utils::Transpose: Epetra::RowMatrixTransposer returned value of " + buf.str();
        throw(Exceptions::RuntimeError(msg));
      }

      RCP<Epetra_CrsMatrix> rcpA(A);
      //RCP<Epetra_CrsMatrix> rcpA = Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::simple_EpetraTranspose(epetraOp);
      RCP<EpetraCrsMatrix> AA = rcp(new EpetraCrsMatrix(rcpA) );
      RCP<Xpetra::CrsMatrix<SC> > AAA = rcp_implicit_cast<Xpetra::CrsMatrix<SC> >(AA);
      RCP<Xpetra::CrsMatrixWrap<SC> > AAAA = rcp( new Xpetra::CrsMatrixWrap<SC>(AAA) );
      AAAA->fillComplete(Op->getRangeMap(),Op->getDomainMap());
      return AAAA;
#else
      throw(Exceptions::RuntimeError("Epetra (Err. 2)"));
#endif
    }

  } //Transpose

  // -- ------------------------------------------------------- --

  void Utils2<double,int,int>::MyOldScaleMatrix_Epetra(RCP<Matrix> &Op, Teuchos::ArrayRCP<SC> const &scalingVector,
                               bool doFillComplete,
                               bool doOptimizeStorage)
  {
#ifdef HAVE_MUELU_EPETRA
    RCP<const Epetra_CrsMatrix> epOp;
    try {
      epOp = Utils<double,int,int>::Op2NonConstEpetraCrs(Op);
    }
    catch (...){
      throw(Exceptions::RuntimeError("Only Epetra_CrsMatrix types can be scaled"));
    }

      Epetra_Map const &rowMap = epOp->RowMap();
      int nnz;
      double *vals;
      int *cols;
      for (int i=0; i<rowMap.NumMyElements(); ++i) {
        epOp->ExtractMyRowView(i,nnz,vals,cols);
        for (int j=0; j<nnz; ++j)
          vals[j] *= scalingVector[i];
      }
#else
    throw(Exceptions::RuntimeError("Matrix scaling is not possible because Epetra has not been enabled."));
#endif // HAVE_MUELU_EPETRAEXT

  } //Utils2::MyOldScaleMatrix_Epetra()

  // -- ------------------------------------------------------- --

  void Utils2<double, int, int>::TwoMatrixAdd(RCP<Matrix> const &A, bool transposeA, SC alpha, RCP<Matrix> &B, SC beta)
  {
    /*typedef double Scalar;
    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef Kokkos::DefaultNode::DefaultNodeType Node;
    typedef typename Kokkos::DefaultKernels<double,int,Node>::SparseOps LocalMatOps;*/
    //typedef Kokkos::DefaultKernels<double,int,NO>::SparseOps LocalMatOps;

    if ( !(A->getRowMap()->isSameAs(*(B->getRowMap()))) ) {
      throw(Exceptions::Incompatible("TwoMatrixAdd: matrix row maps are not the same."));
    }

    if (A->getRowMap()->lib() == Xpetra::UseEpetra) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
      RCP<const Epetra_CrsMatrix> epA = Utils<double,int,int>::Op2EpetraCrs(A);
      RCP<Epetra_CrsMatrix> epB = Utils<double,int,int>::Op2NonConstEpetraCrs(B);

      //FIXME is there a bug if beta=0?
      int i = EpetraExt::MatrixMatrix::Add(*epA,transposeA,alpha,*epB,beta);

      if (i != 0) {
        std::ostringstream buf;
        buf << i;
        std::string msg = "EpetraExt::MatrixMatrix::Add return value of " + buf.str();
        throw(Exceptions::RuntimeError(msg));
      }
      //Xpetra::MatrixMatrix::Add<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(*A,transposeA,alpha,*B,beta);
#else
      throw(Exceptions::RuntimeError("MueLu must be compiled with EpetraExt."));
#endif
    } else if(A->getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      RCP<const Tpetra::CrsMatrix<SC, LO, GO, NO, LMO> > tpA = Utils<double,int,int>::Op2TpetraCrs(A);
      RCP<Tpetra::CrsMatrix<SC, LO, GO, NO, LMO> > tpB = Utils<double,int,int>::Op2NonConstTpetraCrs(B);

      Tpetra::MatrixMatrix::Add(*tpA, transposeA, alpha, *tpB, beta);
      //Xpetra::MatrixMatrix::Add<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(*A,transposeA,alpha,*B,beta);
#else
      throw(Exceptions::RuntimeError("MueLu must be compiled with Tpetra."));
#endif
    }

  } //Utils2::TwoMatrixAdd() (specialization)

  // -- ------------------------------------------------------- --

  void Utils2<double,int,int>::TwoMatrixAdd(RCP<Matrix> const &A, bool const &transposeA, SC const &alpha,
                           RCP<Matrix> const &B, bool const &transposeB, SC const &beta,
                           RCP<Matrix> &C, bool const &AHasFixedNnzPerRow)
  {
    if ( !(A->getRowMap()->isSameAs(*(B->getRowMap()))) ) {
      throw(Exceptions::Incompatible("TwoMatrixAdd: matrix row maps are not the same."));
    }
///////////////////////////////
    if (C==Teuchos::null) {
      if (!A->isFillComplete() || !B->isFillComplete())
        TEUCHOS_TEST_FOR_EXCEPTION(true,Exceptions::RuntimeError,"Global statistics are not available for estimates.");

      size_t maxNzInA = A->getGlobalMaxNumRowEntries();
      size_t maxNzInB = B->getGlobalMaxNumRowEntries();
      size_t numLocalRows = A->getNodeNumRows();

      RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      fos->setOutputToRootOnly(0);

      if (maxNzInA == 1 || maxNzInB == 1 || AHasFixedNnzPerRow) {
        //first check if either A or B has at most 1 nonzero per row
        //the case of both having at most 1 nz per row is handled by the ``else''
        Teuchos::ArrayRCP<size_t> exactNnzPerRow(numLocalRows);
        if ( (maxNzInA == 1 && maxNzInB > 1) || AHasFixedNnzPerRow) {
          for (size_t i=0; i<numLocalRows; ++i)
            exactNnzPerRow[i] = B->getNumEntriesInLocalRow(Teuchos::as<LO>(i)) + maxNzInA;
        } else {
          for (size_t i=0; i<numLocalRows; ++i)
            exactNnzPerRow[i] = A->getNumEntriesInLocalRow(Teuchos::as<LO>(i)) + maxNzInB;
        }
        *fos << "Utils::TwoMatrixAdd : special case detected (one matrix has a fixed nnz per row)"
             << ", using static profiling" << std::endl;
        C = rcp( new Xpetra::CrsMatrixWrap<double,int,int,NO,LMO>(A->getRowMap(), exactNnzPerRow, Xpetra::StaticProfile) );
      }
      else {
        //general case
        double nnzPerRowInA = Teuchos::as<double>(A->getGlobalNumEntries()) / A->getGlobalNumRows();
        double nnzPerRowInB = Teuchos::as<double>(B->getGlobalNumEntries()) / B->getGlobalNumRows();
        LO nnzToAllocate = Teuchos::as<LO>( (nnzPerRowInA + nnzPerRowInB) * 1.5) + Teuchos::as<LO>(1);

        LO maxPossible = A->getGlobalMaxNumRowEntries() + B->getGlobalMaxNumRowEntries();
        //Use static profiling (more efficient) if the estimate is at least as big as the max
        //possible nnz's in any single row of the result.
        Xpetra::ProfileType pft = (maxPossible) > nnzToAllocate ? Xpetra::DynamicProfile : Xpetra::StaticProfile;

        *fos << "nnzPerRowInA = " << nnzPerRowInA << ", nnzPerRowInB = " << nnzPerRowInB << std::endl;
        *fos << "Utils::TwoMatrixAdd : space allocated per row = " << nnzToAllocate
             << ", max possible nnz per row in sum = " << maxPossible
             << ", using " << (pft == Xpetra::DynamicProfile ? "dynamic" : "static" ) << " profiling"
             << std::endl;
        C = rcp( new Xpetra::CrsMatrixWrap<double,int,int,NO,LMO>(A->getRowMap(), nnzToAllocate, pft) );
      }
      if (transposeB)
        *fos << "Utils::TwoMatrixAdd : ** WARNING ** estimate could be badly wrong because second summand is transposed" << std::endl;
    }
///////////////////////////////
    if (C==Teuchos::null) {
      if (!A->isFillComplete() || !B->isFillComplete())
        TEUCHOS_TEST_FOR_EXCEPTION(true,Exceptions::RuntimeError,"Global statistics are not available for estimates.");

      double nnzPerRowInA = Teuchos::as<double>(A->getGlobalNumEntries()) / A->getGlobalNumRows();
      double nnzPerRowInB = Teuchos::as<double>(B->getGlobalNumEntries()) / B->getGlobalNumRows();
      LO nnzToAllocate = Teuchos::as<LO>( (nnzPerRowInA + nnzPerRowInB) * 1.5) + Teuchos::as<LO>(1);

      LO maxPossible = A->getGlobalMaxNumRowEntries() + B->getGlobalMaxNumRowEntries();
      //Use static profiling (more efficient) if the estimate is at least as big as the max possible nnz's in any single row of the result.
      Xpetra::ProfileType pft = (maxPossible) > nnzToAllocate ? Xpetra::DynamicProfile : Xpetra::StaticProfile;

      RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout)); fos->setOutputToRootOnly(0);
      *fos << "nnzPerRowInA = " << nnzPerRowInA << ", nnzPerRowInB = " << nnzPerRowInB << std::endl;
      *fos << "Utils::TwoMatrixAdd : space allocated per row = " << nnzToAllocate
           << ", max possible nnz per row in sum = " << maxPossible
           << ", using " << (pft == Xpetra::DynamicProfile ? "dynamic" : "static" ) << " profiling"
           << std::endl;
      C = rcp( new Xpetra::CrsMatrixWrap<double,int,int,NO,LMO>(A->getRowMap(), nnzToAllocate, pft) );
      if (transposeB)
        *fos << "Utils::TwoMatrixAdd : ** WARNING ** estimate could be badly wrong because second summand is transposed" << std::endl;
    }

    if (C->getRowMap()->lib() == Xpetra::UseEpetra) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
      RCP<const Epetra_CrsMatrix> epA = Utils<double,int,int>::Op2EpetraCrs(A);
      RCP<const Epetra_CrsMatrix> epB = Utils<double,int,int>::Op2EpetraCrs(B);
      RCP<Epetra_CrsMatrix>       epC = Utils<double,int,int>::Op2NonConstEpetraCrs(C);
      Epetra_CrsMatrix* ref2epC = &*epC; //to avoid a compiler error...

      //FIXME is there a bug if beta=0?
      int i = EpetraExt::MatrixMatrix::Add(*epA,transposeA,alpha,*epB,transposeB,beta,ref2epC);

      if (i != 0) {
        std::ostringstream buf;
        buf << i;
        std::string msg = "EpetraExt::MatrixMatrix::Add return value of " + buf.str();
        throw(Exceptions::RuntimeError(msg));
      }
      //Xpetra::MatrixMatrix::Add<double, int, int, NO, LMO>(*A,transposeA,alpha,*B,transposeB,beta,C);
#else
      throw(Exceptions::RuntimeError("MueLu must be compile with EpetraExt."));
#endif
    } else if(C->getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      RCP<const Tpetra::CrsMatrix<SC, LO, GO, NO, LMO> > tpA = Utils<double,int,int>::Op2TpetraCrs(A);
      RCP<const Tpetra::CrsMatrix<SC, LO, GO, NO, LMO> > tpB = Utils<double,int,int>::Op2TpetraCrs(B);
      RCP<Tpetra::CrsMatrix<SC, LO, GO, NO, LMO> >       tpC = Utils<double,int,int>::Op2NonConstTpetraCrs(C);

      Tpetra::MatrixMatrix::Add(*tpA, transposeA, alpha, *tpB, transposeB, beta, tpC);
      //Xpetra::MatrixMatrix::Add<SC,LO,GO,NO,LMO>(*A,transposeA,alpha,*B,transposeB,beta,C);
#else
      throw(Exceptions::RuntimeError("MueLu must be compile with Tpetra."));
#endif
    }

    ///////////////////////// EXPERIMENTAL
    if(A->IsView("stridedMaps")) C->CreateView("stridedMaps", A);
    if(B->IsView("stridedMaps")) C->CreateView("stridedMaps", B);
    ///////////////////////// EXPERIMENTAL
  }
}
