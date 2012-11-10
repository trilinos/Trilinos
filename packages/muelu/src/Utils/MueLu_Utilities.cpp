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
#include <MueLu_Utilities_def.hpp>

namespace MueLu {

  RCP<Xpetra::Matrix<double, int, int> > Utils2<double, int, int>::Transpose(RCP<Xpetra::Matrix<double, int, int> > const &Op, bool const & optimizeTranspose)
  {
   typedef Xpetra::Matrix<double,int,int> Matrix;
   typedef double Scalar;
   typedef int LocalOrdinal;
   typedef int GlobalOrdinal;
   typedef Kokkos::DefaultNode::DefaultNodeType Node;
   typedef Kokkos::DefaultKernels<double,int,NO>::SparseOps LocalMatOps;

#ifdef HAVE_MUELU_EPETRAEXT
    std::string TorE = "epetra";
#else
    std::string TorE = "tpetra";
#endif

#ifdef HAVE_MUELU_EPETRAEXT
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
      //     Tpetra::RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> transposer(*tpetraOp); //more than meets the eye
      //     RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > A = transposer.createTranspose(optimizeTranspose ? Tpetra::DoOptimizeStorage : Tpetra::DoNotOptimizeStorage); //couldn't have just used a bool...
      RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > A = Utils<SC,LO,GO>::simple_Transpose(tpetraOp);
      RCP<Xpetra::TpetraCrsMatrix<SC> > AA = rcp(new Xpetra::TpetraCrsMatrix<SC>(A) );
      RCP<Xpetra::CrsMatrix<SC> > AAA = rcp_implicit_cast<Xpetra::CrsMatrix<SC> >(AA);
      RCP<Xpetra::CrsMatrixWrap<SC> > AAAA = rcp( new Xpetra::CrsMatrixWrap<SC> (AAA) );
      AAAA->fillComplete(Op->getRangeMap(),Op->getDomainMap());
      return AAAA;
#else
      throw(Exceptions::RuntimeError("Tpetra"));
#endif
    } else {
#ifdef HAVE_MUELU_EPETRAEXT
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
#ifdef HAVE_MUELU_EPETRAEXT
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
                           RCP<Matrix> &C)
  {
    if ( !(A->getRowMap()->isSameAs(*(B->getRowMap()))) ) {
      throw(Exceptions::Incompatible("TwoMatrixAdd: matrix row maps are not the same."));
    }
    if (C==Teuchos::null)
      //FIXME 5 is a complete guess as to the #nonzeros per row
      C = rcp( new Xpetra::CrsMatrixWrap<double,int,int,NO,LMO>(A->getRowMap(), 5) );

    if (C->getRowMap()->lib() == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRAEXT
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
