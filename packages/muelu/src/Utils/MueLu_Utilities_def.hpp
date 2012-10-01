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
#ifndef MUELU_UTILITIES_DEF_HPP
#define MUELU_UTILITIES_DEF_HPP

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ParameterList.hpp>

#include "MueLu_ConfigDefs.hpp"

#ifdef HAVE_MUELU_EPETRA
# ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
# endif
#endif

#ifdef HAVE_MUELU_EPETRAEXT
#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_RowMatrixOut.h>
#include <Epetra_RowMatrixTransposer.h>
#endif // HAVE_MUELU_EPETRAEXT

#ifdef HAVE_MUELU_TPETRA
#include <TpetraExt_MatrixMatrix.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Xpetra_TpetraMultiVector.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#endif // HAVE_MUELU_TPETRA

#include <Xpetra_Map.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MatrixFactory.hpp>

#include <MueLu_Utilities_decl.hpp>

#ifdef HAVE_MUELU_ML
#include <ml_operator.h>
#include <ml_epetra_utils.h>
#endif

#define scan(rcpComm, in, out)                                        \
  Teuchos::scan(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out));

namespace MueLu {

#ifdef HAVE_MUELU_EPETRA
  using Xpetra::EpetraCrsMatrix;   // TODO: mv in Xpetra_UseShortNamesScalar
  using Xpetra::EpetraMultiVector;
#endif

#ifdef HAVE_MUELU_EPETRA
  //defined after Utils class
  template<typename SC,typename LO,typename GO,typename NO, typename LMO>
  RCP<Xpetra::CrsMatrixWrap<SC,LO,GO,NO,LMO> > Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap(RCP<Epetra_CrsMatrix> &epAB);
#endif

#ifdef HAVE_MUELU_EPETRA
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const Epetra_MultiVector> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2EpetraMV(RCP<MultiVector> const Vec) {
    //rcp<const EpetraMultiVector> tmpVec = rcp_dynamic_cast<EpetraMultiVector>(Vec);
    RCP<const EpetraMultiVector > tmpVec;
    tmpVec = rcp_dynamic_cast<EpetraMultiVector>(Vec);
    if (tmpVec == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::EpetraMultiVector failed"));
    RCP<const Epetra_MultiVector> epVec = tmpVec->getEpetra_MultiVector();
    return epVec;
  } //MV2EpetraMV

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Epetra_MultiVector> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2NonConstEpetraMV(RCP<MultiVector> Vec) {
    RCP<const EpetraMultiVector> tmpVec = rcp_dynamic_cast<EpetraMultiVector>(Vec);
    if (tmpVec == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::EpetraMultiVector failed"));
    RCP<Epetra_MultiVector> epVec = tmpVec->getEpetra_MultiVector();
    return epVec;
  } //MV2EpetraMV

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Epetra_MultiVector& Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2NonConstEpetraMV(MultiVector &Vec) {
    EpetraMultiVector const &tmpVec = dynamic_cast<EpetraMultiVector const&>(Vec);
    RCP<Epetra_MultiVector> epVec = tmpVec.getEpetra_MultiVector();
    return *epVec;
  } //MV2EpetraMV


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Epetra_MultiVector const& Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2EpetraMV(MultiVector const &Vec) {
    EpetraMultiVector const &tmpVec = dynamic_cast<EpetraMultiVector const&>(Vec);
    RCP<Epetra_MultiVector const> epVec = tmpVec.getEpetra_MultiVector();
    return *epVec;
  } //MV2EpetraMV


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const Epetra_CrsMatrix> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2EpetraCrs(RCP<Matrix> Op) {
    RCP<const Epetra_CrsMatrix> A;
    // Get the underlying Epetra Mtx
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
    RCP<const CrsMatrix> tmp_CrsMtx = crsOp->getCrsMatrix();
    const RCP<const EpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const EpetraCrsMatrix>(tmp_CrsMtx);
    if (tmp_ECrsMtx == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed"));
    A = tmp_ECrsMtx->getEpetra_CrsMatrix();
    return A;
  } //Op2EpetraCrs


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Epetra_CrsMatrix> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2NonConstEpetraCrs(RCP<Matrix> Op) {
    RCP<Epetra_CrsMatrix> A;
    // Get the underlying Epetra Mtx
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
    RCP<const CrsMatrix> tmp_CrsMtx = crsOp->getCrsMatrix();
    const RCP<const EpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const EpetraCrsMatrix>(tmp_CrsMtx);
    if (tmp_ECrsMtx == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed"));
    A = tmp_ECrsMtx->getEpetra_CrsMatrixNonConst();
    return A;
  } //Op2NonConstEpetraCrs
#endif

#ifdef HAVE_MUELU_TPETRA

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps> 
  RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2TpetraMV(RCP<MultiVector> const Vec) {
    //rcp<const TpetraMultiVector> tmpVec = rcp_dynamic_cast<TpetraMultiVector>(Vec);
    RCP<const TpetraMultiVector > tmpVec;
    tmpVec = rcp_dynamic_cast<TpetraMultiVector>(Vec);
    if (tmpVec == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::TpetraMultiVector failed"));
    RCP<const Tpetra::MultiVector<SC,LO,GO,NO> > tpVec = tmpVec->getTpetra_MultiVector();
    return tpVec;
  } //MV2TpetraMV


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2NonConstTpetraMV(RCP<MultiVector> Vec) {
    RCP<const TpetraMultiVector> tmpVec = rcp_dynamic_cast<TpetraMultiVector>(Vec);
    if (tmpVec == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::TpetraMultiVector failed"));
    RCP<Tpetra::MultiVector<SC,LO,GO,NO> > tpVec = tmpVec->getTpetra_MultiVector();
    return tpVec;
  } //MV2TpetraMV


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2NonConstTpetraMV(MultiVector &Vec) {
    TpetraMultiVector const &tmpVec = dynamic_cast<TpetraMultiVector const&>(Vec);
    RCP<Tpetra::MultiVector<SC,LO,GO,NO> > tpVec = tmpVec.getTpetra_MultiVector();
    return *tpVec;
  } //MV2TpetraMV


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2NonConstTpetraMV2(MultiVector &Vec) {
    TpetraMultiVector const &tmpVec = dynamic_cast<TpetraMultiVector const&>(Vec);
    RCP<Tpetra::MultiVector<SC,LO,GO,NO> > tpVec = tmpVec.getTpetra_MultiVector();
    return tpVec;
  } //MV2TpetraMV


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>  const& Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2TpetraMV(MultiVector const &Vec) {
    TpetraMultiVector const &tmpVec = dynamic_cast<TpetraMultiVector const&>(Vec);
    RCP<Tpetra::MultiVector<SC,LO,GO,NO>  const> tpVec = tmpVec.getTpetra_MultiVector();
    return *tpVec;
  } //MV2TpetraMV


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2TpetraCrs(RCP<Matrix> Op) {
    RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > A;
    // Get the underlying Tpetra Mtx
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
    RCP<const CrsMatrix> tmp_CrsMtx = crsOp->getCrsMatrix();
    const RCP<const TpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const TpetraCrsMatrix>(tmp_CrsMtx);
    if (tmp_ECrsMtx == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed"));
    A = tmp_ECrsMtx->getTpetra_CrsMatrix();
    return A;
  } //Op2TpetraCrs


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2NonConstTpetraCrs(RCP<Matrix> Op) {
    RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > A;
    // Get the underlying Tpetra Mtx
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
    RCP<const CrsMatrix> tmp_CrsMtx = crsOp->getCrsMatrix();
    const RCP<const TpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const TpetraCrsMatrix>(tmp_CrsMtx);
    if (tmp_ECrsMtx == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed"));
    A = tmp_ECrsMtx->getTpetra_CrsMatrixNonConst();
    return A;
  } //Op2NonConstTpetraCrs

#endif


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TwoMatrixMultiply(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > const &A, bool transposeA,
                                         RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > const &B, bool transposeB,
                                         bool doFillComplete,
                                         bool doOptimizeStorage)
  {
#ifdef HAVE_MPI
int mypid=-1;
static double t0=0,t1=0;
t0 = MPI_Wtime();
#endif

    RCP<Matrix> C;
    //TODO Can we come up with an estimate for nnz-per-row for result C?
    if(transposeA) C = MatrixFactory::Build(A->getDomainMap(), 1);
    else C = MatrixFactory::Build(A->getRowMap(), 1);

    if (!A->isFillComplete())
      throw(Exceptions::RuntimeError("A is not fill-completed"));
    if (!B->isFillComplete())
      throw(Exceptions::RuntimeError("B is not fill-completed"));

    if (C->getRowMap()->lib() == Xpetra::UseEpetra)
      {
#       ifndef HAVE_MUELU_EPETRAEXT
        throw(Exceptions::RuntimeError("MueLu::TwoMatrixMultiply requires EpetraExt to be compiled."));
#       else
        RCP<Epetra_CrsMatrix> epA = Op2NonConstEpetraCrs(A);
        RCP<Epetra_CrsMatrix> epB = Op2NonConstEpetraCrs(B);
        RCP<Epetra_CrsMatrix> epC = Op2NonConstEpetraCrs(C);

        //ML's multiply cannot implicitly tranpose either matrix.
        int canUseML=1;
        if (transposeA || transposeB)
          canUseML=0;

        switch (canUseML) {

        case true:
          //if ML is not enabled, this case falls through to the EpetraExt multiply.
#ifdef HAVE_MUELU_ML
          {
//#define USE_JHU_ML_MULTIPLY
#ifdef USE_JHU_ML_MULTIPLY // Jonathan's ML-MULTIPLY

            //ML matrix multiply wrap that uses ML_Matrix_WrapEpetraCrsMatrix
            ML_Comm* comm;
            ML_Comm_Create(&comm);
            if (comm->ML_mypid == 0)
              std::cout << "****** USING ML's MATRIX MATRIX MULTIPLY ******" << std::endl;
#ifdef HAVE_MPI
            mypid = comm->ML_mypid;
            // ML_Comm uses MPI_COMM_WORLD, so try to use the same communicator as epA.
            const Epetra_MpiComm * Mcomm=dynamic_cast<const Epetra_MpiComm*>(&(epA->Comm()));
            if(Mcomm) ML_Comm_Set_UsrComm(comm,Mcomm->GetMpiComm());
#endif // HAVE_MPI
            //in order to use ML, there must be no indices missing from the matrix column maps.
            EpetraExt::CrsMatrix_SolverMap AcolMapTransform;
            Epetra_CrsMatrix *transA = &(AcolMapTransform(*epA));
            ML_Matrix *mlA = ML_Matrix_Create(comm);
            ML_Matrix_WrapEpetraCrsMatrix(transA,mlA);

            EpetraExt::CrsMatrix_SolverMap BcolMapTransform;
            Epetra_CrsMatrix *transB = &(BcolMapTransform(*epB));
            ML_Matrix *mlB = ML_Matrix_Create(comm);
            ML_Matrix_WrapEpetraCrsMatrix(transB,mlB);

            ML_Matrix *mlAB = ML_Matrix_Create(comm);
            ML_2matmult(mlA,mlB,mlAB,ML_CSR_MATRIX);
    
            /* Wrap back */
            int nnz;
            double time;
            Epetra_CrsMatrix *result;
            ML_Matrix2EpetraCrsMatrix(mlAB,result,nnz,false,time,0,false);
            result->OptimizeStorage();
            ML_Matrix_Destroy(&mlA);
            ML_Matrix_Destroy(&mlB);
            ML_Matrix_Destroy(&mlAB);
            ML_Comm_Destroy(&comm);

            RCP<Epetra_CrsMatrix> epAB(result);
#else // Michael's MLMULTIPLY

            RCP<Epetra_CrsMatrix> epAB = MLTwoMatrixMultiply(*epA, *epB);
#endif
            C = Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(epAB);
          }
          break;
#endif // HAVE_MUELU_ML

        case false:
          {
            int i = EpetraExt::MatrixMatrix::Multiply(*epA,transposeA,*epB,transposeB,*epC,false);
            if (i != 0) {
              std::ostringstream buf;
              buf << i;
              std::string msg = "EpetraExt::MatrixMatrix::Multiply return value of " + buf.str();
              throw(Exceptions::RuntimeError(msg));
            }
          }
          break;

        } // switch (canUseML)

#endif // HAVE_MUELU_EPETRAEXT

      } else if(C->getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tpA = Op2TpetraCrs(A);
      RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tpB = Op2TpetraCrs(B);
      RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >       tpC = Op2NonConstTpetraCrs(C);

      Tpetra::MatrixMatrix::Multiply(*tpA,transposeA,*tpB,transposeB,*tpC,false);        
#else
      throw(Exceptions::RuntimeError("MueLu must be compiled with Tpetra."));
#endif
    }

    if(doFillComplete) {
      RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList());
      params->set("Optimize Storage",doOptimizeStorage);
      C->fillComplete((transposeB) ? B->getRangeMap() : B->getDomainMap(),
		      (transposeA) ? A->getDomainMap() : A->getRangeMap(),
		      params);
    }

    ///////////////////////// EXPERIMENTAL
    C->CreateView("stridedMaps", A, transposeA, B, transposeB);
    ///////////////////////// EXPERIMENTAL
    
#ifdef HAVE_MPI
t1 += MPI_Wtime() - t0;
if (mypid == 0)
  std::cout << "cumulative MM time = " << t1 << std::endl;
#endif

    return C;
  } //TwoMatrixMultiply()

#ifdef HAVE_MUELU_EPETRAEXT

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Epetra_CrsMatrix> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MLTwoMatrixMultiply(const Epetra_CrsMatrix& epA,
      const Epetra_CrsMatrix& epB)
  {
#if defined(HAVE_MUELU_ML)
    ML_Comm* comm;
    ML_Comm_Create(&comm);
    if (comm->ML_mypid == 0)
      std::cout << "****** USING ML's MATRIX MATRIX MULTIPLY (LNM version) ******" << std::endl;
#           ifdef HAVE_MPI
    // ML_Comm uses MPI_COMM_WORLD, so try to use the same communicator as epA.
    const Epetra_MpiComm * Mcomm=dynamic_cast<const Epetra_MpiComm*>(&(epA.Comm()));
    if(Mcomm) ML_Comm_Set_UsrComm(comm,Mcomm->GetMpiComm());
#           endif
    //in order to use ML, there must be no indices missing from the matrix column maps.
    EpetraExt::CrsMatrix_SolverMap Atransform;
    EpetraExt::CrsMatrix_SolverMap Btransform;
    const Epetra_CrsMatrix& A = Atransform(const_cast<Epetra_CrsMatrix&>(epA));
    const Epetra_CrsMatrix& B = Btransform(const_cast<Epetra_CrsMatrix&>(epB));

    if (!A.Filled())    throw(Exceptions::RuntimeError("A has to be FillCompeleted"));
    if (!B.Filled())    throw(Exceptions::RuntimeError("B has to be FillCompeleted"));

    // create ML operators from EpetraCrsMatrix
    ML_Operator* ml_As = ML_Operator_Create(comm);
    ML_Operator* ml_Bs = ML_Operator_Create(comm);
    ML_Operator_WrapEpetraCrsMatrix(const_cast<Epetra_CrsMatrix*>(&A),ml_As); // Should we test if the lightweight wrapper is actually used or if WrapEpetraCrsMatrix fall backs to the heavy one?
    ML_Operator_WrapEpetraCrsMatrix(const_cast<Epetra_CrsMatrix*>(&B),ml_Bs);
    ML_Operator* ml_AtimesB = ML_Operator_Create(comm);
    ML_2matmult(ml_As,ml_Bs,ml_AtimesB,ML_CSR_MATRIX); // do NOT use ML_EpetraCRS_MATRIX!!!
    ML_Operator_Destroy(&ml_As);
    ML_Operator_Destroy(&ml_Bs);

    // For ml_AtimesB we have to reconstruct the column map in global indexing,
    // The following is going down to the salt-mines of ML ...
    // note: we use integers, since ML only works for Epetra...
    int N_local = ml_AtimesB->invec_leng;
    ML_CommInfoOP* getrow_comm = ml_AtimesB->getrow->pre_comm;
    if (!getrow_comm)   throw(Exceptions::RuntimeError("ML_Operator does not have a CommInfo"));
    ML_Comm* comm_AB = ml_AtimesB->comm;   // get comm object
    if (N_local != B.DomainMap().NumMyElements())
      throw(Exceptions::RuntimeError("Mismatch in local row dimension between ML and Epetra"));
    int N_rcvd = 0;
    int N_send = 0;
    int flag   = 0;
    for (int i=0; i<getrow_comm->N_neighbors; i++)
      {
        N_rcvd += (getrow_comm->neighbors)[i].N_rcv;
        N_send += (getrow_comm->neighbors)[i].N_send;
        if ( ((getrow_comm->neighbors)[i].N_rcv !=0) &&
             ((getrow_comm->neighbors)[i].rcv_list != NULL) ) flag = 1;
      }
    // For some unknown reason, ML likes to have stuff one larger than
    // neccessary...
    std::vector<double> dtemp(N_local + N_rcvd + 1); // "double" vector for comm function
    std::vector<int>    cmap (N_local + N_rcvd + 1); // vector for gids
    for (int i=0; i<N_local; ++i)
      {
        cmap[i] = B.DomainMap().GID(i);
        dtemp[i] = (double) cmap[i];
      }
    ML_cheap_exchange_bdry(&dtemp[0],getrow_comm,N_local,N_send,comm_AB); // do communication
    if (flag) // process received data
      {
        int count = N_local;
        const int neighbors = getrow_comm->N_neighbors;
        for (int i=0; i< neighbors; i++)
          {
            const int nrcv = getrow_comm->neighbors[i].N_rcv;
            for (int j=0; j<nrcv; j++)
              cmap[getrow_comm->neighbors[i].rcv_list[j]] = (int) dtemp[count++];
          }
      }
    else
      for (int i=0; i<N_local+N_rcvd; ++i) cmap[i] = (int)dtemp[i];
    dtemp.clear();  // free double array

    // we can now determine a matching column map for the result
    Epetra_Map gcmap(-1,N_local+N_rcvd,&cmap[0],0,A.Comm());

    int allocated=0;
    int rowlength;
    double* val=NULL;
    int* bindx=NULL;
    const int myrowlength = A.RowMap().NumMyElements();
    const Epetra_Map& rowmap = A.RowMap();

    // determine the maximum bandwith for the result matrix.
    // replaces the old, very(!) memory-consuming guess:
    // int guessnpr = A.MaxNumEntries()*B.MaxNumEntries();
    int educatedguess = 0;
    for (int i=0; i<myrowlength; ++i)
      {
        // get local row
        ML_get_matrix_row(ml_AtimesB,1,&i,&allocated,&bindx,&val,&rowlength,0);
        if (rowlength>educatedguess) educatedguess = rowlength;
      }

    // allocate our result matrix and fill it
    RCP<Epetra_CrsMatrix> result
      = rcp(new Epetra_CrsMatrix(::Copy,A.RangeMap(),gcmap,educatedguess,false));

    std::vector<int> gcid(educatedguess);
    for (int i=0; i<myrowlength; ++i)
      {
        const int grid = rowmap.GID(i);
        // get local row
        ML_get_matrix_row(ml_AtimesB,1,&i,&allocated,&bindx,&val,&rowlength,0);
        if (!rowlength) continue;
        if ((int)gcid.size() < rowlength) gcid.resize(rowlength);
        for (int j=0; j<rowlength; ++j)
          {
            gcid[j] = gcmap.GID(bindx[j]);
            if (gcid[j]<0) throw(Exceptions::RuntimeError("Error: cannot find gcid!"));
          }
        int err = result->InsertGlobalValues(grid,rowlength,val,&gcid[0]);
        if (err!=0 && err!=1) throw(Exceptions::RuntimeError("Epetra_CrsMatrix::InsertGlobalValues returned err="+err));
      }
    // free memory
    if (bindx) ML_free(bindx);
    if (val) ML_free(val);
    ML_Operator_Destroy(&ml_AtimesB);
    ML_Comm_Destroy(&comm);

    return result;
#else // no MUELU_ML
    TEUCHOS_TEST_FOR_EXCEPTION( true, Xpetra::Exceptions::RuntimeError,
                                "HAVE_MUELU_ML compiler flag not set. no ML multiply available." );
    return Teuchos::null;
#endif
  }
#endif //ifdef HAVE_MUELU_EPETRAEXT

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TwoMatrixMultiplyBlock(RCP<BlockedCrsMatrix> const &A, bool transposeA,
                                                        RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > const &B, bool transposeB,
                                                        bool doFillComplete,
                                                        bool doOptimizeStorage)
  {
    if(transposeA || transposeB)
      throw(Exceptions::RuntimeError("TwoMatrixMultiply for BlockedCrsMatrix not implemented for transposeA==true or transposeB==true"));

    // todo make sure that A and B are filled and completed

    RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rgmapextractor = A->getRangeMapExtractor();
    RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > domapextractor = B->getDomainMapExtractor();

    RCP<BlockedCrsMatrix> C = rcp(new BlockedCrsMatrix(rgmapextractor,
                                                           domapextractor,
                                                           33 /* TODO fix me */));

    // loop over all block rows of A
    for(size_t i=0; i<A->Rows(); ++i)
      {
        // loop over all block columns of B
        for(size_t j=0; j<B->Cols(); ++j)
          {
            // empty CrsMatrixWrap
            RCP<Matrix> Cij = MatrixFactory::Build(A->getRangeMap(i), 33 /* TODO fix me */);

            // loop for calculating entry C_{ij}
            for(size_t l=0; l<B->Rows(); ++l)
              {
                RCP<CrsMatrix> crmat1 = A->getMatrix(i,l);
                RCP<CrsMatrix> crmat2 = B->getMatrix(l,j);
                RCP<CrsMatrixWrap> crop1 = rcp(new CrsMatrixWrap(crmat1));
                RCP<CrsMatrixWrap> crop2 = rcp(new CrsMatrixWrap(crmat2));

                RCP<Matrix> temp = MueLu::Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::TwoMatrixMultiply(crop1, false, crop2, false);

                // sum up
                MueLu::Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TwoMatrixAdd(temp, false, 1.0, Cij, 1.0);
              }

            Cij->fillComplete(B->getDomainMap(j), A->getRangeMap(i));

            RCP<CrsMatrixWrap> crsCij = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(Cij);
            TEUCHOS_TEST_FOR_EXCEPTION( Cij==Teuchos::null, Xpetra::Exceptions::BadCast,
                                        "MatrixFactory failed in generating a CrsMatrixWrap." );

            RCP<CrsMatrix> crsMatCij = crsCij->getCrsMatrix();
            C->setMatrix(i,j,crsMatCij);

          }
      }

    if(doFillComplete)
      C->fillComplete();  // call default fillComplete for BlockCrsMatrixWrap objects

    return C;
  } // TwoMatrixMultiplyBlock


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MatrixPrint(RCP<Matrix> const &Op) {
    std::string label = "unlabeled operator";
    MatrixPrint(Op, label);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MatrixPrint(RCP<Matrix> const &Op, std::string const &label) {
#ifdef HAVE_MUELU_EPETRAEXT 
    RCP<const Epetra_CrsMatrix> epOp = Op2EpetraCrs(Op);
    int mypid = epOp->RowMap().Comm().MyPID();
    if (mypid == 0)
      std::cout << "\n===============\n" << label << "\n==============" << std::endl;

    if (mypid == 0) std::cout << "\n   -- row map -- \n" << std::endl;
    std::cout << epOp->RowMap() << std::endl;
    sleep(1);
    epOp->RowMap().Comm().Barrier();

    if (mypid == 0) std::cout << "\n   -- column map -- \n" << std::endl;
    std::cout << epOp->ColMap() << std::endl;
    sleep(1);
    epOp->RowMap().Comm().Barrier();

    std::cout << *epOp << std::endl;
#endif
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildMatrixDiagonal(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > const &A)
  {
    const RCP<const Map> rowmap = A->getRowMap();
    //std::vector<SC> diag(A->getNodeNumRows());
    std::vector<SC> diag(rowmap->getNodeNumElements());
    Teuchos::ArrayView<const LO> cols;
    Teuchos::ArrayView<const SC> vals;
    //for (size_t i=0; i<A->getNodeNumRows(); ++i)
    for (size_t i=0; i<rowmap->getNodeNumElements(); ++i) {
      A->getLocalRowView(i,cols,vals);
      //for (Teuchos::ArrayView<const LO>::size_type j=0; j<cols.size(); j++)
      for (size_t j=0; j<Teuchos::as<size_t>(cols.size()); j++) { // TODO: cleanup LO vs size_t
        //TODO this will break down if diagonal entry is not present
        //if (!(cols[j] > i)) //JG says this will work ... maybe
        if (cols[j] == Teuchos::as<LO>(i)) {  // TODO: cleanup LO vs size_t
          diag[i] = vals[j];
          break;
        }
      }
    }

    RCP< Matrix > D = rcp( new CrsMatrixWrap(rowmap, 1) );
    std::vector<GO> diagInd(1);
    Teuchos::ArrayView<GO> iv(&diagInd[0],1);
    //for (size_t i=0; i< A->getNodeNumRows(); ++i)
    for (size_t i=0; i< rowmap->getNodeNumElements(); ++i) {
      Teuchos::ArrayView<SC> av(&diag[i],1);
      diagInd[0] = rowmap->getGlobalElement(i);
      D->insertGlobalValues(i,iv,av);
    }
    D->fillComplete();
    //MatrixPrint(D);

    return D;

  } //BuildMatrixDiagonal()

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::ArrayRCP<Scalar> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetMatrixDiagonal(RCP<Matrix> const &A)
  {
    const RCP<const Map> rowmap = A->getRowMap();
    size_t locSize = rowmap->getNodeNumElements();
    Teuchos::ArrayRCP<SC> diag(locSize);
    Teuchos::ArrayView<const LO> cols;
    Teuchos::ArrayView<const SC> vals;
    for (size_t i=0; i<locSize; ++i) {
      A->getLocalRowView(i,cols,vals);
      for (LO j=0; j<cols.size(); ++j) {
        //TODO this will break down if diagonal entry is not present
        //if (!(cols[j] > i))   //JG says this will work ... maybe
        if (Teuchos::as<size_t>(cols[j]) == i) {
          diag[i] = vals[j];
          break;
        }
      }
    }
    //for (int i=0; i<locSize; ++i) std::cout << "diag[" << i << "] = " << diag[i] << std::endl;
    return diag;
  } //GetMatrixDiagonal


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ScaleMatrix(RCP<Matrix> &Op, Teuchos::ArrayRCP<SC> const &scalingVector, bool doInverse)
  {
#ifdef HAVE_MUELU_TPETRA
    RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tpOp;
    try {
      tpOp = Op2NonConstTpetraCrs(Op);
    }
    catch(...) {
      throw(Exceptions::RuntimeError("Sorry, haven't implemented matrix scaling for epetra"));
    }
    Tpetra::Vector<SC,LO,GO,NO> x(tpOp->getRowMap(),scalingVector());
    if(doInverse){
      Tpetra::Vector<SC,LO,GO,NO> xi(tpOp->getRowMap());
      xi.reciprocal(x);
      tpOp->leftScale(xi);
    }
    else
      tpOp->leftScale(x);
#else
    throw(Exceptions::RuntimeError("Sorry, haven't implemented matrix scaling for epetra"));
#endif // HAVE_MUELU_TPETRA
  } //ScaleMatrix()

#ifdef UNUSED // and does not work with SC=complex
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildMatrixInverseDiagonal(RCP<Matrix> const &A)
  {
    const RCP<const Map> rowmap = A->getRowMap();
    //std::vector<SC> diag(A->getNodeNumRows());
    std::vector<SC> diag(rowmap->getNodeNumElements());
    Teuchos::ArrayView<const LO> cols;
    Teuchos::ArrayView<const SC> vals;
    //for (size_t i=0; i<A->getNodeNumRows(); ++i)
    LO rowmapLocalSize = (LO) rowmap->getNodeNumElements();
    for (LO i=0; i<rowmapLocalSize; ++i) {
      A->getLocalRowView(i,cols,vals);
      for (LO j=0; j<cols.size(); ++j) {
        //TODO this will break down if diagonal entry is not present
        if (cols[j] == i) {
          diag[i] = 1 / vals[j];
          break;
        }
      }
    }

    RCP< Matrix > D = rcp( new CrsMatrixWrap(rowmap, 1) );
    std::vector<GO> diagInd(1);
    Teuchos::ArrayView<GO> iv(&diagInd[0],1);
    //for (size_t i=0; i< A->getNodeNumRows(); ++i)

    for (size_t i=0; i< rowmap->getNodeNumElements(); ++i) {
      Teuchos::ArrayView<SC> av(&diag[i],1);
      diagInd[0] = rowmap->getGlobalElement(i);
      D->insertGlobalValues(rowmap->getGlobalElement(i),iv,av); //TODO is this expensive?
    }
    D->fillComplete();

    return D;

  } //BuildMatrixInverseDiagonal()
#endif

  //   typedef typename Teuchos::ScalarTraits<SC>::magnitudeType Magnitude;


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>
  Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ResidualNorm(Matrix const &Op, MultiVector const &X, MultiVector const &RHS)
  {
    //if (X.getNumVectors() != RHS.getNumVectors())
    //  throw(Exceptions::RuntimeError("Number of solution vectors != number of right-hand sides"));
    //const size_t numVecs = X.getNumVectors();
    const size_t numVecs = 1;
    RCP<MultiVector> RES = Residual(Op,X,RHS);
    Teuchos::Array<Magnitude> norms(numVecs);
    RES->norm2(norms);
    return norms;
  }
    

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Residual(Matrix const &Op, MultiVector const &X, MultiVector const &RHS)
  {
    SC one = 1.0;
    SC negone = -1.0;
    //if (X.getNumVectors() != RHS.getNumVectors())
    //  throw(Exceptions::RuntimeError("Number of solution vectors != number of right-hand sides"));
    //const size_t numVecs = X.getNumVectors();
    const size_t numVecs = 1;
    RCP<MultiVector> RES = MultiVectorFactory::Build(Op.getRangeMap(),numVecs);
    Op.apply(X,*RES,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
    RES->update(one,RHS,negone);
    return RES;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Write(std::string const & fileName, Matrix const & Op) {
    CrsMatrixWrap const & crsOp = dynamic_cast<CrsMatrixWrap const &>(Op);
    RCP<const CrsMatrix> tmp_CrsMtx = crsOp.getCrsMatrix();
#ifdef HAVE_MUELU_EPETRAEXT
    const RCP<const EpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const EpetraCrsMatrix>(tmp_CrsMtx);
    if (tmp_ECrsMtx != Teuchos::null) {
#ifdef HAVE_MUELU_EPETRAEXT
      RCP<const Epetra_CrsMatrix> A = tmp_ECrsMtx->getEpetra_CrsMatrix();
      int rv = EpetraExt::RowMatrixToMatrixMarketFile(fileName.c_str(), *A);
      if (rv != 0) {
        std::ostringstream buf;
        buf << rv;
        std::string msg = "EpetraExt::RowMatrixToMatrixMarketFile return value of " + buf.str();
        throw(Exceptions::RuntimeError(msg));
      }
#else
      throw(Exceptions::RuntimeError("Compiled without EpetraExt"));
#endif
      return;
    }
#endif // HAVE_MUELU_EPETRAEXT

#ifdef HAVE_MUELU_TPETRA
    const RCP<const TpetraCrsMatrix> &tmp_TCrsMtx = rcp_dynamic_cast<const TpetraCrsMatrix>(tmp_CrsMtx);    
    if (tmp_TCrsMtx != Teuchos::null) {
      RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > A = tmp_TCrsMtx->getTpetra_CrsMatrix();
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >::writeSparseFile(fileName,A);
      return;
    }
#endif // HAVE_MUELU_TPETRA

    throw(Exceptions::BadCast("Could not cast to EpetraCrsMatrix or TpetraCrsMatrix in matrix writing"));

  } //Write

#include <unistd.h>


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::PauseForDebugger()
  {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
   
    int mypid = comm->getRank();
   
    for (int i = 0; i <comm->getSize(); i++) {
      if (i == mypid ) {
        char buf[80];
        char hostname[80];
        gethostname(hostname, sizeof(hostname));
        int pid = getpid();
        sprintf(buf, "Host: %s\tMPI rank: %d,\tPID: %d\n\tattach %d\n\tcontinue\n",
                hostname, mypid, pid, pid);
        printf("%s\n",buf);
        fflush(stdout);
        sleep(1);
      }
    }
   
    if (mypid == 0) {
      printf( "** Enter a character to continue > "); fflush(stdout);
      char go = ' ';
      scanf("%c",&go);
    }
    comm->barrier();
  } //PauseForDebugger

#ifdef HAVE_MUELU_TPETRA

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::simple_Transpose(RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > const &A)
  {
    LocalOrdinal N=A->getNodeNumRows();
    RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > AT=rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(A->getDomainMap(),0));
    const RCP<const Tpetra::Map<LO,GO,NO> > & rowMap=A->getRowMap();
    const RCP<const Tpetra::Map<LO,GO,NO> > & colMap=A->getColMap();

    for(LO i=0;i<N;i++){
      GO grid= rowMap->getGlobalElement(i);
      Teuchos::ArrayView<const LO> indices;
      Teuchos::ArrayView<const SC> vals;
      A->getLocalRowView(i,indices,vals);
      for(LO j=0;j<indices.size();j++){
        GO gcid=colMap->getGlobalElement(indices[j]);
        AT->insertGlobalValues(gcid,Teuchos::tuple(grid),Teuchos::tuple(vals[j]));
      }
    }
    //AT->fillComplete(A->getRangeMap(),A->getDomainMap());
      
    return AT;
  } //simple_Transpose
#endif // HAVE_MUELU_TPETRA

#ifdef HAVE_MUELU_EPETRAEXT

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Epetra_CrsMatrix> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::simple_EpetraTranspose(RCP<const Epetra_CrsMatrix> const &A)
  {
    int N=A->NumMyRows();
    RCP<Epetra_CrsMatrix> AT=rcp(new Epetra_CrsMatrix(Copy,A->DomainMap(),0));
    const Epetra_Map& rowMap=A->RowMap();
    const Epetra_Map& colMap=A->ColMap();

    for(int i=0;i<N;++i){
      int grid= rowMap.GID(i);
      int *indices,nnz;
      double *vals;
      A->ExtractMyRowView(i,nnz,vals,indices);
      for(int j=0;j<nnz;++j){
        int gcid=colMap.GID(indices[j]);

        //if(AT->RowMap().MyGID(gcid)) // no communication between procs!
        //{
        int rv = AT->InsertGlobalValues(gcid,1,vals+j,&grid);
        if (rv < 0) { // throw on errors, not on warnings...
          std::ostringstream buf;
          buf << rv;
          std::string msg = "Utils::simple_EpetraTranspose: Epetra_CrsMatrix::InsertGlobalValues() returned value of " + buf.str();
          throw(Exceptions::RuntimeError(msg));
        }
        //}
      }
    }
     
    return AT;

  } //simple_Transpose
#endif


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Scalar Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::PowerMethod(Matrix const &A, bool scaleByDiag,
                                                                                    LO niters, Magnitude tolerance, bool verbose, unsigned int seed)
  {
    if ( !(A.getRangeMap()->isSameAs(*(A.getDomainMap()))) ) {
      throw(Exceptions::Incompatible("Utils::PowerMethod: operator must have domain and range maps that are equivalent."));
    }
    // create three vectors, fill z with random numbers
    RCP<MultiVector> q = MultiVectorFactory::Build(A.getRangeMap(),1);
    RCP<MultiVector> r = MultiVectorFactory::Build(A.getRangeMap(),1);
    RCP<MultiVector> z = MultiVectorFactory::Build(A.getRangeMap(),1);
    z->setSeed(seed);  // seed random number generator
    z->randomize(true);// use Xpetra implementation: -> same results for Epetra and Tpetra
      
    Teuchos::Array<Magnitude> norms(1);
  
    //std::vector<Scalar> lambda(1);
    //lambda[0] = 0.0;
    Scalar lambda=0.0;
    Magnitude residual = 0.0;
    // power iteration
    Teuchos::ArrayView<Scalar> avLambda(&lambda,1);
    RCP<Vector> diagVec,oneOverDiagonal;
    if (scaleByDiag) {
      diagVec = VectorFactory::Build(A.getRowMap());
      A.getLocalDiagCopy(*diagVec);
      oneOverDiagonal = VectorFactory::Build(A.getRowMap());
      oneOverDiagonal->reciprocal(*diagVec);
    }
    for (int iter = 0; iter < niters; ++iter) {
      z->norm2(norms);                               // Compute 2-norm of z
      q->update(1.0/norms[0],*z,0.);                 // Set q = z / normz
      A.apply(*q, *z);                               // Compute z = A*q
      if (scaleByDiag) z->elementWiseMultiply(1.0, *oneOverDiagonal, *z, 0.0);
      q->dot(*z,avLambda);                            // Approximate maximum eigenvalue: lamba = dot(q,z)
      if ( iter % 100 == 0 || iter + 1 == niters ) {
        r->update(1.0, *z, -lambda, *q, 0.0);         // Compute A*q - lambda*q
        r->norm2(norms);
        residual = Teuchos::ScalarTraits<Scalar>::magnitude(norms[0] / lambda);
        if (verbose) {
          std::cout << "Iter = " << iter
                    << "  Lambda = " << lambda
                    << "  Residual of A*q - lambda*q = " << residual
                    << std::endl;
        }
      }
      if (residual < tolerance) {
        break;
      }
    }
    return lambda;
  } //PowerMethod


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MyOldScaleMatrix(RCP<Matrix> &Op, Teuchos::ArrayRCP<SC> const &scalingVector, bool doInverse,
                               bool doFillComplete,
                               bool doOptimizeStorage)
  {

    Teuchos::ArrayRCP<SC> sv(scalingVector.size());
    if (doInverse) {
      for (int i=0; i<scalingVector.size(); ++i)
        sv[i] = 1.0 / scalingVector[i];
    } else {
      for (int i=0; i<scalingVector.size(); ++i)
        sv[i] = scalingVector[i];
    }

    switch (Op->getRowMap()->lib()) {
      case Xpetra::UseTpetra:
        MyOldScaleMatrix_Tpetra(Op, sv, doFillComplete, doOptimizeStorage);
        break;
      case Xpetra::UseEpetra:
        Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MyOldScaleMatrix_Epetra(Op, sv, doFillComplete, doOptimizeStorage);
        break;
      default:
        throw(Exceptions::RuntimeError("Only Epetra and Tpetra matrices can be scaled."));
        break;
    } //switch
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MyOldScaleMatrix_Tpetra(RCP<Matrix> &Op, Teuchos::ArrayRCP<SC> const &scalingVector,
                               bool doFillComplete,
                               bool doOptimizeStorage)
  {
#ifdef HAVE_MUELU_TPETRA
    RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tpOp;
    try {
      tpOp = Op2NonConstTpetraCrs(Op);
    }
    catch(...) {
      throw(Exceptions::RuntimeError("Only Tpetra::CrsMatrix types can be scaled (Err.1)"));
    }

      const RCP<const Tpetra::Map<LO,GO,NO> > rowMap = tpOp->getRowMap();
      const RCP<const Tpetra::Map<LO,GO,NO> > domainMap = tpOp->getDomainMap();
      const RCP<const Tpetra::Map<LO,GO,NO> > rangeMap = tpOp->getRangeMap();
      size_t maxRowSize = tpOp->getNodeMaxNumRowEntries();
      if (maxRowSize==(size_t)-1) //hasn't been determined yet
        maxRowSize=20;
      std::vector<SC> scaledVals(maxRowSize);
      if (tpOp->isFillComplete()) {
        tpOp->resumeFill();
      }

      if (Op->isLocallyIndexed() == true) {
        Teuchos::ArrayView<const LO> cols;
        Teuchos::ArrayView<const SC> vals;
        for (size_t i=0; i<rowMap->getNodeNumElements(); ++i) {
          tpOp->getLocalRowView(i,cols,vals);
          size_t nnz = tpOp->getNumEntriesInLocalRow(i);
          if (nnz>maxRowSize) {
            maxRowSize=nnz;
            scaledVals.resize(maxRowSize);
          }
          for (size_t j=0; j<nnz; ++j) {
            scaledVals[j] = vals[j]*scalingVector[i];
          }
          if (nnz>0) {
            Teuchos::ArrayView<const SC> valview(&scaledVals[0],nnz);
            tpOp->replaceLocalValues(i,cols,valview);
          }
        } //for (size_t i=0; ...
      } else {
        Teuchos::ArrayView<const GO> cols;
        Teuchos::ArrayView<const SC> vals;
        for (size_t i=0; i<rowMap->getNodeNumElements(); ++i) {
          GO gid = rowMap->getGlobalElement(i);
          tpOp->getGlobalRowView(gid,cols,vals);
          size_t nnz = tpOp->getNumEntriesInGlobalRow(gid);
          if (nnz>maxRowSize) {
            maxRowSize=nnz;
            scaledVals.resize(maxRowSize);
          }
          for (size_t j=0; j<nnz; ++j) {
            scaledVals[j] = vals[j]*scalingVector[i]; //FIXME i or gid?
          }
          if (nnz>0) {
            Teuchos::ArrayView<const SC> valview(&scaledVals[0],nnz);
            tpOp->replaceGlobalValues(gid,cols,valview);
          }
        } //for (size_t i=0; ...
      }

      if (doFillComplete) {
        if (domainMap == Teuchos::null || rangeMap == Teuchos::null)
          throw(Exceptions::RuntimeError("In Utils::Scaling: cannot fillComplete because the domain and/or range map hasn't been defined"));
        RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList());
        params->set("Optimize Storage",doOptimizeStorage);
        Op->fillComplete(Op->getDomainMap(),Op->getRangeMap(),params);
      }
#else
      throw(Exceptions::RuntimeError("Matrix scaling is not possible because Tpetra has not been enabled."));
#endif
  } //MyOldScaleMatrix_Tpetra()

  typedef Kokkos::DefaultNode::DefaultNodeType KDNT;
  typedef Kokkos::DefaultKernels<void,int,Kokkos::DefaultNode::DefaultNodeType>::SparseOps KDKSO;


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Teuchos::FancyOStream> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MakeFancy(std::ostream & os) {
    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(os));
    return fancy;
  }

#ifdef HAVE_MUELU_EPETRA
//   template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
//   RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap(RCP<Epetra_CrsMatrix> &epAB) {
//     TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap cannot be used with Scalar != double, LocalOrdinal != int, GlobalOrdinal != int");
//     return Teuchos::null;
//   }


//   template<>
//   inline RCP<Xpetra::CrsMatrixWrap<double,int,int,KDNT,KDKSO> > Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap<double,int,int,KDNT,KDKSO > (RCP<Epetra_CrsMatrix> &epAB) {
//     RCP<Xpetra::EpetraCrsMatrix> tmpC1 = rcp(new Xpetra::EpetraCrsMatrix(epAB));
//     RCP<Xpetra::CrsMatrix<double,int,int,KDNT,KDKSO> > tmpC2 = rcp_implicit_cast<Xpetra::CrsMatrix<double,int,int,KDNT,KDKSO> >(tmpC1);
//     RCP<Xpetra::CrsMatrixWrap<double,int,int,KDNT,KDKSO> > tmpC3 = rcp(new Xpetra::CrsMatrixWrap<double,int,int,KDNT,KDKSO>(tmpC2));
//     return tmpC3;
//   }
#endif

//   template<class T>
//   std::string toString(T const &what) {
//     std::ostringstream buf; buf << what;
//     return buf.str();
//   }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Transpose(RCP<Matrix> const &Op, bool const & optimizeTranspose)
  {
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
      RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > A=Utils<Scalar, LocalOrdinal, GlobalOrdinal>::simple_Transpose(tpetraOp);
      RCP<TpetraCrsMatrix> AA = rcp(new TpetraCrsMatrix(A) );
      RCP<CrsMatrix> AAA = rcp_implicit_cast<CrsMatrix>(AA);
      RCP<Matrix> AAAA = rcp( new CrsMatrixWrap(AAA) );
      AAAA->fillComplete(Op->getRangeMap(),Op->getDomainMap());
      return AAAA;
#else
      throw(Exceptions::RuntimeError("Tpetra"));
#endif
    } 

    //epetra case
    std::cout << "Utilities::Transpose() not implemented for Epetra" << std::endl;
    return Teuchos::null;
     
  } //Transpose

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MyOldScaleMatrix_Epetra(RCP<Matrix> &Op, Teuchos::ArrayRCP<SC> const &scalingVector,
                               bool doFillComplete,
                               bool doOptimizeStorage)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MyOldScalematrix and Epetra cannot be used with Scalar != double, LocalOrdinal != int, GlobalOrdinal != int");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TwoMatrixAdd(RCP<Matrix> const &A, bool transposeA, SC alpha, RCP<Matrix> &B, SC beta)
  {
    if ( !(A->getRowMap()->isSameAs(*(B->getRowMap()))) ) {
      throw(Exceptions::Incompatible("TwoMatrixAdd: matrix row maps are not the same."));
    }

    if (A->getRowMap()->lib() == Xpetra::UseEpetra) {
      throw(Exceptions::RuntimeError("You cannot use Epetra::MatrixMatrix::Add with Scalar!=double or Ordinal!=int"));
    } else if(A->getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tpA = Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::Op2TpetraCrs(A);
      RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tpB = Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::Op2NonConstTpetraCrs(B);
        
      Tpetra::MatrixMatrix::Add(*tpA, transposeA, alpha, *tpB, beta);
#else
      throw(Exceptions::RuntimeError("MueLu must be compiled with Tpetra."));
#endif
    }

  } //Utils2::TwoMatrixAdd()


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TwoMatrixAdd(RCP<Matrix> const &A, bool const &transposeA, SC const &alpha,
                           RCP<Matrix> const &B, bool const &transposeB, SC const &beta,
                           RCP<Matrix> &C)
  {
    if ( !(A->getRowMap()->isSameAs(*(B->getRowMap()))) ) {
      throw(Exceptions::Incompatible("TwoMatrixAdd: matrix row maps are not the same."));
    }
    if (C==Teuchos::null)
      //FIXME 5 is a complete guess as to the #nonzeros per row
      C = rcp( new CrsMatrixWrap(A->getRowMap(), 5) );

    if (C->getRowMap()->lib() == Xpetra::UseEpetra) {
      throw(Exceptions::RuntimeError("You cannot use Epetra::MatrixMatrix::Add with Scalar!=double or Ordinal!=int"));
    } else if(C->getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tpA = Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::Op2TpetraCrs(A);
      RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tpB = Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::Op2TpetraCrs(B);
      RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >       tpC = Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::Op2NonConstTpetraCrs(C);

      Tpetra::MatrixMatrix::Add(*tpA, transposeA, alpha, *tpB, transposeB, beta, tpC);
#else
      throw(Exceptions::RuntimeError("MueLu must be compile with Tpetra."));
#endif
    }

    ///////////////////////// EXPERIMENTAL
    if(A->IsView("stridedMaps")) C->CreateView("stridedMaps", A);
    if(B->IsView("stridedMaps")) C->CreateView("stridedMaps", B);
    ///////////////////////// EXPERIMENTAL

  } //Utils2::TwoMatrixAdd()


} //namespace MueLu

#define MUELU_UTILITIES_SHORT
#endif // MUELU_UTILITIES_DEF_HPP

//  LocalWords:  LocalOrdinal
