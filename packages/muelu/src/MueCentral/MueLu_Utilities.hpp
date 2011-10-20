#ifndef MUELU_UTILITIES_HPP
#define MUELU_UTILITIES_HPP

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_Utils.hpp>

#include "MueLu_ConfigDefs.hpp"

#include <Xpetra_Map.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsOperator.hpp>
#include <Xpetra_OperatorFactory.hpp>
#include <Xpetra_BlockedCrsOperator.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#ifdef HAVE_MUELU_EPETRA
#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_EpetraVector.hpp>
#include <Xpetra_EpetraMultiVector.hpp>
#include "Epetra_RowMatrixTransposer.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#endif

#include "MueLu_MatrixFactory.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Memory.hpp"

#ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#endif

#ifdef HAVE_MUELU_TPETRA
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Xpetra_TpetraVector.hpp>
#include <Xpetra_TpetraMultiVector.hpp>
#include "TpetraExt_MatrixMatrix.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_RowMatrixTransposer_decl.hpp"
#include "Tpetra_RowMatrixTransposer_def.hpp"
#endif

// MPI helper
#define sumAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out));
#define minAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MIN, in, Teuchos::outArg(out));
#define maxAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MAX, in, Teuchos::outArg(out));

#ifdef HAVE_MUELU_ML
#include "ml_operator.h"
#include "ml_epetra_utils.h"
#endif

namespace MueLu {

#ifdef HAVE_MUELU_EPETRA
  using Xpetra::EpetraCrsMatrix;   // TODO: mv in Xpetra_UseShortNamesScalar
  using Xpetra::EpetraMultiVector;
#endif

#ifdef HAVE_MUELU_EPETRA
//defined after Utils class
template<typename SC,typename LO,typename GO,typename NO, typename LMO>
RCP<Xpetra::CrsOperator<SC,LO,GO,NO,LMO> > Convert_Epetra_CrsMatrix_ToXpetra_CrsOperator(RCP<Epetra_CrsMatrix> &epAB);
#endif

/*!
  @class Utils
  @brief MueLu utility class.

  This class provides a number of static helper methods.  Some are temporary and will eventually
  go away, while others should be moved to Xpetra.
*/
  template <class Scalar, 
            class LocalOrdinal  = int, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node          = Kokkos::DefaultNode::DefaultNodeType, 
            class LocalMatOps   = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps > //TODO: or BlockSparseOp ?
  class Utils {
    
#include "MueLu_UseShortNames.hpp"

  public:
#ifdef HAVE_MUELU_EPETRA
    //! @brief Helper utility to pull out the underlying Epetra_MultiVector from an Xpetra::MultiVector.
    static RCP<const Epetra_MultiVector> MV2EpetraMV(RCP<MultiVector> const Vec) {
      //rcp<const EpetraMultiVector> tmpVec = rcp_dynamic_cast<EpetraMultiVector>(Vec);
      RCP<const EpetraMultiVector > tmpVec;
      tmpVec = rcp_dynamic_cast<EpetraMultiVector>(Vec);
      if (tmpVec == Teuchos::null)
        throw(Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::EpetraMultiVector failed"));
      RCP<const Epetra_MultiVector> epVec = tmpVec->getEpetra_MultiVector();
      return epVec;
    } //MV2EpetraMV

    //! @brief Helper utility to pull out the underlying Epetra_MultiVector from an Xpetra::MultiVector.
    static RCP<Epetra_MultiVector> MV2NonConstEpetraMV(RCP<MultiVector> Vec) {
      RCP<const EpetraMultiVector> tmpVec = rcp_dynamic_cast<EpetraMultiVector>(Vec);
      if (tmpVec == Teuchos::null)
        throw(Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::EpetraMultiVector failed"));
      RCP<Epetra_MultiVector> epVec = tmpVec->getEpetra_MultiVector();
      return epVec;
    } //MV2EpetraMV

    //! @brief Helper utility to pull out the underlying Epetra_MultiVector from an Xpetra::MultiVector.
    static Epetra_MultiVector& MV2NonConstEpetraMV(MultiVector &Vec) {
      EpetraMultiVector const &tmpVec = dynamic_cast<EpetraMultiVector const&>(Vec);
      RCP<Epetra_MultiVector> epVec = tmpVec.getEpetra_MultiVector();
      return *epVec;
    } //MV2EpetraMV

    static Epetra_MultiVector const& MV2EpetraMV(MultiVector const &Vec) {
      EpetraMultiVector const &tmpVec = dynamic_cast<EpetraMultiVector const&>(Vec);
      RCP<Epetra_MultiVector const> epVec = tmpVec.getEpetra_MultiVector();
      return *epVec;
    } //MV2EpetraMV

    //! @brief Helper utility to pull out the underlying Epetra_CrsMatrix from an Xpetra::Operator.
   static RCP<const Epetra_CrsMatrix> Op2EpetraCrs(RCP<Operator> Op) {
      RCP<const Epetra_CrsMatrix> A;
      // Get the underlying Epetra Mtx
      RCP<const CrsOperator> crsOp = rcp_dynamic_cast<const CrsOperator>(Op);
      if (crsOp == Teuchos::null)
        throw(Exceptions::BadCast("Cast from Xpetra::Operator to Xpetra::CrsOperator failed"));
      RCP<const CrsMatrix> tmp_CrsMtx = crsOp->getCrsMatrix();
      const RCP<const EpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const EpetraCrsMatrix>(tmp_CrsMtx);
      if (tmp_ECrsMtx == Teuchos::null)
        throw(Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed"));
      A = tmp_ECrsMtx->getEpetra_CrsMatrix();
      return A;
    } //Op2EpetraCrs


    //! @brief Helper utility to pull out the underlying Epetra_CrsMatrix from an Xpetra::Operator.
   static RCP<Epetra_CrsMatrix> Op2NonConstEpetraCrs(RCP<Operator> Op) {
      RCP<Epetra_CrsMatrix> A;
      // Get the underlying Epetra Mtx
      RCP<const CrsOperator> crsOp = rcp_dynamic_cast<const CrsOperator>(Op);
      if (crsOp == Teuchos::null)
        throw(Exceptions::BadCast("Cast from Xpetra::Operator to Xpetra::CrsOperator failed"));
      RCP<const CrsMatrix> tmp_CrsMtx = crsOp->getCrsMatrix();
      const RCP<const EpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const EpetraCrsMatrix>(tmp_CrsMtx);
      if (tmp_ECrsMtx == Teuchos::null)
        throw(Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed"));
      A = tmp_ECrsMtx->getEpetra_CrsMatrixNonConst();
      return A;
    } //Op2NonConstEpetraCrs
#endif

#ifdef HAVE_MUELU_TPETRA
    //! @brief Helper utility to pull out the underlying Tpetra::MultiVector from an Xpetra::MultiVector.
    static RCP<const Tpetra::MultiVector<SC,LO,GO,NO> > MV2TpetraMV(RCP<MultiVector> const Vec) {
      //rcp<const TpetraMultiVector> tmpVec = rcp_dynamic_cast<TpetraMultiVector>(Vec);
      RCP<const TpetraMultiVector > tmpVec;
      tmpVec = rcp_dynamic_cast<TpetraMultiVector>(Vec);
      if (tmpVec == Teuchos::null)
        throw(Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::TpetraMultiVector failed"));
      RCP<const Tpetra::MultiVector<SC,LO,GO,NO> > tpVec = tmpVec->getTpetra_MultiVector();
      return tpVec;
    } //MV2TpetraMV

    //! @brief Helper utility to pull out the underlying Tpetra::MultiVector from an Xpetra::MultiVector.
    static RCP<Tpetra::MultiVector<SC,LO,GO,NO> > MV2NonConstTpetraMV(RCP<MultiVector> Vec) {
      RCP<const TpetraMultiVector> tmpVec = rcp_dynamic_cast<TpetraMultiVector>(Vec);
      if (tmpVec == Teuchos::null)
        throw(Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::TpetraMultiVector failed"));
      RCP<Tpetra::MultiVector<SC,LO,GO,NO> > tpVec = tmpVec->getTpetra_MultiVector();
      return tpVec;
    } //MV2TpetraMV

    //! @brief Helper utility to pull out the underlying Tpetra::MultiVector from an Xpetra::MultiVector.
    static Tpetra::MultiVector<SC,LO,GO,NO> & MV2NonConstTpetraMV(MultiVector &Vec) {
      TpetraMultiVector const &tmpVec = dynamic_cast<TpetraMultiVector const&>(Vec);
      RCP<Tpetra::MultiVector<SC,LO,GO,NO> > tpVec = tmpVec.getTpetra_MultiVector();
      return *tpVec;
    } //MV2TpetraMV

    //! @brief Helper utility to pull out the underlying Tpetra::MultiVector from an Xpetra::MultiVector.
    static RCP<Tpetra::MultiVector<SC,LO,GO,NO> > MV2NonConstTpetraMV2(MultiVector &Vec) {
      TpetraMultiVector const &tmpVec = dynamic_cast<TpetraMultiVector const&>(Vec);
      RCP<Tpetra::MultiVector<SC,LO,GO,NO> > tpVec = tmpVec.getTpetra_MultiVector();
      return tpVec;
    } //MV2TpetraMV

    static Tpetra::MultiVector<SC,LO,GO,NO>  const& MV2TpetraMV(MultiVector const &Vec) {
      TpetraMultiVector const &tmpVec = dynamic_cast<TpetraMultiVector const&>(Vec);
      RCP<Tpetra::MultiVector<SC,LO,GO,NO>  const> tpVec = tmpVec.getTpetra_MultiVector();
      return *tpVec;
    } //MV2TpetraMV
    //! @brief Helper utility to pull out the underlying Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> from an Xpetra::Operator.
    static RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > Op2TpetraCrs(RCP<Operator> Op) {
     RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > A;
      // Get the underlying Tpetra Mtx
      RCP<const CrsOperator> crsOp = rcp_dynamic_cast<const CrsOperator>(Op);
      if (crsOp == Teuchos::null)
        throw(Exceptions::BadCast("Cast from Xpetra::Operator to Xpetra::CrsOperator failed"));
      RCP<const CrsMatrix> tmp_CrsMtx = crsOp->getCrsMatrix();
      const RCP<const TpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const TpetraCrsMatrix>(tmp_CrsMtx);
      if (tmp_ECrsMtx == Teuchos::null)
        throw(Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed"));
      A = tmp_ECrsMtx->getTpetra_CrsMatrix();
      return A;
    } //Op2TpetraCrs

    //! @brief Helper utility to pull out the underlying Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> from an Xpetra::Operator.
   static RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > Op2NonConstTpetraCrs(RCP<Operator> Op) {
      RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > A;
      // Get the underlying Tpetra Mtx
      RCP<const CrsOperator> crsOp = rcp_dynamic_cast<const CrsOperator>(Op);
      if (crsOp == Teuchos::null)
        throw(Exceptions::BadCast("Cast from Xpetra::Operator to Xpetra::CrsOperator failed"));
      RCP<const CrsMatrix> tmp_CrsMtx = crsOp->getCrsMatrix();
      const RCP<const TpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const TpetraCrsMatrix>(tmp_CrsMtx);
      if (tmp_ECrsMtx == Teuchos::null)
        throw(Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed"));
      A = tmp_ECrsMtx->getTpetra_CrsMatrixNonConst();
      return A;
    } //Op2NonConstTpetraCrs

#endif

    /*! @brief Helper function to do matrix-matrix multiply "in-place"

      Returns RCP to non-constant Xpetra::Operator.

      @param A left matrix
      @param transposeA if true, use the transpose of A
      @param B right matrix
      @param transposeB if true, use the transpose of B
      @param callFillCompleteOnResult if true, the resulting matrix should be fillComplete'd
    */
   static RCP<Operator> TwoMatrixMultiply(RCP<Operator> const &A, bool transposeA,
                                          RCP<Operator> const &B, bool transposeB,
                                          bool doFillComplete=true,
                                          bool doOptimizeStorage=true)
    {
      RCP<Operator> C;
      //TODO Can we come up with an estimate for nnz-per-row for result C?
      if(transposeA) C = OperatorFactory::Build(A->getDomainMap(), 1);
      else C = OperatorFactory::Build(A->getRowMap(), 1);

      if (!A->isFillComplete())
        throw(Exceptions::RuntimeError("A is not fill-completed"));
      if (!B->isFillComplete())
        throw(Exceptions::RuntimeError("B is not fill-completed"));

      if (C->getRowMap()->lib() == Xpetra::UseEpetra)
      {
#       ifndef HAVE_MUELU_EPETRA_AND_EPETRAEXT
        throw(Exceptions::RuntimeError("MueLu::TwoMatrixMultiply requires EpetraExt to be compiled."));
#       else
        RCP<Epetra_CrsMatrix> epA = Op2NonConstEpetraCrs(A);
        RCP<Epetra_CrsMatrix> epB = Op2NonConstEpetraCrs(B);
        RCP<Epetra_CrsMatrix>       epC = Op2NonConstEpetraCrs(C);

        //ML's multiply cannot implicitly tranpose either matrix.
        bool canUseML=true;
        if (transposeA || transposeB)
          canUseML=false;

        switch (canUseML) {

          case true:
#if 0 // Jonathan's ML-MULTIPLY
            //if ML is not enabled, this case falls through to the EpetraExt multiply.
#           if defined(HAVE_MUELU_ML)
            {
            //ML matrix multiply wrap that uses ML_Operator_WrapEpetraCrsMatrix
            ML_Comm* comm;
            ML_Comm_Create(&comm);
            if (comm->ML_mypid == 0)
              std::cout << "****** USING ML's MATRIX MATRIX MULTIPLY ******" << std::endl;
#           ifdef HAVE_MPI
            // ML_Comm uses MPI_COMM_WORLD, so try to use the same communicator as epA.
            const Epetra_MpiComm * Mcomm=dynamic_cast<const Epetra_MpiComm*>(&(epA->Comm()));
            if(Mcomm) ML_Comm_Set_UsrComm(comm,Mcomm->GetMpiComm());
#           endif
            //in order to use ML, there must be no indices missing from the matrix column maps.
            EpetraExt::CrsMatrix_SolverMap AcolMapTransform;
            Epetra_CrsMatrix *transA = &(AcolMapTransform(*epA));
            ML_Operator *mlA = ML_Operator_Create(comm);
            ML_Operator_WrapEpetraCrsMatrix(transA,mlA);

            EpetraExt::CrsMatrix_SolverMap BcolMapTransform;
            Epetra_CrsMatrix *transB = &(BcolMapTransform(*epB));
            ML_Operator *mlB = ML_Operator_Create(comm);
            ML_Operator_WrapEpetraCrsMatrix(transB,mlB);

            ML_Operator *mlAB = ML_Operator_Create(comm);
            ML_2matmult(mlA,mlB,mlAB,ML_CSR_MATRIX);
    
            /* Wrap back */
            int nnz;
            double time;
            Epetra_CrsMatrix *result;
            ML_Operator2EpetraCrsMatrix(mlAB,result,nnz,false,time,0,false);
            result->OptimizeStorage();
            ML_Operator_Destroy(&mlA);
            ML_Operator_Destroy(&mlB);
            ML_Operator_Destroy(&mlAB);
            ML_Comm_Destroy(&comm);

            RCP<Epetra_CrsMatrix> epAB(result);
#else // Michael's MLMULTIPLY
            if (comm->ML_mypid == 0)
              std::cout << "****** USING ML's MATRIX MATRIX MULTIPLY (LNM)******" << std::endl;

            RCP<Epetra_CrsMatrix> epAB = MLTwoMatrixMultiply(
                        RCP<Epetra_CrsMatrix> epA,
                        RCP<Epetra_CrsMatrix> epB);
#endif
            RCP<CrsOperator> tmpC3 = Convert_Epetra_CrsMatrix_ToXpetra_CrsOperator<SC,LO,GO,NO,LMO>(epAB);
            C = tmpC3;
            }
            break;
#           endif //if HAVE_MUELU_ML

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

        } //switch (canUseML)

#       endif //ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT

      } else if(C->getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
        RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > tpA = Op2TpetraCrs(A);
        RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > tpB = Op2TpetraCrs(B);
        RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> >       tpC = Op2NonConstTpetraCrs(C);

        Tpetra::MatrixMatrix::Multiply(*tpA,transposeA,*tpB,transposeB,*tpC,false);        
#else
        throw(Exceptions::RuntimeError("MueLu must be compiled with Tpetra."));
#endif
      }

      if (!doOptimizeStorage) {
        C->fillComplete((transposeB) ? B->getRangeMap() : B->getDomainMap(),
                        (transposeA) ? A->getDomainMap() : A->getRangeMap(),
                        Xpetra::DoNotOptimizeStorage);
      } else {
        C->fillComplete((transposeB) ? B->getRangeMap() : B->getDomainMap(),
                        (transposeA) ? A->getDomainMap() : A->getRangeMap(),
                        Xpetra::DoOptimizeStorage);
      }

      return C;
    } //TwoMatrixMultiply()

   // Michael Gee's MLMultiply
   static RCP<Epetra_CrsMatrix> MLTwoMatrixMultiply(RCP<Epetra_CrsMatrix> epA,
            RCP<Epetra_CrsMatrix> epB)
    {
#if defined(HAVE_MUELU_ML)
        ML_Comm* comm;
        ML_Comm_Create(&comm);
        if (comm->ML_mypid == 0)
          std::cout << "****** USING ML's MATRIX MATRIX MULTIPLY (LNM version) ******" << std::endl;
#           ifdef HAVE_MPI
        // ML_Comm uses MPI_COMM_WORLD, so try to use the same communicator as epA.
        const Epetra_MpiComm * Mcomm=dynamic_cast<const Epetra_MpiComm*>(&(epA->Comm()));
        if(Mcomm) ML_Comm_Set_UsrComm(comm,Mcomm->GetMpiComm());
#           endif
        //in order to use ML, there must be no indices missing from the matrix column maps.
        EpetraExt::CrsMatrix_SolverMap Atransform;
        EpetraExt::CrsMatrix_SolverMap Btransform;
        const Epetra_CrsMatrix& A = Atransform(*epA);
        const Epetra_CrsMatrix& B = Btransform(*epB);

        if (!A.Filled())    throw(Exceptions::RuntimeError("A has to be FillCompeleted"));
        if (!B.Filled())    throw(Exceptions::RuntimeError("B has to be FillCompeleted"));

        // create ML operators from EpetraCrsMatrix
        ML_Operator* ml_As = ML_Operator_Create(comm);
        ML_Operator* ml_Bs = ML_Operator_Create(comm);
        ML_Operator_WrapEpetraCrsMatrix(const_cast<Epetra_CrsMatrix*>(&A),ml_As);
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

        return result;
#else // no MUELU_ML
        TEUCHOS_TEST_FOR_EXCEPTION( true, Xpetra::Exceptions::RuntimeError,
                         "HAVE_MUELU_ML compiler flag not set. no ML multiply available." );
        return Teuchos::null;
#endif
    }

   /*! @brief Helper function to do matrix-matrix multiply "in-place"

     Returns RCP to non-constant Xpetra::BlockedCrsOperator.

     @param A left matrix
     @param transposeA if true, use the transpose of A
     @param B right matrix
     @param transposeB if true, use the transpose of B
     @param callFillCompleteOnResult if true, the resulting matrix should be fillComplete'd
   */
  static RCP<BlockedCrsOperator> TwoMatrixMultiplyBlock(RCP<BlockedCrsOperator> const &A, bool transposeA,
                                         RCP<BlockedCrsOperator> const &B, bool transposeB,
                                         bool doFillComplete=true,
                                         bool doOptimizeStorage=true)
  {
    if(transposeA || transposeB)
      throw(Exceptions::RuntimeError("TwoMatrixMultiply for BlockedCrsOperator not implemented for transposeA==true or transposeB==true"));

    // todo make sure that A and B are filled and completed

    RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rgmapextractor = A->getRangeMapExtractor();
    RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > domapextractor = B->getDomainMapExtractor();

    RCP<BlockedCrsOperator> C = rcp(new BlockedCrsOperator(rgmapextractor,
                     domapextractor,
                     33 /* TODO fix me */));

    // loop over all block rows of A
    for(size_t i=0; i<A->Rows(); ++i)
    {
      // loop over all block columns of B
      for(size_t j=0; j<B->Cols(); ++j)
      {
        // empty CrsOperator
        RCP<Operator> Cij = OperatorFactory::Build(A->getRangeMap(i), 33 /* TODO fix me */);

        // loop for calculating entry C_{ij}
        for(size_t l=0; l<B->Rows(); ++l)
        {
          RCP<CrsMatrix> crmat1 = A->getMatrix(i,l);
          RCP<CrsMatrix> crmat2 = B->getMatrix(l,j);
          RCP<CrsOperator> crop1 = rcp(new CrsOperator(crmat1));
          RCP<CrsOperator> crop2 = rcp(new CrsOperator(crmat2));

          RCP<Operator> temp = MueLu::Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::TwoMatrixMultiply(crop1, false, crop2, false);

          // sum up
          MueLu::Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TwoMatrixAdd(temp, false, 1.0, Cij, 1.0);
        }

        Cij->fillComplete(B->getDomainMap(j), A->getRangeMap(i));

        RCP<CrsOperator> crsCij = Teuchos::rcp_dynamic_cast<CrsOperator>(Cij);
        TEUCHOS_TEST_FOR_EXCEPTION( Cij==Teuchos::null, Xpetra::Exceptions::BadCast,
                 "OperatorFactory failed in generating a CrsOperator." );

        RCP<CrsMatrix> crsMatCij = crsCij->getCrsMatrix();
        C->setMatrix(i,j,crsMatCij);

      }
    }

    if(doFillComplete)
      C->fillComplete();  // call default fillComplete for BlockCrsOperator objects

    return C;
  } // TwoMatrixMultiplyBlock

    /*! @brief Helper function to calculate B = alpha*A + beta*B.

      @param A      left matrix operand
      @param transposeA indicate whether to use transpose of A
      @param alpha  scalar multiplier for A
      @param B      right matrix operand
      @param beta   scalar multiplier for B

      @return sum in B.

      Note that B does not have to be fill-completed.
    */
   static void TwoMatrixAdd(RCP<Operator> const &A, bool transposeA, SC alpha, RCP<Operator> &B, SC beta)
   {
      if ( !(A->getRowMap()->isSameAs(*(B->getRowMap()))) ) {
        throw(Exceptions::Incompatible("TwoMatrixAdd: matrix row maps are not the same."));
      }

      if (A->getRowMap()->lib() == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT
        RCP<const Epetra_CrsMatrix> epA = Op2EpetraCrs(A);
        RCP<Epetra_CrsMatrix> epB = Op2NonConstEpetraCrs(B);
        
        //FIXME is there a bug if beta=0?
        int i = EpetraExt::MatrixMatrix::Add(*epA,transposeA,(double)alpha,*epB,(double)beta);

        if (i != 0) {
          std::ostringstream buf;
          buf << i;
          std::string msg = "EpetraExt::MatrixMatrix::Add return value of " + buf.str();
          throw(Exceptions::RuntimeError(msg));
        }
#else
      throw(Exceptions::RuntimeError("MueLu must be compile with EpetraExt."));
#endif
      } else if(A->getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
        RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > tpA = Op2TpetraCrs(A);
        RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > tpB = Op2NonConstTpetraCrs(B);
        
        Tpetra::MatrixMatrix::Add(*tpA, transposeA, alpha, *tpB, beta);
#else
        throw(Exceptions::RuntimeError("MueLu must be compiled with Tpetra."));
#endif
      }

   } //TwoMatrixAdd()

    /*! @brief Helper function to calculate C = alpha*A + beta*B.

      @param A          left matrix operand
      @param transposeA indicate whether to use transpose of A
      @param alpha      scalar multiplier for A, defaults to 1.0
      @param B          right matrix operand
      @param transposeB indicate whether to use transpose of B
      @param beta       scalar multiplier for B, defaults to 1.0
      @param C          resulting sum

      It is up to the caller to ensure that the resulting matrix sum is fillComplete'd.
    */
   static void TwoMatrixAdd(RCP<Operator> const &A, bool const &transposeA, SC const &alpha,
                                     RCP<Operator> const &B, bool const &transposeB, SC const &beta,
                                     RCP<Operator> &C)
   {
      if ( !(A->getRowMap()->isSameAs(*(B->getRowMap()))) ) {
        throw(Exceptions::Incompatible("TwoMatrixAdd: matrix row maps are not the same."));
      }
      if (C==Teuchos::null)
        //FIXME 5 is a complete guess as to the #nonzeros per row
        C = rcp( new CrsOperator(A->getRowMap(), 5) );

      if (C->getRowMap()->lib() == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT
        RCP<const Epetra_CrsMatrix> epA = Op2EpetraCrs(A);
        RCP<const Epetra_CrsMatrix> epB = Op2EpetraCrs(B);
        RCP<Epetra_CrsMatrix>       epC = Op2NonConstEpetraCrs(C);
        Epetra_CrsMatrix* ref2epC = &*epC; //to avoid a compiler error...

        //FIXME is there a bug if beta=0?
        int i = EpetraExt::MatrixMatrix::Add(*epA,transposeA,(double)alpha,*epB,transposeB,(double)beta,ref2epC);

        if (i != 0) {
          std::ostringstream buf;
          buf << i;
          std::string msg = "EpetraExt::MatrixMatrix::Add return value of " + buf.str();
          throw(Exceptions::RuntimeError(msg));
        }
#else
      throw(Exceptions::RuntimeError("MueLu must be compile with EpetraExt."));
#endif
      } else if(C->getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
        RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > tpA = Op2TpetraCrs(A);
        RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > tpB = Op2TpetraCrs(B);
        RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> >       tpC = Op2NonConstTpetraCrs(C);

        Tpetra::MatrixMatrix::Add(*tpA, transposeA, alpha, *tpB, transposeB, beta, tpC);
#else
        throw(Exceptions::RuntimeError("MueLu must be compile with Tpetra."));
#endif
      }

   } //TwoMatrixAdd()

    static void MatrixPrint(RCP<Operator> const &Op) {
      std::string label = "unlabeled operator";
      MatrixPrint(Op, label);
    }

    static void MatrixPrint(RCP<Operator> const &Op, std::string const &label) {
#ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT 
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

    /*! @brief Get Operator Diagonal
     */
   static RCP<Operator> BuildMatrixDiagonal(RCP<Operator> const &A)
    {
      const RCP<const Map> rowmap = A->getRowMap();
      //std::vector<SC> diag(A->getNodeNumRows());
      std::vector<SC> diag(rowmap->getNodeNumElements());
      Teuchos::ArrayView<const LO> cols;
      Teuchos::ArrayView<const SC> vals;
      //for (size_t i=0; i<A->getNodeNumRows(); ++i) {
      for (size_t i=0; i<rowmap->getNodeNumElements(); ++i) {
        A->getLocalRowView(i,cols,vals);
        //for (Teuchos::ArrayView<const LO>::size_type j=0; j<cols.size(); j++) {
        for (size_t j=0; j<cols.size(); j++) {
          //TODO this will break down if diagonal entry is not present
          //if (!(cols[j] > i)) {  //JG says this will work ... maybe
          if (cols[j] == i) {
            diag[i] = vals[j];
            break;
          }
        }
      }

      RCP< Operator > D = rcp( new CrsOperator(rowmap, 1) );
      std::vector<LO> diagInd(1);
      Teuchos::ArrayView<GO> iv(&diagInd[0],1);
      //for (size_t i=0; i< A->getNodeNumRows(); ++i) {
      for (size_t i=0; i< rowmap->getNodeNumElements(); ++i) {
        Teuchos::ArrayView<SC> av(&diag[i],1);
        diagInd[0] = rowmap->getGlobalElement(i);
        D->insertGlobalValues(i,iv,av);
      }
      D->fillComplete();
      //MatrixPrint(D);

      return D;

    } //BuildMatrixDiagonal()

    /*! @brief Extract Operator Diagonal

        Returns Operator diagonal in ArrayRCP.

        Note -- it's assumed that A has been fillComplete'd.
    */
    static Teuchos::ArrayRCP<SC> GetMatrixDiagonal(RCP<Operator> const &A)
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

    /*! @brief Left scale matrix by an arbitrary vector.

       Algorithmically, this left scales a matrix by a diagonal matrix.
       The inverse of a diagonal matrix can also be applied.

       @param Op matrix to be scaled
       @param scalingVector vector that represents diagonal matrix
       @doInverse Indicates whether the inverse of the diagonal matrix should be applied.  (Default is to use inverse.)
     */
   static void ScaleMatrix(RCP<Operator> &Op, Teuchos::ArrayRCP<SC> const &scalingVector, bool doInverse=true)
   {
#ifdef HAVE_MUELU_TPETRA
      RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > tpOp;
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

    /*! @brief Get reciprocal of Operator diagonal
     */

   static RCP<Operator> BuildMatrixInverseDiagonal(RCP<Operator> const &A)
    {
      const RCP<const Map> rowmap = A->getRowMap();
      //std::vector<SC> diag(A->getNodeNumRows());
      std::vector<SC> diag(rowmap->getNodeNumElements());
      Teuchos::ArrayView<const LO> cols;
      Teuchos::ArrayView<const SC> vals;
      //for (size_t i=0; i<A->getNodeNumRows(); ++i) {
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

      RCP< Operator > D = rcp( new CrsOperator(rowmap, 1) );
      std::vector<GO> diagInd(1);
      Teuchos::ArrayView<GO> iv(&diagInd[0],1);
      //for (size_t i=0; i< A->getNodeNumRows(); ++i) {

      for (size_t i=0; i< rowmap->getNodeNumElements(); ++i) {
        Teuchos::ArrayView<SC> av(&diag[i],1);
        diagInd[0] = rowmap->getGlobalElement(i);
        D->insertGlobalValues(rowmap->getGlobalElement(i),iv,av); //TODO is this expensive?
      }
      D->fillComplete();

      return D;

    } //BuildMatrixInverseDiagonal()

   typedef typename Teuchos::ScalarTraits<SC>::magnitudeType Magnitude;

    // TODO: should NOT return an Array. Definition must be changed to:
    // - ArrayRCP<> ResidualNorm(Operator const &Op, MultiVector const &X, MultiVector const &RHS)
    // or
    // - void ResidualNorm(Operator const &Op, MultiVector const &X, MultiVector const &RHS, Array &)
   static Teuchos::Array<Magnitude>
   ResidualNorm(Operator const &Op, MultiVector const &X, MultiVector const &RHS)
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
    
    static RCP<MultiVector> Residual(Operator const &Op, MultiVector const &X, MultiVector const &RHS)
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

   /*! @brief Save matrix to file in Matrix Market format.

     TODO Move this to Xpetra?
   */
   static void Write(std::string const & fileName, Operator const & Op) {
    CrsOperator const & crsOp = dynamic_cast<CrsOperator const &>(Op);
    RCP<const CrsMatrix> tmp_CrsMtx = crsOp.getCrsMatrix();
#ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT
    const RCP<const EpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const EpetraCrsMatrix>(tmp_CrsMtx);
    if (tmp_ECrsMtx != Teuchos::null) {
#ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT
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
#endif // HAVE_MUELU_EPETRA_AND_EPETRAEXT

#ifdef HAVE_MUELU_TPETRA
    const RCP<const TpetraCrsMatrix> &tmp_TCrsMtx = rcp_dynamic_cast<const TpetraCrsMatrix>(tmp_CrsMtx);    
    if (tmp_TCrsMtx != Teuchos::null) {
      RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > A = tmp_TCrsMtx->getTpetra_CrsMatrix();
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> >::writeSparseFile(fileName,A);
      return;
    }
#endif // HAVE_MUELU_TPETRA

    throw(Exceptions::BadCast("Could not cast to EpetraCrsMatrix or TpetraCrsMatrix in matrix writing"));

   } //Write

#include <unistd.h>

   static void PauseForDebugger()
   {
     RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
   
     int mypid = comm->getRank();
   
     for (int i = 0; i <comm->getSize() ; i++) {
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


   /*! @brief Simple transpose for Tpetra::CrsMatrix types

      Note:  This is very inefficient, as it inserts one entry at a time.
   */
#ifdef HAVE_MUELU_TPETRA
   static RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > simple_Transpose(RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > const &A)
   {
      LocalOrdinal N=A->getNodeNumRows();
      RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > AT=rcp(new Tpetra::CrsMatrix<SC,LO,GO,NO,LMO>(A->getDomainMap(),0));
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

#ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT
   /*! @brief Simple transpose for Epetra_CrsMatrix types

      Note:  This is very inefficient, as it inserts one entry at a time.
   */
   static RCP<Epetra_CrsMatrix> simple_EpetraTranspose(RCP<const Epetra_CrsMatrix> const &A)
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

    /*! @brief Power method.

      @param A matrix
      @param scaleByDiag if true, estimate the largest eigenvalue of \f$ D^{-1} A \f$.
      @param niters maximum number of iterations
      @param tolerance stopping tolerance
      @verbose if true, print iteration information
      
      (Shamelessly grabbed from tpetra/examples.)
    */
    static Scalar PowerMethod(Operator const &A, bool scaleByDiag=true,
                              LO niters=10, Magnitude tolerance=1e-2, bool verbose=false, unsigned int seed = 123)
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

   static void MyOldScaleMatrix(RCP<Operator> &Op, Teuchos::ArrayRCP<SC> const &scalingVector, bool doInverse=true,
                                bool doFillComplete=true,
                                bool doOptimizeStorage=true)
   {
   //Note: Epetra and Tpetra could be enabled simultaneously.
#ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT
     std::string TorE = "epetra";
#else
     std::string TorE = "tpetra";
#endif

#ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT
      RCP<const Epetra_CrsMatrix> epOp;
      try {
        epOp = Op2NonConstEpetraCrs(Op);
      }
      catch (...){
        TorE = "tpetra";
      }
#endif

#ifdef HAVE_MUELU_TPETRA
      RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > tpOp;
      if (TorE=="tpetra") {
        try {
          tpOp = Op2NonConstTpetraCrs(Op);
        }
        catch(...) {
          throw(Exceptions::RuntimeError("Only Epetra_CrsMatrix or Tpetra::CrsMatrix types can be scaled (Err.1)"));
        }
      } //if
#endif

      Teuchos::ArrayRCP<SC> sv(scalingVector.size());
      if (doInverse) {
        for (int i=0; i<scalingVector.size(); ++i)
          sv[i] = 1.0 / scalingVector[i];
      } else {
        for (int i=0; i<scalingVector.size(); ++i)
          sv[i] = scalingVector[i];
      }

      if (TorE == "tpetra") {
#ifdef HAVE_MUELU_TPETRA
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
              scaledVals[j] = vals[j]*sv[i];
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
              scaledVals[j] = vals[j]*sv[i]; //FIXME i or gid?
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
          if (doOptimizeStorage)
            Op->fillComplete(Op->getDomainMap(),Op->getRangeMap(),Xpetra::DoOptimizeStorage);
          else
            Op->fillComplete(Op->getDomainMap(),Op->getRangeMap(),Xpetra::DoNotOptimizeStorage);
        }
#else
        throw(Exceptions::RuntimeError("Tpetra"));   
#endif // HAVE_MUELU_TPETRA
      } 

      if (TorE == "epetra") {
#ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT
        Epetra_Map const &rowMap = epOp->RowMap();
        int nnz;
        double *vals;
        int *cols;
        for (int i=0; i<rowMap.NumMyElements(); ++i) {
          epOp->ExtractMyRowView(i,nnz,vals,cols);
          for (int j=0; j<nnz; ++j)
            vals[j] *= sv[i];
        }
#else
        throw(Exceptions::RuntimeError("Epetra (Err. 1)"));   
#endif // HAVE_MUELU_EPETRA_AND_EPETRAEXT
      }

      if (TorE != "epetra" && TorE != "tpetra")
        //throw should already have occured, thus should never get here
        throw(Exceptions::RuntimeError("Only Epetra_CrsMatrix or Tpetra::CrsMatrix types can be scaled (Err. 2)"));

   } //ScaleMatrix()

   static RCP<Teuchos::FancyOStream> MakeFancy(std::ostream & os) {
     RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(os));
     return fancy;
   }

  }; // class Utils

#ifdef HAVE_MUELU_EPETRA
//This non-member templated function exists so that the matrix-matrix multiply will compile if Epetra, Tpetra, and ML are enabled.
template<class SC,class LO,class GO,class NO, class LMO>
RCP<Xpetra::CrsOperator<SC,LO,GO,NO,LMO> > Convert_Epetra_CrsMatrix_ToXpetra_CrsOperator(RCP<Epetra_CrsMatrix> &epAB) {
   TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Convert_Epetra_CrsMatrix_ToXpetra_CrsOperator cannot be used with Scalar != double, LocalOrdinal != int, GlobalOrdinal != int");
   return Teuchos::null;
}

typedef Kokkos::DefaultNode::DefaultNodeType KDNT;
typedef Kokkos::DefaultKernels<void,int,Kokkos::DefaultNode::DefaultNodeType>::SparseOps KDKSO;

//specialization for the case of ScalarType=double and LocalOrdinal=GlobalOrdinal=int
template<>
inline RCP<Xpetra::CrsOperator<double,int,int,KDNT,KDKSO> > Convert_Epetra_CrsMatrix_ToXpetra_CrsOperator<double,int,int,KDNT,KDKSO > (RCP<Epetra_CrsMatrix> &epAB) {
  RCP<Xpetra::EpetraCrsMatrix> tmpC1 = rcp(new Xpetra::EpetraCrsMatrix(epAB));
   RCP<Xpetra::CrsMatrix<double,int,int,KDNT,KDKSO> > tmpC2 = rcp_implicit_cast<Xpetra::CrsMatrix<double,int,int,KDNT,KDKSO> >(tmpC1);
   RCP<Xpetra::CrsOperator<double,int,int,KDNT,KDKSO> > tmpC3 = rcp(new Xpetra::CrsOperator<double,int,int,KDNT,KDKSO>(tmpC2));
   return tmpC3;
}
#endif

//! Little helper function to convert non-string types to strings
template<class T>
std::string toString(T const &what) {
  std::ostringstream buf; buf << what;
  return buf.str();
}

//RCP<Xpetra::CrsOperator<double,int,int,KDNT,KDKSO> > Convert_Epetra_CrsMatrix_ToXpetra_CrsOperator<double,int,int,KDNT,KDKSO > (RCP<Epetra_CrsMatrix> epAB) {

/*
  Separate class for Utilities that need a specialization for Epetra.
*/
  template <class Scalar, 
            class LocalOrdinal  = int,
            class GlobalOrdinal = LocalOrdinal,
            class Node          = Kokkos::DefaultNode::DefaultNodeType,
            class LocalMatOps   = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps > //TODO: or BlockSparseOp ?
  class Utils2 {

#include "MueLu_UseShortNames.hpp"

public:

    /*! @brief Transpose a Xpetra::Operator

      Note: Currently, an error is thrown if the matrix isn't a Tpetra::CrsMatrix or Epetra_CrsMatrix.
      In principle, however, we could allow any Epetra_RowMatrix because the Epetra transposer does.
    */

   static RCP<Operator> Transpose(RCP<Operator> const &Op, bool const & optimizeTranspose=false)
   {
#ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT
     std::string TorE = "epetra";
#else
     std::string TorE = "tpetra";
#endif

#ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT
     RCP<Epetra_CrsMatrix> epetraOp;
     try {
       epetraOp = Utils<SC,LO,GO,NO,LMO>::Op2NonConstEpetraCrs(Op);
     }
     catch (...) {
       TorE = "tpetra";
     }
#endif

#ifdef HAVE_MUELU_TPETRA
     RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > tpetraOp;
     if (TorE=="tpetra") {
       try {
         tpetraOp = Utils<SC,LO,GO,NO,LMO>::Op2TpetraCrs(Op);
       }
       catch (...) {
         throw(Exceptions::RuntimeError("Utils::Transpose: Can only transpose Crs matrices"));
       }
     } //if
#endif

     if (TorE == "tpetra") {
#ifdef HAVE_MUELU_TPETRA
       //     Tpetra::RowMatrixTransposer<SC,LO,GO,NO,LMO> transposer(*tpetraOp); //more than meets the eye
       //     RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > A = transposer.createTranspose(optimizeTranspose ? Tpetra::DoOptimizeStorage : Tpetra::DoNotOptimizeStorage); //couldn't have just used a bool...
       RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > A=Utils<SC,LO,GO>::simple_Transpose(tpetraOp);
       RCP<TpetraCrsMatrix> AA = rcp(new TpetraCrsMatrix(A) );
       RCP<CrsMatrix> AAA = rcp_implicit_cast<CrsMatrix>(AA);
       RCP<Operator> AAAA = rcp( new CrsOperator(AAA) );
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
  }; // class Utils2

  // specialization Utils2 for SC=double
  template<>
  class Utils2<double,int,int>//, Kokkos::DefaultNode::DefaultNodeType,
               //Kokkos::DefaultKernels<double,int,Kokkos::DefaultNode::DefaultNodeType>::SparseOps >
  {
   typedef Xpetra::Operator<double,int,int> Operator;
   typedef double SC;
   typedef int LO;
   typedef int GO;
   typedef Kokkos::DefaultNode::DefaultNodeType NO;
   typedef Kokkos::DefaultKernels<double,int,NO>::SparseOps LMO;

public:

   static RCP<Operator> Transpose(RCP<Operator> const &Op, bool const & optimizeTranspose=false)
   {
#ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT
     std::string TorE = "epetra";
#else
     std::string TorE = "tpetra";
#endif

#ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT
     RCP<Epetra_CrsMatrix> epetraOp;
     try {
       epetraOp = Utils<SC,LO,GO,NO,LMO>::Op2NonConstEpetraCrs(Op);
     }
     catch (...) {
       TorE = "tpetra";
     }
#endif

#ifdef HAVE_MUELU_TPETRA
     RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > tpetraOp;
     if (TorE=="tpetra") {
       try {
         tpetraOp = Utils<SC,LO,GO,NO,LMO>::Op2TpetraCrs(Op);
       }
       catch (...) {
         throw(Exceptions::RuntimeError("Utils::Transpose: Can only transpose Crs matrices"));
       }
     } //if
#endif

     if (TorE == "tpetra") {
#ifdef HAVE_MUELU_TPETRA
       //     Tpetra::RowMatrixTransposer<SC,LO,GO,NO,LMO> transposer(*tpetraOp); //more than meets the eye
       //     RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > A = transposer.createTranspose(optimizeTranspose ? Tpetra::DoOptimizeStorage : Tpetra::DoNotOptimizeStorage); //couldn't have just used a bool...
       RCP<Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > A=Utils<SC,LO,GO>::simple_Transpose(tpetraOp);
       RCP<Xpetra::TpetraCrsMatrix<SC> > AA = rcp(new Xpetra::TpetraCrsMatrix<SC>(A) );
       RCP<Xpetra::CrsMatrix<SC> > AAA = rcp_implicit_cast<Xpetra::CrsMatrix<SC> >(AA);
       RCP<Xpetra::CrsOperator<SC> > AAAA = rcp( new Xpetra::CrsOperator<SC> (AAA) );
       AAAA->fillComplete(Op->getRangeMap(),Op->getDomainMap());
       return AAAA;
#else
       throw(Exceptions::RuntimeError("Tpetra"));
#endif
     } else {
#ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT
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
       //RCP<Epetra_CrsMatrix> rcpA = Utils<SC,LO,GO,NO,LMO>::simple_EpetraTranspose(epetraOp);
       RCP<EpetraCrsMatrix> AA = rcp(new EpetraCrsMatrix(rcpA) );
       RCP<Xpetra::CrsMatrix<SC> > AAA = rcp_implicit_cast<Xpetra::CrsMatrix<SC> >(AA);
       RCP<Xpetra::CrsOperator<SC> > AAAA = rcp( new Xpetra::CrsOperator<SC>(AAA) );
       AAAA->fillComplete(Op->getRangeMap(),Op->getDomainMap());
       return AAAA;
#else
       throw(Exceptions::RuntimeError("Epetra (Err. 2)"));
#endif
     }
     
   } //Transpose
  }; //specialization to Scalar=double


} //namespace MueLu
#define MUELU_UTILITIES_SHORT

#endif //ifndef MUELU_UTILITIES_HPP
