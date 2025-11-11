// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PACKAGES_XPETRA_SUP_UTILS_XPETRA_MATRIXMATRIX_DECL_HPP_
#define PACKAGES_XPETRA_SUP_UTILS_XPETRA_MATRIXMATRIX_DECL_HPP_

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_BlockedCrsMatrix.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_MapExtractor.hpp"
#include "Xpetra_Map.hpp"
#include "Xpetra_MatrixFactory.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_StridedMapFactory.hpp"
#include "Xpetra_StridedMap.hpp"

#include "Xpetra_Helpers.hpp"

#ifdef HAVE_XPETRA_EPETRA
#include <Xpetra_EpetraCrsMatrix_fwd.hpp>
#endif

#ifdef HAVE_XPETRA_EPETRAEXT
#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_RowMatrixOut.h>
#include <Epetra_RowMatrixTransposer.h>
#endif  // HAVE_XPETRA_EPETRAEXT

#ifdef HAVE_XPETRA_TPETRA
#include <TpetraExt_MatrixMatrix.hpp>
#include <Tpetra_RowMatrixTransposer.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Xpetra_TpetraBlockCrsMatrix.hpp>
#include <Tpetra_BlockCrsMatrix_Helpers.hpp>
#include <Xpetra_TpetraMultiVector.hpp>
#include <Xpetra_TpetraVector.hpp>
#endif  // HAVE_XPETRA_TPETRA

namespace Xpetra {

template <class Scalar,
          class LocalOrdinal /*= int*/,
          class GlobalOrdinal /*= LocalOrdinal*/,
          class Node /*= Tpetra::KokkosClassic::DefaultNode::DefaultNodeType*/>
class MatrixMatrix {
#undef XPETRA_MATRIXMATRIX_SHORT
#include "Xpetra_UseShortNames.hpp"

 public:
  /** Given CrsMatrix objects A, B and C, form the product C = A*B.
    In a parallel setting, A and B need not have matching distributions,
    but C needs to have the same row-map as A (if transposeA is false).
    At this time C=AT*B and C=A*BT are known to not work. However,
    C=A*B and C=AT*BT are known to work, Kurtis Nusbaum 03/24/2011

    @param A Input, must already have had 'FillComplete()' called.
    @param transposeA Input, whether to use transpose of matrix A.
    @param B Input, must already have had 'FillComplete()' called.
    @param transposeB Input, whether to use transpose of matrix B.
    @param C Result. On entry to this method, it doesn't matter whether
    FillComplete() has already been called on C or not. If it has,
    then C's graph must already contain all nonzero locations that
    will be produced when forming the product A*B. On exit,
    C.FillComplete() will have been called, unless the last argument
    to this function is specified to be false.
    @param call_FillComplete_on_result Optional argument, defaults to true.
    Power users may specify this argument to be false if they *DON'T*
    want this function to call C.FillComplete. (It is often useful
    to allow this function to call C.FillComplete, in cases where
    one or both of the input matrices are rectangular and it is not
    trivial to know which maps to use for the domain- and range-maps.)

*/
  static void Multiply(const Matrix& A, bool transposeA,
                       const Matrix& B, bool transposeB,
                       Matrix& C,
                       bool call_FillComplete_on_result = true,
                       bool doOptimizeStorage           = true,
                       const std::string& label         = std::string(),
                       const RCP<ParameterList>& params = null);

  /**
    @brief Helper function to do matrix-matrix multiply

    Given CrsMatrix objects A, B and C, form the product C = A*B.
    In a parallel setting, A and B need not have matching distributions,
    but C needs to have the same row-map as A (if transposeA is false).
    At this time C=AT*B and C=A*BT are known to not work. However,
    C=A*B and C=AT*BT are known to work, Kurtis Nusbaum 03/24/2011

    @param A Input, must already have had 'FillComplete()' called.
    @param transposeA Input, whether to use transpose of matrix A.
    @param B Input, must already have had 'FillComplete()' called.
    @param transposeB Input, whether to use transpose of matrix B.
    @param C Result. If Teuchos::null, a new CrsMatrix is created with optimal number of nnz per row.
    @param call_FillComplete_on_result Optional argument, defaults to true.
    Power users may specify this argument to be false if they *DON'T*
    want this function to call C.FillComplete. (It is often useful
    to allow this function to call C.FillComplete, in cases where
    one or both of the input matrices are rectangular and it is not
    trivial to know which maps to use for the domain- and range-maps.)

*/
  static RCP<Matrix> Multiply(const Matrix& A, bool transposeA, const Matrix& B, bool transposeB, RCP<Matrix> C_in,
                              Teuchos::FancyOStream& fos,
                              bool doFillComplete              = true,
                              bool doOptimizeStorage           = true,
                              const std::string& label         = std::string(),
                              const RCP<ParameterList>& params = null);

  /*! @brief Helper function to do matrix-matrix multiply

    Returns C = AB.

    @param A left matrix
    @param transposeA if true, use the transpose of A
    @param B right matrix
    @param transposeB if true, use the transpose of B
    @param callFillCompleteOnResult if true, the resulting matrix should be fillComplete'd
    */
  static RCP<Matrix> Multiply(const Matrix& A, bool transposeA, const Matrix& B, bool transposeB, Teuchos::FancyOStream& fos,
                              bool callFillCompleteOnResult = true, bool doOptimizeStorage = true, const std::string& label = std::string(),
                              const RCP<ParameterList>& params = null);

#ifdef HAVE_XPETRA_EPETRAEXT
  // Michael Gee's MLMultiply
  static RCP<Epetra_CrsMatrix> MLTwoMatrixMultiply(const Epetra_CrsMatrix& epA,
                                                   const Epetra_CrsMatrix& epB,
                                                   Teuchos::FancyOStream& fos);
#endif  // ifdef HAVE_XPETRA_EPETRAEXT

  /*! @brief Helper function to do matrix-matrix multiply "in-place"

    Returns RCP to non-constant Xpetra::BlockedCrsMatrix.

    @param A left matrix
    @param transposeA if true, use the transpose of A
    @param B right matrix
    @param transposeB if true, use the transpose of B
    @param doOptimizeStorage if true, the resulting matrix should be fillComplete'd
    */
  static RCP<BlockedCrsMatrix> TwoMatrixMultiplyBlock(const BlockedCrsMatrix& A, bool transposeA,
                                                      const BlockedCrsMatrix& B, bool transposeB,
                                                      Teuchos::FancyOStream& fos,
                                                      bool doFillComplete    = true,
                                                      bool doOptimizeStorage = true);

  /*! @brief Helper function to calculate B = alpha*A + beta*B.

    @param A      left matrix operand
    @param transposeA indicate whether to use transpose of A
    @param alpha  scalar multiplier for A
    @param B      right matrix operand
    @param beta   scalar multiplier for B

    @return sum in B.

    Note that B does not have to be fill-completed.
    */
  static void TwoMatrixAdd(const Matrix& A, bool transposeA, SC alpha, Matrix& B, SC beta);

  /*! @brief Helper function to calculate C = alpha*A + beta*B.

    @param A          left matrix operand
    @param transposeA indicate whether to use transpose of A
    @param alpha      scalar multiplier for A, defaults to 1.0
    @param B          right matrix operand
    @param transposeB indicate whether to use transpose of B
    @param beta       scalar multiplier for B, defaults to 1.0
    @param C          resulting sum
    @param fos        output stream for printing to screen
    @param AHasFixedNnzPerRow

    It is up to the caller to ensure that the resulting matrix sum is fillComplete'd.
    */
  static void TwoMatrixAdd(const Matrix& A, bool transposeA, const SC& alpha,
                           const Matrix& B, bool transposeB, const SC& beta,
                           RCP<Matrix>& C, Teuchos::FancyOStream& fos, bool AHasFixedNnzPerRow = false);

};  // class MatrixMatrix

#ifdef HAVE_XPETRA_EPETRA
// specialization MatrixMatrix for SC=double, LO=GO=int
template <>
class MatrixMatrix<double, int, int, EpetraNode> {
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef EpetraNode Node;
#include "Xpetra_UseShortNames.hpp"

 public:
  /** Given CrsMatrix objects A, B and C, form the product C = A*B.
    In a parallel setting, A and B need not have matching distributions,
    but C needs to have the same row-map as A (if transposeA is false).
    At this time C=AT*B and C=A*BT are known to not work. However,
    C=A*B and C=AT*BT are known to work, Kurtis Nusbaum 03/24/2011

    @param A Input, must already have had 'FillComplete()' called.
    @param transposeA Input, whether to use transpose of matrix A.
    @param B Input, must already have had 'FillComplete()' called.
    @param transposeB Input, whether to use transpose of matrix B.
    @param C Result. On entry to this method, it doesn't matter whether
    FillComplete() has already been called on C or not. If it has,
    then C's graph must already contain all nonzero locations that
    will be produced when forming the product A*B. On exit,
    C.FillComplete() will have been called, unless the last argument
    to this function is specified to be false.
    @param call_FillComplete_on_result Optional argument, defaults to true.
    Power users may specify this argument to be false if they *DON'T*
    want this function to call C.FillComplete. (It is often useful
    to allow this function to call C.FillComplete, in cases where
    one or both of the input matrices are rectangular and it is not
    trivial to know which maps to use for the domain- and range-maps.)

*/
  static void Multiply(const Matrix& A, bool transposeA,
                       const Matrix& B, bool transposeB,
                       Matrix& C,
                       bool call_FillComplete_on_result = true,
                       bool doOptimizeStorage           = true,
                       const std::string& label         = std::string(),
                       const RCP<ParameterList>& params = null) {
    TEUCHOS_TEST_FOR_EXCEPTION(transposeA == false && C.getRowMap()->isSameAs(*A.getRowMap()) == false,
                               Xpetra::Exceptions::RuntimeError, "XpetraExt::MatrixMatrix::Multiply: row map of C is not same as row map of A");
    TEUCHOS_TEST_FOR_EXCEPTION(transposeA == true && C.getRowMap()->isSameAs(*A.getDomainMap()) == false,
                               Xpetra::Exceptions::RuntimeError, "XpetraExt::MatrixMatrix::Multiply: row map of C is not same as domain map of A");

    TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), Xpetra::Exceptions::RuntimeError, "A is not fill-completed");
    TEUCHOS_TEST_FOR_EXCEPTION(!B.isFillComplete(), Xpetra::Exceptions::RuntimeError, "B is not fill-completed");

    bool haveMultiplyDoFillComplete = call_FillComplete_on_result && doOptimizeStorage;

    using helpers = Xpetra::Helpers<SC, LO, GO, NO>;

    if (C.getRowMap()->lib() == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT)
      helpers::epetraExtMult(A, transposeA, B, transposeB, C, haveMultiplyDoFillComplete);
#else
      throw(Xpetra::Exceptions::RuntimeError("Xpetra::MatrixMatrix::Multiply requires EpetraExt to be compiled."));
#endif
    } else if (C.getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
      throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra <double,int,int> ETI enabled."));
#else
      if (helpers::isTpetraCrs(A) && helpers::isTpetraCrs(B) && helpers::isTpetraCrs(C)) {
        // All matrices are Crs
        const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpA = helpers::Op2TpetraCrs(A);
        const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpB = helpers::Op2TpetraCrs(B);
        Tpetra::CrsMatrix<SC, LO, GO, NO>& tpC       = helpers::Op2NonConstTpetraCrs(C);

        // 18Feb2013 JJH I'm reenabling the code that allows the matrix matrix multiply to do the fillComplete.
        // Previously, Tpetra's matrix matrix multiply did not support fillComplete.
        Tpetra::MatrixMatrix::Multiply(tpA, transposeA, tpB, transposeB, tpC, haveMultiplyDoFillComplete, label, params);
      } else if (helpers::isTpetraBlockCrs(A) && helpers::isTpetraBlockCrs(B)) {
        // All matrices are BlockCrs (except maybe Ac)
        // FIXME: For the moment we're just going to clobber the innards of Ac, so no reuse. Once we have a reuse kernel,
        // we'll need to think about refactoring BlockCrs so we can do something smarter here.

        if (!A.getRowMap()->getComm()->getRank())
          std::cout << "WARNING: Using inefficient BlockCrs Multiply Placeholder" << std::endl;

        const Tpetra::BlockCrsMatrix<SC, LO, GO, NO>& tpA = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraBlockCrs(A);
        const Tpetra::BlockCrsMatrix<SC, LO, GO, NO>& tpB = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraBlockCrs(B);
        using CRS                                         = Tpetra::CrsMatrix<SC, LO, GO, NO>;
        RCP<const CRS> Acrs                               = Tpetra::convertToCrsMatrix(tpA);
        RCP<const CRS> Bcrs                               = Tpetra::convertToCrsMatrix(tpB);

        // We need the global constants to do the copy back to BlockCrs
        RCP<ParameterList> new_params;
        if (!params.is_null()) {
          new_params = rcp(new Teuchos::ParameterList(*params));
          new_params->set("compute global constants", true);
        }

        // FIXME: The lines below only works because we're assuming Ac is Point
        RCP<CRS> tempAc = Teuchos::rcp(new CRS(Acrs->getRowMap(), 0));
        Tpetra::MatrixMatrix::Multiply(*Acrs, transposeA, *Bcrs, transposeB, *tempAc, haveMultiplyDoFillComplete, label, new_params);

        // Temporary output matrix
        RCP<Tpetra::BlockCrsMatrix<SC, LO, GO, NO> > Ac_t       = Tpetra::convertToBlockCrsMatrix(*tempAc, A.GetStorageBlockSize());
        RCP<Xpetra::TpetraBlockCrsMatrix<SC, LO, GO, NO> > Ac_x = Teuchos::rcp(new Xpetra::TpetraBlockCrsMatrix<SC, LO, GO, NO>(Ac_t));
        RCP<Xpetra::CrsMatrix<SC, LO, GO, NO> > Ac_p            = Ac_x;

        // We can now cheat and replace the innards of Ac
        RCP<Xpetra::CrsMatrixWrap<SC, LO, GO, NO> > Ac_w = Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<SC, LO, GO, NO> >(Teuchos::rcpFromRef(C));
        Ac_w->replaceCrsMatrix(Ac_p);
      } else {
        // Mix and match
        TEUCHOS_TEST_FOR_EXCEPTION(1, Exceptions::RuntimeError, "Mix-and-match Crs/BlockCrs Multiply not currently supported");
      }
#endif
#else
      throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra."));
#endif
    }

    if (call_FillComplete_on_result && !haveMultiplyDoFillComplete) {
      RCP<Teuchos::ParameterList> fillParams = rcp(new Teuchos::ParameterList());
      fillParams->set("Optimize Storage", doOptimizeStorage);
      C.fillComplete((transposeB) ? B.getRangeMap() : B.getDomainMap(),
                     (transposeA) ? A.getDomainMap() : A.getRangeMap(),
                     fillParams);
    }

    // transfer striding information
    RCP<Matrix> rcpA = Teuchos::rcp_const_cast<Matrix>(Teuchos::rcpFromRef(A));
    RCP<Matrix> rcpB = Teuchos::rcp_const_cast<Matrix>(Teuchos::rcpFromRef(B));
    C.CreateView("stridedMaps", rcpA, transposeA, rcpB, transposeB);  // TODO use references instead of RCPs
  }                                                                   // end Multiply

  /**
    @brief Helper function to do matrix-matrix multiply

    Given CrsMatrix objects A, B and C, form the product C = A*B.
    In a parallel setting, A and B need not have matching distributions,
    but C needs to have the same row-map as A (if transposeA is false).
    At this time C=AT*B and C=A*BT are known to not work. However,
    C=A*B and C=AT*BT are known to work, Kurtis Nusbaum 03/24/2011

    @param A Input, must already have had 'FillComplete()' called.
    @param transposeA Input, whether to use transpose of matrix A.
    @param B Input, must already have had 'FillComplete()' called.
    @param transposeB Input, whether to use transpose of matrix B.
    @param C Result. If Teuchos::null, a new CrsMatrix is created with optimal number of nnz per row.
    @param call_FillComplete_on_result Optional argument, defaults to true.
    Power users may specify this argument to be false if they *DON'T*
    want this function to call C.FillComplete. (It is often useful
    to allow this function to call C.FillComplete, in cases where
    one or both of the input matrices are rectangular and it is not
    trivial to know which maps to use for the domain- and range-maps.)

*/
  static RCP<Matrix> Multiply(const Matrix& A, bool transposeA,
                              const Matrix& B, bool transposeB,
                              RCP<Matrix> C_in,
                              Teuchos::FancyOStream& fos,
                              bool doFillComplete              = true,
                              bool doOptimizeStorage           = true,
                              const std::string& label         = std::string(),
                              const RCP<ParameterList>& params = null) {
    TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), Exceptions::RuntimeError, "A is not fill-completed");
    TEUCHOS_TEST_FOR_EXCEPTION(!B.isFillComplete(), Exceptions::RuntimeError, "B is not fill-completed");

    // Optimization using ML Multiply when available and requested
    // This feature is currently not supported. We would have to introduce the HAVE_XPETRA_ML_MMM flag
#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT) && defined(HAVE_XPETRA_ML_MMM)
    if (B.getDomainMap()->lib() == Xpetra::UseEpetra && !transposeA && !transposeB) {
      RCP<const Epetra_CrsMatrix> epA = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(rcpFromRef(A));
      RCP<const Epetra_CrsMatrix> epB = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(rcpFromRef(B));
      RCP<Epetra_CrsMatrix> epC       = MLTwoMatrixMultiply(*epA, *epB, fos);

      RCP<Matrix> C = Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap<SC, LO, GO, NO>(epC);
      if (doFillComplete) {
        RCP<Teuchos::ParameterList> fillParams = rcp(new Teuchos::ParameterList());
        fillParams->set("Optimize Storage", doOptimizeStorage);
        C->fillComplete(B.getDomainMap(), A.getRangeMap(), fillParams);
      }

      // Fill strided maps information
      // This is necessary since the ML matrix matrix multiplication routine has no handling for this
      // TODO: move this call to MLMultiply...
      C->CreateView("stridedMaps", rcpFromRef(A), transposeA, rcpFromRef(B), transposeB);

      return C;
    }
#endif  // EPETRA + EPETRAEXT + ML

    // Default case: Xpetra Multiply
    RCP<Matrix> C = C_in;

    if (C == Teuchos::null) {
      double nnzPerRow = Teuchos::as<double>(0);

#if 0
        if (A.getDomainMap()->lib() == Xpetra::UseTpetra) {
          // For now, follow what ML and Epetra do.
          GO numRowsA = A.getGlobalNumRows();
          GO numRowsB = B.getGlobalNumRows();
          nnzPerRow = sqrt(Teuchos::as<double>(A.getGlobalNumEntries())/numRowsA) +
              sqrt(Teuchos::as<double>(B.getGlobalNumEntries())/numRowsB) - 1;
          nnzPerRow *=  nnzPerRow;
          double totalNnz = nnzPerRow * A.getGlobalNumRows() * 0.75 + 100;
          double minNnz = Teuchos::as<double>(1.2 * A.getGlobalNumEntries());
          if (totalNnz < minNnz)
            totalNnz = minNnz;
          nnzPerRow = totalNnz / A.getGlobalNumRows();

          fos << "Matrix product nnz per row estimate = " << Teuchos::as<LO>(nnzPerRow) << std::endl;
        }
#endif

      if (transposeA)
        C = MatrixFactory::Build(A.getDomainMap(), Teuchos::as<LO>(nnzPerRow));
      else
        C = MatrixFactory::Build(A.getRowMap(), Teuchos::as<LO>(nnzPerRow));

    } else {
      C->resumeFill();  // why this is not done inside of Tpetra MxM?
      fos << "Reuse C pattern" << std::endl;
    }

    Multiply(A, transposeA, B, transposeB, *C, doFillComplete, doOptimizeStorage, label, params);  // call Multiply routine from above

    return C;
  }

  /*! @brief Helper function to do matrix-matrix multiply

    Returns C = AB.

    @param A left matrix
    @param transposeA if true, use the transpose of A
    @param B right matrix
    @param transposeB if true, use the transpose of B
    @param callFillCompleteOnResult if true, the resulting matrix should be fillComplete'd
    */
  static RCP<Matrix> Multiply(const Matrix& A, bool transposeA,
                              const Matrix& B, bool transposeB,
                              Teuchos::FancyOStream& fos,
                              bool callFillCompleteOnResult    = true,
                              bool doOptimizeStorage           = true,
                              const std::string& label         = std::string(),
                              const RCP<ParameterList>& params = null) {
    return Multiply(A, transposeA, B, transposeB, Teuchos::null, fos, callFillCompleteOnResult, doOptimizeStorage, label, params);
  }

#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT)
  // Michael Gee's MLMultiply
  static RCP<Epetra_CrsMatrix> MLTwoMatrixMultiply(const Epetra_CrsMatrix& epA,
                                                   const Epetra_CrsMatrix& epB,
                                                   Teuchos::FancyOStream& fos) {
#if defined(HAVE_XPETRA_ML_MMM)  // Note: this is currently not supported
    ML_Comm* comm;
    ML_Comm_Create(&comm);
    fos << "****** USING ML's MATRIX MATRIX MULTIPLY ******" << std::endl;
#ifdef HAVE_MPI
    // ML_Comm uses MPI_COMM_WORLD, so try to use the same communicator as epA.
    const Epetra_MpiComm* Mcomm = dynamic_cast<const Epetra_MpiComm*>(&(epA.Comm()));
    if (Mcomm)
      ML_Comm_Set_UsrComm(comm, Mcomm->GetMpiComm());
#endif
    // in order to use ML, there must be no indices missing from the matrix column maps.
    EpetraExt::CrsMatrix_SolverMap Atransform;
    EpetraExt::CrsMatrix_SolverMap Btransform;
    const Epetra_CrsMatrix& A = Atransform(const_cast<Epetra_CrsMatrix&>(epA));
    const Epetra_CrsMatrix& B = Btransform(const_cast<Epetra_CrsMatrix&>(epB));

    if (!A.Filled()) throw Exceptions::RuntimeError("A has to be FillCompleted");
    if (!B.Filled()) throw Exceptions::RuntimeError("B has to be FillCompleted");

    // create ML operators from EpetraCrsMatrix
    ML_Operator* ml_As = ML_Operator_Create(comm);
    ML_Operator* ml_Bs = ML_Operator_Create(comm);
    ML_Operator_WrapEpetraCrsMatrix(const_cast<Epetra_CrsMatrix*>(&A), ml_As);  // Should we test if the lightweight wrapper is actually used or if WrapEpetraCrsMatrix fall backs to the heavy one?
    ML_Operator_WrapEpetraCrsMatrix(const_cast<Epetra_CrsMatrix*>(&B), ml_Bs);
    ML_Operator* ml_AtimesB = ML_Operator_Create(comm);
    {
      Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("ML_2matmult kernel"));
      ML_2matmult(ml_As, ml_Bs, ml_AtimesB, ML_CSR_MATRIX);  // do NOT use ML_EpetraCRS_MATRIX!!!
    }
    ML_Operator_Destroy(&ml_As);
    ML_Operator_Destroy(&ml_Bs);

    // For ml_AtimesB we have to reconstruct the column map in global indexing,
    // The following is going down to the salt-mines of ML ...
    // note: we use integers, since ML only works for Epetra...
    int N_local                = ml_AtimesB->invec_leng;
    ML_CommInfoOP* getrow_comm = ml_AtimesB->getrow->pre_comm;
    if (!getrow_comm) throw(Exceptions::RuntimeError("ML_Operator does not have a CommInfo"));
    ML_Comm* comm_AB = ml_AtimesB->comm;  // get comm object
    if (N_local != B.DomainMap().NumMyElements())
      throw(Exceptions::RuntimeError("Mismatch in local row dimension between ML and Epetra"));
    int N_rcvd = 0;
    int N_send = 0;
    int flag   = 0;
    for (int i = 0; i < getrow_comm->N_neighbors; i++) {
      N_rcvd += (getrow_comm->neighbors)[i].N_rcv;
      N_send += (getrow_comm->neighbors)[i].N_send;
      if (((getrow_comm->neighbors)[i].N_rcv != 0) &&
          ((getrow_comm->neighbors)[i].rcv_list != NULL)) flag = 1;
    }
    // For some unknown reason, ML likes to have stuff one larger than
    // neccessary...
    std::vector<double> dtemp(N_local + N_rcvd + 1);  // "double" vector for comm function
    std::vector<int> cmap(N_local + N_rcvd + 1);      // vector for gids
    for (int i = 0; i < N_local; ++i) {
      cmap[i]  = B.DomainMap().GID(i);
      dtemp[i] = (double)cmap[i];
    }
    ML_cheap_exchange_bdry(&dtemp[0], getrow_comm, N_local, N_send, comm_AB);  // do communication
    if (flag) {                                                                // process received data
      int count           = N_local;
      const int neighbors = getrow_comm->N_neighbors;
      for (int i = 0; i < neighbors; i++) {
        const int nrcv = getrow_comm->neighbors[i].N_rcv;
        for (int j = 0; j < nrcv; j++)
          cmap[getrow_comm->neighbors[i].rcv_list[j]] = (int)dtemp[count++];
      }
    } else {
      for (int i = 0; i < N_local + N_rcvd; ++i)
        cmap[i] = (int)dtemp[i];
    }
    dtemp.clear();  // free double array

    // we can now determine a matching column map for the result
    Epetra_Map gcmap(-1, N_local + N_rcvd, &cmap[0], B.ColMap().IndexBase(), A.Comm());

    int allocated = 0;
    int rowlength;
    double* val = NULL;
    int* bindx  = NULL;

    const int myrowlength    = A.RowMap().NumMyElements();
    const Epetra_Map& rowmap = A.RowMap();

    // Determine the maximum bandwith for the result matrix.
    // replaces the old, very(!) memory-consuming guess:
    // int guessnpr = A.MaxNumEntries()*B.MaxNumEntries();
    int educatedguess = 0;
    for (int i = 0; i < myrowlength; ++i) {
      // get local row
      ML_get_matrix_row(ml_AtimesB, 1, &i, &allocated, &bindx, &val, &rowlength, 0);
      if (rowlength > educatedguess)
        educatedguess = rowlength;
    }

    // allocate our result matrix and fill it
    RCP<Epetra_CrsMatrix> result = rcp(new Epetra_CrsMatrix(::Copy, A.RangeMap(), gcmap, educatedguess, false));

    std::vector<int> gcid(educatedguess);
    for (int i = 0; i < myrowlength; ++i) {
      const int grid = rowmap.GID(i);
      // get local row
      ML_get_matrix_row(ml_AtimesB, 1, &i, &allocated, &bindx, &val, &rowlength, 0);
      if (!rowlength) continue;
      if ((int)gcid.size() < rowlength) gcid.resize(rowlength);
      for (int j = 0; j < rowlength; ++j) {
        gcid[j] = gcmap.GID(bindx[j]);
        if (gcid[j] < 0)
          throw Exceptions::RuntimeError("Error: cannot find gcid!");
      }
      int err = result->InsertGlobalValues(grid, rowlength, val, &gcid[0]);
      if (err != 0 && err != 1) {
        std::ostringstream errStr;
        errStr << "Epetra_CrsMatrix::InsertGlobalValues returned err=" << err;
        throw Exceptions::RuntimeError(errStr.str());
      }
    }
    // free memory
    if (bindx) ML_free(bindx);
    if (val) ML_free(val);
    ML_Operator_Destroy(&ml_AtimesB);
    ML_Comm_Destroy(&comm);

    return result;
#else  // no MUELU_ML
    (void)epA;
    (void)epB;
    (void)fos;
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
                               "No ML multiplication available. This feature is currently not supported by Xpetra.");
    TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
#endif
  }
#endif  // ifdef HAVE_XPETRA_EPETRAEXT

  /*! @brief Helper function to do matrix-matrix multiply "in-place"

    Returns RCP to non-constant Xpetra::BlockedCrsMatrix.

    @param A left matrix
    @param transposeA if true, use the transpose of A
    @param B right matrix
    @param transposeB if true, use the transpose of B
    @param doOptimizeStorage if true, the resulting matrix should be fillComplete'd
    */
  static RCP<BlockedCrsMatrix> TwoMatrixMultiplyBlock(const BlockedCrsMatrix& A, bool transposeA,
                                                      const BlockedCrsMatrix& B, bool transposeB,
                                                      Teuchos::FancyOStream& fos,
                                                      bool doFillComplete    = true,
                                                      bool doOptimizeStorage = true) {
    TEUCHOS_TEST_FOR_EXCEPTION(transposeA || transposeB, Exceptions::RuntimeError,
                               "TwoMatrixMultiply for BlockedCrsMatrix not implemented for transposeA==true or transposeB==true");

    // Preconditions
    TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), Exceptions::RuntimeError, "A is not fill-completed");
    TEUCHOS_TEST_FOR_EXCEPTION(!B.isFillComplete(), Exceptions::RuntimeError, "B is not fill-completed");

    RCP<const MapExtractor> rgmapextractor = A.getRangeMapExtractor();
    RCP<const MapExtractor> domapextractor = B.getDomainMapExtractor();

    RCP<BlockedCrsMatrix> C = rcp(new BlockedCrsMatrix(rgmapextractor, domapextractor, 33 /* TODO fix me */));

    for (size_t i = 0; i < A.Rows(); ++i) {    // loop over all block rows of A
      for (size_t j = 0; j < B.Cols(); ++j) {  // loop over all block columns of B
        RCP<Matrix> Cij;

        for (size_t l = 0; l < B.Rows(); ++l) {  // loop for calculating entry C_{ij}
          RCP<Matrix> crmat1 = A.getMatrix(i, l);
          RCP<Matrix> crmat2 = B.getMatrix(l, j);

          if (crmat1.is_null() || crmat2.is_null())
            continue;

          // try unwrapping 1x1 blocked matrix
          {
            auto unwrap1 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(crmat1);
            auto unwrap2 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(crmat2);

            if (unwrap1.is_null() != unwrap2.is_null()) {
              if (unwrap1 != Teuchos::null && unwrap1->Rows() == 1 && unwrap1->Cols() == 1)
                crmat1 = unwrap1->getCrsMatrix();
              if (unwrap2 != Teuchos::null && unwrap2->Rows() == 1 && unwrap2->Cols() == 1)
                crmat2 = unwrap2->getCrsMatrix();
            }
          }

          RCP<CrsMatrixWrap> crop1 = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(crmat1);
          RCP<CrsMatrixWrap> crop2 = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(crmat2);
          TEUCHOS_TEST_FOR_EXCEPTION(crop1.is_null() != crop2.is_null(), Xpetra::Exceptions::RuntimeError,
                                     "A and B must be either both (compatible) BlockedCrsMatrix objects or both CrsMatrixWrap objects.");

          // Forcibly compute the global constants if we don't have them (only works for real CrsMatrices, not nested blocks)
          if (!crop1.is_null())
            Teuchos::rcp_const_cast<CrsGraph>(crmat1->getCrsGraph())->computeGlobalConstants();
          if (!crop2.is_null())
            Teuchos::rcp_const_cast<CrsGraph>(crmat2->getCrsGraph())->computeGlobalConstants();

          TEUCHOS_TEST_FOR_EXCEPTION(!crmat1->haveGlobalConstants(), Exceptions::RuntimeError,
                                     "crmat1 does not have global constants");
          TEUCHOS_TEST_FOR_EXCEPTION(!crmat2->haveGlobalConstants(), Exceptions::RuntimeError,
                                     "crmat2 does not have global constants");

          if (crmat1->getGlobalNumEntries() == 0 || crmat2->getGlobalNumEntries() == 0)
            continue;

          // temporary matrix containing result of local block multiplication
          RCP<Matrix> temp = Teuchos::null;

          if (crop1 != Teuchos::null && crop2 != Teuchos::null)
            temp = Multiply(*crop1, false, *crop2, false, fos);
          else {
            RCP<BlockedCrsMatrix> bop1 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(crmat1);
            RCP<BlockedCrsMatrix> bop2 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(crmat2);
            TEUCHOS_TEST_FOR_EXCEPTION(bop1.is_null() == true, Xpetra::Exceptions::BadCast, "A is not a BlockedCrsMatrix. (TwoMatrixMultiplyBlock)");
            TEUCHOS_TEST_FOR_EXCEPTION(bop2.is_null() == true, Xpetra::Exceptions::BadCast, "B is not a BlockedCrsMatrix. (TwoMatrixMultiplyBlock)");
            TEUCHOS_TEST_FOR_EXCEPTION(bop1->Cols() != bop2->Rows(), Xpetra::Exceptions::RuntimeError, "A has " << bop1->Cols() << " columns and B has " << bop2->Rows() << " rows. Matrices are not compatible! (TwoMatrixMultiplyBlock)");
            TEUCHOS_TEST_FOR_EXCEPTION(bop1->getDomainMap()->isSameAs(*(bop2->getRangeMap())) == false, Xpetra::Exceptions::RuntimeError, "Domain map of A is not the same as range map of B. Matrices are not compatible! (TwoMatrixMultiplyBlock)");

            // recursive multiplication call
            temp = TwoMatrixMultiplyBlock(*bop1, transposeA, *bop2, transposeB, fos, doFillComplete, doOptimizeStorage);

            RCP<BlockedCrsMatrix> btemp = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(temp);
            TEUCHOS_TEST_FOR_EXCEPTION(btemp->Rows() != bop1->Rows(), Xpetra::Exceptions::RuntimeError, "Number of block rows of local blocked operator is " << btemp->Rows() << " but should be " << bop1->Rows() << ". (TwoMatrixMultiplyBlock)");
            TEUCHOS_TEST_FOR_EXCEPTION(btemp->Cols() != bop2->Cols(), Xpetra::Exceptions::RuntimeError, "Number of block cols of local blocked operator is " << btemp->Cols() << " but should be " << bop2->Cols() << ". (TwoMatrixMultiplyBlock)");
            TEUCHOS_TEST_FOR_EXCEPTION(btemp->getRangeMapExtractor()->getFullMap()->isSameAs(*(bop1->getRangeMapExtractor()->getFullMap())) == false, Xpetra::Exceptions::RuntimeError, "Range map of local blocked operator should be same as first operator. (TwoMatrixMultiplyBlock)");
            TEUCHOS_TEST_FOR_EXCEPTION(btemp->getDomainMapExtractor()->getFullMap()->isSameAs(*(bop2->getDomainMapExtractor()->getFullMap())) == false, Xpetra::Exceptions::RuntimeError, "Domain map of local blocked operator should be same as second operator. (TwoMatrixMultiplyBlock)");
            TEUCHOS_TEST_FOR_EXCEPTION(btemp->getRangeMapExtractor()->getThyraMode() != bop1->getRangeMapExtractor()->getThyraMode(), Xpetra::Exceptions::RuntimeError, "Thyra mode of local range map extractor incompatible with range map extractor of A (TwoMatrixMultiplyBlock)");
            TEUCHOS_TEST_FOR_EXCEPTION(btemp->getDomainMapExtractor()->getThyraMode() != bop2->getDomainMapExtractor()->getThyraMode(), Xpetra::Exceptions::RuntimeError, "Thyra mode of local domain map extractor incompatible with domain map extractor of B (TwoMatrixMultiplyBlock)");
          }

          TEUCHOS_TEST_FOR_EXCEPTION(temp->isFillComplete() == false, Xpetra::Exceptions::RuntimeError, "Local block is not filled. (TwoMatrixMultiplyBlock)");

          RCP<Matrix> addRes = null;
          if (Cij.is_null())
            Cij = temp;
          else {
            MatrixMatrix::TwoMatrixAdd(*temp, false, 1.0, *Cij, false, 1.0, addRes, fos);
            Cij = addRes;
          }
        }

        if (!Cij.is_null()) {
          if (Cij->isFillComplete())
            Cij->resumeFill();
          Cij->fillComplete(B.getDomainMap(j), A.getRangeMap(i));
          C->setMatrix(i, j, Cij);
        } else {
          C->setMatrix(i, j, Teuchos::null);
        }
      }
    }

    if (doFillComplete)
      C->fillComplete();  // call default fillComplete for BlockCrsMatrixWrap objects

    return C;
  }  // TwoMatrixMultiplyBlock

  /*! @brief Helper function to calculate B = alpha*A + beta*B.

    @param A      left matrix operand
    @param transposeA indicate whether to use transpose of A
    @param alpha  scalar multiplier for A
    @param B      right matrix operand
    @param beta   scalar multiplier for B

    @return sum in B.

    Note that B does not have to be fill-completed.
    */
  static void TwoMatrixAdd(const Matrix& A, bool transposeA, SC alpha, Matrix& B, SC beta) {
    TEUCHOS_TEST_FOR_EXCEPTION(!(A.getRowMap()->isSameAs(*B.getRowMap())), Exceptions::Incompatible,
                               "TwoMatrixAdd: matrix row maps are not the same.");

    if (A.getRowMap()->lib() == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT)
      const Epetra_CrsMatrix& epA = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(A);
      Epetra_CrsMatrix& epB       = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstEpetraCrs(B);

      // FIXME is there a bug if beta=0?
      int rv = EpetraExt::MatrixMatrix::Add(epA, transposeA, alpha, epB, beta);

      if (rv != 0)
        throw Exceptions::RuntimeError("EpetraExt::MatrixMatrix::Add return value " + Teuchos::toString(rv));
      std::ostringstream buf;
#else
      throw Exceptions::RuntimeError("Xpetra must be compiled with EpetraExt.");
#endif
    } else if (A.getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
      throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra GO=int enabled."));
#else
      const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpA = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraCrs(A);
      Tpetra::CrsMatrix<SC, LO, GO, NO>& tpB = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstTpetraCrs(B);

      Tpetra::MatrixMatrix::Add(tpA, transposeA, alpha, tpB, beta);
#endif
#else
      throw Exceptions::RuntimeError("Xpetra must be compiled with Tpetra.");
#endif
    }
  }  // MatrixMatrix::TwoMatrixAdd()

  /*! @brief Helper function to calculate C = alpha*A + beta*B.

    @param A          left matrix operand
    @param transposeA indicate whether to use transpose of A
    @param alpha      scalar multiplier for A, defaults to 1.0
    @param B          right matrix operand
    @param transposeB indicate whether to use transpose of B
    @param beta       scalar multiplier for B, defaults to 1.0
    @param C          resulting sum
    @param fos        output stream for printing to screen
    @param AHasFixedNnzPerRow

    It is up to the caller to ensure that the resulting matrix sum is fillComplete'd.
    */
  static void TwoMatrixAdd(const Matrix& A, bool transposeA, const SC& alpha,
                           const Matrix& B, bool transposeB, const SC& beta,
                           RCP<Matrix>& C, Teuchos::FancyOStream& fos, bool AHasFixedNnzPerRow = false) {
    using helpers                       = Xpetra::Helpers<SC, LO, GO, NO>;
    RCP<const Matrix> rcpA              = Teuchos::rcpFromRef(A);
    RCP<const Matrix> rcpB              = Teuchos::rcpFromRef(B);
    RCP<const BlockedCrsMatrix> rcpBopA = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrix>(rcpA);
    RCP<const BlockedCrsMatrix> rcpBopB = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrix>(rcpB);

    if (rcpBopA == Teuchos::null && rcpBopB == Teuchos::null) {
      if (!(A.getRowMap()->isSameAs(*(B.getRowMap()))))
        throw Exceptions::Incompatible("TwoMatrixAdd: matrix row maps are not the same.");

      auto lib = A.getRowMap()->lib();
      if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT)
        const Epetra_CrsMatrix& epA = helpers::Op2EpetraCrs(A);
        const Epetra_CrsMatrix& epB = helpers::Op2EpetraCrs(B);
        if (C.is_null()) {
          size_t maxNzInA     = 0;
          size_t maxNzInB     = 0;
          size_t numLocalRows = 0;
          if (A.isFillComplete() && B.isFillComplete()) {
            maxNzInA     = A.getLocalMaxNumRowEntries();
            maxNzInB     = B.getLocalMaxNumRowEntries();
            numLocalRows = A.getLocalNumRows();
          }

          if (maxNzInA == 1 || maxNzInB == 1 || AHasFixedNnzPerRow) {
            // first check if either A or B has at most 1 nonzero per row
            // the case of both having at most 1 nz per row is handled by the ``else''
            Teuchos::ArrayRCP<size_t> exactNnzPerRow(numLocalRows);

            if ((maxNzInA == 1 && maxNzInB > 1) || AHasFixedNnzPerRow) {
              for (size_t i = 0; i < numLocalRows; ++i)
                exactNnzPerRow[i] = B.getNumEntriesInLocalRow(Teuchos::as<LO>(i)) + maxNzInA;

            } else {
              for (size_t i = 0; i < numLocalRows; ++i)
                exactNnzPerRow[i] = A.getNumEntriesInLocalRow(Teuchos::as<LO>(i)) + maxNzInB;
            }

            fos << "MatrixMatrix::TwoMatrixAdd : special case detected (one matrix has a fixed nnz per row)"
                << ", using static profiling" << std::endl;
            C = rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(A.getRowMap(), exactNnzPerRow));

          } else {
            // general case
            LO maxPossibleNNZ = A.getLocalMaxNumRowEntries() + B.getLocalMaxNumRowEntries();
            C                 = rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(A.getRowMap(), maxPossibleNNZ));
          }
          if (transposeB)
            fos << "MatrixMatrix::TwoMatrixAdd : ** WARNING ** estimate could be badly wrong because second summand is transposed" << std::endl;
        }
        RCP<Epetra_CrsMatrix> epC = helpers::Op2NonConstEpetraCrs(C);
        Epetra_CrsMatrix* ref2epC = &*epC;  // to avoid a compiler error...

        // FIXME is there a bug if beta=0?
        int rv = EpetraExt::MatrixMatrix::Add(epA, transposeA, alpha, epB, transposeB, beta, ref2epC);

        if (rv != 0)
          throw Exceptions::RuntimeError("EpetraExt::MatrixMatrix::Add return value of " + Teuchos::toString(rv));
#else
        throw Exceptions::RuntimeError("MueLu must be compile with EpetraExt.");
#endif
      } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
        throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra GO=int enabled."));
#else
        using tcrs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NO>;
        const tcrs_matrix_type& tpA = helpers::Op2TpetraCrs(A);
        const tcrs_matrix_type& tpB = helpers::Op2TpetraCrs(B);
        C = helpers::tpetraAdd(tpA, transposeA, alpha, tpB, transposeB, beta);
#endif
#else
        throw Exceptions::RuntimeError("Xpetra must be compile with Tpetra.");
#endif
      }

      ///////////////////////// EXPERIMENTAL
      if (A.IsView("stridedMaps")) C->CreateView("stridedMaps", rcpFromRef(A));
      if (B.IsView("stridedMaps")) C->CreateView("stridedMaps", rcpFromRef(B));
      ///////////////////////// EXPERIMENTAL
    }
    // the first matrix is of type CrsMatrixWrap, the second is a blocked operator
    else if (rcpBopA == Teuchos::null && rcpBopB != Teuchos::null) {
      RCP<const MapExtractor> rgmapextractor = rcpBopB->getRangeMapExtractor();
      RCP<const MapExtractor> domapextractor = rcpBopB->getDomainMapExtractor();

      C                        = rcp(new BlockedCrsMatrix(rgmapextractor, domapextractor, 33 /* TODO fix me */));
      RCP<BlockedCrsMatrix> bC = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(C);

      size_t i = 0;
      for (size_t j = 0; j < rcpBopB->Cols(); ++j) {  // loop over all block columns of B
        RCP<Matrix> Cij = Teuchos::null;
        if (rcpA != Teuchos::null &&
            rcpBopB->getMatrix(i, j) != Teuchos::null) {
          // recursive call
          TwoMatrixAdd(*rcpA, transposeA, alpha,
                       *(rcpBopB->getMatrix(i, j)), transposeB, beta,
                       Cij, fos, AHasFixedNnzPerRow);
        } else if (rcpA == Teuchos::null &&
                   rcpBopB->getMatrix(i, j) != Teuchos::null) {
          Cij = rcpBopB->getMatrix(i, j);
        } else if (rcpA != Teuchos::null &&
                   rcpBopB->getMatrix(i, j) == Teuchos::null) {
          Cij = Teuchos::rcp_const_cast<Matrix>(rcpA);
        } else {
          Cij = Teuchos::null;
        }

        if (!Cij.is_null()) {
          if (Cij->isFillComplete())
            Cij->resumeFill();
          Cij->fillComplete();
          bC->setMatrix(i, j, Cij);
        } else {
          bC->setMatrix(i, j, Teuchos::null);
        }
      }  // loop over columns j
    }
    // the second matrix is of type CrsMatrixWrap, the first is a blocked operator
    else if (rcpBopA != Teuchos::null && rcpBopB == Teuchos::null) {
      RCP<const MapExtractor> rgmapextractor = rcpBopA->getRangeMapExtractor();
      RCP<const MapExtractor> domapextractor = rcpBopA->getDomainMapExtractor();

      C                        = rcp(new BlockedCrsMatrix(rgmapextractor, domapextractor, 33 /* TODO fix me */));
      RCP<BlockedCrsMatrix> bC = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(C);

      size_t j = 0;
      for (size_t i = 0; i < rcpBopA->Rows(); ++i) {  // loop over all block rows of A
        RCP<Matrix> Cij = Teuchos::null;
        if (rcpBopA->getMatrix(i, j) != Teuchos::null &&
            rcpB != Teuchos::null) {
          // recursive call
          TwoMatrixAdd(*(rcpBopA->getMatrix(i, j)), transposeA, alpha,
                       *rcpB, transposeB, beta,
                       Cij, fos, AHasFixedNnzPerRow);
        } else if (rcpBopA->getMatrix(i, j) == Teuchos::null &&
                   rcpB != Teuchos::null) {
          Cij = Teuchos::rcp_const_cast<Matrix>(rcpB);
        } else if (rcpBopA->getMatrix(i, j) != Teuchos::null &&
                   rcpB == Teuchos::null) {
          Cij = rcpBopA->getMatrix(i, j);
        } else {
          Cij = Teuchos::null;
        }

        if (!Cij.is_null()) {
          if (Cij->isFillComplete())
            Cij->resumeFill();
          Cij->fillComplete();
          bC->setMatrix(i, j, Cij);
        } else {
          bC->setMatrix(i, j, Teuchos::null);
        }
      }  // loop over rows i
    } else {
      // This is the version for blocked matrices

      // check the compatibility of the blocked operators
      TEUCHOS_TEST_FOR_EXCEPTION(rcpBopA.is_null() == true, Xpetra::Exceptions::BadCast, "A is not a BlockedCrsMatrix. (TwoMatrixAdd)");
      TEUCHOS_TEST_FOR_EXCEPTION(rcpBopB.is_null() == true, Xpetra::Exceptions::BadCast, "B is not a BlockedCrsMatrix. (TwoMatrixAdd)");
      TEUCHOS_TEST_FOR_EXCEPTION(rcpBopA->Rows() != rcpBopB->Rows(), Xpetra::Exceptions::RuntimeError, "A has " << rcpBopA->Rows() << " rows and B has " << rcpBopA->Rows() << " rows. Matrices are not compatible! (TwoMatrixAdd)");
      TEUCHOS_TEST_FOR_EXCEPTION(rcpBopA->Rows() != rcpBopB->Rows(), Xpetra::Exceptions::RuntimeError, "A has " << rcpBopA->Cols() << " cols and B has " << rcpBopA->Cols() << " cols. Matrices are not compatible! (TwoMatrixAdd)");
      TEUCHOS_TEST_FOR_EXCEPTION(rcpBopA->getRangeMap()->isSameAs(*(rcpBopB->getRangeMap())) == false, Xpetra::Exceptions::RuntimeError, "Range map of A is not the same as range map of B. Matrices are not compatible! (TwoMatrixAdd)");
      TEUCHOS_TEST_FOR_EXCEPTION(rcpBopA->getDomainMap()->isSameAs(*(rcpBopB->getDomainMap())) == false, Xpetra::Exceptions::RuntimeError, "Domain map of A is not the same as domain map of B. Matrices are not compatible! (TwoMatrixAdd)");
      TEUCHOS_TEST_FOR_EXCEPTION(rcpBopA->getRangeMapExtractor()->getThyraMode() != rcpBopB->getRangeMapExtractor()->getThyraMode(), Xpetra::Exceptions::RuntimeError, "Different Thyra/Xpetra style gids in RangeMapExtractor of A and B. Matrices are not compatible! (TwoMatrixAdd)");
      TEUCHOS_TEST_FOR_EXCEPTION(rcpBopA->getDomainMapExtractor()->getThyraMode() != rcpBopB->getDomainMapExtractor()->getThyraMode(), Xpetra::Exceptions::RuntimeError, "Different Thyra/Xpetra style gids in DomainMapExtractor of A and B. Matrices are not compatible! (TwoMatrixAdd)");

      RCP<const MapExtractor> rgmapextractor = rcpBopA->getRangeMapExtractor();
      RCP<const MapExtractor> domapextractor = rcpBopB->getDomainMapExtractor();

      C                        = rcp(new BlockedCrsMatrix(rgmapextractor, domapextractor, 33 /* TODO fix me */));
      RCP<BlockedCrsMatrix> bC = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(C);

      for (size_t i = 0; i < rcpBopA->Rows(); ++i) {    // loop over all block rows of A
        for (size_t j = 0; j < rcpBopB->Cols(); ++j) {  // loop over all block columns of B

          RCP<Matrix> Cij = Teuchos::null;
          if (rcpBopA->getMatrix(i, j) != Teuchos::null &&
              rcpBopB->getMatrix(i, j) != Teuchos::null) {
            // recursive call

            TwoMatrixAdd(*(rcpBopA->getMatrix(i, j)), transposeA, alpha,
                         *(rcpBopB->getMatrix(i, j)), transposeB, beta,
                         Cij, fos, AHasFixedNnzPerRow);
          } else if (rcpBopA->getMatrix(i, j) == Teuchos::null &&
                     rcpBopB->getMatrix(i, j) != Teuchos::null) {
            Cij = rcpBopB->getMatrix(i, j);
          } else if (rcpBopA->getMatrix(i, j) != Teuchos::null &&
                     rcpBopB->getMatrix(i, j) == Teuchos::null) {
            Cij = rcpBopA->getMatrix(i, j);
          } else {
            Cij = Teuchos::null;
          }

          if (!Cij.is_null()) {
            if (Cij->isFillComplete())
              Cij->resumeFill();
            // Cij->fillComplete(rcpBopA->getDomainMap(j), rcpBopA->getRangeMap(i));
            Cij->fillComplete();
            bC->setMatrix(i, j, Cij);
          } else {
            bC->setMatrix(i, j, Teuchos::null);
          }
        }  // loop over columns j
      }    // loop over rows i

    }  // end blocked recursive algorithm
  }    // MatrixMatrix::TwoMatrixAdd()
};     // end specialization on SC=double, GO=int and NO=EpetraNode

// specialization MatrixMatrix for SC=double, GO=long long and NO=EptraNode
template <>
class MatrixMatrix<double, int, long long, EpetraNode> {
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef long long GlobalOrdinal;
  typedef EpetraNode Node;
#include "Xpetra_UseShortNames.hpp"

 public:
  /** Given CrsMatrix objects A, B and C, form the product C = A*B.
    In a parallel setting, A and B need not have matching distributions,
    but C needs to have the same row-map as A (if transposeA is false).
    At this time C=AT*B and C=A*BT are known to not work. However,
    C=A*B and C=AT*BT are known to work, Kurtis Nusbaum 03/24/2011

    @param A Input, must already have had 'FillComplete()' called.
    @param transposeA Input, whether to use transpose of matrix A.
    @param B Input, must already have had 'FillComplete()' called.
    @param transposeB Input, whether to use transpose of matrix B.
    @param C Result. On entry to this method, it doesn't matter whether
    FillComplete() has already been called on C or not. If it has,
    then C's graph must already contain all nonzero locations that
    will be produced when forming the product A*B. On exit,
    C.FillComplete() will have been called, unless the last argument
    to this function is specified to be false.
    @param call_FillComplete_on_result Optional argument, defaults to true.
    Power users may specify this argument to be false if they *DON'T*
    want this function to call C.FillComplete. (It is often useful
    to allow this function to call C.FillComplete, in cases where
    one or both of the input matrices are rectangular and it is not
    trivial to know which maps to use for the domain- and range-maps.)

*/
  static void Multiply(const Matrix& A, bool transposeA,
                       const Matrix& B, bool transposeB,
                       Matrix& C,
                       bool call_FillComplete_on_result = true,
                       bool doOptimizeStorage           = true,
                       const std::string& label         = std::string(),
                       const RCP<ParameterList>& params = null) {
    using helpers = Xpetra::Helpers<SC, LO, GO, NO>;
    TEUCHOS_TEST_FOR_EXCEPTION(transposeA == false && C.getRowMap()->isSameAs(*A.getRowMap()) == false,
                               Xpetra::Exceptions::RuntimeError, "XpetraExt::MatrixMatrix::Multiply: row map of C is not same as row map of A");
    TEUCHOS_TEST_FOR_EXCEPTION(transposeA == true && C.getRowMap()->isSameAs(*A.getDomainMap()) == false,
                               Xpetra::Exceptions::RuntimeError, "XpetraExt::MatrixMatrix::Multiply: row map of C is not same as domain map of A");

    TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), Xpetra::Exceptions::RuntimeError, "A is not fill-completed");
    TEUCHOS_TEST_FOR_EXCEPTION(!B.isFillComplete(), Xpetra::Exceptions::RuntimeError, "B is not fill-completed");

    bool haveMultiplyDoFillComplete = call_FillComplete_on_result && doOptimizeStorage;

    if (C.getRowMap()->lib() == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT)
      helpers::epetraExtMult(A, transposeA, B, transposeB, C, haveMultiplyDoFillComplete);
#else
      throw(Xpetra::Exceptions::RuntimeError("Xpetra::MatrixMatrix::Multiply requires EpetraExt to be compiled."));
#endif
    } else if (C.getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
      throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra <double,int,long long, EpetraNode> ETI enabled."));
#else
      if (helpers::isTpetraCrs(A) && helpers::isTpetraCrs(B) && helpers::isTpetraCrs(C)) {
        // All matrices are Crs
        const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpA = helpers::Op2TpetraCrs(A);
        const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpB = helpers::Op2TpetraCrs(B);
        Tpetra::CrsMatrix<SC, LO, GO, NO>& tpC = helpers::Op2NonConstTpetraCrs(C);

        // 18Feb2013 JJH I'm reenabling the code that allows the matrix matrix multiply to do the fillComplete.
        // Previously, Tpetra's matrix matrix multiply did not support fillComplete.
        Tpetra::MatrixMatrix::Multiply(tpA, transposeA, tpB, transposeB, tpC, haveMultiplyDoFillComplete, label, params);
      } else if (helpers::isTpetraBlockCrs(A) && helpers::isTpetraBlockCrs(B)) {
        // All matrices are BlockCrs (except maybe Ac)
        // FIXME: For the moment we're just going to clobber the innards of Ac, so no reuse. Once we have a reuse kernel,
        // we'll need to think about refactoring BlockCrs so we can do something smarter here.

        if (!A.getRowMap()->getComm()->getRank())
          std::cout << "WARNING: Using inefficient BlockCrs Multiply Placeholder" << std::endl;

        const Tpetra::BlockCrsMatrix<SC, LO, GO, NO>& tpA = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraBlockCrs(A);
        const Tpetra::BlockCrsMatrix<SC, LO, GO, NO>& tpB = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraBlockCrs(B);
        using CRS = Tpetra::CrsMatrix<SC, LO, GO, NO>;
        RCP<const CRS> Acrs = Tpetra::convertToCrsMatrix(tpA);
        RCP<const CRS> Bcrs = Tpetra::convertToCrsMatrix(tpB);

        // We need the global constants to do the copy back to BlockCrs
        RCP<ParameterList> new_params;
        if (!params.is_null()) {
          new_params = rcp(new Teuchos::ParameterList(*params));
          new_params->set("compute global constants", true);
        }

        // FIXME: The lines below only works because we're assuming Ac is Point
        RCP<CRS> tempAc = Teuchos::rcp(new CRS(Acrs->getRowMap(), 0));
        Tpetra::MatrixMatrix::Multiply(*Acrs, transposeA, *Bcrs, transposeB, *tempAc, haveMultiplyDoFillComplete, label, new_params);

        // Temporary output matrix
        RCP<Tpetra::BlockCrsMatrix<SC, LO, GO, NO> > Ac_t = Tpetra::convertToBlockCrsMatrix(*tempAc, A.GetStorageBlockSize());
        RCP<Xpetra::TpetraBlockCrsMatrix<SC, LO, GO, NO> > Ac_x = Teuchos::rcp(new Xpetra::TpetraBlockCrsMatrix<SC, LO, GO, NO>(Ac_t));
        RCP<Xpetra::CrsMatrix<SC, LO, GO, NO> > Ac_p = Ac_x;

        // We can now cheat and replace the innards of Ac
        RCP<Xpetra::CrsMatrixWrap<SC, LO, GO, NO> > Ac_w = Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<SC, LO, GO, NO> >(Teuchos::rcpFromRef(C));
        Ac_w->replaceCrsMatrix(Ac_p);
      } else {
        // Mix and match
        TEUCHOS_TEST_FOR_EXCEPTION(1, Exceptions::RuntimeError, "Mix-and-match Crs/BlockCrs Multiply not currently supported");
      }

#endif
#else
      throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra."));
#endif
    }

    if (call_FillComplete_on_result && !haveMultiplyDoFillComplete) {
      RCP<Teuchos::ParameterList> fillParams = rcp(new Teuchos::ParameterList());
      fillParams->set("Optimize Storage", doOptimizeStorage);
      C.fillComplete((transposeB) ? B.getRangeMap() : B.getDomainMap(),
                     (transposeA) ? A.getDomainMap() : A.getRangeMap(),
                     fillParams);
    }

    // transfer striding information
    RCP<Matrix> rcpA = Teuchos::rcp_const_cast<Matrix>(Teuchos::rcpFromRef(A));
    RCP<Matrix> rcpB = Teuchos::rcp_const_cast<Matrix>(Teuchos::rcpFromRef(B));
    C.CreateView("stridedMaps", rcpA, transposeA, rcpB, transposeB);  // TODO use references instead of RCPs
  }                                                                   // end Multiply

  /**
    @brief Helper function to do matrix-matrix multiply

    Given CrsMatrix objects A, B and C, form the product C = A*B.
    In a parallel setting, A and B need not have matching distributions,
    but C needs to have the same row-map as A (if transposeA is false).
    At this time C=AT*B and C=A*BT are known to not work. However,
    C=A*B and C=AT*BT are known to work, Kurtis Nusbaum 03/24/2011

    @param A Input, must already have had 'FillComplete()' called.
    @param transposeA Input, whether to use transpose of matrix A.
    @param B Input, must already have had 'FillComplete()' called.
    @param transposeB Input, whether to use transpose of matrix B.
    @param C Result. If Teuchos::null, a new CrsMatrix is created with optimal number of nnz per row.
    @param call_FillComplete_on_result Optional argument, defaults to true.
    Power users may specify this argument to be false if they *DON'T*
    want this function to call C.FillComplete. (It is often useful
    to allow this function to call C.FillComplete, in cases where
    one or both of the input matrices are rectangular and it is not
    trivial to know which maps to use for the domain- and range-maps.)

*/
  static RCP<Matrix> Multiply(const Matrix& A, bool transposeA,
                              const Matrix& B, bool transposeB,
                              RCP<Matrix> C_in,
                              Teuchos::FancyOStream& fos,
                              bool doFillComplete              = true,
                              bool doOptimizeStorage           = true,
                              const std::string& label         = std::string(),
                              const RCP<ParameterList>& params = null) {
    TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), Exceptions::RuntimeError, "A is not fill-completed");
    TEUCHOS_TEST_FOR_EXCEPTION(!B.isFillComplete(), Exceptions::RuntimeError, "B is not fill-completed");

    // Default case: Xpetra Multiply
    RCP<Matrix> C = C_in;

    if (C == Teuchos::null) {
      double nnzPerRow = Teuchos::as<double>(0);

#if 0
        if (A.getDomainMap()->lib() == Xpetra::UseTpetra) {
          // For now, follow what ML and Epetra do.
          GO numRowsA = A.getGlobalNumRows();
          GO numRowsB = B.getGlobalNumRows();
          nnzPerRow = sqrt(Teuchos::as<double>(A.getGlobalNumEntries())/numRowsA) +
              sqrt(Teuchos::as<double>(B.getGlobalNumEntries())/numRowsB) - 1;
          nnzPerRow *=  nnzPerRow;
          double totalNnz = nnzPerRow * A.getGlobalNumRows() * 0.75 + 100;
          double minNnz = Teuchos::as<double>(1.2 * A.getGlobalNumEntries());
          if (totalNnz < minNnz)
            totalNnz = minNnz;
          nnzPerRow = totalNnz / A.getGlobalNumRows();

          fos << "Matrix product nnz per row estimate = " << Teuchos::as<LO>(nnzPerRow) << std::endl;
        }
#endif
      if (transposeA)
        C = MatrixFactory::Build(A.getDomainMap(), Teuchos::as<LO>(nnzPerRow));
      else
        C = MatrixFactory::Build(A.getRowMap(), Teuchos::as<LO>(nnzPerRow));

    } else {
      C->resumeFill();  // why this is not done inside of Tpetra MxM?
      fos << "Reuse C pattern" << std::endl;
    }

    Multiply(A, transposeA, B, transposeB, *C, doFillComplete, doOptimizeStorage, label, params);  // call Multiply routine from above

    return C;
  }

  /*! @brief Helper function to do matrix-matrix multiply

    Returns C = AB.

    @param A left matrix
    @param transposeA if true, use the transpose of A
    @param B right matrix
    @param transposeB if true, use the transpose of B
    @param callFillCompleteOnResult if true, the resulting matrix should be fillComplete'd
    */
  static RCP<Matrix> Multiply(const Matrix& A, bool transposeA,
                              const Matrix& B, bool transposeB,
                              Teuchos::FancyOStream& fos,
                              bool callFillCompleteOnResult    = true,
                              bool doOptimizeStorage           = true,
                              const std::string& label         = std::string(),
                              const RCP<ParameterList>& params = null) {
    return Multiply(A, transposeA, B, transposeB, Teuchos::null, fos, callFillCompleteOnResult, doOptimizeStorage, label, params);
  }

#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT)
  // Michael Gee's MLMultiply
  static RCP<Epetra_CrsMatrix> MLTwoMatrixMultiply(const Epetra_CrsMatrix& /* epA */,
                                                   const Epetra_CrsMatrix& /* epB */,
                                                   Teuchos::FancyOStream& /* fos */) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
                               "No ML multiplication available. This feature is currently not supported by Xpetra.");
    TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
  }
#endif  // ifdef HAVE_XPETRA_EPETRAEXT

  /*! @brief Helper function to do matrix-matrix multiply "in-place"

    Returns RCP to non-constant Xpetra::BlockedCrsMatrix.

    @param A left matrix
    @param transposeA if true, use the transpose of A
    @param B right matrix
    @param transposeB if true, use the transpose of B
    @param doOptimizeStorage if true, the resulting matrix should be fillComplete'd
    */
  static RCP<BlockedCrsMatrix> TwoMatrixMultiplyBlock(const BlockedCrsMatrix& A, bool transposeA,
                                                      const BlockedCrsMatrix& B, bool transposeB,
                                                      Teuchos::FancyOStream& fos,
                                                      bool doFillComplete    = true,
                                                      bool doOptimizeStorage = true) {
    TEUCHOS_TEST_FOR_EXCEPTION(transposeA || transposeB, Exceptions::RuntimeError,
                               "TwoMatrixMultiply for BlockedCrsMatrix not implemented for transposeA==true or transposeB==true");

    // Preconditions
    TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), Exceptions::RuntimeError, "A is not fill-completed");
    TEUCHOS_TEST_FOR_EXCEPTION(!B.isFillComplete(), Exceptions::RuntimeError, "B is not fill-completed");

    RCP<const MapExtractor> rgmapextractor = A.getRangeMapExtractor();
    RCP<const MapExtractor> domapextractor = B.getDomainMapExtractor();

    RCP<BlockedCrsMatrix> C = rcp(new BlockedCrsMatrix(rgmapextractor, domapextractor, 33 /* TODO fix me */));

    for (size_t i = 0; i < A.Rows(); ++i) {    // loop over all block rows of A
      for (size_t j = 0; j < B.Cols(); ++j) {  // loop over all block columns of B
        RCP<Matrix> Cij;

        for (size_t l = 0; l < B.Rows(); ++l) {  // loop for calculating entry C_{ij}
          RCP<Matrix> crmat1 = A.getMatrix(i, l);
          RCP<Matrix> crmat2 = B.getMatrix(l, j);

          if (crmat1.is_null() || crmat2.is_null())
            continue;

          // try unwrapping 1x1 blocked matrix
          {
            auto unwrap1 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(crmat1);
            auto unwrap2 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(crmat2);

            if (unwrap1.is_null() != unwrap2.is_null()) {
              if (unwrap1 != Teuchos::null && unwrap1->Rows() == 1 && unwrap1->Cols() == 1)
                crmat1 = unwrap1->getCrsMatrix();
              if (unwrap2 != Teuchos::null && unwrap2->Rows() == 1 && unwrap2->Cols() == 1)
                crmat2 = unwrap2->getCrsMatrix();
            }
          }

          RCP<CrsMatrixWrap> crop1 = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(crmat1);
          RCP<CrsMatrixWrap> crop2 = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(crmat2);
          TEUCHOS_TEST_FOR_EXCEPTION(crop1.is_null() != crop2.is_null(), Xpetra::Exceptions::RuntimeError,
                                     "A and B must be either both (compatible) BlockedCrsMatrix objects or both CrsMatrixWrap objects.");

          // Forcibly compute the global constants if we don't have them (only works for real CrsMatrices, not nested blocks)
          if (!crop1.is_null())
            Teuchos::rcp_const_cast<CrsGraph>(crmat1->getCrsGraph())->computeGlobalConstants();
          if (!crop2.is_null())
            Teuchos::rcp_const_cast<CrsGraph>(crmat2->getCrsGraph())->computeGlobalConstants();

          TEUCHOS_TEST_FOR_EXCEPTION(!crmat1->haveGlobalConstants(), Exceptions::RuntimeError,
                                     "crmat1 does not have global constants");
          TEUCHOS_TEST_FOR_EXCEPTION(!crmat2->haveGlobalConstants(), Exceptions::RuntimeError,
                                     "crmat2 does not have global constants");

          if (crmat1->getGlobalNumEntries() == 0 || crmat2->getGlobalNumEntries() == 0)
            continue;

          // temporary matrix containing result of local block multiplication
          RCP<Matrix> temp = Teuchos::null;

          if (crop1 != Teuchos::null && crop2 != Teuchos::null)
            temp = Multiply(*crop1, false, *crop2, false, fos);
          else {
            RCP<BlockedCrsMatrix> bop1 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(crmat1);
            RCP<BlockedCrsMatrix> bop2 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(crmat2);
            TEUCHOS_TEST_FOR_EXCEPTION(bop1.is_null() == true, Xpetra::Exceptions::BadCast, "A is not a BlockedCrsMatrix. (TwoMatrixMultiplyBlock)");
            TEUCHOS_TEST_FOR_EXCEPTION(bop2.is_null() == true, Xpetra::Exceptions::BadCast, "B is not a BlockedCrsMatrix. (TwoMatrixMultiplyBlock)");
            TEUCHOS_TEST_FOR_EXCEPTION(bop1->Cols() != bop2->Rows(), Xpetra::Exceptions::RuntimeError, "A has " << bop1->Cols() << " columns and B has " << bop2->Rows() << " rows. Matrices are not compatible! (TwoMatrixMultiplyBlock)");
            TEUCHOS_TEST_FOR_EXCEPTION(bop1->getDomainMap()->isSameAs(*(bop2->getRangeMap())) == false, Xpetra::Exceptions::RuntimeError, "Domain map of A is not the same as range map of B. Matrices are not compatible! (TwoMatrixMultiplyBlock)");

            // recursive multiplication call
            temp = TwoMatrixMultiplyBlock(*bop1, transposeA, *bop2, transposeB, fos, doFillComplete, doOptimizeStorage);

            RCP<BlockedCrsMatrix> btemp = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(temp);
            TEUCHOS_TEST_FOR_EXCEPTION(btemp->Rows() != bop1->Rows(), Xpetra::Exceptions::RuntimeError, "Number of block rows of local blocked operator is " << btemp->Rows() << " but should be " << bop1->Rows() << ". (TwoMatrixMultiplyBlock)");
            TEUCHOS_TEST_FOR_EXCEPTION(btemp->Cols() != bop2->Cols(), Xpetra::Exceptions::RuntimeError, "Number of block cols of local blocked operator is " << btemp->Cols() << " but should be " << bop2->Cols() << ". (TwoMatrixMultiplyBlock)");
            TEUCHOS_TEST_FOR_EXCEPTION(btemp->getRangeMapExtractor()->getFullMap()->isSameAs(*(bop1->getRangeMapExtractor()->getFullMap())) == false, Xpetra::Exceptions::RuntimeError, "Range map of local blocked operator should be same as first operator. (TwoMatrixMultiplyBlock)");
            TEUCHOS_TEST_FOR_EXCEPTION(btemp->getDomainMapExtractor()->getFullMap()->isSameAs(*(bop2->getDomainMapExtractor()->getFullMap())) == false, Xpetra::Exceptions::RuntimeError, "Domain map of local blocked operator should be same as second operator. (TwoMatrixMultiplyBlock)");
            TEUCHOS_TEST_FOR_EXCEPTION(btemp->getRangeMapExtractor()->getThyraMode() != bop1->getRangeMapExtractor()->getThyraMode(), Xpetra::Exceptions::RuntimeError, "Thyra mode of local range map extractor incompatible with range map extractor of A (TwoMatrixMultiplyBlock)");
            TEUCHOS_TEST_FOR_EXCEPTION(btemp->getDomainMapExtractor()->getThyraMode() != bop2->getDomainMapExtractor()->getThyraMode(), Xpetra::Exceptions::RuntimeError, "Thyra mode of local domain map extractor incompatible with domain map extractor of B (TwoMatrixMultiplyBlock)");
          }

          TEUCHOS_TEST_FOR_EXCEPTION(temp->isFillComplete() == false, Xpetra::Exceptions::RuntimeError, "Local block is not filled. (TwoMatrixMultiplyBlock)");

          RCP<Matrix> addRes = null;
          if (Cij.is_null())
            Cij = temp;
          else {
            MatrixMatrix::TwoMatrixAdd(*temp, false, 1.0, *Cij, false, 1.0, addRes, fos);
            Cij = addRes;
          }
        }

        if (!Cij.is_null()) {
          if (Cij->isFillComplete())
            Cij->resumeFill();
          Cij->fillComplete(B.getDomainMap(j), A.getRangeMap(i));
          C->setMatrix(i, j, Cij);
        } else {
          C->setMatrix(i, j, Teuchos::null);
        }
      }
    }

    if (doFillComplete)
      C->fillComplete();  // call default fillComplete for BlockCrsMatrixWrap objects

    return C;
  }  // TwoMatrixMultiplyBlock

  /*! @brief Helper function to calculate B = alpha*A + beta*B.

    @param A      left matrix operand
    @param transposeA indicate whether to use transpose of A
    @param alpha  scalar multiplier for A
    @param B      right matrix operand
    @param beta   scalar multiplier for B

    @return sum in B.

    Note that B does not have to be fill-completed.
    */
  static void TwoMatrixAdd(const Matrix& A, bool transposeA, SC alpha, Matrix& B, SC beta) {
    TEUCHOS_TEST_FOR_EXCEPTION(!(A.getRowMap()->isSameAs(*(B.getRowMap()))), Exceptions::Incompatible,
                               "TwoMatrixAdd: matrix row maps are not the same.");

    if (A.getRowMap()->lib() == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT)
      const Epetra_CrsMatrix& epA = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(A);
      Epetra_CrsMatrix& epB       = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstEpetraCrs(B);

      // FIXME is there a bug if beta=0?
      int rv = EpetraExt::MatrixMatrix::Add(epA, transposeA, alpha, epB, beta);

      if (rv != 0)
        throw Exceptions::RuntimeError("EpetraExt::MatrixMatrix::Add return value " + Teuchos::toString(rv));
      std::ostringstream buf;
#else
      throw Exceptions::RuntimeError("Xpetra must be compiled with EpetraExt.");
#endif
    } else if (A.getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
      throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra GO=long long enabled."));
#else
      const Tpetra::CrsMatrix<SC, LO, GO, NO>& tpA = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraCrs(A);
      Tpetra::CrsMatrix<SC, LO, GO, NO>& tpB = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstTpetraCrs(B);

      Tpetra::MatrixMatrix::Add(tpA, transposeA, alpha, tpB, beta);
#endif
#else
      throw Exceptions::RuntimeError("Xpetra must be compiled with Tpetra.");
#endif
    }
  }  // MatrixMatrix::TwoMatrixAdd()

  /*! @brief Helper function to calculate C = alpha*A + beta*B.

    @param A          left matrix operand
    @param transposeA indicate whether to use transpose of A
    @param alpha      scalar multiplier for A, defaults to 1.0
    @param B          right matrix operand
    @param transposeB indicate whether to use transpose of B
    @param beta       scalar multiplier for B, defaults to 1.0
    @param C          resulting sum
    @param fos        output stream for printing to screen
    @param AHasFixedNnzPerRow

    It is up to the caller to ensure that the resulting matrix sum is fillComplete'd.
    */
  static void TwoMatrixAdd(const Matrix& A, bool transposeA, const SC& alpha,
                           const Matrix& B, bool transposeB, const SC& beta,
                           RCP<Matrix>& C, Teuchos::FancyOStream& fos, bool AHasFixedNnzPerRow = false) {
    RCP<const Matrix> rcpA              = Teuchos::rcpFromRef(A);
    RCP<const Matrix> rcpB              = Teuchos::rcpFromRef(B);
    RCP<const BlockedCrsMatrix> rcpBopA = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrix>(rcpA);
    RCP<const BlockedCrsMatrix> rcpBopB = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrix>(rcpB);

    if (rcpBopA == Teuchos::null && rcpBopB == Teuchos::null) {
      TEUCHOS_TEST_FOR_EXCEPTION(!(A.getRowMap()->isSameAs(*(B.getRowMap()))), Exceptions::Incompatible,
                                 "TwoMatrixAdd: matrix row maps are not the same.");
      auto lib = A.getRowMap()->lib();
      if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT)
        const Epetra_CrsMatrix& epA = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(A);
        const Epetra_CrsMatrix& epB = Xpetra::Helpers<SC, LO, GO, NO>::Op2EpetraCrs(B);
        if (C.is_null()) {
          size_t maxNzInA     = 0;
          size_t maxNzInB     = 0;
          size_t numLocalRows = 0;
          if (A.isFillComplete() && B.isFillComplete()) {
            maxNzInA     = A.getLocalMaxNumRowEntries();
            maxNzInB     = B.getLocalMaxNumRowEntries();
            numLocalRows = A.getLocalNumRows();
          }

          if (maxNzInA == 1 || maxNzInB == 1 || AHasFixedNnzPerRow) {
            // first check if either A or B has at most 1 nonzero per row
            // the case of both having at most 1 nz per row is handled by the ``else''
            Teuchos::ArrayRCP<size_t> exactNnzPerRow(numLocalRows);

            if ((maxNzInA == 1 && maxNzInB > 1) || AHasFixedNnzPerRow) {
              for (size_t i = 0; i < numLocalRows; ++i)
                exactNnzPerRow[i] = B.getNumEntriesInLocalRow(Teuchos::as<LO>(i)) + maxNzInA;

            } else {
              for (size_t i = 0; i < numLocalRows; ++i)
                exactNnzPerRow[i] = A.getNumEntriesInLocalRow(Teuchos::as<LO>(i)) + maxNzInB;
            }

            fos << "MatrixMatrix::TwoMatrixAdd : special case detected (one matrix has a fixed nnz per row)"
                << ", using static profiling" << std::endl;
            C = rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(A.getRowMap(), exactNnzPerRow));

          } else {
            // general case
            LO maxPossibleNNZ = A.getLocalMaxNumRowEntries() + B.getLocalMaxNumRowEntries();
            C                 = rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(A.getRowMap(), maxPossibleNNZ));
          }
          if (transposeB)
            fos << "MatrixMatrix::TwoMatrixAdd : ** WARNING ** estimate could be badly wrong because second summand is transposed" << std::endl;
        }
        RCP<Epetra_CrsMatrix> epC = Xpetra::Helpers<SC, LO, GO, NO>::Op2NonConstEpetraCrs(C);
        Epetra_CrsMatrix* ref2epC = &*epC;  // to avoid a compiler error...

        // FIXME is there a bug if beta=0?
        int rv = EpetraExt::MatrixMatrix::Add(epA, transposeA, alpha, epB, transposeB, beta, ref2epC);

        if (rv != 0)
          throw Exceptions::RuntimeError("EpetraExt::MatrixMatrix::Add return value of " + Teuchos::toString(rv));
#else
        throw Exceptions::RuntimeError("MueLu must be compile with EpetraExt.");
#endif
      } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        throw(Xpetra::Exceptions::RuntimeError("Xpetra must be compiled with Tpetra GO=long long enabled."));
#else
        using helpers = Xpetra::Helpers<SC, LO, GO, NO>;
        using tcrs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NO>;
        const tcrs_matrix_type& tpA = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraCrs(A);
        const tcrs_matrix_type& tpB = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraCrs(B);
        C = helpers::tpetraAdd(tpA, transposeA, alpha, tpB, transposeB, beta);
#endif
#else
        throw Exceptions::RuntimeError("Xpetra must be compile with Tpetra.");
#endif
      }

      ///////////////////////// EXPERIMENTAL
      if (A.IsView("stridedMaps")) C->CreateView("stridedMaps", rcpFromRef(A));
      if (B.IsView("stridedMaps")) C->CreateView("stridedMaps", rcpFromRef(B));
      ///////////////////////// EXPERIMENTAL
    }
    // the first matrix is of type CrsMatrixWrap, the second is a blocked operator
    else if (rcpBopA == Teuchos::null && rcpBopB != Teuchos::null) {
      RCP<const MapExtractor> rgmapextractor = rcpBopB->getRangeMapExtractor();
      RCP<const MapExtractor> domapextractor = rcpBopB->getDomainMapExtractor();

      C                        = rcp(new BlockedCrsMatrix(rgmapextractor, domapextractor, 33 /* TODO fix me */));
      RCP<BlockedCrsMatrix> bC = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(C);

      size_t i = 0;
      for (size_t j = 0; j < rcpBopB->Cols(); ++j) {  // loop over all block columns of B
        RCP<Matrix> Cij = Teuchos::null;
        if (rcpA != Teuchos::null &&
            rcpBopB->getMatrix(i, j) != Teuchos::null) {
          // recursive call
          TwoMatrixAdd(*rcpA, transposeA, alpha,
                       *(rcpBopB->getMatrix(i, j)), transposeB, beta,
                       Cij, fos, AHasFixedNnzPerRow);
        } else if (rcpA == Teuchos::null &&
                   rcpBopB->getMatrix(i, j) != Teuchos::null) {
          Cij = rcpBopB->getMatrix(i, j);
        } else if (rcpA != Teuchos::null &&
                   rcpBopB->getMatrix(i, j) == Teuchos::null) {
          Cij = Teuchos::rcp_const_cast<Matrix>(rcpA);
        } else {
          Cij = Teuchos::null;
        }

        if (!Cij.is_null()) {
          if (Cij->isFillComplete())
            Cij->resumeFill();
          Cij->fillComplete();
          bC->setMatrix(i, j, Cij);
        } else {
          bC->setMatrix(i, j, Teuchos::null);
        }
      }  // loop over columns j
    }
    // the second matrix is of type CrsMatrixWrap, the first is a blocked operator
    else if (rcpBopA != Teuchos::null && rcpBopB == Teuchos::null) {
      RCP<const MapExtractor> rgmapextractor = rcpBopA->getRangeMapExtractor();
      RCP<const MapExtractor> domapextractor = rcpBopA->getDomainMapExtractor();

      C                        = rcp(new BlockedCrsMatrix(rgmapextractor, domapextractor, 33 /* TODO fix me */));
      RCP<BlockedCrsMatrix> bC = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(C);

      size_t j = 0;
      for (size_t i = 0; i < rcpBopA->Rows(); ++i) {  // loop over all block rows of A
        RCP<Matrix> Cij = Teuchos::null;
        if (rcpBopA->getMatrix(i, j) != Teuchos::null &&
            rcpB != Teuchos::null) {
          // recursive call
          TwoMatrixAdd(*(rcpBopA->getMatrix(i, j)), transposeA, alpha,
                       *rcpB, transposeB, beta,
                       Cij, fos, AHasFixedNnzPerRow);
        } else if (rcpBopA->getMatrix(i, j) == Teuchos::null &&
                   rcpB != Teuchos::null) {
          Cij = Teuchos::rcp_const_cast<Matrix>(rcpB);
        } else if (rcpBopA->getMatrix(i, j) != Teuchos::null &&
                   rcpB == Teuchos::null) {
          Cij = rcpBopA->getMatrix(i, j);
        } else {
          Cij = Teuchos::null;
        }

        if (!Cij.is_null()) {
          if (Cij->isFillComplete())
            Cij->resumeFill();
          Cij->fillComplete();
          bC->setMatrix(i, j, Cij);
        } else {
          bC->setMatrix(i, j, Teuchos::null);
        }
      }  // loop over rows i
    } else {
      // This is the version for blocked matrices

      // check the compatibility of the blocked operators
      TEUCHOS_TEST_FOR_EXCEPTION(rcpBopA.is_null() == true, Xpetra::Exceptions::BadCast, "A is not a BlockedCrsMatrix. (TwoMatrixAdd)");
      TEUCHOS_TEST_FOR_EXCEPTION(rcpBopB.is_null() == true, Xpetra::Exceptions::BadCast, "B is not a BlockedCrsMatrix. (TwoMatrixAdd)");
      TEUCHOS_TEST_FOR_EXCEPTION(rcpBopA->Rows() != rcpBopB->Rows(), Xpetra::Exceptions::RuntimeError, "A has " << rcpBopA->Rows() << " rows and B has " << rcpBopA->Rows() << " rows. Matrices are not compatible! (TwoMatrixAdd)");
      TEUCHOS_TEST_FOR_EXCEPTION(rcpBopA->Rows() != rcpBopB->Rows(), Xpetra::Exceptions::RuntimeError, "A has " << rcpBopA->Cols() << " cols and B has " << rcpBopA->Cols() << " cols. Matrices are not compatible! (TwoMatrixAdd)");
      TEUCHOS_TEST_FOR_EXCEPTION(rcpBopA->getRangeMap()->isSameAs(*(rcpBopB->getRangeMap())) == false, Xpetra::Exceptions::RuntimeError, "Range map of A is not the same as range map of B. Matrices are not compatible! (TwoMatrixAdd)");
      TEUCHOS_TEST_FOR_EXCEPTION(rcpBopA->getDomainMap()->isSameAs(*(rcpBopB->getDomainMap())) == false, Xpetra::Exceptions::RuntimeError, "Domain map of A is not the same as domain map of B. Matrices are not compatible! (TwoMatrixAdd)");
      TEUCHOS_TEST_FOR_EXCEPTION(rcpBopA->getRangeMapExtractor()->getThyraMode() != rcpBopB->getRangeMapExtractor()->getThyraMode(), Xpetra::Exceptions::RuntimeError, "Different Thyra/Xpetra style gids in RangeMapExtractor of A and B. Matrices are not compatible! (TwoMatrixAdd)");
      TEUCHOS_TEST_FOR_EXCEPTION(rcpBopA->getDomainMapExtractor()->getThyraMode() != rcpBopB->getDomainMapExtractor()->getThyraMode(), Xpetra::Exceptions::RuntimeError, "Different Thyra/Xpetra style gids in DomainMapExtractor of A and B. Matrices are not compatible! (TwoMatrixAdd)");

      RCP<const MapExtractor> rgmapextractor = rcpBopA->getRangeMapExtractor();
      RCP<const MapExtractor> domapextractor = rcpBopB->getDomainMapExtractor();

      C                        = rcp(new BlockedCrsMatrix(rgmapextractor, domapextractor, 33 /* TODO fix me */));
      RCP<BlockedCrsMatrix> bC = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(C);

      for (size_t i = 0; i < rcpBopA->Rows(); ++i) {    // loop over all block rows of A
        for (size_t j = 0; j < rcpBopB->Cols(); ++j) {  // loop over all block columns of B

          RCP<Matrix> Cij = Teuchos::null;
          if (rcpBopA->getMatrix(i, j) != Teuchos::null &&
              rcpBopB->getMatrix(i, j) != Teuchos::null) {
            // recursive call

            TwoMatrixAdd(*(rcpBopA->getMatrix(i, j)), transposeA, alpha,
                         *(rcpBopB->getMatrix(i, j)), transposeB, beta,
                         Cij, fos, AHasFixedNnzPerRow);
          } else if (rcpBopA->getMatrix(i, j) == Teuchos::null &&
                     rcpBopB->getMatrix(i, j) != Teuchos::null) {
            Cij = rcpBopB->getMatrix(i, j);
          } else if (rcpBopA->getMatrix(i, j) != Teuchos::null &&
                     rcpBopB->getMatrix(i, j) == Teuchos::null) {
            Cij = rcpBopA->getMatrix(i, j);
          } else {
            Cij = Teuchos::null;
          }

          if (!Cij.is_null()) {
            if (Cij->isFillComplete())
              Cij->resumeFill();
            // Cij->fillComplete(rcpBopA->getDomainMap(j), rcpBopA->getRangeMap(i));
            Cij->fillComplete();
            bC->setMatrix(i, j, Cij);
          } else {
            bC->setMatrix(i, j, Teuchos::null);
          }
        }  // loop over columns j
      }    // loop over rows i

    }  // end blocked recursive algorithm
  }    // MatrixMatrix::TwoMatrixAdd()
};     // end specialization on GO=long long and NO=EpetraNode

#endif  // HAVE_XPETRA_EPETRA

}  // end namespace Xpetra

#define XPETRA_MATRIXMATRIX_SHORT

#endif /* PACKAGES_XPETRA_SUP_UTILS_XPETRA_MATRIXMATRIX_DECL_HPP_ */
