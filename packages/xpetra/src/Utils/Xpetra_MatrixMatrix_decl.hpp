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

#include <TpetraExt_MatrixMatrix.hpp>
#include <Tpetra_RowMatrixTransposer.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Xpetra_TpetraBlockCrsMatrix.hpp>
#include <Tpetra_BlockCrsMatrix_Helpers.hpp>
#include <Xpetra_TpetraMultiVector.hpp>
#include <Xpetra_TpetraVector.hpp>

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

}  // end namespace Xpetra

#define XPETRA_MATRIXMATRIX_SHORT

#endif /* PACKAGES_XPETRA_SUP_UTILS_XPETRA_MATRIXMATRIX_DECL_HPP_ */
