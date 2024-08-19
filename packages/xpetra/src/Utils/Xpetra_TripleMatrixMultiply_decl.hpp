// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PACKAGES_XPETRA_SUP_UTILS_XPETRA_TRIPLEMATRIXMULTIPLY_DECL_HPP_
#define PACKAGES_XPETRA_SUP_UTILS_XPETRA_TRIPLEMATRIXMULTIPLY_DECL_HPP_

#include "Xpetra_ConfigDefs.hpp"

// #include "Xpetra_BlockedCrsMatrix.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_MapExtractor.hpp"
#include "Xpetra_Map.hpp"
#include "Xpetra_MatrixFactory.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_StridedMapFactory.hpp"
#include "Xpetra_StridedMap.hpp"
#include "Xpetra_IO.hpp"

#include <TpetraExt_TripleMatrixMultiply.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Tpetra_BlockCrsMatrix.hpp>
#include <Tpetra_BlockCrsMatrix_Helpers.hpp>
// #include <Xpetra_TpetraMultiVector.hpp>
// #include <Xpetra_TpetraVector.hpp>

namespace Xpetra {

template <class Scalar,
          class LocalOrdinal /*= int*/,
          class GlobalOrdinal /*= LocalOrdinal*/,
          class Node /*= Tpetra::KokkosClassic::DefaultNode::DefaultNodeType*/>
class TripleMatrixMultiply {
#undef XPETRA_TRIPLEMATRIXMULTIPLY_SHORT
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
  static void MultiplyRAP(const Matrix& R, bool transposeR,
                          const Matrix& A, bool transposeA,
                          const Matrix& P, bool transposeP,
                          Matrix& Ac,
                          bool call_FillComplete_on_result = true,
                          bool doOptimizeStorage           = true,
                          const std::string& label         = std::string(),
                          const RCP<ParameterList>& params = null);

};  // class TripleMatrixMultiply

}  // end namespace Xpetra

#define XPETRA_TRIPLEMATRIXMULTIPLY_SHORT

#endif /* PACKAGES_XPETRA_SUP_UTILS_XPETRA_TRIPLEMATRIXMULTIPLY_DECL_HPP_ */
