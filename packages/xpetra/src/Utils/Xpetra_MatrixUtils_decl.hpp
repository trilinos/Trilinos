// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PACKAGES_XPETRA_SUP_MATRIX_UTILS_DECL_HPP_
#define PACKAGES_XPETRA_SUP_MATRIX_UTILS_DECL_HPP_

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_Map.hpp"
#include "Xpetra_MapUtils.hpp"
#include "Xpetra_StridedMap.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_MapExtractor.hpp"
#include "Xpetra_MapExtractorFactory.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_MatrixFactory.hpp"
#include "Xpetra_BlockedCrsMatrix.hpp"
#include "Xpetra_MatrixMatrix.hpp"
#include "Xpetra_Helpers.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraMultiVector.hpp"
#include <Tpetra_RowMatrixTransposer.hpp>
#include <Tpetra_Details_extractBlockDiagonal.hpp>
#include <Tpetra_Details_scaleBlockDiagonal.hpp>
#endif

namespace Xpetra {

/*!
  @class MatrixUtils
  @brief Xpetra utility class for common matrix-related routines

  The routines should be independent from Epetra/Tpetra and be purely implemented in Xpetra.
  Other matrix-related routines are out-sourced into other helper classes (e.g. MatrixMatrix for
  MM multiplication and addition).

*/
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class MatrixUtils {
#undef XPETRA_MATRIXUTILS_SHORT
#include "Xpetra_UseShortNames.hpp"

 public:
  static Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpetraGidNumbering2ThyraGidNumbering(
      const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& input);

  static Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> findColumnSubMap(
      const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& input,
      const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& domainMap);

  /** Given a matrix A split it into a nxm blocked matrix using the map extractors.

    @param input Input matrix, must already have had 'FillComplete()' called.
    @param rangeMapExtractor MapExtractor object describing the splitting of rows of the output block matrix
    @param domainMapExtractor MapExtractor object describing the splitting of columns of the output block matrix
    @param columnMapExtractor (not fully clear whether we need that. is always Teuchos::null)
    @param bThyraMode If true, build a n x n blocked operator using Thyra GIDs

    @return Fill-completed block version of intput matrix
  */
  static Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> SplitMatrix(
      const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& input,
      Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>> rangeMapExtractor,
      Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>> domainMapExtractor,
      Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>> columnMapExtractor = Teuchos::null,
      bool bThyraMode                                                                                        = false);

  /** Given a matrix A, detect too small diagonals and replace any found with ones. */

  static void CheckRepairMainDiagonal(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& Ac,
                                      bool const& repairZeroDiagonals, Teuchos::FancyOStream& fos,
                                      const typename Teuchos::ScalarTraits<Scalar>::magnitudeType threshold = Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero(),
                                      const Scalar replacementValue                                         = Teuchos::ScalarTraits<Scalar>::one());

  /** Given a matrix A, boost the diagonal to a relative floor.  Multiple PDEs will be scaled differently.
      Each PDE can be given its own relative threshold, or a single threshold can be used for all PDEs
      NOTE: This is not Kokkos-ized
      Precondition: A->GetFixedBlockSize() == relativeThreshold.size()  OR relativeThreshold.size() == 1
  **/

  static void RelativeDiagonalBoost(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
                                    const Teuchos::ArrayView<const double>& relativeThreshold, Teuchos::FancyOStream& fos);

  // Extracting the block diagonal of a matrix
  static void extractBlockDiagonal(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                                   Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& diagonal);

  // Inverse scaling by a block-diagonal matrix
  static void inverseScaleBlockDiagonal(Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& blockDiagonal,
                                        bool doTranspose,
                                        Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& toBeScaled);

  static void checkLocalRowMapMatchesColMap(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A);

  /*!
    \@brief Convert matrix to strided row and column maps

    @param matrix Matrix to be converted
    @param rangeStridingInfo Striding information for row/range map
    @param domainStridingInfo Striding information for column/domain map
  */
  static void convertMatrixToStridedMaps(
      Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> matrix,
      std::vector<size_t>& rangeStridingInfo, std::vector<size_t>& domainStridingInfo);
};

}  // end namespace Xpetra

#define XPETRA_MATRIXUTILS_SHORT

#endif  // PACKAGES_XPETRA_SUP_MATRIX_UTILS_DECL_HPP_
