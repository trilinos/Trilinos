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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
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

#include "Xpetra_IO.hpp"

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
