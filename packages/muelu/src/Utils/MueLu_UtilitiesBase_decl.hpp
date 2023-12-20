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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_UTILITIESBASE_DECL_HPP
#define MUELU_UTILITIESBASE_DECL_HPP

#include <string>

#include "MueLu_ConfigDefs.hpp"

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_ParameterList.hpp>

#include "Kokkos_ArithTraits.hpp"

#include <Xpetra_BlockedCrsMatrix_fwd.hpp>
#include <Xpetra_BlockedMap_fwd.hpp>
#include <Xpetra_BlockedVector_fwd.hpp>
#include <Xpetra_CrsGraphFactory_fwd.hpp>
#include <Xpetra_CrsGraph_fwd.hpp>
#include <Xpetra_CrsMatrix_fwd.hpp>
#include <Xpetra_CrsMatrixWrap_fwd.hpp>
#include <Xpetra_Import_fwd.hpp>
#include <Xpetra_ImportFactory_fwd.hpp>
#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>
#include <Xpetra_Operator_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>

namespace MueLu {

// MPI helpers
#define MueLu_sumAll(rcpComm, in, out) \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out))
#define MueLu_minAll(rcpComm, in, out) \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MIN, in, Teuchos::outArg(out))
#define MueLu_maxAll(rcpComm, in, out) \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MAX, in, Teuchos::outArg(out))

/*!
  @class UtilitiesBase
  @brief MueLu utility class.

  This class provides a number of static helper methods. Some are temporary and will eventually
  go away, while others should be moved to Xpetra.
*/
template <class Scalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class UtilitiesBase {
 public:
#undef MUELU_UTILITIESBASE_SHORT
#include "MueLu_UseShortNames.hpp"
 public:
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;

  static RCP<Matrix> Crs2Op(RCP<CrsMatrix> Op);

  /*! @brief Threshold a matrix

    Returns matrix filtered with a threshold value.

    NOTE -- it's assumed that A has been fillComplete'd.
  */
  static RCP<CrsMatrixWrap> GetThresholdedMatrix(const RCP<Matrix>& Ain, const Scalar threshold, const bool keepDiagonal = true, const GlobalOrdinal expectedNNZperRow = -1);

  /*! @brief Threshold a graph

    Returns graph filtered with a threshold value.

    NOTE -- it's assumed that A has been fillComplete'd.
  */
  static RCP<Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>> GetThresholdedGraph(const RCP<Matrix>& A, const Magnitude threshold, const GlobalOrdinal expectedNNZperRow = -1);

  /*! @brief Extract Matrix Diagonal

    Returns Matrix diagonal in ArrayRCP.

    NOTE -- it's assumed that A has been fillComplete'd.
  */
  static Teuchos::ArrayRCP<Scalar> GetMatrixDiagonal_arcp(const Matrix& A);

  /*! @brief Extract Matrix Diagonal

 Returns Matrix diagonal in RCP<Vector>.

 NOTE -- it's assumed that A has been fillComplete'd.
 */
  static RCP<Vector> GetMatrixDiagonal(const Matrix& A);

  /*! @brief Extract Matrix Diagonal

    Returns inverse of the Matrix diagonal in ArrayRCP.

    NOTE -- it's assumed that A has been fillComplete'd.
  */
  // static RCP<Vector> GetMatrixDiagonalInverse(const Matrix& A, Magnitude tol = Teuchos::ScalarTraits<Scalar>::eps()*100, Scalar valReplacement = Teuchos::ScalarTraits<Scalar>::zero());
  static RCP<Vector> GetMatrixDiagonalInverse(const Matrix& A, Magnitude tol = Teuchos::ScalarTraits<Scalar>::eps() * 100, Scalar valReplacement = Teuchos::ScalarTraits<Scalar>::zero(), const bool doLumped = false);

  /*! @brief Extract Matrix Diagonal of lumped matrix

    Returns Matrix diagonal of lumped matrix in RCP<Vector>.

    NOTE -- it's assumed that A has been fillComplete'd.
  */
  static Teuchos::RCP<Vector> GetLumpedMatrixDiagonal(Matrix const& A, const bool doReciprocal = false,
                                                      Magnitude tol                            = Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero()),
                                                      Scalar valReplacement                    = Teuchos::ScalarTraits<Scalar>::zero(),
                                                      const bool replaceSingleEntryRowWithZero = false,
                                                      const bool useAverageAbsDiagVal          = false);

  /*! @brief Return vector containing: max_{i\not=k}(-a_ik), for each for i in the matrix
   *
   * @param[in] A: input matrix
   * @ret: vector containing max_{i\not=k}(-a_ik)
   */

  static Teuchos::ArrayRCP<Magnitude> GetMatrixMaxMinusOffDiagonal(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A);

  static Teuchos::ArrayRCP<Magnitude> GetMatrixMaxMinusOffDiagonal(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A, const Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>& BlockNumber);

  /*! @brief Return vector containing inverse of input vector
   *
   * @param[in] v: input vector
   * @param[in] tol: tolerance. If entries of input vector are smaller than tolerance they are replaced by valReplacement (see below). The default value for tol is 100*eps (machine precision)
   * @param[in] valReplacement: Value put in for undefined entries in output vector (default: 0.0)
   * @ret: vector containing inverse values of input vector v
   */
  static Teuchos::RCP<Vector> GetInverse(Teuchos::RCP<const Vector> v, Magnitude tol = Teuchos::ScalarTraits<Scalar>::eps() * 100, Scalar valReplacement = Teuchos::ScalarTraits<Scalar>::zero());

  /*! @brief Extract Overlapped Matrix Diagonal

    Returns overlapped Matrix diagonal in ArrayRCP.

    The local overlapped diagonal has an entry for each index in A's column map.
    NOTE -- it's assumed that A has been fillComplete'd.
  */
  static RCP<Vector> GetMatrixOverlappedDiagonal(const Matrix& A);

  /*! @brief Extract Overlapped Matrix Deleted Rowsum

    Returns overlapped Matrix deleted Rowsum in ArrayRCP.

    The local overlapped deleted Rowsum has an entry for each index in A's column map.
    NOTE -- it's assumed that A has been fillComplete'd.
  */

  static RCP<Vector> GetMatrixOverlappedDeletedRowsum(const Matrix& A);

  static RCP<Xpetra::Vector<Magnitude, LocalOrdinal, GlobalOrdinal, Node>>
  GetMatrixOverlappedAbsDeletedRowsum(const Matrix& A);

  // TODO: should NOT return an Array. Definition must be changed to:
  // - ArrayRCP<> ResidualNorm(Matrix const &Op, MultiVector const &X, MultiVector const &RHS)
  // or
  // - void ResidualNorm(Matrix const &Op, MultiVector const &X, MultiVector const &RHS, Array &)
  static Teuchos::Array<Magnitude> ResidualNorm(const Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op, const MultiVector& X, const MultiVector& RHS);

  static Teuchos::Array<Magnitude> ResidualNorm(const Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op, const MultiVector& X, const MultiVector& RHS, MultiVector& Resid);

  static RCP<MultiVector> Residual(const Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op, const MultiVector& X, const MultiVector& RHS);

  static void Residual(const Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op, const MultiVector& X, const MultiVector& RHS, MultiVector& Resid);

  /*! @brief Power method.

    @param A matrix
    @param scaleByDiag if true, estimate the largest eigenvalue of \f$ D^; A \f$.
    @param niters maximum number of iterations
    @param tolerance stopping tolerance
    @verbose if true, print iteration information
    @seed  seed for randomizing initial guess

    (Shamelessly grabbed from tpetra/examples.)
  */
  static Scalar PowerMethod(const Matrix& A, bool scaleByDiag = true,
                            LocalOrdinal niters = 10, Magnitude tolerance = 1e-2, bool verbose = false, unsigned int seed = 123);

  /*! @brief Power method.

    @param A matrix
    @param diagInvVec reciprocal of matrix diagonal
    @param niters maximum number of iterations
    @param tolerance stopping tolerance
    @verbose if true, print iteration information
    @seed  seed for randomizing initial guess

    (Shamelessly grabbed from tpetra/examples.)
  */
  static Scalar PowerMethod(const Matrix& A, const RCP<Vector>& diagInvVec,
                            LocalOrdinal niters = 10, Magnitude tolerance = 1e-2, bool verbose = false, unsigned int seed = 123);

  static RCP<Teuchos::FancyOStream> MakeFancy(std::ostream& os);

  /*! @brief Squared distance between two rows in a multivector

    Used for coordinate vectors.
  */
  static typename Teuchos::ScalarTraits<Scalar>::magnitudeType Distance2(const Teuchos::Array<Teuchos::ArrayRCP<const Scalar>>& v, LocalOrdinal i0, LocalOrdinal i1);

  /*! @brief Weighted squared distance between two rows in a multivector

    Used for coordinate vectors.
  */
  static typename Teuchos::ScalarTraits<Scalar>::magnitudeType Distance2(const Teuchos::ArrayView<double>& weight, const Teuchos::Array<Teuchos::ArrayRCP<const Scalar>>& v, LocalOrdinal i0, LocalOrdinal i1);

  /*! @brief Detect Dirichlet rows

    The routine assumes, that if there is only one nonzero per row, it is on the diagonal and therefore a DBC.
    This is safe for most of our applications, but one should be aware of that.

    There is an alternative routine (see DetectDirichletRowsExt)

    @param[in] A matrix
    @param[in] tol If a row entry's magnitude is less than or equal to this tolerance, the entry is treated as zero.

    @return boolean array.  The ith entry is true iff row i is a Dirichlet row.
  */
  static Teuchos::ArrayRCP<const bool> DetectDirichletRows(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A, const Magnitude& tol = Teuchos::ScalarTraits<Magnitude>::zero(), bool count_twos_as_dirichlet = false);

  /*! @brief Detect Dirichlet rows

    @param[in] A matrix
    @param[in] tol If a row entry's magnitude is less than or equal to this tolerance, the entry is treated as zero.

    @return boolean array.  The ith entry is true iff row i is a Dirichlet row.
  */
  static Kokkos::View<bool*, typename NO::device_type> DetectDirichletRows_kokkos(const Matrix& A, const Magnitude& tol = Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<SC>::magnitudeType>::zero(), const bool count_twos_as_dirichlet = false);

  /*! @brief Detect Dirichlet rows (extended version)

    Look at each matrix row and mark it as Dirichlet if there is only one
    "not small" nonzero on the diagonal. In determining whether a nonzero
    is "not small" use
    \f abs(A(i,j)) / sqrt(abs(diag[i]*diag[j])) > tol

    @param[in] A matrix
    @param[in/out] bHasZeroDiagonal Reference to boolean variable. Returns true if there is a zero on the diagonal in the local part of the Matrix. Otherwise it is false. Different processors might return a different value. There is no global reduction!
    @param[in] tol If a row entry's magnitude is less than or equal to this tolerance, the entry is treated as zero.
    @return boolean array.  The ith entry is true iff row i is a Dirichlet row.
  */
  static Teuchos::ArrayRCP<const bool> DetectDirichletRowsExt(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A, bool& bHasZeroDiagonal, const Magnitude& tol = Teuchos::ScalarTraits<Scalar>::zero());

  /*! @brief Find non-zero values in an ArrayRCP
    Compares the value to 2 * machine epsilon

    @param[in]  vals - ArrayRCP<const Scalar> of values to be tested
    @param[out] nonzeros - ArrayRCP<bool> of true/false values for whether each entry in vals is nonzero
  */

  static void FindNonZeros(const Teuchos::ArrayRCP<const Scalar> vals,
                           Teuchos::ArrayRCP<bool> nonzeros);

  /*! @brief Find non-zero values in an ArrayRCP
    Compares the value to 2 * machine epsilon

    @param[in]  vals - ArrayRCP<const Scalar> of values to be tested
    @param[out] nonzeros - ArrayRCP<bool> of true/false values for whether each entry in vals is nonzero
  */
  static void FindNonZeros(const typename Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dual_view_type::t_dev_const_um vals,
                           Kokkos::View<bool*, typename Node::device_type> nonzeros);

  /*! @brief Detects Dirichlet columns & domains from a list of Dirichlet rows

    @param[in] A - Matrix on which to apply Dirichlet column detection
    @param[in] dirichletRows - ArrayRCP<bool> of indicators as to which rows are Dirichlet
    @param[out] dirichletCols - ArrayRCP<bool> of indicators as to which cols are Dirichlet
    @param[out] dirichletDomain - ArrayRCP<bool> of indicators as to which domains are Dirichlet
  */

  static void DetectDirichletColsAndDomains(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                                            const Teuchos::ArrayRCP<bool>& dirichletRows,
                                            Teuchos::ArrayRCP<bool> dirichletCols,
                                            Teuchos::ArrayRCP<bool> dirichletDomain);

  /*! @brief Detects Dirichlet columns & domains from a list of Dirichlet rows

    @param[in] A - Matrix on which to apply Dirichlet column detection
    @param[in] dirichletRows - View<bool> of indicators as to which rows are Dirichlet
    @param[out] dirichletCols - View<bool> of indicators as to which cols are Dirichlet
    @param[out] dirichletDomain - View<bool> of indicators as to which domains are Dirichlet
  */
  static void DetectDirichletColsAndDomains(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                                            const Kokkos::View<bool*, typename Node::device_type>& dirichletRows,
                                            Kokkos::View<bool*, typename Node::device_type> dirichletCols,
                                            Kokkos::View<bool*, typename Node::device_type> dirichletDomain);

  /*! @brief Apply Rowsum Criterion

    Flags a row i as dirichlet if:

    \sum_{j\not=i} A_ij > A_ii * tol

    @param[in] A matrix
    @param[in] rowSumTol See above
    @param[in/out] dirichletRows boolean array.  The ith entry is true if the above criterion is satisfied (or if it was already set to true)

  */
  static void ApplyRowSumCriterion(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A, const Magnitude rowSumTol, Teuchos::ArrayRCP<bool>& dirichletRows);

  static void ApplyRowSumCriterion(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A, const Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>& BlockNumber, const Magnitude rowSumTol, Teuchos::ArrayRCP<bool>& dirichletRows);

  static void ApplyRowSumCriterion(const Matrix& A,
                                   const typename Teuchos::ScalarTraits<Scalar>::magnitudeType rowSumTol,
                                   Kokkos::View<bool*, typename NO::device_type>& dirichletRows);

  /*! @brief Detect Dirichlet columns based on Dirichlet rows

    The routine finds all column indices that are in Dirichlet rows, where Dirichlet rows are described by dirichletRows,
    as returned by DetectDirichletRows.

    @param[in] A matrix
    @param[in] dirichletRows array of Dirichlet rows.

    @return boolean array.  The ith entry is true iff row i is a Dirichlet column.
  */
  static Teuchos::ArrayRCP<const bool> DetectDirichletCols(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                                                           const Teuchos::ArrayRCP<const bool>& dirichletRows);

  /*! @brief Detect Dirichlet columns based on Dirichlet rows

    The routine finds all column indices that are in Dirichlet rows, where Dirichlet rows are described by dirichletRows,
    as returned by DetectDirichletRows.

    @param[in] A matrix
    @param[in] dirichletRows array of Dirichlet rows.

    @return boolean array.  The ith entry is true iff row i is a Dirichlet column.
  */
  static Kokkos::View<bool*, typename NO::device_type> DetectDirichletCols(const Matrix& A, const Kokkos::View<const bool*, typename NO::device_type>& dirichletRows);

  /*! @brief Frobenius inner product of two matrices

    Used in energy minimization algorithms
  */
  static Scalar Frobenius(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A, const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B);

  /*! @brief Set seed for random number generator.

    Distribute the seeds evenly in [1,INT_MAX-1].  This guarantees nothing
    about where random number streams on difference processes will intersect.
    This does avoid overflow situations in parallel when multiplying by a PID.
    It also avoids the pathological case of having the *same* random number stream
    on each process.
  */

  static void SetRandomSeed(const Teuchos::Comm<int>& comm);

  // Finds the OAZ Dirichlet rows for this matrix
  // so far only used in IntrepidPCoarsenFactory
  // TODO check whether we can use DetectDirichletRows instead
  static void FindDirichletRows(Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
                                std::vector<LocalOrdinal>& dirichletRows, bool count_twos_as_dirichlet = false);

  // Applies Ones-and-Zeros to matrix rows
  // Takes a vector of row indices
  static void ApplyOAZToMatrixRows(Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
                                   const std::vector<LocalOrdinal>& dirichletRows);

  // Applies Ones-and-Zeros to matrix rows
  // Takes a Boolean array.
  static void ApplyOAZToMatrixRows(Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
                                   const Teuchos::ArrayRCP<const bool>& dirichletRows);

  static void ApplyOAZToMatrixRows(RCP<Matrix>& A, const Kokkos::View<const bool*, typename Node::device_type>& dirichletRows);

  // Zeros out rows
  // Takes a vector containg Dirichlet row indices
  static void ZeroDirichletRows(Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
                                const std::vector<LocalOrdinal>& dirichletRows,
                                Scalar replaceWith = Teuchos::ScalarTraits<Scalar>::zero());

  // Zeros out rows
  // Takes a Boolean ArrayRCP
  static void ZeroDirichletRows(Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
                                const Teuchos::ArrayRCP<const bool>& dirichletRows,
                                Scalar replaceWith = Teuchos::ScalarTraits<Scalar>::zero());

  // Zeros out rows
  // Takes a Boolean ArrayRCP
  static void ZeroDirichletRows(Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& X,
                                const Teuchos::ArrayRCP<const bool>& dirichletRows,
                                Scalar replaceWith = Teuchos::ScalarTraits<Scalar>::zero());

  static void ZeroDirichletRows(RCP<Matrix>& A, const Kokkos::View<const bool*, typename NO::device_type>& dirichletRows, SC replaceWith = Teuchos::ScalarTraits<SC>::zero());

  static void ZeroDirichletRows(RCP<MultiVector>& X, const Kokkos::View<const bool*, typename NO::device_type>& dirichletRows, SC replaceWith = Teuchos::ScalarTraits<SC>::zero());

  // Zeros out columns
  // Takes a Boolean vector
  static void ZeroDirichletCols(Teuchos::RCP<Matrix>& A,
                                const Teuchos::ArrayRCP<const bool>& dirichletCols,
                                Scalar replaceWith = Teuchos::ScalarTraits<Scalar>::zero());

  static void ZeroDirichletCols(RCP<Matrix>& A, const Kokkos::View<const bool*, typename NO::device_type>& dirichletCols, SC replaceWith = Teuchos::ScalarTraits<SC>::zero());

  // Finds the OAZ Dirichlet rows for this matrix
  static void FindDirichletRowsAndPropagateToCols(Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
                                                  Teuchos::RCP<Xpetra::Vector<int, LocalOrdinal, GlobalOrdinal, Node>>& isDirichletRow,
                                                  Teuchos::RCP<Xpetra::Vector<int, LocalOrdinal, GlobalOrdinal, Node>>& isDirichletCol);

  //! Creates a copy of a matrix where the non-zero entries are replaced by ones.
  // You can use this to de-normalize a tenative prolongator, for instance
  static RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> ReplaceNonZerosWithOnes(const RCP<Matrix>& original);

  // This routine takes a BlockedMap and an Importer (assuming that the BlockedMap matches the source of the importer) and generates a BlockedMap corresponding
  // to the Importer's target map.  We assume that the targetMap is unique (which, is not a strict requirement of an Importer, but is here and no, we don't check)
  // This is largely intended to be used in repartitioning of blocked matrices
  static RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>> GeneratedBlockedTargetMap(const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>& sourceBlockedMap,
                                                                                                    const Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>& Importer);

  // Checks to see if the first chunk of the colMap is also the row map.  This simiplifies a bunch of
  // operation in coarsening
  static bool MapsAreNested(const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& rowMap, const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& colMap);

  /*! Perform a Cuthill-McKee (CM) or Reverse Cuthill-McKee (RCM) ordering of the local component of the matrix.
    Kokkos-Kernels has an RCM implementation, so we reverse that here if we call CM.
  */
  static RCP<Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>> CuthillMcKee(const Matrix& Op);

  /*! Perform a Reverse Cuthill-McKee (RCM) ordering of the local component of the matrix.
   */
  static RCP<Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>> ReverseCuthillMcKee(const Matrix& Op);

};  // class UtilitiesBase

// Useful Kokkos conversions
template <class View, unsigned AppendValue>
struct AppendTrait {
  // static_assert(false, "Error: NOT a Kokkos::View");
};

template <class MT, unsigned T>
struct CombineMemoryTraits {
  // empty
};

template <unsigned U, unsigned T>
struct CombineMemoryTraits<Kokkos::MemoryTraits<U>, T> {
  typedef Kokkos::MemoryTraits<U | T> type;
};

template <class DataType, unsigned T, class... Pack>
struct AppendTrait<Kokkos::View<DataType, Pack...>, T> {
  typedef Kokkos::View<DataType, Pack...> view_type;
  using type = Kokkos::View<DataType, typename view_type::array_layout, typename view_type::device_type, typename CombineMemoryTraits<typename view_type::memory_traits, T>::type>;
};

///////////////////////////////////////////

}  // namespace MueLu

#define MUELU_UTILITIESBASE_SHORT
#endif  // MUELU_UTILITIESBASE_DECL_HPP
