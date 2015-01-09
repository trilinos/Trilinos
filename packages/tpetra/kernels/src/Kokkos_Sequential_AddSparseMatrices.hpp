//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef KOKKOS_SEQUENTIAL_ADDSPARSEMATRICES_HPP
#define KOKKOS_SEQUENTIAL_ADDSPARSEMATRICES_HPP

#include <Kokkos_ConfigDefs.hpp>
#include <functional> // std::plus

namespace Kokkos {
namespace Sequential {

/// \fn countMergedRowEntries
/// \brief Compute number of entries in row i of the sum of the two
///   sparse matrices.
///
/// Each sequence of indices must be sorted.  It need not be unique,
/// but if not unique, then the resulting output sequence will not be
/// unique.
///
/// The sum does not consider numerical cancellation, that is, the sum
/// of the values making the value of an entry exactly zero.
///
/// \param ptr1 [in] Row offsets of sparse matrix A.
/// \param ind1 [in] Column indices of sparse matrix A.
/// \param ptr2 [in] Row offsets of sparse matrix B.
/// \param ind2 [in] Column indices of sparse matrix B.
/// \param i [in] Row index.
///
/// \return Number of entries in row i of the sum of the two sparse
///   matrices.
template<class OffsetType, class OrdinalType>
OrdinalType
countMergedRowEntries (const OffsetType ptr1[],
                       const OrdinalType ind1[],
                       const OffsetType ptr2[],
                       const OrdinalType ind2[],
                       const OrdinalType i);

/// \fn mergeRows
/// \brief Compute C(i,:) = f( A(i,:), B(i,:) ).
///
/// Merge row i of the sparse matrices A and B, putting the result in
/// C.  All three matrices are stored in three-array compressed sparse
/// row format.  We assume that row i of C contains the correct number
/// of entries.  The merge still stores values that are numerically
/// zero.
///
/// Each sequence of indices must be sorted.  It need not be unique,
/// but if not unique, then the resulting output sequence will not be
/// unique.
///
/// \param ptrOut [in] Previously computed row offsets of sparse matrix C.
/// \param indOut [out] Column indices of sparse matrix C.
/// \param valOut [out] Values of sparse matrix C.
/// \param ptr1 [in] Row offsets of sparse matrix A.
/// \param ind1 [in] Column indices of sparse matrix A.
/// \param val1 [in] Values of sparse matrix A.
/// \param ptr2 [in] Row offsets of sparse matrix B.
/// \param ind2 [in] Column indices of sparse matrix B.
/// \param val2 [in] Values of sparse matrix B.
/// \param i [in] Row index.
///
/// \return One plus the last position in indOut and valOut to which
///   this function wrote.
template<class OffsetType,
         class OrdinalType,
         class ScalarType,
         class BinaryFunctionType>
OrdinalType
mergeRows (const OffsetType ptrOut[],
           OrdinalType indOut[],
           ScalarType valOut[],
           const OffsetType ptr1[],
           const OrdinalType ind1[],
           const ScalarType val1[],
           const OffsetType ptr2[],
           const OrdinalType ind2[],
           const ScalarType val2[],
           const OrdinalType i,
           BinaryFunctionType f);

/// \class AddMatrixEntries
/// \brief Compute C_ij = alpha * A_ij + beta * B_ij.
///
/// You can use this function object in mergeRows() to compute the sum
/// C(i,:) = alpha*A(i,:) + beta*B(i,:).
template<class ScalarType>
class AddMatrixEntries {
public:
  AddMatrixEntries (const ScalarType& alpha, const ScalarType& beta) :
    alpha_ (alpha), beta_ (beta) {}

  inline ScalarType
  operator () (const ScalarType& A_ij, const ScalarType& B_ij) const {
    return alpha_ * A_ij + beta_ * B_ij;
  }
private:
  const ScalarType alpha_;
  const ScalarType beta_;
};

/// \fn addSparseMatrices
/// \brief Compute the sparse matrix sum C = alpha*A + beta*B.
///
/// \param ptrResult [out] Row offsets of sparse matrix C.
/// \param indResult [out] Column indices of sparse matrix C.
/// \param valResult [out] Values of sparse matrix C.
/// \param alpha [in] Coefficient of A in the sum.
/// \param ptr1 [in] Row offsets of sparse matrix A.
/// \param ind1 [in] Column indices of sparse matrix A.
/// \param val1 [in] Values of sparse matrix A.
/// \param beta [in] Coefficient of B in the sum.
/// \param ptr2 [in] Row offsets of sparse matrix B.
/// \param ind2 [in] Column indices of sparse matrix B.
/// \param val2 [in] Values of sparse matrix B.
/// \param numRows [in] Number of rows in A, B, and C.
///
/// The returned arrays are allocated using \c new and must be
/// deallocated using <tt>delete []</tt>.  This function satisfies the
/// strong exception guarantee.
template<class OffsetType, class OrdinalType, class ScalarType>
void
addSparseMatrices (OffsetType*& ptrResult,
                   OrdinalType*& indResult,
                   ScalarType*& valResult,
                   const ScalarType alpha,
                   const OffsetType ptr1[],
                   const OrdinalType ind1[],
                   const ScalarType val1[],
                   const ScalarType beta,
                   const OffsetType ptr2[],
                   const OrdinalType ind2[],
                   const ScalarType val2[],
                   const OrdinalType numRows);

template<class OffsetType, class OrdinalType>
OrdinalType
countMergedRowEntries (const OffsetType ptr1[],
                       const OrdinalType ind1[],
                       const OffsetType ptr2[],
                       const OrdinalType ind2[],
                       const OrdinalType i)
{
  const OffsetType start1 = ptr1[i];
  const OffsetType end1 = ptr1[i+1];
  const OffsetType start2 = ptr2[i];
  const OffsetType end2 = ptr2[i+1];

  OffsetType mark1 = start1, mark2 = start2;
  OrdinalType count = 0;
  while (mark1 < end1 && mark2 < end2) {
    if (ind1[mark1] == ind2[mark2]) {
      ++mark1;
      ++mark2;
    } else if (ind1[mark1] < ind2[mark2]) {
      ++mark1;
    } else { // ind1[mark1] > ind2[mark2]
      ++mark2;
    }
    ++count;
  }
  // Include any remaining entries.
  count += end1 - mark1;
  count += end2 - mark2;
  return count;
}

template<class OffsetType,
         class OrdinalType,
         class ScalarType,
         class BinaryFunctionType>
OrdinalType
mergeRows (const OffsetType ptrOut[],
           OrdinalType indOut[],
           ScalarType valOut[],
           const OffsetType ptr1[],
           const OrdinalType ind1[],
           const ScalarType val1[],
           const OffsetType ptr2[],
           const OrdinalType ind2[],
           const ScalarType val2[],
           const OrdinalType i,
           BinaryFunctionType f)
{
  const OffsetType startOut = ptrOut[i];
  const OffsetType endOut = ptrOut[i+1];
  const OffsetType start1 = ptr1[i];
  const OffsetType end1 = ptr1[i+1];
  const OffsetType start2 = ptr2[i];
  const OffsetType end2 = ptr2[i+1];
  const ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero ();

  OffsetType mark1 = start1, mark2 = start2, markOut = startOut;
  while (mark1 < end1 && mark2 < end2 && markOut < endOut) {
    if (ind1[mark1] == ind2[mark2]) {
      indOut[markOut] = ind1[mark1];
      valOut[markOut] = f (val1[mark1], val2[mark2]);
      ++mark1;
      ++mark2;
    } else if (ind1[mark1] < ind2[mark2]) {
      indOut[markOut] = ind1[mark1];
      // This makes sense if f(x,y) is alpha*x + beta*y.
      valOut[markOut] = f (val1[mark1], zero);
      ++mark1;
    } else { // ind1[mark1] > ind2[mark2]
      indOut[markOut] = ind2[mark2];
      // This makes sense if f(x,y) is alpha*x + beta*y.
      valOut[markOut] = f (zero, val2[mark2]);
      ++mark2;
    }
    ++markOut;
  }
  // Include any remaining entries.
  while (mark1 < end1 && markOut < endOut) {
    indOut[markOut] = ind1[mark1];
    // This makes sense if f(x,y) is alpha*x + beta*y.
    valOut[markOut] = f (val1[mark1], zero);
    ++mark1;
    ++markOut;
  }
  while (mark2 < end2 && markOut < endOut) {
    indOut[markOut] = ind2[mark2];
    // This makes sense if f(x,y) is alpha*x + beta*y.
    valOut[markOut] = f (zero, val2[mark2]);
    ++mark2;
    ++markOut;
  }
  // This is a logic error, because it means either that
  // countMergedRowEntries didn't work, or that it was called
  // incorrectly for this row.
  TEUCHOS_TEST_FOR_EXCEPTION(
    markOut >= endOut && (mark1 < end1 || mark2 < end2),
    std::logic_error,
    "Kokkos::Sequential::mergeRows: Row " << i << " of the output array has "
    << (end1 - mark1) << " + " << (end2 - mark2) << " too few entries.");
  return markOut;
}

template<class OffsetType, class OrdinalType, class ScalarType>
void
addSparseMatrices (OffsetType*& ptrResult,
                   OrdinalType*& indResult,
                   ScalarType*& valResult,
                   const ScalarType alpha,
                   const OffsetType ptr1[],
                   const OrdinalType ind1[],
                   const ScalarType val1[],
                   const ScalarType beta,
                   const OffsetType ptr2[],
                   const OrdinalType ind2[],
                   const ScalarType val2[],
                   const OrdinalType numRows)
{
  typedef Teuchos::ScalarTraits<ScalarType> STS;

  // We don't allocate using ArrayRCP's constructor (that takes the
  // array size), because it initializes the arrays itself.  We want
  // to initialize ourselves, to ensure first-touch allocation.

  OffsetType* ptrOut = NULL;
  OrdinalType* indOut = NULL;
  ScalarType* valOut = NULL;

  // Allocate and fill ptrOut.
  try {
    // We should be careful not to initialize the array here, so that
    // the parallelizable loop that assigns to it below will first-touch
    // initialize it.
    ptrOut = new OffsetType [numRows + 1];

    // Expect that the column indices are sorted and merged.  Compute
    // the number of entries in each row.  We'll do this in place in the
    // row offsets array.  This may be done safely in parallel.  Doing
    // so would first-touch initialize ptrOut.

    if (alpha == STS::zero ()) {
      if (beta == STS::zero ()) {
        // Special case: alpha == 0, beta == 0.
        // The resulting sum has zero entries.
        std::fill (ptrOut, ptrOut + numRows + 1, Teuchos::as<OffsetType> (0));
        ptrResult = ptrOut;
        indResult = NULL;
        valResult = NULL;
        return;
      } else { // alpha == 0, beta != 0
        // Just copy the second matrix, suitably scaled, into the output.
        memcpy (ptrOut, ptr2, (numRows+1) * sizeof (OffsetType));
        const OffsetType numEntries = ptr2[numRows];
        indOut = new OrdinalType [numEntries];
        memcpy (indOut, ind2, numEntries * sizeof (OrdinalType));
        valOut = new ScalarType [numEntries];
        if (beta == STS::one ()) {
          memcpy (valOut, val2, numEntries * sizeof (ScalarType));
        } else { // beta != 1
          for (OrdinalType k = 0; k < numEntries; ++k) {
            valOut[k] = beta * val2[k];
          }
        }
        ptrResult = ptrOut;
        indResult = indOut;
        valResult = valOut;
        return;
      }
    }
    else if (beta == STS::zero ()) { // alpha != 0, beta == 0
      // Just copy the first matrix into the output.
      memcpy (ptrOut, ptr1, (numRows+1) * sizeof (OffsetType));
      const OffsetType numEntries = ptr1[numRows];
      indOut = new OrdinalType [numEntries];
      memcpy (indOut, ind1, numEntries * sizeof (OrdinalType));
      valOut = new ScalarType [numEntries];
      if (alpha == STS::one ()) {
        memcpy (valOut, val1, numEntries * sizeof (ScalarType));
      } else { // alpha != 1
        for (OrdinalType k = 0; k < numEntries; ++k) {
          valOut[k] = alpha * val1[k];
        }
      }
      ptrResult = ptrOut;
      indResult = indOut;
      valResult = valOut;
      return;
    }
    else { // alpha != 0 and beta != 0
      ptrOut[0] = 0;
      for (OrdinalType i = 0; i < numRows; ++i) {
        ptrOut[i+1] = countMergedRowEntries (ptr1, ind1, ptr2, ind2, i);
      }
      // Sum-scan to compute row offsets.
      // This may be parallelized via e.g., TBB's parallel_scan().
      for (OrdinalType i = 1; i < numRows; ++i) {
        ptrOut[i+1] += ptrOut[i];
      }
    }
  } catch (...) {
    if (ptrOut != NULL) {
      delete [] ptrOut;
    }
    if (indOut != NULL) {
      delete [] indOut;
    }
    if (valOut != NULL) {
      delete [] valOut;
    }
    throw;
  }
  //
  // Allocate storage for indices and values.
  //
  const OffsetType numEntries = ptrOut[numRows];

  // Allocate and fill indOut and valOut.
  try {
    indOut = new OrdinalType [numEntries];
    valOut = new ScalarType [numEntries];

    // Merge and add the matrices.  This may be done safely in
    // parallel, since all the arrays have correct sizes and writes to
    // different rows are independent.  We've also already tested for
    // the special cases alpha == 0 and/or beta == 0.
    if (alpha == STS::one () && beta == STS::one ()) {
      for (OrdinalType i = 0; i < numRows; ++i) {
        (void) mergeRows (ptrOut, indOut, valOut,
                          ptr1, ind1, val1,
                          ptr2, ind2, val2, i,
                          std::plus<ScalarType> ());
      }
    } else {
      AddMatrixEntries<ScalarType> f (alpha, beta);
      for (OrdinalType i = 0; i < numRows; ++i) {
        (void) mergeRows (ptrOut, indOut, valOut,
                          ptr1, ind1, val1,
                          ptr2, ind2, val2,
                          i, f);
      }
    }
  } catch (...) {
    if (indOut != NULL) {
      delete [] indOut;
    }
    if (valOut != NULL) {
      delete [] valOut;
    }
    if (ptrOut != NULL) { // we know it's not, but doesn't hurt to test
      delete [] ptrOut;
    }
    throw;
  }

  // "Commit" the output arguments.
  ptrResult = ptrOut;
  indResult = indOut;
  valResult = valOut;
}

} // namespace Sequential
} // namespace Kokkos

#endif // KOKKOS_SEQUENTIAL_ADDSPARSEMATRICES_HPP

