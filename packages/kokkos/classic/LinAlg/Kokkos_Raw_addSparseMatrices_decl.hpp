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

#ifndef __Kokkos_Raw_addSparseMatrices_decl_hpp
#define __Kokkos_Raw_addSparseMatrices_decl_hpp

#include <Kokkos_ConfigDefs.hpp>

namespace Kokkos {
namespace Raw {

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

} // namespace Raw
} // namespace Kokkos

#endif // __Kokkos_Raw_addSparseMatrices_decl_hpp

