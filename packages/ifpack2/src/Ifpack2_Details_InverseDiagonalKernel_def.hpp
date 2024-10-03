// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_DETAILS_INVERSEDIAGONALKERNEL_DEF_HPP
#define IFPACK2_DETAILS_INVERSEDIAGONALKERNEL_DEF_HPP

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Export_decl.hpp"
#include "Tpetra_Import_decl.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "Teuchos_Assert.hpp"
#include <type_traits>
#include "KokkosSparse_spmv_impl.hpp"

namespace Ifpack2 {
namespace Details {
namespace Impl {

/// \brief Functor for extracting the inverse diagonal of a matrix
///
/// This is an implementation detail of the relaxation smoothers.
template<class DVector,
         class AMatrix,
         class DiagOffsetType,
         bool do_L1,
         bool fix_tiny>
struct InverseDiagonalWithExtraction {

  using execution_space = typename AMatrix::execution_space;
  using LO = typename AMatrix::non_const_ordinal_type;
  using value_type = typename AMatrix::non_const_value_type;
  using team_policy = typename Kokkos::TeamPolicy<execution_space>;
  using team_member = typename team_policy::member_type;
  using ATV = Kokkos::ArithTraits<value_type>;
  // using IST = typename vector_type::impl_scalar_type;
  using magnitude_type = typename ATV::mag_type;
  using MATV = Kokkos::ArithTraits<magnitude_type>;

  DVector m_d;
  AMatrix m_A;
  DiagOffsetType m_offsets;
  magnitude_type m_L1Eta;
  magnitude_type m_MinDiagonalValue;

  InverseDiagonalWithExtraction (const DVector& m_d_,
                                 const AMatrix& m_A_,
                                 const DiagOffsetType& m_offsets_,
                                 const magnitude_type m_L1Eta_,
                                 const magnitude_type m_MinDiagonalValue_) :
    m_d (m_d_),
    m_A (m_A_),
    m_offsets (m_offsets_),
    m_L1Eta (m_L1Eta_),
    m_MinDiagonalValue (m_MinDiagonalValue_)
  {
    const size_t numRows = m_A.numRows ();

    TEUCHOS_ASSERT( numRows == size_t (m_d.extent (0)) );
    TEUCHOS_ASSERT( numRows == size_t (m_offsets.extent (0)) );
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const LO lclRow) const
  {
    const size_t INV = Tpetra::Details::OrdinalTraits<size_t>::invalid ();
    const value_type one = ATV::one();

    // In case the row has no entries
    m_d(lclRow,0) = ATV::zero();

    if (m_offsets(lclRow) != INV) {
      auto curRow = m_A.rowConst (lclRow);
      value_type d = curRow.value(m_offsets(lclRow));

      if (do_L1) {
        // Setup for L1 Methods.
        // Here we add half the value of the off-processor entries in the row,
        // but only if diagonal isn't sufficiently large.
        //
        // This follows from Equation (6.5) in: Baker, Falgout, Kolev and
        // Yang.  "Multigrid Smoothers for Ultraparallel Computing."  SIAM
        // J. Sci. Comput., Vol. 33, No. 5. (2011), pp. 2864-2887.

        const magnitude_type half = MATV::one () / (MATV::one () + MATV::one ());
        const LO numRows = static_cast<LO> (m_A.numRows ());
        const LO row_length = static_cast<LO> (curRow.length);
        magnitude_type diagonal_boost = MATV::zero();
        for (LO iEntry = 0; iEntry < row_length; iEntry++) {
          if (curRow.colidx(iEntry) >= numRows)
            diagonal_boost += ATV::magnitude(curRow.value(iEntry));
        }
        diagonal_boost *= half;
        if (ATV::magnitude(d) < m_L1Eta * diagonal_boost)
          d += diagonal_boost;
      }

      if (fix_tiny) {
        // Replace diagonal entries that are too small.

        if (ATV::magnitude(d) <= m_MinDiagonalValue)
          d = m_MinDiagonalValue;
      }

      // invert diagonal entries
      m_d(lclRow,0) = one / d;
    }
  }

};

} // namespace Impl


template<class TpetraOperatorType>
InverseDiagonalKernel<TpetraOperatorType>::
InverseDiagonalKernel (const Teuchos::RCP<const operator_type>& A)
{
  setMatrix (A);
}

template<class TpetraOperatorType>
void
InverseDiagonalKernel<TpetraOperatorType>::
setMatrix (const Teuchos::RCP<const operator_type>& A)
{
  if (A_op_.get () != A.get ()) {
    A_op_ = A;

    using Teuchos::rcp_dynamic_cast;
    A_crs_ = rcp_dynamic_cast<const crs_matrix_type> (A);

    TEUCHOS_TEST_FOR_EXCEPTION
      (A_crs_.is_null(), std::logic_error,
       "Ifpack2::Details::InverseDiagonalKernel: operator A must be a Tpetra::CrsMatrix.");

    const size_t lclNumRows = A_crs_->getRowMap ()->getLocalNumElements ();

    if (offsets_.extent (0) < lclNumRows) {
      using Kokkos::view_alloc;
      using Kokkos::WithoutInitializing;
      using offsets_view_type = Kokkos::View<size_t*, device_type>;

      offsets_ = offsets_view_type (); // clear 1st to save mem
      auto howAlloc = view_alloc ("offsets", WithoutInitializing);
      offsets_ = offsets_view_type (howAlloc, lclNumRows);
    }

    A_crs_->getCrsGraph ()->getLocalDiagOffsets (offsets_);
  }
}

template<class TpetraOperatorType>
void
InverseDiagonalKernel<TpetraOperatorType>::
compute (vector_type& D_inv,
         bool do_l1, magnitude_type L1Eta,
         bool fixTinyDiagEntries, magnitude_type MinDiagonalValue)
{


    // Canonicalize template arguments to avoid redundant instantiations.
  using d_type = typename vector_type::dual_view_type::t_dev;
  //  using h_matrix_type = typename crs_matrix_type::local_matrix_host_type;
  using d_matrix_type = typename crs_matrix_type::local_matrix_device_type;

  const char kernel_label[] = "inverse_diagonal_kernel";
  using execution_space = typename NT::execution_space;
  using range_type = Kokkos::RangePolicy<execution_space, LO>;
  const size_t lclNumRows = A_crs_->getRowMap ()->getLocalNumElements ();
  auto policy = range_type(0, lclNumRows);

  d_type d = D_inv.getLocalViewDevice(Tpetra::Access::OverwriteAll);
  d_matrix_type a = A_crs_->getLocalMatrixDevice();

  if (do_l1) {
    constexpr bool do_l1_template = true;
    if (fixTinyDiagEntries) {
      constexpr bool fix_tiny_template = true;
      using functor_type =
        Impl::InverseDiagonalWithExtraction<d_type,
                               d_matrix_type,
                               offset_type,
                               do_l1_template,
                               fix_tiny_template>;
      functor_type func (d, a, offsets_, L1Eta, MinDiagonalValue);
      Kokkos::parallel_for (kernel_label, policy, func);
    } else {
      constexpr bool fix_tiny_template = false;
      using functor_type =
        Impl::InverseDiagonalWithExtraction<d_type,
                               d_matrix_type,
                               offset_type,
                               do_l1_template,
                               fix_tiny_template>;
      functor_type func (d, a, offsets_, L1Eta, MinDiagonalValue);
      Kokkos::parallel_for (kernel_label, policy, func);
    }
  } else {
    constexpr bool do_l1_template = false;
    if (fixTinyDiagEntries) {
      constexpr bool fix_tiny_template = true;
      using functor_type =
        Impl::InverseDiagonalWithExtraction<d_type,
                               d_matrix_type,
                               offset_type,
                               do_l1_template,
                               fix_tiny_template>;
      functor_type func (d, a, offsets_, L1Eta, MinDiagonalValue);
      Kokkos::parallel_for (kernel_label, policy, func);
    } else {
      constexpr bool fix_tiny_template = false;
      using functor_type =
        Impl::InverseDiagonalWithExtraction<d_type,
                               d_matrix_type,
                               offset_type,
                               do_l1_template,
                               fix_tiny_template>;
      functor_type func (d, a, offsets_, L1Eta, MinDiagonalValue);
      Kokkos::parallel_for (kernel_label, policy, func);
    }
  }
}

} // namespace Details
} // namespace Ifpack2

#define IFPACK2_DETAILS_INVERSEDIAGONALKERNEL_INSTANT(SC,LO,GO,NT) \
  template class Ifpack2::Details::InverseDiagonalKernel<Tpetra::Operator<SC, LO, GO, NT> >;

#endif // IFPACK2_DETAILS_INVERSEDIAGONALKERNEL_DEF_HPP
