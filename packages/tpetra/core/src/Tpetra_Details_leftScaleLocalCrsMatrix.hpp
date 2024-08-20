// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_LEFTSCALELOCALCRSMATRIX_HPP
#define TPETRA_DETAILS_LEFTSCALELOCALCRSMATRIX_HPP

/// \file Tpetra_Details_leftScaleLocalCrsMatrix.hpp
/// \brief Declaration and definition of
///   Tpetra::Details::leftScaleLocalCrsMatrix
///
/// \warning This file, and its contents, are implementation details
///   of Tpetra.  They may disappear or change at any time.

#include "TpetraCore_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"
#include <type_traits>

namespace Tpetra {
namespace Details {

/// \brief Kokkos::parallel_for functor that left-scales a
///   KokkosSparse::CrsMatrix.
///
/// \tparam LocalSparseMatrixType KokkosSparse::CrsMatrix specialization.
/// \tparam ScalingFactorsViewType Kokkos::View specialization storing
///   scaling factors by which to divide the rows of the local sparse
///   matrix.
/// \tparam divide If true, divide, else multiply.
template<class LocalSparseMatrixType,
         class ScalingFactorsViewType,
         const bool divide>
class LeftScaleLocalCrsMatrix {
public:
  using val_type =
    typename std::remove_const<typename LocalSparseMatrixType::value_type>::type;
  using mag_type = typename ScalingFactorsViewType::non_const_value_type;
  static_assert (ScalingFactorsViewType::rank == 1,
                 "scalingFactors must be a rank-1 Kokkos::View.");
  using device_type = typename LocalSparseMatrixType::device_type;
  using LO = typename LocalSparseMatrixType::ordinal_type;
  using policy_type = Kokkos::TeamPolicy<typename device_type::execution_space, LO>;

  /// \param A_lcl [in/out] The local sparse matrix.
  ///
  /// \param scalingFactors [in] rowNorms from EquilibrationInfo.
  ///
  /// \param assumeSymmetric [in] Whether to assume that the (global)
  ///   sparse matrix is symmetric.  If true, divde matrix entries by
  ///   square roots of scaling factors; else, divide by the scaling
  ///   factors themselves.
  LeftScaleLocalCrsMatrix (const LocalSparseMatrixType& A_lcl,
                           const ScalingFactorsViewType& scalingFactors,
                           const bool assumeSymmetric) :
    A_lcl_ (A_lcl),
    scalingFactors_ (scalingFactors),
    assumeSymmetric_ (assumeSymmetric)
  {}

  KOKKOS_INLINE_FUNCTION void
  operator () (const typename policy_type::member_type & team) const
  {
    using KAM = Kokkos::ArithTraits<mag_type>;

    const LO lclRow = team.league_rank();
    const mag_type curRowNorm = scalingFactors_(lclRow);
    // Users are responsible for any divisions or multiplications by
    // zero.
    const mag_type scalingFactor = assumeSymmetric_ ?
      KAM::sqrt (curRowNorm) : curRowNorm;
    auto curRow = A_lcl_.row (lclRow);
    const LO numEnt = curRow.length;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, numEnt), [&](const LO k) {
      if (divide) { // constexpr, so should get compiled out
        curRow.value (k) = curRow.value(k) / scalingFactor;
      }
      else {
        curRow.value (k) = curRow.value(k) * scalingFactor;
      }
    });
  }

private:
  LocalSparseMatrixType A_lcl_;
  typename ScalingFactorsViewType::const_type scalingFactors_;
  bool assumeSymmetric_;
};

/// \brief Left-scale a KokkosSparse::CrsMatrix.
///
/// \tparam LocalSparseMatrixType KokkosSparse::CrsMatrix specialization.
/// \tparam ScalingFactorsViewType Kokkos::View specialization storing
///   scaling factors by which to divide the rows of the local sparse
///   matrix.
///
/// \param A_lcl [in/out] The local sparse matrix.
/// \param scalingFactors [in] Row scaling factors.
/// \param assumeSymmetric [in] If true, divide matrix entries by
///   square roots of scaling factors; else, divide by the scaling
///   factors themselves.
/// \param divide [in] If true, divide; else multiply.
template<class LocalSparseMatrixType, class ScalingFactorsViewType>
void
leftScaleLocalCrsMatrix (const LocalSparseMatrixType& A_lcl,
                         const ScalingFactorsViewType& scalingFactors,
                         const bool assumeSymmetric,
                         const bool divide = true)
{
  using device_type = typename LocalSparseMatrixType::device_type;
  using execution_space = typename device_type::execution_space;
  using LO = typename LocalSparseMatrixType::ordinal_type;
  using policy_type = Kokkos::TeamPolicy<execution_space, LO>;

  const LO lclNumRows = A_lcl.numRows ();
  if (divide) {
    using functor_type =
      LeftScaleLocalCrsMatrix<LocalSparseMatrixType,
        typename ScalingFactorsViewType::const_type, true>;
    functor_type functor (A_lcl, scalingFactors, assumeSymmetric);
    Kokkos::parallel_for ("leftScaleLocalCrsMatrix",
                          policy_type (lclNumRows, Kokkos::AUTO), functor);
  }
  else {
    using functor_type =
      LeftScaleLocalCrsMatrix<LocalSparseMatrixType,
        typename ScalingFactorsViewType::const_type, false>;
    functor_type functor (A_lcl, scalingFactors, assumeSymmetric);
    Kokkos::parallel_for ("leftScaleLocalCrsMatrix",
                          policy_type (lclNumRows, Kokkos::AUTO), functor);
  }
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_LEFTSCALELOCALCRSMATRIX_HPP
