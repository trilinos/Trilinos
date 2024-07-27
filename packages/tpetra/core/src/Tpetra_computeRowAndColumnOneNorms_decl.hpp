// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_COMPUTEROWANDCOLUMNONENORMS_DECL_HPP
#define TPETRA_COMPUTEROWANDCOLUMNONENORMS_DECL_HPP

/// \file Tpetra_computeRowAndColumnOneNorms_decl.hpp
/// \brief Declaration of Tpetra::computeRowAndColumnOneNorms

#include "TpetraCore_config.h"
#include "Kokkos_ArithTraits.hpp"
#include "Tpetra_Details_EquilibrationInfo.hpp"
#include "Tpetra_RowMatrix_fwd.hpp"

namespace Tpetra {

/// \brief Compute global row one-norms ("row sums") of the input
///   sparse matrix A, in a way suitable for one-sided (left only)
///   equilibration.
///
/// \note This function is collective over A's communicator, and may
///   need to communicate, depending on A's Maps.
///
/// \param A [in] The input sparse matrix A.
///
/// \return Input to leftAndOrRightScaleCrsMatrix (which see).  The
///   result is only safe to use for left scaling, not for right
///   scaling.
template<class SC, class LO, class GO, class NT>
Details::EquilibrationInfo<typename Kokkos::ArithTraits<SC>::val_type,
                           typename NT::device_type>
computeRowOneNorms (const Tpetra::RowMatrix<SC, LO, GO, NT>& A);

/// \brief Compute global row and column one-norms ("row sums" and
///   "column sums") of the input sparse matrix A, in a way suitable
///   for two-sided (left and right) equilibration.
///
/// \note For AztecOO users: If you set assumeSymmetric=true, this
///   function should behave like setting the <tt>AZ_scaling</tt>
///   option to <tt>AZ_sym_row_sum</tt>.
///
/// \note This function is collective over A's communicator, and may
///   need to communicate, depending on A's Maps.
///
/// For the nonsymmetric case, this function works like a sparse
/// version of LAPACK's DGEEQU routine, except that it uses one norms
/// (sums of absolute values) instead of infinity norms (maximum
/// absolute value).  The resulting row and column scaling is NOT
/// symmetric.  For the symmetric case, this function computes the row
/// norms and uses those for the column norms.  The resulting scaling
/// is symmetric IF you take square roots.
///
/// \param A [in] The input sparse matrix A.
///
/// \param assumeSymmetric [in] Whether to assume that the matrix A is
///   (globally) symmetric.  If so, don't compute row-scaled column
///   norms separately from row norms.
///
/// \return Input to leftAndOrRightScaleCrsMatrix (which see).
template<class SC, class LO, class GO, class NT>
Details::EquilibrationInfo<typename Kokkos::ArithTraits<SC>::val_type,
                           typename NT::device_type>
computeRowAndColumnOneNorms (const Tpetra::RowMatrix<SC, LO, GO, NT>& A,
                             const bool assumeSymmetric);

} // namespace Tpetra

#endif // TPETRA_COMPUTEROWANDCOLUMNONENORMS_DECL_HPP
