// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_LEFTANDORRIGHTSCALECRSMATRIX_DECL_HPP
#define TPETRA_LEFTANDORRIGHTSCALECRSMATRIX_DECL_HPP

/// \file Tpetra_leftAndOrRightScaleCrsMatrix_decl.hpp
/// \brief Declaration of Tpetra::leftAndOrRightScaleCrsMatrix

#include "TpetraCore_config.h"
#include "Kokkos_ArithTraits.hpp"
#include "Kokkos_Core.hpp"
#include "Tpetra_CrsMatrix_fwd.hpp"
#include "Tpetra_Vector_fwd.hpp"

namespace Tpetra {

/// \brief Whether "scaling" a matrix means multiplying or dividing
///   its entries.
///
/// leftAndOrRightScaleCrsMatrix (see below) takes this as input.
enum EScaling {
  SCALING_MULTIPLY,
  SCALING_DIVIDE
};

/// \brief Left-scale and/or right-scale (in that order) the entries
///   of the input Tpetra::CrsMatrix A.
///
/// \param A [in/out] The sparse matrix A to scale.  It must have a
///   valid KokkosSparse::CrsMatrix.  This is true if fillComplete has
///   been called on it at least once, or if the matrix was created
///   with a local sparse matrix.
///
/// \param rowScalingFactors [in] Use
///   Details::EquilibrationInfo::rowNorms.
///
/// \param colScalingFactors [in] If assumeSymmetric, use
///   Details::EquilibrationInfo::colNorms, else use
///   Details::EquilibrationInfo::rowScaledColNorms.
///
/// \param leftScale [in] Whether to left-scale A.  Left scaling
///   happens first.
///
/// \param rightScale [in] Whether to right-scale A.  Right scaling
///   happens last.
///
/// \param Whether to assume symmetric scaling, that is, whether to
///   take square roots of scaling factors before scaling.
///
/// \param scaling [in] If SCALING_DIVIDE, "scale" means "divide by";
///   if SCALING_MULTIPLY, it means "multiply by."
template<class SC, class LO, class GO, class NT>
void
leftAndOrRightScaleCrsMatrix (Tpetra::CrsMatrix<SC, LO, GO, NT>& A,
                              const Kokkos::View<
                                const typename Kokkos::ArithTraits<SC>::mag_type*,
                                typename NT::device_type>& rowScalingFactors,
                              const Kokkos::View<
                                const typename Kokkos::ArithTraits<SC>::mag_type*,
                                typename NT::device_type>& colScalingFactors,
                              const bool leftScale,
                              const bool rightScale,
                              const bool assumeSymmetric,
                              const EScaling scaling);

/// \brief Left-scale and/or right-scale (in that order) the entries
///   of the input Tpetra::CrsMatrix A.
///
/// \param A [in/out] The sparse matrix A to scale.  It must have a
///   valid KokkosSparse::CrsMatrix.  This is true if fillComplete has
///   been called on it at least once, or if the matrix was created
///   with a local sparse matrix.
///
/// \param rowScalingFactors [in] The row scaling factors; must be
///   distributed according to the matrix's row Map.  This function
///   does NOT redistribute the Vector for you.
///
/// \param colScalingFactors [in] The column scaling factors; must be
///   distributed according to the matrix's column Map.  This function
///   does NOT redistribute the Vector for you.
///
/// \param leftScale [in] Whether to left-scale A.  Left scaling
///   happens first.
///
/// \param rightScale [in] Whether to right-scale A.  Right scaling
///   happens last.
///
/// \param Whether to assume symmetric scaling, that is, whether to
///   take square roots of scaling factors before scaling.
///
/// \param scaling [in] If SCALING_DIVIDE, "scale" means "divide by";
///   if SCALING_MULTIPLY, it means "multiply by."
template<class SC, class LO, class GO, class NT>
void
leftAndOrRightScaleCrsMatrix (Tpetra::CrsMatrix<SC, LO, GO, NT>& A,
                              const Tpetra::Vector<typename Kokkos::ArithTraits<SC>::mag_type,
                                LO, GO, NT>& rowScalingFactors,
                              const Tpetra::Vector<typename Kokkos::ArithTraits<SC>::mag_type,
                                LO, GO, NT>& colScalingFactors,
                              const bool leftScale,
                              const bool rightScale,
                              const bool assumeSymmetric,
                              const EScaling scaling);

} // namespace Tpetra

#endif // TPETRA_LEFTANDORRIGHTSCALECRSMATRIX_DECL_HPP
