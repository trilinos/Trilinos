// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

#ifndef TPETRA_LEFTANDORRIGHTSCALECRSMATRIX_DECL_HPP
#define TPETRA_LEFTANDORRIGHTSCALECRSMATRIX_DECL_HPP

/// \file Tpetra_leftAndOrRightScaleCrsMatrix_decl.hpp
/// \brief Declaration of Tpetra::leftAndOrRightScaleCrsMatrix

#include "TpetraCore_config.h"
#include "Kokkos_ArithTraits.hpp"
#include "Kokkos_Core.hpp"

namespace Tpetra {

//
// Dear users: These are just forward declarations.  Please skip over
// them and go down to the leftAndOrRightScaleCrsMatrix function
// declaration.  Thank you.
//
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class SC, class LO, class GO, class N>
class CrsMatrix;

template<class SC, class LO, class GO, class N>
class Vector;
#endif // DOXYGEN_SHOULD_SKIP_THIS

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
