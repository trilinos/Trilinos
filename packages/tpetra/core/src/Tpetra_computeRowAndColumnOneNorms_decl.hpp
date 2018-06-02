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

#ifndef TPETRA_COMPUTEROWANDCOLUMNONENORMS_DECL_HPP
#define TPETRA_COMPUTEROWANDCOLUMNONENORMS_DECL_HPP

/// \file Tpetra_computeRowAndColumnOneNorms_decl.hpp
/// \brief Declaration of Tpetra::computeRowAndColumnOneNorms

#include "TpetraCore_config.h"
#include "Kokkos_ArithTraits.hpp"
#include "Tpetra_Details_EquilibrationInfo.hpp"

namespace Tpetra {

//
// Dear users: These are just forward declarations.  Please skip over
// them and go down to the leftAndOrRightScaleCrsMatrix function
// declaration.  Thank you.
//
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class SC, class LO, class GO, class N>
class RowMatrix;

namespace Details {
template<class ValueType, class DeviceType>
class EquilibrationInfo;
} // namespace Details
#endif // DOXYGEN_SHOULD_SKIP_THIS

/// \brief Compute global row and column one-norms ("row sums" and
///   "column sums") of the input sparse matrix A, in a way suitable
///   for equilibration.
///
/// \note USERS: This is a function you want.
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
