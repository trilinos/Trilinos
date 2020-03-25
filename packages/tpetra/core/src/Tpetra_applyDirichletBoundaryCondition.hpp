/*
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
*/

#ifndef TPETRA_APPLYDIRICHLETBOUNDARYCONDITION_HPP
#define TPETRA_APPLYDIRICHLETBOUNDARYCONDITION_HPP

/// \file Tpetra_applyDirichletBoundaryCondition.hpp
/// \brief Declare and define
///   Tpetra::applyDirichletBoundaryConditionToLocalMatrixRows.

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "Kokkos_ArithTraits.hpp"

namespace Tpetra {

/// \brief For all k in <tt>[0, lclRowInds.extent(0))</tt>, set local
///   row <tt>lclRowInds[k]</tt> of A to have 1 on the diagonal and 0
///   elsewhere.  Run on device.
///
/// \tparam CrsMatrixType Specialization of Tpetra::CrsMatrix.
///
/// \param execSpace [in] Execution space instance on which to run.
/// \param A [in/out] Sparse matrix to modify.
/// \param lclRowInds [in] Local indices of the rows of A to modify.
template<class CrsMatrixType>
void
applyDirichletBoundaryConditionToLocalMatrixRows
(const typename CrsMatrixType::execution_space& execSpace,
 CrsMatrixType& A,
 const Kokkos::View<
 typename CrsMatrixType::local_ordinal_type*,
 typename CrsMatrixType::device_type> & lclRowInds);

/// \brief For all k in <tt>[0, lclRowInds.extent(0))</tt>, set local
///   row <tt>lclRowInds[k]</tt> of A to have 1 on the diagonal and 0
///   elsewhere.  Run on device, using the default execution space
///   instance.
///
/// \tparam CrsMatrixType Specialization of Tpetra::CrsMatrix.
///
/// \param A [in/out] Sparse matrix to modify.
/// \param lclRowInds [in] Local indices of the rows of A to modify.
template<class CrsMatrixType>
void
applyDirichletBoundaryConditionToLocalMatrixRows
(CrsMatrixType& A,
 const Kokkos::View<
 typename CrsMatrixType::local_ordinal_type*,
 typename CrsMatrixType::device_type>& lclRowInds);
  
/// \brief For all k in <tt>[0, lclRowInds.extent(0))</tt>, set local
///   row <tt>lclRowInds[k]</tt> of A to have 1 on the diagonal and 0
///   elsewhere.  Run on host, using the default host execution space
///   instance.
///
/// \tparam CrsMatrixType Specialization of Tpetra::CrsMatrix.
///
/// \param A [in/out] Sparse matrix to modify.
/// \param lclRowInds [in] Local indices of the rows of A to modify.
template<class CrsMatrixType>
void
applyDirichletBoundaryConditionToLocalMatrixRows
(CrsMatrixType& A,
 const Kokkos::View<
 typename CrsMatrixType::local_ordinal_type*,
 Kokkos::HostSpace> & lclRowInds);

namespace Details {

template<class SC, class LO, class GO, class NT>
struct ApplyDirichletBoundaryConditionToLocalMatrixRows {
  using crs_matrix_type = ::Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using execution_space = typename crs_matrix_type::execution_space;
  using local_row_indices_type =
    Kokkos::View<const LO*, Kokkos::AnonymousSpace>;

  static void
  run (const execution_space& execSpace,
       crs_matrix_type& A,
       const local_row_indices_type& lclRowInds,
       const bool runOnHost)
  {
    using IST = typename crs_matrix_type::impl_scalar_type;
    using KAT = Kokkos::ArithTraits<IST>;

    const auto rowMap = A.getRowMap ();
    TEUCHOS_TEST_FOR_EXCEPTION
      (rowMap.get () == nullptr, std::invalid_argument,
       "The matrix must have a row Map.");
    const auto colMap = A.getColMap ();
    TEUCHOS_TEST_FOR_EXCEPTION
      (colMap.get () == nullptr, std::invalid_argument,
       "The matrix must have a column Map.");
    auto A_lcl = A.getLocalMatrix ();

    const LO lclNumRows = static_cast<LO> (rowMap->getNodeNumElements ());
    TEUCHOS_TEST_FOR_EXCEPTION
      (lclNumRows != 0 && static_cast<LO>(A_lcl.graph.numRows ()) != lclNumRows,
       std::invalid_argument, "The matrix must have been either created "
       "with a KokkosSparse::CrsMatrix, or must have been fill-completed "
       "at least once.");
    
    auto lclRowMap = A.getRowMap ()->getLocalMap ();
    auto lclColMap = A.getColMap ()->getLocalMap ();
    auto rowptr = A_lcl.graph.row_map;
    auto colind = A_lcl.graph.entries;
    auto values = A_lcl.values;

    const bool wasFillComplete = A.isFillComplete ();
    if (wasFillComplete) {
      A.resumeFill ();
    }

    const LO numInputRows = lclRowInds.extent (0);
    if (! runOnHost) {
      using range_type = Kokkos::RangePolicy<execution_space, LO>;
      Kokkos::parallel_for
        ("Tpetra::CrsMatrix apply Dirichlet: Device",
         range_type (execSpace, 0, numInputRows),
         KOKKOS_LAMBDA (const LO i) {
          const GO row_gid = lclRowMap.getGlobalElement (lclRowInds(i));
          for (auto j = rowptr(i); j < rowptr(i+1); ++j) {
            const bool diagEnt =
              lclColMap.getGlobalElement (colind(j)) == row_gid;
            values(j) = diagEnt ? KAT::one () : KAT::zero ();
          }
        });
    }
    else {
      using range_type =
        Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace, LO>;
      Kokkos::parallel_for
        ("Tpetra::CrsMatrix apply Dirichlet: Host",
         range_type (0, numInputRows),
         [&] (const LO i) {
          const GO row_gid = lclRowMap.getGlobalElement (lclRowInds(i));
          for (auto j = rowptr(i); j < rowptr(i+1); ++j) {
            const bool diagEnt =
              lclColMap.getGlobalElement (colind(j)) == row_gid;
            values(j) = diagEnt ? KAT::one () : KAT::zero ();
          }
        });
    }
    if (wasFillComplete) {
      A.fillComplete (A.getDomainMap (), A.getRangeMap ());
    }
  }
};

} // namespace Details

template<class CrsMatrixType>
void
applyDirichletBoundaryConditionToLocalMatrixRows
(const typename CrsMatrixType::execution_space& execSpace,
 CrsMatrixType& A,
 const Kokkos::View<
 typename CrsMatrixType::local_ordinal_type*,
 typename CrsMatrixType::device_type> & lclRowInds)
{
  using SC = typename CrsMatrixType::scalar_type;
  using LO = typename CrsMatrixType::local_ordinal_type;
  using GO = typename CrsMatrixType::global_ordinal_type;
  using NT = typename CrsMatrixType::node_type;

  using local_row_indices_type =
    Kokkos::View<const LO*, Kokkos::AnonymousSpace>;
  const local_row_indices_type lclRowInds_a (lclRowInds);

  using Details::ApplyDirichletBoundaryConditionToLocalMatrixRows;
  using impl_type =
    ApplyDirichletBoundaryConditionToLocalMatrixRows<SC, LO, GO, NT>;
  const bool runOnHost = false;
  impl_type::run (execSpace, A, lclRowInds_a, runOnHost);
}

template<class CrsMatrixType>
void
applyDirichletBoundaryConditionToLocalMatrixRows
(CrsMatrixType& A,
 const Kokkos::View<
 typename CrsMatrixType::local_ordinal_type*,
 typename CrsMatrixType::device_type> & lclRowInds)
{
  using execution_space = typename CrsMatrixType::execution_space;
  applyDirichletBoundaryConditionToLocalMatrixRows (execution_space (), A, lclRowInds);
}

template<class CrsMatrixType>
void
applyDirichletBoundaryConditionToLocalMatrixRows
(CrsMatrixType& A,
 const Kokkos::View<
 typename CrsMatrixType::local_ordinal_type*,
 Kokkos::HostSpace> & lclRowInds)
{
  using SC = typename CrsMatrixType::scalar_type;
  using LO = typename CrsMatrixType::local_ordinal_type;
  using GO = typename CrsMatrixType::global_ordinal_type;
  using NT = typename CrsMatrixType::node_type;
  using crs_matrix_type = ::Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using execution_space = typename crs_matrix_type::execution_space;

  using local_row_indices_type =
    Kokkos::View<const LO*, Kokkos::AnonymousSpace>;
  const local_row_indices_type lclRowInds_a (lclRowInds);

  using Details::ApplyDirichletBoundaryConditionToLocalMatrixRows;
  using impl_type =
    ApplyDirichletBoundaryConditionToLocalMatrixRows<SC, LO, GO, NT>;
  const bool runOnHost = true;
  impl_type::run (execution_space (), A, lclRowInds_a, runOnHost);
}

} // namespace Tpetra

#endif // TPETRA_APPLYDIRICHLETBOUNDARYCONDITION_HPP
