// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_APPLYDIRICHLETBOUNDARYCONDITION_HPP
#define TPETRA_APPLYDIRICHLETBOUNDARYCONDITION_HPP

/// \file Tpetra_applyDirichletBoundaryCondition.hpp
/// \brief Declare and define
///   Tpetra::applyDirichletBoundaryConditionToLocalMatrixRows.

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
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


/// \brief For all k in <tt>[0, lclRowInds.extent(0))</tt>, set local
///   row and column <tt>lclRowInds[k]</tt> of A to have 1 on the diagonal and 0
///   elsewhere.  Run on host, using the default host execution space
///   instance.
///
/// \tparam CrsMatrixType Specialization of Tpetra::CrsMatrix.
///
/// \param A [in/out] Sparse matrix to modify.
/// \param lclRowInds [in] Local indices of the rows of A to modify.
template<class CrsMatrixType>
void
applyDirichletBoundaryConditionToLocalMatrixRowsAndColumns
(const typename CrsMatrixType::execution_space& execSpace,
 CrsMatrixType& A,
 const Kokkos::View<
 typename CrsMatrixType::local_ordinal_type*,
 typename CrsMatrixType::device_type> & lclRowInds);




/// \brief For all k in <tt>[0, lclRowInds.extent(0))</tt>, set local
///   row and column <tt>lclRowInds[k]</tt> of A to have 1 on the diagonal and 0
///   elsewhere.  Run on host, using the default host execution space
///   instance.
///
/// \tparam CrsMatrixType Specialization of Tpetra::CrsMatrix.
///
/// \param A [in/out] Sparse matrix to modify.

/// \param lclRowInds [in] Local indices of the rows of A to modify.
template<class CrsMatrixType>
void
applyDirichletBoundaryConditionToLocalMatrixRowsAndColumns
(CrsMatrixType& A,
 const Kokkos::View<
 typename CrsMatrixType::local_ordinal_type*,
 Kokkos::HostSpace> & lclRowInds);

/// \brief For all k in <tt>[0, lclRowInds.extent(0))</tt>, set local
///   row and column <tt>lclRowInds[k]</tt> of A to have 1 on the diagonal and 0
///   elsewhere.  Run on host, using the default host execution space
///   instance.
///
/// \tparam CrsMatrixType Specialization of Tpetra::CrsMatrix.
///
/// \param A [in/out] Sparse matrix to modify.
/// \param lclRowInds [in] Local indices of the rows of A to modify.
template<class CrsMatrixType>
void
applyDirichletBoundaryConditionToLocalMatrixRowsAndColumns
(CrsMatrixType& A,
 const Kokkos::View<
 typename CrsMatrixType::local_ordinal_type*,
 typename CrsMatrixType::device_type> & lclRowInds);


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
    // Notes for future refactoring:  This routine seems to have one more layer 
    // of options than it probably needs.  For instance, if you passed a Kokkos::Serial 
    // execution_space instance as the first argument you probably wound't need the runOnHost
    // option and then the code below could be collapsed out removing one of the parallel_for's

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
    auto A_lcl = A.getLocalMatrixDevice ();

    const LO lclNumRows = static_cast<LO> (rowMap->getLocalNumElements ());
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
          LO row = lclRowInds(i);
          const GO row_gid = lclRowMap.getGlobalElement(row);
          for (auto j = rowptr(row); j < rowptr(row+1); ++j) {
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
          LO row = lclRowInds(i);
          const GO row_gid = lclRowMap.getGlobalElement(row);
          for (auto j = rowptr(row); j < rowptr(row+1); ++j) {
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



template<class SC, class LO, class GO, class NT>
struct ApplyDirichletBoundaryConditionToLocalMatrixColumns {
  using crs_matrix_type = ::Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using execution_space = typename crs_matrix_type::execution_space;
  using local_col_flag_type =
    Kokkos::View<bool*, Kokkos::AnonymousSpace>;


  static void
  run (const execution_space& execSpace,
       crs_matrix_type& A,
       const local_col_flag_type& lclColFlags,
       const bool runOnHost)
  {
    // Notes for future refactoring:  This routine seems to have one more layer 
    // of options than it probably needs.  For instance, if you passed a Kokkos::Serial 
    // execution_space instance as the first argument you probably wound't need the runOnHost
    // option and then the code below could be collapsed out removing one of the parallel_for's

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
    auto A_lcl = A.getLocalMatrixDevice ();

    const LO lclNumRows = static_cast<LO> (rowMap->getLocalNumElements ());
    TEUCHOS_TEST_FOR_EXCEPTION
      (lclNumRows != 0 && static_cast<LO>(A_lcl.graph.numRows ()) != lclNumRows,
       std::invalid_argument, "The matrix must have been either created "
       "with a KokkosSparse::CrsMatrix, or must have been fill-completed "
       "at least once.");
    
    auto lclRowMap = A.getRowMap()->getLocalMap ();
    auto lclColMap = A.getColMap()->getLocalMap ();
    auto rowptr = A_lcl.graph.row_map;
    auto colind = A_lcl.graph.entries;
    auto values = A_lcl.values;

    const bool wasFillComplete = A.isFillComplete ();
    if (wasFillComplete) {
      A.resumeFill ();
    }

    const LO numRows = (LO) A.getRowMap()->getLocalNumElements();
    if (! runOnHost) {
      using range_type = Kokkos::RangePolicy<execution_space, LO>;
      Kokkos::parallel_for
        ("Tpetra::CrsMatrix apply Dirichlet cols: Device",
         range_type (execSpace, 0, numRows),
         KOKKOS_LAMBDA (const LO i) {
          for (auto j = rowptr(i); j < rowptr(i+1); ++j) {
            if(lclColFlags[colind[j]])
              values(j) = KAT::zero();
          }
        });
    }
    else {
      using range_type =
        Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace, LO>;
      Kokkos::parallel_for
        ("Tpetra::CrsMatrix apply Dirichlet cols: Host",
         range_type (0, numRows),
         KOKKOS_LAMBDA (const LO i) {
          for (auto j = rowptr(i); j < rowptr(i+1); ++j) {
            if(lclColFlags[colind[j]])
              values(j) = KAT::zero();
          }
        });
    }
    if (wasFillComplete) {
      A.fillComplete (A.getDomainMap (), A.getRangeMap ());
    }
  }
};



template<class SC, class LO, class GO, class NT>
void localRowsToColumns(const typename ::Tpetra::CrsMatrix<SC, LO, GO, NT>::execution_space& execSpace, const ::Tpetra::CrsMatrix<SC, LO, GO, NT>& A,const Kokkos::View<const LO*,Kokkos::AnonymousSpace> & dirichletRowIds, Kokkos::View<bool*,Kokkos::AnonymousSpace> & dirichletColFlags) {
  using crs_matrix_type = ::Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using execution_space = typename crs_matrix_type::execution_space;
  using memory_space = typename crs_matrix_type::device_type::memory_space;

  // Need a colMap
  TEUCHOS_TEST_FOR_EXCEPTION(A.getColMap().get () == nullptr, std::invalid_argument,"The matrix must have a column Map.");
  
  // NOTE: We assume that the RowMap and DomainMap of the matrix match.
  // This could get relaxed at a later date, if we need that functionality
  TEUCHOS_TEST_FOR_EXCEPTION(!A.getRowMap()->isSameAs(*A.getDomainMap()),std::invalid_argument, "localRowsToColumns: Row/Domain maps do not match");

  // Assume that the dirichletColFlags array is the correct size
  TEUCHOS_TEST_FOR_EXCEPTION(A.getColMap()->getLocalNumElements() != dirichletColFlags.size(), std::invalid_argument,"localRowsToColumns: dirichletColFlags must be the correct size");

  LO numDirichletRows = (LO) dirichletRowIds.size();
  LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();

  if(A.getCrsGraph()->getImporter().is_null()) {
    // Serial case:  If A doesn't have an importer, just set the flags from the dirichletRowIds
    using range_type = Kokkos::RangePolicy<execution_space, LO>;
    auto lclRowMap = A.getRowMap()->getLocalMap();
    auto lclColMap = A.getColMap()->getLocalMap();

    Kokkos::deep_copy(execSpace,dirichletColFlags,false);
    using range_type = Kokkos::RangePolicy<execution_space, LO>;
    Kokkos::parallel_for
        ("Tpetra::CrsMatrix flag Dirichlet cols",
         range_type (execSpace, 0, numDirichletRows),
         KOKKOS_LAMBDA (const LO i) {          
          GO row_gid = lclRowMap.getGlobalElement(dirichletRowIds[i]);
          LO col_lid = lclColMap.getLocalElement(row_gid);
          if(col_lid != LO_INVALID)
            dirichletColFlags[col_lid] = true;          
        });  
  }
  else {
    // Parallel case: Use A's importer to halo-exchange Dirichlet information
    auto Importer  = A.getCrsGraph()->getImporter();
    auto lclRowMap = A.getRowMap()->getLocalMap();
    auto lclColMap = A.getColMap()->getLocalMap();
    ::Tpetra::Vector<LO,LO,GO,NT> domainDirichlet(A.getDomainMap());
    ::Tpetra::Vector<LO,LO,GO,NT> colDirichlet(A.getColMap());
    const LO one = Teuchos::OrdinalTraits<LO>::one();
    using range_type = Kokkos::RangePolicy<execution_space, LO>;
    {
      auto domain_data = domainDirichlet.template getLocalView<memory_space>(Access::ReadWrite);
      Kokkos::parallel_for
        ("Tpetra::CrsMatrix flag Dirichlet domains",
         range_type (execSpace, 0, numDirichletRows),
         KOKKOS_LAMBDA (const LO i) {          
          GO row_gid = lclRowMap.getGlobalElement(dirichletRowIds[i]);
          LO col_lid = lclColMap.getLocalElement(row_gid);
          if(col_lid != LO_INVALID)
            domain_data(col_lid,0) = one;          
        });  
    }
    colDirichlet.doImport(domainDirichlet,*Importer,::Tpetra::INSERT);
    LO numCols = (LO) A.getColMap()->getLocalNumElements();
    {
      auto col_data = colDirichlet.template getLocalView<memory_space>(Access::ReadOnly);
      Kokkos::parallel_for
        ("Tpetra::CrsMatrix flag Dirichlet cols",
         range_type (execSpace, 0, numCols),
         KOKKOS_LAMBDA (const LO i) {          
          dirichletColFlags[i] = (col_data(i,0) == one) ? true : false;
        });  
    }
  }
}

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
  using memory_space = typename crs_matrix_type::device_type::memory_space;

  using Details::ApplyDirichletBoundaryConditionToLocalMatrixRows;
  using impl_type =
    ApplyDirichletBoundaryConditionToLocalMatrixRows<SC, LO, GO, NT>;

  // Only run on host if we can access the data
  const bool runOnHost = Kokkos::SpaceAccessibility<Kokkos::Serial,memory_space>::accessible;
  if(runOnHost) {
    using local_row_indices_type = Kokkos::View<const LO*, Kokkos::AnonymousSpace>;
    const local_row_indices_type lclRowInds_a (lclRowInds);
    impl_type::run (execution_space (), A, lclRowInds_a, true);
  }
  else {
    auto  lclRowInds_a = Kokkos::create_mirror_view_and_copy(execution_space(),lclRowInds);
    impl_type::run (execution_space (), A, lclRowInds_a, false);
  }
}


template<class CrsMatrixType>
void
applyDirichletBoundaryConditionToLocalMatrixRowsAndColumns
(CrsMatrixType& A,
 const Kokkos::View<
 typename CrsMatrixType::local_ordinal_type*,
 Kokkos::HostSpace> & lclRowInds) {
  using SC = typename CrsMatrixType::scalar_type;
  using LO = typename CrsMatrixType::local_ordinal_type;
  using GO = typename CrsMatrixType::global_ordinal_type;
  using NT = typename CrsMatrixType::node_type;
  using crs_matrix_type = ::Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using execution_space = typename crs_matrix_type::execution_space;
  using memory_space = typename crs_matrix_type::device_type::memory_space;

  TEUCHOS_TEST_FOR_EXCEPTION(A.getColMap().get () == nullptr, std::invalid_argument,"The matrix must have a column Map.");
  
  // Copy the Host array to device
  auto lclRowInds_d = Kokkos::create_mirror_view_and_copy(execution_space(),lclRowInds);

  Kokkos::View<bool*,memory_space> dirichletColFlags("dirichletColFlags",A.getColMap()->getLocalNumElements());
  Kokkos::View<bool*, Kokkos::AnonymousSpace> dirichletColFlags_a(dirichletColFlags);
  Details::localRowsToColumns<SC,LO,GO,NT>(execution_space(),A,lclRowInds_d,dirichletColFlags_a);

  Details::ApplyDirichletBoundaryConditionToLocalMatrixColumns<SC, LO, GO, NT>::run(execution_space(),A,dirichletColFlags,false);
  Details::ApplyDirichletBoundaryConditionToLocalMatrixRows<SC, LO, GO, NT>::run(execution_space(),A,lclRowInds_d,false);
}


template<class CrsMatrixType>
void
applyDirichletBoundaryConditionToLocalMatrixRowsAndColumns
(CrsMatrixType& A,
 const Kokkos::View<
 typename CrsMatrixType::local_ordinal_type*,
 typename CrsMatrixType::device_type> & lclRowInds_d) {
  using SC = typename CrsMatrixType::scalar_type;
  using LO = typename CrsMatrixType::local_ordinal_type;
  using GO = typename CrsMatrixType::global_ordinal_type;
  using NT = typename CrsMatrixType::node_type;
  using crs_matrix_type = ::Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using execution_space = typename crs_matrix_type::execution_space;
  using memory_space = typename crs_matrix_type::device_type::memory_space;

  TEUCHOS_TEST_FOR_EXCEPTION(A.getColMap().get () == nullptr, std::invalid_argument,"The matrix must have a column Map.");
  
  Kokkos::View<bool*,memory_space> dirichletColFlags("dirichletColFlags",A.getColMap()->getLocalNumElements());
  Kokkos::View<bool*, Kokkos::AnonymousSpace> dirichletColFlags_a(dirichletColFlags);
  Details::localRowsToColumns<SC,LO,GO,NT>(execution_space(),A,lclRowInds_d,dirichletColFlags_a);

  Details::ApplyDirichletBoundaryConditionToLocalMatrixColumns<SC, LO, GO, NT>::run(execution_space(),A,dirichletColFlags,false);
  Details::ApplyDirichletBoundaryConditionToLocalMatrixRows<SC, LO, GO, NT>::run(execution_space(),A,lclRowInds_d,false);

}


template<class CrsMatrixType>
void
applyDirichletBoundaryConditionToLocalMatrixRowsAndColumns
(const typename CrsMatrixType::execution_space& execSpace,
 CrsMatrixType& A,
 const Kokkos::View<
 typename CrsMatrixType::local_ordinal_type*,
 typename CrsMatrixType::device_type> & lclRowInds_d) {
  using SC = typename CrsMatrixType::scalar_type;
  using LO = typename CrsMatrixType::local_ordinal_type;
  using GO = typename CrsMatrixType::global_ordinal_type;
  using NT = typename CrsMatrixType::node_type;
  using crs_matrix_type = ::Tpetra::CrsMatrix<SC, LO, GO, NT>;
  //  using execution_space = typename crs_matrix_type::execution_space;
  using memory_space = typename crs_matrix_type::device_type::memory_space;

  TEUCHOS_TEST_FOR_EXCEPTION(A.getColMap().get () == nullptr, std::invalid_argument,"The matrix must have a column Map.");
  
  Kokkos::View<bool*,memory_space> dirichletColFlags("dirichletColFlags",A.getColMap()->getLocalNumElements());
  Kokkos::View<bool*, Kokkos::AnonymousSpace> dirichletColFlags_a(dirichletColFlags);
  Details::localRowsToColumns<SC,LO,GO,NT>(execSpace,A,lclRowInds_d,dirichletColFlags_a);

  Details::ApplyDirichletBoundaryConditionToLocalMatrixColumns<SC, LO, GO, NT>::run(execSpace,A,dirichletColFlags,false);
  Details::ApplyDirichletBoundaryConditionToLocalMatrixRows<SC, LO, GO, NT>::run(execSpace,A,lclRowInds_d,false);

}


} // namespace Tpetra

#endif // TPETRA_APPLYDIRICHLETBOUNDARYCONDITION_HPP
