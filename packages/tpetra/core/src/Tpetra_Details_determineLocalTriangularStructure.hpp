// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_DETERMINELOCALTRIANGULARSTRUCTURE_HPP
#define TPETRA_DETAILS_DETERMINELOCALTRIANGULARSTRUCTURE_HPP

/// \file Tpetra_Details_determineLocalTriangularStructure.hpp
/// \brief Declaration and definition of
///   Tpetra::Details::determineLocalTriangularStructure.
///
/// \warning This file, and its contents, are an implementation detail
///   of Tpetra.  They may change or disappear at any time.

#include "Kokkos_Core.hpp"
#include "Tpetra_Details_OrdinalTraits.hpp"

namespace Tpetra {
namespace Details {

/// \brief Return value of determineLocalTriangularStructure.
///
/// \tparam LO Same as Tpetra classes' \c LocalOrdinal template parameter.
template<class LO>
struct LocalTriangularStructureResult {
  //! (Local) number of populated diagonal entries.
  LO diagCount;
  //! Maximum number of entries over all local rows.
  LO maxNumRowEnt;
  //! Whether the graph is locally structurally lower triangular.
  bool couldBeLowerTriangular;
  //! Whether the graph is locally structurally upper triangular.
  bool couldBeUpperTriangular;
};

namespace Impl {
  /// \brief Implementation of
  ///   Tpetra::Details::determineLocalTriangularStructure (which see
  ///   below).
  ///
  /// \warning This is an implementation detail of an implementation
  ///   detail.  It may change or disappear at any time.
  ///
  /// Kokkos::parallel_reduce functor for counting the local
  /// number of diagonal entries in a sparse graph, and determining
  /// whether the graph is lower or upper triangular (or neither).
  ///
  /// \tparam LocalGraphType Kokkos::StaticCrsGraph specialization
  /// \tparam LocalMapType Result of Tpetra::Map::getLocalGraph()
  template<class LocalGraphType, class LocalMapType>
  class DetermineLocalTriangularStructure {
  public:
    // Result can't be more than the number of local rows, so
    // local_ordinal_type is appropriate.
    using result_type =
      LocalTriangularStructureResult<typename LocalMapType::local_ordinal_type>;

    /// \brief Constructor.
    ///
    /// \param G [in] The local sparse graph
    /// \param rowMap [in] The graph's local row Map
    /// \param colMap [in] The graph's local column Map
    /// \param ignoreMapsForTriangularStructure [in] If true, ignore the
    ///   Maps when determining whether the graph is structurally lower or
    ///   upper triangular (or neither).  See GitHub Issue #2658.
    ///   Regardless, use the Maps to count diagonal entries.
    DetermineLocalTriangularStructure (const LocalGraphType& G,
                                       const LocalMapType& rowMap,
                                       const LocalMapType& colMap,
                                       const bool ignoreMapsForTriangularStructure) :
      G_ (G),
      rowMap_ (rowMap),
      colMap_ (colMap),
      ignoreMapsForTriangularStructure_ (ignoreMapsForTriangularStructure)
    {}

    //! Set the initial value of the reduction result.
    KOKKOS_INLINE_FUNCTION void init (result_type& dst) const
    {
      dst.diagCount = 0;
      dst.maxNumRowEnt = 0;
      dst.couldBeLowerTriangular = true; // well, we don't know yet, do we?
      dst.couldBeUpperTriangular = true; // ditto
    }

    KOKKOS_INLINE_FUNCTION void
    join (result_type& dst,
          const result_type& src) const
    {
      dst.diagCount += src.diagCount;
      dst.maxNumRowEnt = (src.maxNumRowEnt > dst.maxNumRowEnt) ?
        src.maxNumRowEnt : dst.maxNumRowEnt;
      dst.couldBeLowerTriangular &= src.couldBeLowerTriangular;
      dst.couldBeUpperTriangular &= src.couldBeUpperTriangular;
    }

    //! Reduction operator: result is (diagonal count, error count).
    KOKKOS_INLINE_FUNCTION void
    operator () (const typename LocalMapType::local_ordinal_type lclRow,
                 result_type& result) const
    {
      using LO = typename LocalMapType::local_ordinal_type;
      using GO = typename LocalMapType::global_ordinal_type;
      using LOT = typename ::Tpetra::Details::OrdinalTraits<LO>;

      auto G_row = G_.rowConst (lclRow);
      const LO numEnt = G_row.length;
      if (numEnt != 0) {
        result.maxNumRowEnt = (numEnt > result.maxNumRowEnt) ?
          numEnt : result.maxNumRowEnt;
        // Use global row and column indices to find the diagonal
        // entry.  Caller promises that local row index is in the row
        // Map on the calling process.
        const GO gblDiagCol = rowMap_.getGlobalElement (lclRow);
        const LO lclDiagCol = colMap_.getLocalElement (gblDiagCol);
        // If it's not in the column Map, then there's no diagonal entry.
        if (lclDiagCol != LOT::invalid ()) {
          // TODO (mfh 25 Apr 2018) Use findRelOffset to optimize for
          // the sorted case, but note that it requires operator[].
          bool foundDiag = false; // don't count duplicates

          if (ignoreMapsForTriangularStructure_) {
            for (LO k = 0; k < numEnt && ! foundDiag; ++k) {
              const LO lclCol = G_row(k);
              if (lclCol == lclDiagCol) {
                foundDiag = true;
              }
            }
            // mfh 30 Apr 2018: See GitHub Issue #2658.  Per
            // current Tpetra::CrsGraph::computeLocalConstants
            // behavior, assume that local column indices are
            // sorted in each row.
            if (numEnt > LO (0)) {
              const LO smallestLclCol = G_row(0);
              const LO largestLclCol = G_row(numEnt-1); // could be same

              if (smallestLclCol < lclRow) {
                result.couldBeUpperTriangular = false;
              }
              if (lclRow < largestLclCol) {
                result.couldBeLowerTriangular = false;
              }
            }
          }
          else {
            for (LO k = 0; k < numEnt &&
                   ((! foundDiag) ||
                    result.couldBeLowerTriangular ||
                    result.couldBeUpperTriangular);
                 ++k) {
              const LO lclCol = G_row(k);
              if (lclCol == lclDiagCol) {
                foundDiag = true;
              }
              else {
                const GO gblCol = colMap_.getGlobalElement (lclCol);
                if (gblCol < gblDiagCol) {
                  result.couldBeUpperTriangular = false;
                }
                if (gblDiagCol < gblCol) {
                  result.couldBeLowerTriangular = false;
                }
              }
            } // for each entry in lclRow
          } // if-else ignoreMapsForTriangularStructure

          if (foundDiag) {
            ++(result.diagCount);
          }
        }
      }
    }

  private:
    LocalGraphType G_;
    LocalMapType rowMap_;
    LocalMapType colMap_;
    bool ignoreMapsForTriangularStructure_;
  };

} // namespace Impl

/// \brief Count the local number of diagonal entries in a local
///   sparse graph, and determine whether the local part of the graph
///   is structurally lower or upper triangular (or neither).
///
/// \warning This is an implementation detail of Tpetra.  It may
///   change or disappear at any time.
///
/// \tparam LocalGraphType Kokkos::StaticCrsGraph specialization
/// \tparam LocalMapType Result of Tpetra::Map::getLocalGraph()
///
/// \param G [in] The local sparse graph
/// \param rowMap [in] The graph's local row Map
/// \param colMap [in] The graph's local column Map
/// \param ignoreMapsForTriangularStructure [in] If true, ignore the
///   Maps when determining whether the graph is structurally lower or
///   upper triangular (or neither).  See GitHub Issue #2658.
///   Regardless, use the Maps to count diagonal entries.
template<class LocalGraphType, class LocalMapType>
LocalTriangularStructureResult<typename LocalMapType::local_ordinal_type>
determineLocalTriangularStructure (const LocalGraphType& G,
                                   const LocalMapType& rowMap,
                                   const LocalMapType& colMap,
                                   const bool ignoreMapsForTriangularStructure)
{
  using LO = typename LocalMapType::local_ordinal_type;
  using execution_space = typename LocalGraphType::device_type::execution_space;
  using range_type = Kokkos::RangePolicy<execution_space, LO>;
  using functor_type =
    Impl::DetermineLocalTriangularStructure<LocalGraphType, LocalMapType>;

  LocalTriangularStructureResult<LO> result {0, 0, true, true};
  Kokkos::parallel_reduce ("Tpetra::Details::determineLocalTriangularStructure",
                           range_type (0, G.numRows ()),
                           functor_type (G, rowMap, colMap,
                                         ignoreMapsForTriangularStructure),
                           result);
  return result;
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_DETERMINELOCALTRIANGULARSTRUCTURE_HPP
