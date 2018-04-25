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

#ifndef TPETRA_DETAILS_GETNUMDIAGS_HPP
#define TPETRA_DETAILS_GETNUMDIAGS_HPP

/// \file Tpetra_Details_getNumDiags.hpp
/// \brief Declaration and definition of getLocalNumDiags and getGlobalNumDiags
///
/// These two functions are meant to help Tpetra developers deprecate
/// and remove the getNodeNumDiags and getGlobalNumDiags methods from
/// various Tpetra classes.  See Trilinos GitHub issue #2630.

#include "Tpetra_CrsGraph.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Tpetra_Details_OrdinalTraits.hpp"

namespace Tpetra {
namespace Details {

namespace Impl {
  /// \brief Kokkos::parallel_reduce functor for counting the local
  ///   number of diagonal entries in a sparse graph.
  ///
  /// \tparam LocalGraphType Kokkos::StaticCrsGraph specialization
  /// \tparam LocalMapType Result of Tpetra::Map::getLocalGraph()
  template<class LocalGraphType, class LocalMapType>
  class CountLocalNumDiags {
  public:
    CountLocalNumDiags (const LocalGraphType& G,
                        const LocalMapType& rowMap,
                        const LocalMapType& colMap) :
      G_ (G), rowMap_ (rowMap), colMap_ (colMap)
    {}

    // Result can't be more than the number of local rows, so
    // local_ordinal_type is appropriate.
    using result_type = typename LocalMapType::local_ordinal_type;

    //! Reduction function: result is (diagonal count, error count).
    KOKKOS_INLINE_FUNCTION void
    operator () (const typename LocalMapType::local_ordinal_type lclRow,
                 result_type& diagCount) const
    {
      using LO = typename LocalMapType::local_ordinal_type;
      using LOT = typename ::Tpetra::Details::OrdinalTraits<LO>;

      auto G_row = G_.rowConst (lclRow);
      const LO numEnt = G_row.length;
      if (numEnt != 0) {
        // Use global row and column indices to find the diagonal
        // entry.  Caller promises that local row index is in the row
        // Map on the calling process.
        const LO lclDiagCol = colMap_.getLocalElement (rowMap_.getGlobalElement (lclRow));
        // If it's not in the column Map, then there's no diagonal entry.
        if (lclDiagCol != LOT::invalid ()) {
          // TODO (mfh 25 Apr 2018) Use findRelOffset to optimize for
          // the sorted case, but note that it requires operator[].
          for (LO k = 0; k < numEnt; ++k) {
            if (lclDiagCol == G_row(k)) {
              ++diagCount;
              break; // don't count duplicates
            }
          }
        }
      }
    }

  private:
    LocalGraphType G_;
    LocalMapType rowMap_;
    LocalMapType colMap_;
  };

  template<class LO, class GO, class NT>
  typename ::Tpetra::CrsGraph<LO, GO, NT>::local_ordinal_type
  countLocalNumDiagsInFillCompleteGraph (const ::Tpetra::CrsGraph<LO, GO, NT>& G)
  {
    using crs_graph_type = ::Tpetra::CrsGraph<LO, GO, NT>;
    using local_map_type = typename crs_graph_type::map_type::local_map_type;
    using local_graph_type = typename crs_graph_type::local_graph_type;
    using functor_type = CountLocalNumDiags<local_graph_type, local_map_type>;
    using execution_space = typename crs_graph_type::device_type::execution_space;
    using policy_type = Kokkos::RangePolicy<execution_space, LO>;

    const auto rowMap = G.getRowMap ();
    const auto colMap = G.getColMap ();
    if (rowMap.get () == nullptr || colMap.get () == nullptr) {
      return 0; // this process does not participate
    }
    else {
      LO lclNumDiags {0};
      functor_type f (G.getLocalGraph (), rowMap->getLocalMap (), colMap->getLocalMap ());
      Kokkos::parallel_reduce (policy_type (0, G.getNodeNumRows ()), f, lclNumDiags);
      return lclNumDiags;
    }
  }

  /// \brief Local columm index of diagonal entry.
  ///
  /// Use global row and column indices to identify the diagonal
  /// entry.  Caller promises that local row index is in the row Map
  /// on the calling process.  Return
  /// Tpetra::Details::OrdinalTraits<LO>::invalid() if there's no
  /// diagonal entry on the calling process, which can happen if the
  /// global row index doesn't live in the column Map on the calling
  /// process.
  template<class MapType>
  typename MapType::local_ordinal_type
  getLocalDiagonalColumnIndex (const typename MapType::local_ordinal_type lclRow,
                               const MapType& rowMap,
                               const MapType& colMap)
  {
    return colMap.getLocalElement (rowMap.getGlobalElement (lclRow));
  }

  //! Return local number of diagonal entries.
  template<class LO, class GO, class NT>
  typename ::Tpetra::RowGraph<LO, GO, NT>::local_ordinal_type
  countLocalNumDiagsInNonFillCompleteLocallyIndexedGraphWithRowViews (const ::Tpetra::RowGraph<LO, GO, NT>& G)
  {
    using LOT = typename ::Tpetra::Details::OrdinalTraits<LO>;

    const auto rowMap = G.getRowMap ();
    const auto colMap = G.getColMap ();
    if (rowMap.get () == nullptr || colMap.get () == nullptr) {
      return 0; // this process does not participate
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION
        (! G.supportsRowViews (), std::logic_error, "Not implemented!");

      Teuchos::ArrayView<const LO> lclColInds;
      const LO lclNumRows = static_cast<LO> (G.getNodeNumRows ());

      LO diagCount = 0;
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        G.getLocalRowView (lclRow, lclColInds);
        const LO numEnt = static_cast<LO> (lclColInds.size ());
        if (numEnt != 0) {
          const LO lclDiagCol = colMap->getLocalElement (rowMap->getGlobalElement (lclRow));
          // If it's not in the column Map, then there's no diagonal entry.
          if (lclDiagCol != LOT::invalid ()) {
            // TODO (mfh 25 Apr 2018) Use findRelOffset to optimize
            // for the sorted case.
            for (LO k = 0; k < numEnt; ++k) {
              if (lclDiagCol == lclColInds[k]) {
                ++diagCount;
                break; // don't count duplicate entries
              }
            } // for each columm index in lclRow
          } // if lclDiagCol is valid
        } // numEnt != 0
      } // for each lclRow

      return diagCount;
    } // if-else
  }

  //! Return local number of diagonal entries.
  template<class LO, class GO, class NT>
  typename ::Tpetra::RowGraph<LO, GO, NT>::local_ordinal_type
  countLocalNumDiagsInNonFillCompleteLocallyIndexedGraphWithoutRowViews (const ::Tpetra::RowGraph<LO, GO, NT>& G)
  {
    using LOT = typename ::Tpetra::Details::OrdinalTraits<LO>;

    const auto rowMap = G.getRowMap ();
    const auto colMap = G.getColMap ();
    if (rowMap.get () == nullptr || colMap.get () == nullptr) {
      return 0; // this process does not participate
    }
    else {
      Teuchos::Array<LO> lclColIndsBuf;
      const LO lclNumRows = static_cast<LO> (G.getNodeNumRows ());

      LO diagCount = 0;
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        size_t numEntSizeT = G.getNumEntriesInLocalRow (lclRow);
        const LO numEnt = static_cast<LO> (numEntSizeT);
        if (static_cast<LO> (lclColIndsBuf.size ()) < numEnt) {
          lclColIndsBuf.resize (numEnt);
        }
        Teuchos::ArrayView<LO> lclColInds = lclColIndsBuf (0, numEnt);
        G.getLocalRowCopy (lclRow, lclColInds, numEntSizeT);

        if (numEnt != 0) {
          const LO lclDiagCol =
            colMap->getLocalElement (rowMap->getGlobalElement (lclRow));
          // If it's not in the column Map, then there's no diagonal entry.
          if (lclDiagCol != LOT::invalid ()) {
            // TODO (mfh 25 Apr 2018) Use findRelOffset to optimize
            // for the sorted case.
            for (LO k = 0; k < numEnt; ++k) {
              if (lclDiagCol == lclColInds[k]) {
                ++diagCount;
                break; // don't count duplicate entries
              }
            } // for each columm index in lclRow
          } // if lclDiagCol is valid
        } // numEnt != 0
      } // for each lclRow

      return diagCount;
    } // if-else
  }

  //! Return local number of diagonal entries.
  template<class LO, class GO, class NT>
  typename ::Tpetra::RowGraph<LO, GO, NT>::local_ordinal_type
  countLocalNumDiagsInNonFillCompleteGloballyIndexedGraphWithRowViews (const ::Tpetra::RowGraph<LO, GO, NT>& G)
  {
    const auto rowMap = G.getRowMap ();
    if (rowMap.get () == nullptr) {
      return 0; // this process does not participate
    }
    else {
      Teuchos::ArrayView<const GO> gblColInds;
      const LO lclNumRows = static_cast<LO> (G.getNodeNumRows ());

      LO diagCount = 0;
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        const GO gblRow = rowMap->getGlobalElement (lclRow);
        G.getGlobalRowView (gblRow, gblColInds);
        const LO numEnt = static_cast<LO> (gblColInds.size ());
        if (numEnt != 0) {
          // TODO (mfh 25 Apr 2018) Use findRelOffset to optimize for
          // the sorted case.
          for (LO k = 0; k < numEnt; ++k) {
            if (gblRow == gblColInds[k]) {
              ++diagCount;
              break; // don't count duplicate entries
            }
          } // for each column index in lclRow
        } // if numEnt != 0
      } // for each lclRow

      return diagCount;
    } // if-else
  }

  //! Return local number of diagonal entries.
  template<class LO, class GO, class NT>
  typename ::Tpetra::RowGraph<LO, GO, NT>::local_ordinal_type
  countLocalNumDiagsInNonFillCompleteGloballyIndexedGraphWithoutRowViews (const ::Tpetra::RowGraph<LO, GO, NT>& G)
  {
    const auto rowMap = G.getRowMap ();
    if (rowMap.get () == nullptr) {
      return 0; // this process does not participate
    }
    else {
      Teuchos::Array<GO> gblColIndsBuf;
      const LO lclNumRows = static_cast<LO> (G.getNodeNumRows ());

      LO diagCount = 0;
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        size_t numEntSizeT = G.getNumEntriesInLocalRow (lclRow);
        const LO numEnt = static_cast<LO> (numEntSizeT);
        if (static_cast<LO> (gblColIndsBuf.size ()) < numEnt) {
          gblColIndsBuf.resize (numEnt);
        }
        Teuchos::ArrayView<GO> gblColInds = gblColIndsBuf (0, numEnt);
        const GO gblRow = rowMap->getGlobalElement (lclRow);
        G.getGlobalRowCopy (gblRow, gblColInds, numEntSizeT);

        if (numEnt != 0) {
          // TODO (mfh 25 Apr 2018) Use findRelOffset to optimize for
          // the sorted case.
          for (LO k = 0; k < numEnt; ++k) {
            if (gblRow == gblColInds[k]) {
              ++diagCount;
              break; // don't count duplicate entries
            }
          } // for each column index in lclRow
        } // if numEnt != 0
      } // for each lclRow

      return diagCount;
    } // if-else
  }

  /// \brief Implementation of Tpetra::Details::getLocalNumDiags (see
  ///   below).
  ///
  /// This is the most generic version.  It expects a RowMatrix,
  /// CrsMatrix, or BlockCrsMatrix -- something that has a getGraph()
  /// method that returns RCP<const RowGraph>.
  template<class MatrixType>
  struct GetLocalNumDiags {
    static typename MatrixType::local_ordinal_type
    getLocalNumDiags (const MatrixType& A)
    {
      using LO = typename MatrixType::local_ordinal_type;
      using GO = typename MatrixType::global_ordinal_type;
      using NT = typename MatrixType::node_type;
      using row_graph_type = ::Tpetra::RowGraph<LO, GO, NT>;

      auto G = A.getGraph ();
      if (G.get () == nullptr) {
        return 0;
      }
      else {
        return GetLocalNumDiags<row_graph_type>::getLocalNumDiags (*G);
      }
    }
  };

  /// \brief Specialization of GetLocalNumDiags for RowGraph.
  template<class LO, class GO, class NT>
  struct GetLocalNumDiags< ::Tpetra::RowGraph<LO, GO, NT> > {
    static LO
    getLocalNumDiags (const ::Tpetra::RowGraph<LO, GO, NT>& G)
    {
      using crs_graph_type = ::Tpetra::CrsGraph<LO, GO, NT>;

      const crs_graph_type* G_crs = dynamic_cast<const crs_graph_type*> (&G);
      if (G_crs != nullptr && G_crs->isFillComplete ()) {
        return countLocalNumDiagsInFillCompleteGraph (*G_crs);
      }
      else {
        if (G.isLocallyIndexed ()) {
          if (G.supportsRowViews ()) {
            return countLocalNumDiagsInNonFillCompleteLocallyIndexedGraphWithRowViews (G);
          }
          else {
            return countLocalNumDiagsInNonFillCompleteLocallyIndexedGraphWithoutRowViews (G);
          }
        }
        else if (G.isGloballyIndexed ()) {
          if (G.supportsRowViews ()) {
            return countLocalNumDiagsInNonFillCompleteGloballyIndexedGraphWithRowViews (G);
          }
          else {
            return countLocalNumDiagsInNonFillCompleteGloballyIndexedGraphWithoutRowViews (G);
          }
        }
        else { // G is empty
          return 0;
        }
      }
    }
  };

  /// \brief Specialization of GetLocalNumDiags for CrsGraph.
  template<class LO, class GO, class NT>
  struct GetLocalNumDiags< ::Tpetra::CrsGraph<LO, GO, NT> > {
    static LO
    getLocalNumDiags (const ::Tpetra::CrsGraph<LO, GO, NT>& G)
    {
      using row_graph_type = ::Tpetra::RowGraph<LO, GO, NT>;
      return GetLocalNumDiags<row_graph_type>::getLocalNumDiags (G);
    }
  };
} // namespace Impl

/// \brief Number of populated diagonal entries in the given sparse
///   graph, on the calling (MPI) process.
template<class CrsGraphType>
typename CrsGraphType::local_ordinal_type
getLocalNumDiags (const CrsGraphType& G)
{
  return Impl::GetLocalNumDiags<CrsGraphType>::getLocalNumDiags (G);
}

/// \brief Number of populated diagonal entries in the given sparse
///   graph, over all processes in the graph's (MPI) communicator.
template<class CrsGraphType>
typename CrsGraphType::global_ordinal_type
getGlobalNumDiags (const CrsGraphType& G)
{
  using GO = typename CrsGraphType::global_ordinal_type;

  const auto map = G.getRowMap ();
  if (map.get () == nullptr) {
    return GO (0); // this process should not participate
  }
  else {
    const auto comm = map->getComm ();
    if (comm.get () == nullptr) {
      return GO (0); // this process should not participate
    }
    else {
      const GO lclNumDiags = static_cast<GO> (getLocalNumDiags (G));

      using Teuchos::REDUCE_SUM;
      using Teuchos::reduceAll;
      using Teuchos::outArg;

      GO gblNumDiags {0};
      reduceAll<int, GO> (*comm, REDUCE_SUM, lclNumDiags, outArg (gblNumDiags));
      return gblNumDiags;
    }
  }
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_GETNUMDIAGS_HPP

