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
// ************************************************************************
// @HEADER

#ifndef TPETRA_CRSGRAPH_DEF_HPP
#define TPETRA_CRSGRAPH_DEF_HPP

/// \file Tpetra_CrsGraph_def.hpp
/// \brief Definition of the Tpetra::CrsGraph class
///
/// If you want to use Tpetra::CrsGraph, include "Tpetra_CrsGraph.hpp"
/// (a file which CMake generates and installs for you).  If you only
/// want the declaration of Tpetra::CrsGraph, include
/// "Tpetra_CrsGraph_decl.hpp".

#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_computeOffsets.hpp"
#include "Tpetra_Details_copyOffsets.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_Details_getGraphDiagOffsets.hpp"
#include "Tpetra_Details_makeColMap.hpp"
#include "Tpetra_Details_Profiling.hpp"
#include "Tpetra_Details_getEntryOnHost.hpp"
#include "Tpetra_Details_crsUtils.hpp"
#include "Tpetra_Distributor.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Import_Util.hpp"
#include "Tpetra_Import_Util2.hpp"
#include "Tpetra_Details_packCrsGraph.hpp"
#include "Tpetra_Details_unpackCrsGraphAndCombine.hpp"
#include "Tpetra_Details_CrsPadding.hpp"
#include "Tpetra_Util.hpp"
#include <algorithm>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace Tpetra {
  namespace Details {
    namespace Impl {

      template<class MapIter>
      void
      verbosePrintMap(std::ostream& out,
                      MapIter beg,
                      MapIter end,
                      const size_t numEnt,
                      const char mapName[])
      {
        using ::Tpetra::Details::Behavior;
        using ::Tpetra::Details::verbosePrintArray;

        out << mapName << ": {";
        const size_t maxNumToPrint =
          Behavior::verbosePrintCountThreshold();
        if (maxNumToPrint == 0) {
          if (numEnt != 0) {
            out << "...";
          }
        }
        else {
          const size_t numToPrint = numEnt > maxNumToPrint ?
            maxNumToPrint : numEnt;
          size_t count = 0;
          for (MapIter it = beg; it != end; ++it) {
            out << "(" << (*it).first << ", ";
            verbosePrintArray(out, (*it).second, "gblColInds",
                              maxNumToPrint);
            out << ")";
            if (count + size_t(1) < numToPrint) {
              out << ", ";
            }
            ++count;
          }
          if (count < numEnt) {
            out << ", ...";
          }
        }
        out << "}";
      }

      template<class LO, class GO, class Node>
      Teuchos::ArrayView<GO>
      getRowGraphGlobalRow(
        std::vector<GO>& gblColIndsStorage,
        const RowGraph<LO, GO, Node>& graph,
        const GO gblRowInd)
      {
        size_t origNumEnt = graph.getNumEntriesInGlobalRow(gblRowInd);
        if (gblColIndsStorage.size() < origNumEnt) {
          gblColIndsStorage.resize(origNumEnt);
        }
        Teuchos::ArrayView<GO> gblColInds(gblColIndsStorage.data(),
                                          origNumEnt);
        graph.getGlobalRowCopy(gblRowInd, gblColInds, origNumEnt);
        return gblColInds;
      }

      template<class LO, class GO, class DT, class OffsetType, class NumEntType>
      class ConvertColumnIndicesFromGlobalToLocal {
      public:
        ConvertColumnIndicesFromGlobalToLocal (const ::Kokkos::View<LO*, DT>& lclColInds,
                                               const ::Kokkos::View<const GO*, DT>& gblColInds,
                                               const ::Kokkos::View<const OffsetType*, DT>& ptr,
                                               const ::Tpetra::Details::LocalMap<LO, GO, DT>& lclColMap,
                                               const ::Kokkos::View<const NumEntType*, DT>& numRowEnt) :
          lclColInds_ (lclColInds),
          gblColInds_ (gblColInds),
          ptr_ (ptr),
          lclColMap_ (lclColMap),
          numRowEnt_ (numRowEnt)
        {}

        KOKKOS_FUNCTION void
        operator () (const LO& lclRow, OffsetType& curNumBad) const
        {
          const OffsetType offset = ptr_(lclRow);
          // NOTE (mfh 26 Jun 2016) It's always legal to cast the number
          // of entries in a row to LO, as long as the row doesn't have
          // too many duplicate entries.
          const LO numEnt = static_cast<LO> (numRowEnt_(lclRow));
          for (LO j = 0; j < numEnt; ++j) {
            const GO gid = gblColInds_(offset + j);
            const LO lid = lclColMap_.getLocalElement (gid);
            lclColInds_(offset + j) = lid;
            if (lid == ::Tpetra::Details::OrdinalTraits<LO>::invalid ()) {
              ++curNumBad;
            }
          }
        }

        static OffsetType
        run (const ::Kokkos::View<LO*, DT>& lclColInds,
             const ::Kokkos::View<const GO*, DT>& gblColInds,
             const ::Kokkos::View<const OffsetType*, DT>& ptr,
             const ::Tpetra::Details::LocalMap<LO, GO, DT>& lclColMap,
             const ::Kokkos::View<const NumEntType*, DT>& numRowEnt)
        {
          typedef ::Kokkos::RangePolicy<typename DT::execution_space, LO> range_type;
          typedef ConvertColumnIndicesFromGlobalToLocal<LO, GO, DT, OffsetType, NumEntType> functor_type;

          const LO lclNumRows = ptr.extent (0) == 0 ?
            static_cast<LO> (0) : static_cast<LO> (ptr.extent (0) - 1);
          OffsetType numBad = 0;
          // Count of "bad" column indices is a reduction over rows.
          ::Kokkos::parallel_reduce (range_type (0, lclNumRows),
                                     functor_type (lclColInds, gblColInds, ptr,
                                                   lclColMap, numRowEnt),
                                     numBad);
          return numBad;
        }

      private:
        ::Kokkos::View<LO*, DT> lclColInds_;
        ::Kokkos::View<const GO*, DT> gblColInds_;
        ::Kokkos::View<const OffsetType*, DT> ptr_;
        ::Tpetra::Details::LocalMap<LO, GO, DT> lclColMap_;
        ::Kokkos::View<const NumEntType*, DT> numRowEnt_;
      };

    } // namespace Impl

    /// \brief Convert a CrsGraph's global column indices into local
    ///   column indices.
    ///
    /// \param lclColInds [out] On output: The graph's local column
    ///   indices.  This may alias gblColInds, if LO == GO.
    /// \param gblColInds [in] On input: The graph's global column
    ///   indices.  This may alias lclColInds, if LO == GO.
    /// \param ptr [in] The graph's row offsets.
    /// \param lclColMap [in] "Local" (threaded-kernel-worthy) version
    ///   of the column Map.
    /// \param numRowEnt [in] Array with number of entries in each row.
    ///
    /// \return the number of "bad" global column indices (that don't
    ///   live in the column Map on the calling process).
    template<class LO, class GO, class DT, class OffsetType, class NumEntType>
    OffsetType
    convertColumnIndicesFromGlobalToLocal (const Kokkos::View<LO*, DT>& lclColInds,
                                           const Kokkos::View<const GO*, DT>& gblColInds,
                                           const Kokkos::View<const OffsetType*, DT>& ptr,
                                           const LocalMap<LO, GO, DT>& lclColMap,
                                           const Kokkos::View<const NumEntType*, DT>& numRowEnt)
    {
      using Impl::ConvertColumnIndicesFromGlobalToLocal;
      typedef ConvertColumnIndicesFromGlobalToLocal<LO, GO, DT, OffsetType, NumEntType> impl_type;
      return impl_type::run (lclColInds, gblColInds, ptr, lclColMap, numRowEnt);
    }

    template<class ViewType, class LO>
    class MaxDifference {
    public:
      MaxDifference (const ViewType& ptr) : ptr_ (ptr) {}

      KOKKOS_INLINE_FUNCTION void init (LO& dst) const {
        dst = 0;
      }

      KOKKOS_INLINE_FUNCTION void
      join (volatile LO& dst, const volatile LO& src) const
      {
        dst = (src > dst) ? src : dst;
      }

      KOKKOS_INLINE_FUNCTION void
      operator () (const LO lclRow, LO& maxNumEnt) const
      {
        const LO numEnt = static_cast<LO> (ptr_(lclRow+1) - ptr_(lclRow));
        maxNumEnt = (numEnt > maxNumEnt) ? numEnt : maxNumEnt;
      }
    private:
      typename ViewType::const_type ptr_;
    };

    template<class ViewType, class LO>
    typename ViewType::non_const_value_type
    maxDifference (const char kernelLabel[],
                   const ViewType& ptr,
                   const LO lclNumRows)
    {
      if (lclNumRows == 0) {
        // mfh 07 May 2018: Weirdly, I need this special case,
        // otherwise I get the wrong answer.
        return static_cast<LO> (0);
      }
      else {
        using execution_space = typename ViewType::execution_space;
        using range_type = Kokkos::RangePolicy<execution_space, LO>;
        LO theMaxNumEnt {0};
        Kokkos::parallel_reduce (kernelLabel,
                                 range_type (0, lclNumRows),
                                 MaxDifference<ViewType, LO> (ptr),
                                 theMaxNumEnt);
        return theMaxNumEnt;
      }
    }

  } // namespace Details

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getDebug() {
    return Details::Behavior::debug("CrsGraph");
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getVerbose() {
    return Details::Behavior::verbose("CrsGraph");
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
            const size_t maxNumEntriesPerRow,
            const ProfileType /* pftype */,
            const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
    , rowMap_ (rowMap)
    , numAllocForAllRows_ (maxNumEntriesPerRow)
  {
    const char tfecfFuncName[] =
      "CrsGraph(rowMap,maxNumEntriesPerRow,pftype,params): ";
    staticAssertions ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (maxNumEntriesPerRow == Teuchos::OrdinalTraits<size_t>::invalid (),
       std::invalid_argument, "The allocation hint maxNumEntriesPerRow must be "
       "a valid size_t value, which in this case means it must not be "
       "Teuchos::OrdinalTraits<size_t>::invalid().");
    resumeFill (params);
    checkInternalState ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
            const Teuchos::RCP<const map_type>& colMap,
            const size_t maxNumEntriesPerRow,
            const ProfileType /* pftype */,
            const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
    , rowMap_ (rowMap)
    , colMap_ (colMap)
    , numAllocForAllRows_ (maxNumEntriesPerRow)
  {
    const char tfecfFuncName[] =
      "CrsGraph(rowMap,colMap,maxNumEntriesPerRow,pftype,params): ";
    staticAssertions ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      maxNumEntriesPerRow == Teuchos::OrdinalTraits<size_t>::invalid (),
      std::invalid_argument, "The allocation hint maxNumEntriesPerRow must be "
      "a valid size_t value, which in this case means it must not be "
      "Teuchos::OrdinalTraits<size_t>::invalid().");
    resumeFill (params);
    checkInternalState ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
            const Teuchos::ArrayView<const size_t>& numEntPerRow,
            const ProfileType /* pftype */,
            const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
    , rowMap_ (rowMap)
    , numAllocForAllRows_ (0)
  {
    const char tfecfFuncName[] =
      "CrsGraph(rowMap,numEntPerRow,pftype,params): ";
    staticAssertions ();

    const size_t lclNumRows = rowMap.is_null () ?
      static_cast<size_t> (0) : rowMap->getNodeNumElements ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (numEntPerRow.size ()) != lclNumRows,
      std::invalid_argument, "numEntPerRow has length " << numEntPerRow.size ()
      << " != the local number of rows " << lclNumRows << " as specified by "
      "the input row Map.");

    if (debug_) {
      for (size_t r = 0; r < lclNumRows; ++r) {
        const size_t curRowCount = numEntPerRow[r];
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (curRowCount == Teuchos::OrdinalTraits<size_t>::invalid (),
           std::invalid_argument, "numEntPerRow(" << r << ") "
           "specifies an invalid number of entries "
           "(Teuchos::OrdinalTraits<size_t>::invalid()).");
      }
    }

    // Deep-copy the (host-accessible) input into k_numAllocPerRow_.
    // The latter is a const View, so we have to copy into a nonconst
    // View first, then assign.
    typedef decltype (k_numAllocPerRow_) out_view_type;
    typedef typename out_view_type::non_const_type nc_view_type;
    typedef Kokkos::View<const size_t*,
                         typename nc_view_type::array_layout,
                         Kokkos::HostSpace,
                         Kokkos::MemoryUnmanaged> in_view_type;
    in_view_type numAllocPerRowIn (numEntPerRow.getRawPtr (), lclNumRows);
    nc_view_type numAllocPerRowOut ("Tpetra::CrsGraph::numAllocPerRow",
                                    lclNumRows);
    Kokkos::deep_copy (numAllocPerRowOut, numAllocPerRowIn);
    k_numAllocPerRow_ = numAllocPerRowOut;

    resumeFill (params);
    checkInternalState ();
  }



  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
            const Kokkos::DualView<const size_t*, execution_space>& numEntPerRow,
            const ProfileType /* pftype */,
            const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
    , rowMap_ (rowMap)
    , k_numAllocPerRow_ (numEntPerRow.h_view)
    , numAllocForAllRows_ (0)
  {
    const char tfecfFuncName[] =
      "CrsGraph(rowMap,numEntPerRow,pftype,params): ";
    staticAssertions ();

    const size_t lclNumRows = rowMap.is_null () ?
      static_cast<size_t> (0) : rowMap->getNodeNumElements ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (numEntPerRow.extent (0)) != lclNumRows,
      std::invalid_argument, "numEntPerRow has length " <<
      numEntPerRow.extent (0) << " != the local number of rows " <<
      lclNumRows << " as specified by " "the input row Map.");

    if (debug_) {
      for (size_t r = 0; r < lclNumRows; ++r) {
        const size_t curRowCount = numEntPerRow.h_view(r);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (curRowCount == Teuchos::OrdinalTraits<size_t>::invalid (),
           std::invalid_argument, "numEntPerRow(" << r << ") "
           "specifies an invalid number of entries "
           "(Teuchos::OrdinalTraits<size_t>::invalid()).");
      }
    }

    resumeFill (params);
    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
            const Teuchos::RCP<const map_type>& colMap,
            const Kokkos::DualView<const size_t*, execution_space>& numEntPerRow,
            const ProfileType /* pftype */,
            const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
    , rowMap_ (rowMap)
    , colMap_ (colMap)
    , k_numAllocPerRow_ (numEntPerRow.h_view)
    , numAllocForAllRows_ (0)
  {
    const char tfecfFuncName[] =
      "CrsGraph(rowMap,colMap,numEntPerRow,pftype,params): ";
    staticAssertions ();

    const size_t lclNumRows = rowMap.is_null () ?
      static_cast<size_t> (0) : rowMap->getNodeNumElements ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (numEntPerRow.extent (0)) != lclNumRows,
      std::invalid_argument, "numEntPerRow has length " <<
      numEntPerRow.extent (0) << " != the local number of rows " <<
      lclNumRows << " as specified by " "the input row Map.");

    if (debug_) {
      for (size_t r = 0; r < lclNumRows; ++r) {
        const size_t curRowCount = numEntPerRow.h_view(r);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (curRowCount == Teuchos::OrdinalTraits<size_t>::invalid (),
           std::invalid_argument, "numEntPerRow(" << r << ") "
           "specifies an invalid number of entries "
           "(Teuchos::OrdinalTraits<size_t>::invalid()).");
      }
    }

    resumeFill (params);
    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
            const Teuchos::RCP<const map_type>& colMap,
            const Teuchos::ArrayView<const size_t>& numEntPerRow,
            const ProfileType /* pftype */,
            const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
    , rowMap_ (rowMap)
    , colMap_ (colMap)
    , numAllocForAllRows_ (0)
  {
    const char tfecfFuncName[] =
      "CrsGraph(rowMap,colMap,numEntPerRow,pftype,params): ";
    staticAssertions ();

    const size_t lclNumRows = rowMap.is_null () ?
      static_cast<size_t> (0) : rowMap->getNodeNumElements ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (numEntPerRow.size ()) != lclNumRows,
      std::invalid_argument, "numEntPerRow has length " << numEntPerRow.size ()
      << " != the local number of rows " << lclNumRows << " as specified by "
      "the input row Map.");

    if (debug_) {
      for (size_t r = 0; r < lclNumRows; ++r) {
        const size_t curRowCount = numEntPerRow[r];
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (curRowCount == Teuchos::OrdinalTraits<size_t>::invalid (),
           std::invalid_argument, "numEntPerRow(" << r << ") "
           "specifies an invalid number of entries "
           "(Teuchos::OrdinalTraits<size_t>::invalid()).");
      }
    }

    // Deep-copy the (host-accessible) input into k_numAllocPerRow_.
    // The latter is a const View, so we have to copy into a nonconst
    // View first, then assign.
    typedef decltype (k_numAllocPerRow_) out_view_type;
    typedef typename out_view_type::non_const_type nc_view_type;
    typedef Kokkos::View<const size_t*,
                         typename nc_view_type::array_layout,
                         Kokkos::HostSpace,
                         Kokkos::MemoryUnmanaged> in_view_type;
    in_view_type numAllocPerRowIn (numEntPerRow.getRawPtr (), lclNumRows);
    nc_view_type numAllocPerRowOut ("Tpetra::CrsGraph::numAllocPerRow",
                                    lclNumRows);
    Kokkos::deep_copy (numAllocPerRowOut, numAllocPerRowIn);
    k_numAllocPerRow_ = numAllocPerRowOut;

    resumeFill (params);
    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
            const Teuchos::RCP<const map_type>& colMap,
            const typename local_graph_type::row_map_type& rowPointers,
            const typename local_graph_type::entries_type::non_const_type& columnIndices,
            const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
    , rowMap_(rowMap)
    , colMap_(colMap)
    , numAllocForAllRows_(0)
    , storageStatus_(Details::STORAGE_1D_PACKED)
    , indicesAreAllocated_(true)
    , indicesAreLocal_(true)
  {
    staticAssertions ();
    if (! params.is_null() && params->isParameter("sorted") &&
        ! params->get<bool>("sorted")) {
      indicesAreSorted_ = false;
    }
    else {
      indicesAreSorted_ = true;
    }
    setAllIndices (rowPointers, columnIndices);
    checkInternalState ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
            const Teuchos::RCP<const map_type>& colMap,
            const Teuchos::ArrayRCP<size_t>& rowPointers,
            const Teuchos::ArrayRCP<LocalOrdinal> & columnIndices,
            const Teuchos::RCP<Teuchos::ParameterList>& params) :
    dist_object_type (rowMap)
    , rowMap_ (rowMap)
    , colMap_ (colMap)
    , numAllocForAllRows_ (0)
    , storageStatus_ (Details::STORAGE_1D_PACKED)
    , indicesAreAllocated_ (true)
    , indicesAreLocal_ (true)
  {
    staticAssertions ();
    if (! params.is_null() && params->isParameter("sorted") &&
        ! params->get<bool>("sorted")) {
      indicesAreSorted_ = false;
    }
    else {
      indicesAreSorted_ = true;
    }
    setAllIndices (rowPointers, columnIndices);
    checkInternalState ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
            const Teuchos::RCP<const map_type>& colMap,
            const local_graph_type& k_local_graph_,
            const Teuchos::RCP<Teuchos::ParameterList>& params)
    : CrsGraph (k_local_graph_,
                rowMap,
                colMap,
                Teuchos::null,
                Teuchos::null,
                params)
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  CrsGraph (const local_graph_type& k_local_graph_,
            const Teuchos::RCP<const map_type>& rowMap,
            const Teuchos::RCP<const map_type>& colMap,
            const Teuchos::RCP<const map_type>& domainMap,
            const Teuchos::RCP<const map_type>& rangeMap,
            const Teuchos::RCP<Teuchos::ParameterList>& params)
    : DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, node_type> (rowMap)
    , rowMap_ (rowMap)
    , colMap_ (colMap)
    , lclGraph_ (k_local_graph_)
    , numAllocForAllRows_ (0)
    , storageStatus_ (Details::STORAGE_1D_PACKED)
    , indicesAreAllocated_ (true)
    , indicesAreLocal_ (true)
  {
    staticAssertions();
    const char tfecfFuncName[] = "CrsGraph(Kokkos::LocalStaticCrsGraph,Map,Map,Map,Map)";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      colMap.is_null (), std::runtime_error,
      ": The input column Map must be nonnull.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      k_local_graph_.numRows () != rowMap->getNodeNumElements (),
      std::runtime_error,
      ": The input row Map and the input local graph need to have the same "
      "number of rows.  The row Map claims " << rowMap->getNodeNumElements ()
      << " row(s), but the local graph claims " << k_local_graph_.numRows ()
      << " row(s).");

    // NOTE (mfh 17 Mar 2014) getNodeNumRows() returns
    // rowMap_->getNodeNumElements(), but it doesn't have to.
    // TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    //   k_local_graph_.numRows () != getNodeNumRows (), std::runtime_error,
    //   ": The input row Map and the input local graph need to have the same "
    //   "number of rows.  The row Map claims " << getNodeNumRows () << " row(s), "
    //   "but the local graph claims " << k_local_graph_.numRows () << " row(s).");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      k_lclInds1D_.extent (0) != 0 || k_gblInds1D_.extent (0) != 0, std::logic_error,
      ": cannot have 1D data structures allocated.");

    if(! params.is_null() && params->isParameter("sorted") &&
       ! params->get<bool>("sorted")) {
      indicesAreSorted_ = false;
    }
    else {
      indicesAreSorted_ = true;
    }

    setDomainRangeMaps (domainMap.is_null() ? rowMap_ : domainMap,
                        rangeMap .is_null() ? rowMap_ : rangeMap);
    Teuchos::Array<int> remotePIDs (0); // unused output argument
    this->makeImportExport (remotePIDs, false);

    k_lclInds1D_ = lclGraph_.entries;
    k_rowPtrs_ = lclGraph_.row_map;

    const bool callComputeGlobalConstants = params.get () == nullptr ||
      params->get ("compute global constants", true);

    if (callComputeGlobalConstants) {
      this->computeGlobalConstants ();
    }
    this->fillComplete_ = true;
    this->checkInternalState ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  CrsGraph (const local_graph_type& lclGraph,
            const Teuchos::RCP<const map_type>& rowMap,
            const Teuchos::RCP<const map_type>& colMap,
            const Teuchos::RCP<const map_type>& domainMap,
            const Teuchos::RCP<const map_type>& rangeMap,
            const Teuchos::RCP<const import_type>& importer,
            const Teuchos::RCP<const export_type>& exporter,
            const Teuchos::RCP<Teuchos::ParameterList>& params) :
    DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, node_type> (rowMap),
    rowMap_ (rowMap),
    colMap_ (colMap),
    rangeMap_ (rangeMap.is_null () ? rowMap : rangeMap),
    domainMap_ (domainMap.is_null () ? rowMap : domainMap),
    importer_ (importer),
    exporter_ (exporter),
    lclGraph_ (lclGraph),
    numAllocForAllRows_ (0),
    storageStatus_ (Details::STORAGE_1D_PACKED),
    indicesAreAllocated_ (true),
    indicesAreLocal_ (true)
  {
    staticAssertions();
    const char tfecfFuncName[] = "Tpetra::CrsGraph(local_graph_type,"
      "Map,Map,Map,Map,Import,Export,params): ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (colMap.is_null (), std::runtime_error,
       "The input column Map must be nonnull.");

    k_lclInds1D_ = lclGraph_.entries;
    k_rowPtrs_ = lclGraph_.row_map;

    if (! params.is_null() && params->isParameter("sorted") &&
        ! params->get<bool>("sorted")) {
      indicesAreSorted_ = false;
    }
    else {
      indicesAreSorted_ = true;
    }

    const bool callComputeGlobalConstants =
      params.get () == nullptr ||
      params->get ("compute global constants", true);
    if (callComputeGlobalConstants) {
      this->computeGlobalConstants ();
    }
    fillComplete_ = true;
    checkInternalState ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Teuchos::ParameterList>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getValidParameters () const
  {
    using Teuchos::RCP;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;

    RCP<ParameterList> params = parameterList ("Tpetra::CrsGraph");

    // Make a sublist for the Import.
    RCP<ParameterList> importSublist = parameterList ("Import");

    // FIXME (mfh 02 Apr 2012) We should really have the Import and
    // Export objects fill in these lists.  However, we don't want to
    // create an Import or Export unless we need them.  For now, we
    // know that the Import and Export just pass the list directly to
    // their Distributor, so we can create a Distributor here
    // (Distributor's constructor is a lightweight operation) and have
    // it fill in the list.

    // Fill in Distributor default parameters by creating a
    // Distributor and asking it to do the work.
    Distributor distributor (rowMap_->getComm (), importSublist);
    params->set ("Import", *importSublist, "How the Import performs communication.");

    // Make a sublist for the Export.  For now, it's a clone of the
    // Import sublist.  It's not a shallow copy, though, since we
    // might like the Import to do communication differently than the
    // Export.
    params->set ("Export", *importSublist, "How the Export performs communication.");

    return params;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    Teuchos::RCP<const Teuchos::ParameterList> validParams =
      getValidParameters ();
    params->validateParametersAndSetDefaults (*validParams);
    this->setMyParamList (params);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  global_size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getGlobalNumRows () const
  {
    return rowMap_->getGlobalNumElements ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  global_size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getGlobalNumCols () const
  {
    const char tfecfFuncName[] = "getGlobalNumCols: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isFillComplete () || getDomainMap ().is_null (), std::runtime_error,
      "The graph does not have a domain Map.  You may not call this method in "
      "that case.");
    return getDomainMap ()->getGlobalNumElements ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getNodeNumRows () const
  {
    return this->rowMap_.is_null () ?
      static_cast<size_t> (0) :
      this->rowMap_->getNodeNumElements ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getNodeNumCols () const
  {
    const char tfecfFuncName[] = "getNodeNumCols: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! hasColMap (), std::runtime_error,
      "The graph does not have a column Map.  You may not call this method "
      "unless the graph has a column Map.  This requires either that a custom "
      "column Map was given to the constructor, or that fillComplete() has "
      "been called.");
    return colMap_.is_null () ? static_cast<size_t> (0) :
      colMap_->getNodeNumElements ();
  }



  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const typename CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::map_type>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getRowMap () const
  {
    return rowMap_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const typename CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::map_type>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getColMap () const
  {
    return colMap_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const typename CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::map_type>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getDomainMap () const
  {
    return domainMap_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const typename CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::map_type>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getRangeMap () const
  {
    return rangeMap_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const typename CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::import_type>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getImporter () const
  {
    return importer_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const typename CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::export_type>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getExporter () const
  {
    return exporter_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  hasColMap () const
  {
    return ! colMap_.is_null ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  isStorageOptimized () const
  {
    // FIXME (mfh 07 Aug 2014) Why wouldn't storage be optimized if
    // getNodeNumRows() is zero?

    const bool isOpt = indicesAreAllocated_ &&
      k_numRowEntries_.extent (0) == 0 &&
      getNodeNumRows () > 0;

    return isOpt;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ProfileType
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getProfileType () const
  {
    return StaticProfile;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  global_size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getGlobalNumEntries () const
  {
    const char tfecfFuncName[] = "getGlobalNumEntries: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! this->haveGlobalConstants_, std::logic_error,
       "The graph does not have global constants computed, "
       "but the user has requested them.");

    return globalNumEntries_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getNodeNumEntries () const
  {
    typedef LocalOrdinal LO;

    if (this->indicesAreAllocated_) {
      const LO lclNumRows = this->getNodeNumRows ();
      if (lclNumRows == 0) {
        return static_cast<size_t> (0);
      }
      else {
        // Avoid the "*this capture" issue by creating a local Kokkos::View.
        auto numEntPerRow = this->k_numRowEntries_;
        const LO numNumEntPerRow = numEntPerRow.extent (0);
        if (numNumEntPerRow == 0) {
          if (static_cast<LO> (this->lclGraph_.row_map.extent (0)) <
              static_cast<LO> (lclNumRows + 1)) {
            return static_cast<size_t> (0);
          }
          else {
            return ::Tpetra::Details::getEntryOnHost (this->lclGraph_.row_map, lclNumRows);
          }
        }
        else { // k_numRowEntries_ is populated
          // k_numRowEntries_ is actually be a host View, so we run
          // the sum in its native execution space.  This also means
          // that we can use explicit capture (which could perhaps
          // improve build time) instead of KOKKOS_LAMBDA, and avoid
          // any CUDA build issues with trying to run a __device__ -
          // only function on host.
          typedef typename num_row_entries_type::execution_space
            host_exec_space;
          typedef Kokkos::RangePolicy<host_exec_space, LO> range_type;

          const LO upperLoopBound = lclNumRows < numNumEntPerRow ?
            lclNumRows :
            numNumEntPerRow;
          size_t nodeNumEnt = 0;
          Kokkos::parallel_reduce ("Tpetra::CrsGraph::getNumNodeEntries",
                                   range_type (0, upperLoopBound),
                                   [=] (const LO& k, size_t& lclSum) {
                                     lclSum += numEntPerRow(k);
                                   }, nodeNumEnt);
          return nodeNumEnt;
        }
      }
    }
    else { // nothing allocated on this process, so no entries
      return static_cast<size_t> (0);
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  global_size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getGlobalMaxNumRowEntries () const
  {
    const char tfecfFuncName[] = "getGlobalMaxNumRowEntries: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! this->haveGlobalConstants_, std::logic_error,
       "The graph does not have global constants computed, "
       "but the user has requested them.");

    return globalMaxNumRowEntries_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getNodeMaxNumRowEntries () const
  {
    return nodeMaxNumRowEntries_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  isFillComplete () const
  {
    return fillComplete_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  isFillActive () const
  {
    return ! fillComplete_;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  isLocallyIndexed () const
  {
    return indicesAreLocal_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  isGloballyIndexed () const
  {
    return indicesAreGlobal_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getNodeAllocationSize () const
  {
    typedef LocalOrdinal LO;

    if (this->indicesAreAllocated_) {
      const LO lclNumRows = this->getNodeNumRows ();
      if (lclNumRows == 0) {
        return static_cast<size_t> (0);
      }
      else if (storageStatus_ == Details::STORAGE_1D_PACKED) {
        if (static_cast<LO> (this->lclGraph_.row_map.extent (0)) <
            static_cast<LO> (lclNumRows + 1)) {
          return static_cast<size_t> (0);
        }
        else {
          return ::Tpetra::Details::getEntryOnHost (this->lclGraph_.row_map, lclNumRows);
        }
      }
      else if (storageStatus_ == Details::STORAGE_1D_UNPACKED) {
        if (this->k_rowPtrs_.extent (0) == 0) {
          return static_cast<size_t> (0);
        }
        else {
          return ::Tpetra::Details::getEntryOnHost (this->k_rowPtrs_, lclNumRows);
        }
      }
      else {
        return static_cast<size_t> (0);
      }
    }
    else {
      return Tpetra::Details::OrdinalTraits<size_t>::invalid ();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Teuchos::Comm<int> >
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getComm () const
  {
    return this->rowMap_.is_null () ? Teuchos::null : this->rowMap_->getComm ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  GlobalOrdinal
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getIndexBase () const
  {
    return rowMap_->getIndexBase ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  indicesAreAllocated () const
  {
    return indicesAreAllocated_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  isSorted () const
  {
    return indicesAreSorted_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  isMerged () const
  {
    return noRedundancies_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  setLocallyModified ()
  {
    // FIXME (mfh 07 May 2013) How do we know that the change
    // introduced a redundancy, or even that it invalidated the sorted
    // order of indices?  CrsGraph has always made this conservative
    // guess.  It could be a bit costly to check at insertion time,
    // though.
    indicesAreSorted_ = false;
    noRedundancies_ = false;

    // We've modified the graph, so we'll have to recompute local
    // constants like the number of diagonal entries on this process.
    haveLocalConstants_ = false;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  allocateIndices (const ELocalGlobal lg, const bool verbose)
  {
    using Details::ProfilingRegion;
    using Teuchos::arcp;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using std::endl;
    typedef Teuchos::ArrayRCP<size_t>::size_type size_type;
    typedef typename local_graph_type::row_map_type::non_const_type
      non_const_row_map_type;
    typedef typename local_graph_type::entries_type::non_const_type
      lcl_col_inds_type;
    typedef Kokkos::View<GlobalOrdinal*,
      typename lcl_col_inds_type::array_layout,
      device_type> gbl_col_inds_type;
    const char tfecfFuncName[] = "allocateIndices: ";
    const char suffix[] =
      "  Please report this bug to the Tpetra developers.";
    ProfilingRegion profRegion("Tpetra::CrsGraph::allocateIndices");

    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsGraph", tfecfFuncName);
      std::ostringstream os;
      os << *prefix << "Start: lg="
         << (lg == GlobalIndices ? "GlobalIndices" : "LocalIndices")
         << ", numRows: " << this->getNodeNumRows() << endl;
      std::cerr << os.str();
    }

    // This is a protected function, only callable by us.  If it was
    // called incorrectly, it is our fault.  That's why the tests
    // below throw std::logic_error instead of std::invalid_argument.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (isLocallyIndexed () && lg == GlobalIndices, std::logic_error,
       ": The graph is locally indexed, but Tpetra code is calling "
       "this method with lg=GlobalIndices." << suffix);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (isGloballyIndexed () && lg == LocalIndices, std::logic_error,
       ": The graph is globally indexed, but Tpetra code is calling "
       "this method with lg=LocalIndices." << suffix);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (indicesAreAllocated (), std::logic_error, ": The graph's "
       "indices are already allocated, but Tpetra is calling "
       "allocateIndices again." << suffix);
    const size_t numRows = this->getNodeNumRows ();

    //
    //  STATIC ALLOCATION PROFILE
    //
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Allocate k_rowPtrs: " << (numRows+1) << endl;
      std::cerr << os.str();
    }
    non_const_row_map_type k_rowPtrs ("Tpetra::CrsGraph::ptr", numRows + 1);

    if (this->k_numAllocPerRow_.extent (0) != 0) {
      // It's OK to throw std::invalid_argument here, because we
      // haven't incurred any side effects yet.  Throwing that
      // exception (and not, say, std::logic_error) implies that the
      // instance can recover.
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (this->k_numAllocPerRow_.extent (0) != numRows,
         std::invalid_argument, "k_numAllocPerRow_ is allocated, that is, "
         "has nonzero length " << this->k_numAllocPerRow_.extent (0)
         << ", but its length != numRows = " << numRows << ".");

      // k_numAllocPerRow_ is a host View, but k_rowPtrs (the thing
      // we want to compute here) lives on device.  That's OK;
      // computeOffsetsFromCounts can handle this case.
      using Details::computeOffsetsFromCounts;

      // FIXME (mfh 27 Jun 2016) Currently, computeOffsetsFromCounts
      // doesn't attempt to check its input for "invalid" flag
      // values.  For now, we omit that feature of the sequential
      // code disabled below.
      computeOffsetsFromCounts (k_rowPtrs, k_numAllocPerRow_);
    }
    else {
      // It's OK to throw std::invalid_argument here, because we
      // haven't incurred any side effects yet.  Throwing that
      // exception (and not, say, std::logic_error) implies that the
      // instance can recover.
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (this->numAllocForAllRows_ ==
         Tpetra::Details::OrdinalTraits<size_t>::invalid (),
         std::invalid_argument, "numAllocForAllRows_ has an invalid value, "
         "namely Tpetra::Details::OrdinalTraits<size_t>::invalid() = " <<
         Tpetra::Details::OrdinalTraits<size_t>::invalid () << ".");

      using Details::computeOffsetsFromConstantCount;
      computeOffsetsFromConstantCount (k_rowPtrs, this->numAllocForAllRows_);
    }

    // "Commit" the resulting row offsets.
    this->k_rowPtrs_ = k_rowPtrs;

    const size_type numInds =
      Details::getEntryOnHost(this->k_rowPtrs_, numRows);
    if (lg == LocalIndices) {
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Allocate local column indices "
          "k_lclInds1D_: " << numInds << endl;
        std::cerr << os.str();
      }
      k_lclInds1D_ = lcl_col_inds_type ("Tpetra::CrsGraph::ind", numInds);
    }
    else {
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Allocate global column indices "
          "k_gblInds1D_: " << numInds << endl;
        std::cerr << os.str();
      }
      k_gblInds1D_ = gbl_col_inds_type ("Tpetra::CrsGraph::ind", numInds);
    }
    storageStatus_ = Details::STORAGE_1D_UNPACKED;

    this->indicesAreLocal_  = (lg == LocalIndices);
    this->indicesAreGlobal_ = (lg == GlobalIndices);

    if (numRows > 0) { // reallocate k_numRowEntries_ & fill w/ 0s
      using Kokkos::ViewAllocateWithoutInitializing;
      typedef decltype (k_numRowEntries_) row_ent_type;
      const char label[] = "Tpetra::CrsGraph::numRowEntries";
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Allocate k_numRowEntries_: " << numRows
           << endl;
        std::cerr << os.str();
      }
      row_ent_type numRowEnt (ViewAllocateWithoutInitializing (label), numRows);
      Kokkos::deep_copy (numRowEnt, static_cast<size_t> (0)); // fill w/ 0s
      this->k_numRowEntries_ = numRowEnt; // "commit" our allocation
    }

    // Once indices are allocated, CrsGraph needs to free this information.
    this->numAllocForAllRows_ = 0;
    this->k_numAllocPerRow_ = decltype (k_numAllocPerRow_) ();
    this->indicesAreAllocated_ = true;

    try {
      this->checkInternalState ();
    }
    catch (std::logic_error& e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::logic_error, "At end of allocateIndices, "
         "checkInternalState threw std::logic_error: "
         << e.what ());
    }
    catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::runtime_error, "At end of allocateIndices, "
         "checkInternalState threw std::exception: "
         << e.what ());
    }
    catch (...) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::runtime_error, "At end of allocateIndices, "
         "checkInternalState threw an exception "
         "not a subclass of std::exception.");
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayView<const LocalOrdinal>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getLocalView (const RowInfo& rowinfo) const
  {
    using Kokkos::subview;
    typedef LocalOrdinal LO;
    typedef Kokkos::View<const LO*, execution_space,
      Kokkos::MemoryUnmanaged> row_view_type;

    if (rowinfo.allocSize == 0) {
      return Teuchos::ArrayView<const LO> ();
    }
    else { // nothing in the row to view
      if (k_lclInds1D_.extent (0) != 0) { // 1-D storage
        const size_t start = rowinfo.offset1D;
        const size_t len = rowinfo.allocSize;
        const std::pair<size_t, size_t> rng (start, start + len);
        // mfh 23 Nov 2015: Don't just create a subview of
        // k_lclInds1D_ directly, because that first creates a
        // _managed_ subview, then returns an unmanaged version of
        // that.  That touches the reference count, which costs
        // performance in a measurable way.
        row_view_type rowView = subview (row_view_type (k_lclInds1D_), rng);
        const LO* const rowViewRaw = (len == 0) ? nullptr : rowView.data ();
        return Teuchos::ArrayView<const LO> (rowViewRaw, len, Teuchos::RCP_DISABLE_NODE_LOOKUP);
      }
      else {
        return Teuchos::ArrayView<const LO> (); // nothing in the row to view
      }
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getLocalViewRawConst (const LocalOrdinal*& lclInds,
                        LocalOrdinal& capacity,
                        const RowInfo& rowInfo) const
  {
    lclInds = nullptr;
    capacity = 0;

    if (rowInfo.allocSize != 0 && k_lclInds1D_.extent (0) != 0) {
      if (debug_) {
        if (rowInfo.offset1D + rowInfo.allocSize >
            static_cast<size_t> (k_lclInds1D_.extent (0))) {
          return static_cast<LocalOrdinal> (-1);
        }
      }
      lclInds = k_lclInds1D_.data () + rowInfo.offset1D;
      capacity = rowInfo.allocSize;
    }
    return static_cast<LocalOrdinal> (0);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayView<LocalOrdinal>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getLocalViewNonConst (const RowInfo& rowinfo)
  {
    using Kokkos::subview;
    typedef LocalOrdinal LO;
    typedef Kokkos::View<LO*, execution_space,
      Kokkos::MemoryUnmanaged> row_view_type;

    if (rowinfo.allocSize == 0) { // nothing in the row to view
      return Teuchos::ArrayView<LO> ();
    }
    else {
      if (k_lclInds1D_.extent (0) != 0) { // 1-D storage
        const size_t start = rowinfo.offset1D;
        const size_t len = rowinfo.allocSize;
        const std::pair<size_t, size_t> rng (start, start + len);
        // mfh 23 Nov 2015: Don't just create a subview of
        // k_lclInds1D_ directly, because that first creates a
        // _managed_ subview, then returns an unmanaged version of
        // that.  That touches the reference count, which costs
        // performance in a measurable way.
        row_view_type rowView = subview (row_view_type (k_lclInds1D_), rng);
        LO* const rowViewRaw = (len == 0) ? nullptr : rowView.data ();
        return Teuchos::ArrayView<LO> (rowViewRaw, len, Teuchos::RCP_DISABLE_NODE_LOOKUP);
      }
      else {
        return Teuchos::ArrayView<LO> (); // nothing in the row to view
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Kokkos::View<const LocalOrdinal*,
               typename CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::execution_space,
               Kokkos::MemoryUnmanaged>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getLocalKokkosRowView (const RowInfo& rowInfo) const
  {
    typedef LocalOrdinal LO;
    typedef Kokkos::View<const LO*, execution_space,
      Kokkos::MemoryUnmanaged> row_view_type;

    if (rowInfo.allocSize == 0) {
      return row_view_type ();
    }
    else { // nothing in the row to view
      if (k_lclInds1D_.extent (0) != 0) { // 1-D storage
        const size_t start = rowInfo.offset1D;
        const size_t len = rowInfo.allocSize;
        const std::pair<size_t, size_t> rng (start, start + len);
        // mfh 23 Nov 2015: Don't just create a subview of
        // k_lclInds1D_ directly, because that first creates a
        // _managed_ subview, then returns an unmanaged version of
        // that.  That touches the reference count, which costs
        // performance in a measurable way.
        return Kokkos::subview (row_view_type (k_lclInds1D_), rng);
      }
      else {
        return row_view_type (); // nothing in the row to view
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Kokkos::View<LocalOrdinal*,
               typename CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::execution_space,
               Kokkos::MemoryUnmanaged>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getLocalKokkosRowViewNonConst (const RowInfo& rowInfo)
  {
    using row_view_type = Kokkos::View<LocalOrdinal*,
      execution_space, Kokkos::MemoryUnmanaged>;

    if (rowInfo.allocSize == 0) {
      return row_view_type ();
    }
    else { // nothing in the row to view
      if (k_lclInds1D_.extent (0) != 0) { // 1-D storage
        const size_t start = rowInfo.offset1D;
        const size_t len = rowInfo.allocSize;
        const std::pair<size_t, size_t> rng (start, start + len);
        // mfh 23 Nov 2015: Don't just create a subview of
        // k_lclInds1D_ directly, because that first creates a
        // _managed_ subview, then returns an unmanaged version of
        // that.  That touches the reference count, which costs
        // performance in a measurable way.
        return Kokkos::subview (row_view_type (this->k_lclInds1D_), rng);
      }
      else {
        return row_view_type (); // nothing in the row to view
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Kokkos::View<const GlobalOrdinal*,
               typename CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::execution_space,
               Kokkos::MemoryUnmanaged>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getGlobalKokkosRowView (const RowInfo& rowinfo) const
  {
    using row_view_type = Kokkos::View<const GlobalOrdinal*,
      execution_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    if (rowinfo.allocSize == 0) {
      return row_view_type ();
    }
    else { // nothing in the row to view
      if (this->k_gblInds1D_.extent (0) != 0) { // 1-D storage
        const size_t start = rowinfo.offset1D;
        const size_t len = rowinfo.allocSize;
        const std::pair<size_t, size_t> rng (start, start + len);
        // mfh 23 Nov 2015: Don't just create a subview of
        // k_gblInds1D_ directly, because that first creates a
        // _managed_ subview, then returns an unmanaged version of
        // that.  That touches the reference count, which costs
        // performance in a measurable way.
        return Kokkos::subview (row_view_type (this->k_gblInds1D_), rng);
      }
      else {
        return row_view_type (); // nothing in the row to view
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayView<const GlobalOrdinal>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getGlobalView (const RowInfo& rowinfo) const
  {
    using GO = global_ordinal_type;

    Teuchos::ArrayView<const GO> view;
    if (rowinfo.allocSize > 0 && k_gblInds1D_.extent (0) != 0) {
      const auto rng =
        std::make_pair (rowinfo.offset1D,
                        rowinfo.offset1D + rowinfo.allocSize);
      // mfh 23 Nov 2015: Don't just create a subview of
      // k_gblInds1D_ directly, because that first creates a
      // _managed_ subview, then returns an unmanaged version of
      // that.  That touches the reference count, which costs
      // performance in a measurable way.
      using row_view_type = Kokkos::View<const GO*,
        execution_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
      row_view_type k_gblInds1D_unmanaged = k_gblInds1D_;
      using Kokkos::Compat::getConstArrayView;
      using Kokkos::subview;
      view = getConstArrayView (subview (k_gblInds1D_unmanaged, rng));
    }
    return view;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getGlobalViewRawConst (const GlobalOrdinal*& gblInds,
                         LocalOrdinal& capacity,
                         const RowInfo& rowInfo) const
  {
    gblInds = nullptr;
    capacity = 0;

    if (rowInfo.allocSize != 0 && k_gblInds1D_.extent (0) != 0) {
      if (debug_) {
        if (rowInfo.offset1D + rowInfo.allocSize >
            static_cast<size_t> (k_gblInds1D_.extent (0))) {
          return static_cast<LocalOrdinal> (-1);
        }
      }
      gblInds = k_gblInds1D_.data () + rowInfo.offset1D;
      capacity = rowInfo.allocSize;
    }
    return static_cast<LocalOrdinal> (0);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayView<GlobalOrdinal>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getGlobalViewNonConst (const RowInfo& rowinfo)
  {
    using GO = global_ordinal_type;

    Teuchos::ArrayView<GO> view;
    if (rowinfo.allocSize > 0 && k_gblInds1D_.extent (0) != 0) {
      const auto rng =
        std::make_pair (rowinfo.offset1D,
                        rowinfo.offset1D + rowinfo.allocSize);
      // mfh 23 Nov 2015: Don't just create a subview of
      // k_gblInds1D_ directly, because that first creates a
      // _managed_ subview, then returns an unmanaged version of
      // that.  That touches the reference count, which costs
      // performance in a measurable way.
      using row_view_type = Kokkos::View<GO*, execution_space,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
      row_view_type k_gblInds1D_unmanaged = k_gblInds1D_;
      using Kokkos::Compat::getArrayView;
      using Kokkos::subview;
      view = getArrayView (subview (k_gblInds1D_unmanaged, rng));
    }
    return view;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  RowInfo
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getRowInfo (const LocalOrdinal myRow) const
  {
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid ();
    RowInfo ret;
    if (this->rowMap_.is_null () || ! this->rowMap_->isNodeLocalElement (myRow)) {
      ret.localRow = STINV;
      ret.allocSize = 0;
      ret.numEntries = 0;
      ret.offset1D = STINV;
      return ret;
    }

    ret.localRow = static_cast<size_t> (myRow);
    if (this->indicesAreAllocated ()) {
      // Offsets tell us the allocation size in this case.
      if (this->k_rowPtrs_.extent (0) == 0) {
        ret.offset1D  = 0;
        ret.allocSize = 0;
      }
      else {
        ret.offset1D  = this->k_rowPtrs_(myRow);
        ret.allocSize = this->k_rowPtrs_(myRow+1) - this->k_rowPtrs_(myRow);
      }

      ret.numEntries = (this->k_numRowEntries_.extent (0) == 0) ?
        ret.allocSize :
        this->k_numRowEntries_(myRow);
    }
    else { // haven't performed allocation yet; probably won't hit this code
      // FIXME (mfh 07 Aug 2014) We want graph's constructors to
      // allocate, rather than doing lazy allocation at first insert.
      // This will make k_numAllocPerRow_ obsolete.
      ret.allocSize = (this->k_numAllocPerRow_.extent (0) != 0) ?
        this->k_numAllocPerRow_(myRow) : // this is a host View
        this->numAllocForAllRows_;
      ret.numEntries = 0;
      ret.offset1D = STINV;
    }

    return ret;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  RowInfo
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getRowInfoFromGlobalRowIndex (const GlobalOrdinal gblRow) const
  {
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid ();
    RowInfo ret;
    if (this->rowMap_.is_null ()) {
      ret.localRow = STINV;
      ret.allocSize = 0;
      ret.numEntries = 0;
      ret.offset1D = STINV;
      return ret;
    }
    const LocalOrdinal myRow = this->rowMap_->getLocalElement (gblRow);
    if (myRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid ()) {
      ret.localRow = STINV;
      ret.allocSize = 0;
      ret.numEntries = 0;
      ret.offset1D = STINV;
      return ret;
    }

    ret.localRow = static_cast<size_t> (myRow);
    if (this->indicesAreAllocated ()) {
      // graph data structures have the info that we need
      //
      // if static graph, offsets tell us the allocation size
      if (this->k_rowPtrs_.extent (0) == 0) {
        ret.offset1D  = 0;
        ret.allocSize = 0;
      }
      else {
        ret.offset1D  = this->k_rowPtrs_(myRow);
        ret.allocSize = this->k_rowPtrs_(myRow+1) - this->k_rowPtrs_(myRow);
      }

      ret.numEntries = (this->k_numRowEntries_.extent (0) == 0) ?
        ret.allocSize :
        this->k_numRowEntries_(myRow);
    }
    else { // haven't performed allocation yet; probably won't hit this code
      // FIXME (mfh 07 Aug 2014) We want graph's constructors to
      // allocate, rather than doing lazy allocation at first insert.
      // This will make k_numAllocPerRow_ obsolete.
      ret.allocSize = (this->k_numAllocPerRow_.extent (0) != 0) ?
        this->k_numAllocPerRow_(myRow) : // this is a host View
        this->numAllocForAllRows_;
      ret.numEntries = 0;
      ret.offset1D = STINV;
    }

    return ret;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  staticAssertions () const
  {
    using Teuchos::OrdinalTraits;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef global_size_t GST;

    // Assumption: sizeof(GlobalOrdinal) >= sizeof(LocalOrdinal):
    //     This is so that we can store local indices in the memory
    //     formerly occupied by global indices.
    static_assert (sizeof (GlobalOrdinal) >= sizeof (LocalOrdinal),
                   "Tpetra::CrsGraph: sizeof(GlobalOrdinal) must be >= sizeof(LocalOrdinal).");
    // Assumption: max(size_t) >= max(LocalOrdinal)
    //     This is so that we can represent any LocalOrdinal as a size_t.
    static_assert (sizeof (size_t) >= sizeof (LocalOrdinal),
                   "Tpetra::CrsGraph: sizeof(size_t) must be >= sizeof(LocalOrdinal).");
    static_assert (sizeof(GST) >= sizeof(size_t),
                   "Tpetra::CrsGraph: sizeof(Tpetra::global_size_t) must be >= sizeof(size_t).");

    // FIXME (mfh 30 Sep 2015) We're not using
    // Teuchos::CompileTimeAssert any more.  Can we do these checks
    // with static_assert?

    // can't call max() with CompileTimeAssert, because it isn't a
    // constant expression; will need to make this a runtime check
    const char msg[] = "Tpetra::CrsGraph: Object cannot be created with the "
      "given template arguments: size assumptions are not valid.";
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t> (Teuchos::OrdinalTraits<LO>::max ()) > Teuchos::OrdinalTraits<size_t>::max (),
      std::runtime_error, msg);
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<GST> (Teuchos::OrdinalTraits<LO>::max ()) > static_cast<GST> (Teuchos::OrdinalTraits<GO>::max ()),
      std::runtime_error, msg);
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t> (Teuchos::OrdinalTraits<GO>::max ()) > Teuchos::OrdinalTraits<GST>::max(),
      std::runtime_error, msg);
    TEUCHOS_TEST_FOR_EXCEPTION(
      Teuchos::OrdinalTraits<size_t>::max () > Teuchos::OrdinalTraits<GST>::max (),
      std::runtime_error, msg);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  insertIndices (RowInfo& rowinfo,
                 const SLocalGlobalViews &newInds,
                 const ELocalGlobal lg,
                 const ELocalGlobal I)
  {
    using Teuchos::ArrayView;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    const char tfecfFuncName[] = "insertIndices: ";

    size_t oldNumEnt = 0;
    if (debug_) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (lg != GlobalIndices && lg != LocalIndices, std::invalid_argument,
         "lg must be either GlobalIndices or LocalIndices.");
      oldNumEnt = this->getNumEntriesInLocalRow (rowinfo.localRow);
    }

    size_t numNewInds = 0;
    if (lg == GlobalIndices) { // input indices are global
      ArrayView<const GO> new_ginds = newInds.ginds;
      numNewInds = new_ginds.size();
      if (I == GlobalIndices) { // store global indices
        ArrayView<GO> gind_view = this->getGlobalViewNonConst (rowinfo);
        if (debug_) {
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (static_cast<size_t> (gind_view.size ()) <
             rowinfo.numEntries + numNewInds, std::logic_error,
             "gind_view.size() = " << gind_view.size ()
             << " < rowinfo.numEntries (= " << rowinfo.numEntries
             << ") + numNewInds (= " << numNewInds << ").");
        }
        GO* const gblColInds_out = gind_view.getRawPtr () + rowinfo.numEntries;
        for (size_t k = 0; k < numNewInds; ++k) {
          gblColInds_out[k] = new_ginds[k];
        }
      }
      else if (I == LocalIndices) { // store local indices
        ArrayView<LO> lind_view = this->getLocalViewNonConst (rowinfo);
        if (debug_) {
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (static_cast<size_t> (lind_view.size ()) <
             rowinfo.numEntries + numNewInds, std::logic_error,
             "lind_view.size() = " << lind_view.size ()
             << " < rowinfo.numEntries (= " << rowinfo.numEntries
             << ") + numNewInds (= " << numNewInds << ").");
        }
        LO* const lclColInds_out = lind_view.getRawPtr () + rowinfo.numEntries;
        for (size_t k = 0; k < numNewInds; ++k) {
          lclColInds_out[k] = colMap_->getLocalElement (new_ginds[k]);
        }
      }
    }
    else if (lg == LocalIndices) { // input indices are local
      ArrayView<const LO> new_linds = newInds.linds;
      numNewInds = new_linds.size();
      if (I == LocalIndices) { // store local indices
        ArrayView<LO> lind_view = this->getLocalViewNonConst (rowinfo);
        if (debug_) {
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (static_cast<size_t> (lind_view.size ()) <
             rowinfo.numEntries + numNewInds, std::logic_error,
             "lind_view.size() = " << lind_view.size ()
             << " < rowinfo.numEntries (= " << rowinfo.numEntries
             << ") + numNewInds (= " << numNewInds << ").");
        }
        LO* const lclColInds_out = lind_view.getRawPtr () + rowinfo.numEntries;
        for (size_t k = 0; k < numNewInds; ++k) {
          lclColInds_out[k] = new_linds[k];
        }
      }
      else if (I == GlobalIndices) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (true, std::logic_error, "The case where the input indices are local "
           "and the indices to write are global (lg=LocalIndices, I="
           "GlobalIndices) is not implemented, because it does not make sense."
           << std::endl << "If you have correct local column indices, that "
           "means the graph has a column Map.  In that case, you should be "
           "storing local indices.");
      }
    }

    rowinfo.numEntries += numNewInds;
    this->k_numRowEntries_(rowinfo.localRow) += numNewInds;
    this->setLocallyModified ();

    if (debug_) {
      const size_t chkNewNumEnt =
        this->getNumEntriesInLocalRow (rowinfo.localRow);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (chkNewNumEnt != oldNumEnt + numNewInds, std::logic_error,
         "chkNewNumEnt = " << chkNewNumEnt
         << " != oldNumEnt (= " << oldNumEnt
         << ") + numNewInds (= " << numNewInds << ").");
    }

    return numNewInds;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  insertGlobalIndicesImpl (const LocalOrdinal lclRow,
                           const GlobalOrdinal inputGblColInds[],
                           const size_t numInputInds)
  {
    return this->insertGlobalIndicesImpl (this->getRowInfo (lclRow),
                                          inputGblColInds, numInputInds);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  insertGlobalIndicesImpl (const RowInfo& rowInfo,
                           const GlobalOrdinal inputGblColInds[],
                           const size_t numInputInds,
                           std::function<void(const size_t, const size_t, const size_t)> fun)
  {
    using Details::verbosePrintArray;
    using Kokkos::View;
    using Kokkos::subview;
    using Kokkos::MemoryUnmanaged;
    using Teuchos::ArrayView;
    using LO = LocalOrdinal;
    using GO = GlobalOrdinal;
    const char tfecfFuncName[] = "insertGlobalIndicesImpl: ";
    const LO lclRow = static_cast<LO> (rowInfo.localRow);

    auto numEntries = rowInfo.numEntries;
    using inp_view_type = View<const GO*, execution_space, MemoryUnmanaged>;
    inp_view_type inputInds(inputGblColInds, numInputInds);
    size_t numInserted = Details::insertCrsIndices(lclRow, k_rowPtrs_,
      this->k_gblInds1D_, numEntries, inputInds, fun);

    const bool insertFailed =
      numInserted == Teuchos::OrdinalTraits<size_t>::invalid();
    if(insertFailed) {
      constexpr size_t ONE (1);
      const int myRank = this->getComm()->getRank();
      std::ostringstream os;

      os << "Proc " << myRank << ": Not enough capacity to insert "
         << numInputInds
         << " ind" << (numInputInds != ONE ? "ices" : "ex")
         << " into local row " << lclRow << ", which currently has "
         << rowInfo.numEntries
         << " entr" << (rowInfo.numEntries != ONE ? "ies" : "y")
         << " and total allocation size " << rowInfo.allocSize
         << ".  ";
      const size_t maxNumToPrint =
        Details::Behavior::verbosePrintCountThreshold();
      ArrayView<const GO> inputGblColIndsView(inputGblColInds,
                                              numInputInds);
      verbosePrintArray(os, inputGblColIndsView, "Input global "
                        "column indices", maxNumToPrint);
      os << ", ";
      const GO* const curGblColInds =
        k_gblInds1D_.data() + rowInfo.offset1D;
      ArrayView<const GO> curGblColIndsView(curGblColInds,
                                            rowInfo.numEntries);
      verbosePrintArray(os, curGblColIndsView, "Current global "
                        "column indices", maxNumToPrint);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::runtime_error, os.str());
    }

    this->k_numRowEntries_(lclRow) += numInserted;
    this->setLocallyModified();
    return numInserted;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  insertLocalIndicesImpl (const LocalOrdinal myRow,
                          const Teuchos::ArrayView<const LocalOrdinal>& indices,
                          std::function<void(const size_t, const size_t, const size_t)> fun)
  {
    using Kokkos::MemoryUnmanaged;
    using Kokkos::subview;
    using Kokkos::View;
    using LO = LocalOrdinal;
    const char tfecfFuncName[] = "insertLocallIndicesImpl: ";

    const RowInfo rowInfo = this->getRowInfo(myRow);

    size_t numNewInds = 0;
    size_t newNumEntries = 0;

    auto numEntries = rowInfo.numEntries;
    // Note: Teuchos::ArrayViews are in HostSpace
    using inp_view_type = View<const LO*, Kokkos::HostSpace, MemoryUnmanaged>;
    inp_view_type inputInds(indices.getRawPtr(), indices.size());
    auto numInserted = Details::insertCrsIndices(myRow, k_rowPtrs_,
      this->k_lclInds1D_, numEntries, inputInds, fun);

    const bool insertFailed =
      numInserted == Teuchos::OrdinalTraits<size_t>::invalid();
    if(insertFailed) {
      constexpr size_t ONE (1);
      const size_t numInputInds(indices.size());
      const int myRank = this->getComm()->getRank();
      std::ostringstream os;
      os << "On MPI Process " << myRank << ": Not enough capacity to "
        "insert " << numInputInds
         << " ind" << (numInputInds != ONE ? "ices" : "ex")
         << " into local row " << myRow << ", which currently has "
         << rowInfo.numEntries
         << " entr" << (rowInfo.numEntries != ONE ? "ies" : "y")
         << " and total allocation size " << rowInfo.allocSize << ".";
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::runtime_error, os.str());
    }
    numNewInds = numInserted;
    newNumEntries = rowInfo.numEntries + numNewInds;

    this->k_numRowEntries_(myRow) += numNewInds;
    this->setLocallyModified ();

    if (debug_) {
      const size_t chkNewNumEntries = this->getNumEntriesInLocalRow (myRow);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (chkNewNumEntries != newNumEntries, std::logic_error,
         "getNumEntriesInLocalRow(" << myRow << ") = " << chkNewNumEntries
         << " != newNumEntries = " << newNumEntries
         << ".  Please report this bug to the Tpetra developers.");
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  findLocalIndices(const RowInfo& rowInfo,
                   const Teuchos::ArrayView<const LocalOrdinal>& indices,
                   std::function<void(const size_t, const size_t, const size_t)> fun) const
  {
    using LO = LocalOrdinal;
    using inp_view_type = Kokkos::View<const LO*, Kokkos::HostSpace,
      Kokkos::MemoryUnmanaged>;
    inp_view_type inputInds(indices.getRawPtr(), indices.size());

    size_t numFound = 0;
    LO lclRow = rowInfo.localRow;
    if (this->isLocallyIndexed())
    {
      numFound = Details::findCrsIndices(lclRow, k_rowPtrs_, rowInfo.numEntries,
        this->k_lclInds1D_, inputInds, fun);
    }
    else if (this->isGloballyIndexed())
    {
      if (this->colMap_.is_null())
        return Teuchos::OrdinalTraits<size_t>::invalid();
      const auto& colMap = *(this->colMap_);
      auto map = [&](LO const lclInd){return colMap.getGlobalElement(lclInd);};
      numFound = Details::findCrsIndices(lclRow, k_rowPtrs_, rowInfo.numEntries,
        this->k_gblInds1D_, inputInds, map, fun);
    }
    return numFound;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  findGlobalIndices(const RowInfo& rowInfo,
                    const Teuchos::ArrayView<const GlobalOrdinal>& indices,
                    std::function<void(const size_t, const size_t, const size_t)> fun) const
  {
    using GO = GlobalOrdinal;
    using Kokkos::View;
    using Kokkos::MemoryUnmanaged;
    auto invalidCount = Teuchos::OrdinalTraits<size_t>::invalid();

    using inp_view_type = View<const GO*, execution_space, MemoryUnmanaged>;
    inp_view_type inputInds(indices.getRawPtr(), indices.size());

    size_t numFound = 0;
    LocalOrdinal lclRow = rowInfo.localRow;
    if (this->isLocallyIndexed())
    {
      if (this->colMap_.is_null())
        return invalidCount;
      const auto& colMap = *(this->colMap_);
      auto map = [&](GO const gblInd){return colMap.getLocalElement(gblInd);};
      numFound = Details::findCrsIndices(lclRow, k_rowPtrs_, rowInfo.numEntries,
        this->k_lclInds1D_, inputInds, map, fun);
    }
    else if (this->isGloballyIndexed())
    {
      numFound = Details::findCrsIndices(lclRow, k_rowPtrs_, rowInfo.numEntries,
        this->k_gblInds1D_, inputInds, fun);
    }
    return numFound;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  sortAndMergeRowIndices (const RowInfo& rowInfo,
                          const bool sorted,
                          const bool merged)
  {
    const size_t origNumEnt = rowInfo.numEntries;
    if (origNumEnt != Tpetra::Details::OrdinalTraits<size_t>::invalid () &&
        origNumEnt != 0) {
      auto lclColInds = this->getLocalKokkosRowViewNonConst (rowInfo);

      LocalOrdinal* const lclColIndsRaw = lclColInds.data ();
      if (! sorted) {
        // FIXME (mfh 08 May 2017) This assumes CUDA UVM.
        std::sort (lclColIndsRaw, lclColIndsRaw + origNumEnt);
      }

      if (! merged) {
        LocalOrdinal* const beg = lclColIndsRaw;
        LocalOrdinal* const end = beg + rowInfo.numEntries;
        // FIXME (mfh 08 May 2017) This assumes CUDA UVM.
        LocalOrdinal* const newend = std::unique (beg, end);
        const size_t newNumEnt = newend - beg;

        // NOTE (mfh 08 May 2017) This is a host View, so it does not assume UVM.
        this->k_numRowEntries_(rowInfo.localRow) = newNumEnt;
        return origNumEnt - newNumEnt; // the number of duplicates in the row
      }
      else {
        return static_cast<size_t> (0); // assume no duplicates
      }
    }
    else {
      return static_cast<size_t> (0); // no entries in the row
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  setDomainRangeMaps (const Teuchos::RCP<const map_type>& domainMap,
                      const Teuchos::RCP<const map_type>& rangeMap)
  {
    // simple pointer comparison for equality
    if (domainMap_ != domainMap) {
      domainMap_ = domainMap;
      importer_ = Teuchos::null;
    }
    if (rangeMap_ != rangeMap) {
      rangeMap_  = rangeMap;
      exporter_ = Teuchos::null;
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  clearGlobalConstants ()
  {
    const auto INV = Teuchos::OrdinalTraits<global_size_t>::invalid();

    globalNumEntries_       = INV;
    globalMaxNumRowEntries_ = INV;
    haveGlobalConstants_    = false;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  checkInternalState () const
  {
    if (debug_) {
      using std::endl;
      const char tfecfFuncName[] = "checkInternalState: ";
      const char suffix[] = "  Please report this bug to the Tpetra developers.";

      std::unique_ptr<std::string> prefix;
      if (verbose_) {
        prefix = this->createPrefix("CrsGraph", "checkInternalState");
        std::ostringstream os;
        os << *prefix << "Start" << endl;
        std::cerr << os.str();
      }

      const global_size_t GSTI = Teuchos::OrdinalTraits<global_size_t>::invalid ();
      //const size_t         STI = Teuchos::OrdinalTraits<size_t>::invalid (); // unused
      // check the internal state of this data structure
      // this is called by numerous state-changing methods, in a debug build, to ensure that the object
      // always remains in a valid state

      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (this->rowMap_.is_null (), std::logic_error,
         "Row Map is null." << suffix);
      // This may access the row Map, so we need to check first (above)
      // whether the row Map is null.
      const LocalOrdinal lclNumRows =
        static_cast<LocalOrdinal> (this->getNodeNumRows ());

      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (this->isFillActive () == this->isFillComplete (), std::logic_error,
         "Graph cannot be both fill active and fill complete." << suffix);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (this->isFillComplete () &&
         (this->colMap_.is_null () ||
          this->rangeMap_.is_null () ||
          this->domainMap_.is_null ()),
         std::logic_error,
         "Graph is full complete, but at least one of {column, range, domain} "
         "Map is null." << suffix);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (this->isStorageOptimized () && ! this->indicesAreAllocated (),
         std::logic_error, "Storage is optimized, but indices are not "
         "allocated, not even trivially." << suffix);

      size_t nodeAllocSize = 0;
      try {
        nodeAllocSize = this->getNodeAllocationSize ();
      }
      catch (std::logic_error& e) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (true, std::runtime_error, "getNodeAllocationSize threw "
           "std::logic_error: " << e.what ());
      }
      catch (std::exception& e) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (true, std::runtime_error, "getNodeAllocationSize threw an "
           "std::exception: " << e.what ());
      }
      catch (...) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (true, std::runtime_error, "getNodeAllocationSize threw an exception "
           "not a subclass of std::exception.");
      }

      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (this->isStorageOptimized () &&
         nodeAllocSize != this->getNodeNumEntries (),
         std::logic_error, "Storage is optimized, but "
         "this->getNodeAllocationSize() = " << nodeAllocSize
         << " != this->getNodeNumEntries() = " << this->getNodeNumEntries ()
         << "." << suffix);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (! this->haveGlobalConstants_ &&
         (this->globalNumEntries_ != GSTI ||
          this->globalMaxNumRowEntries_ != GSTI),
         std::logic_error, "Graph claims not to have global constants, but "
         "some of the global constants are not marked as invalid." << suffix);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (this->haveGlobalConstants_ &&
         (this->globalNumEntries_ == GSTI ||
          this->globalMaxNumRowEntries_ == GSTI),
         std::logic_error, "Graph claims to have global constants, but "
         "some of them are marked as invalid." << suffix);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (this->haveGlobalConstants_ &&
         (this->globalNumEntries_ < this->getNodeNumEntries () ||
          this->globalMaxNumRowEntries_ < this->nodeMaxNumRowEntries_),
         std::logic_error, "Graph claims to have global constants, and "
         "all of the values of the global constants are valid, but "
         "some of the local constants are greater than "
         "their corresponding global constants." << suffix);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (this->indicesAreAllocated () &&
         (this->numAllocForAllRows_ != 0 ||
          this->k_numAllocPerRow_.extent (0) != 0),
         std::logic_error, "The graph claims that its indices are allocated, but "
         "either numAllocForAllRows_ (= " << this->numAllocForAllRows_ << ") is "
         "nonzero, or k_numAllocPerRow_ has nonzero dimension.  In other words, "
         "the graph is supposed to release its \"allocation specifications\" "
         "when it allocates its indices." << suffix);
      if (isGloballyIndexed() && k_rowPtrs_.extent(0) != 0) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (size_t(k_rowPtrs_.extent(0)) != size_t(lclNumRows + 1),
           std::logic_error, "The graph is globally indexed and "
           "k_rowPtrs_ has nonzero size " << k_rowPtrs_.extent(0)
           << ", but that size does not equal lclNumRows+1 = "
           << (lclNumRows+1) << "." << suffix);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (k_rowPtrs_(lclNumRows) != size_t(k_gblInds1D_.extent(0)),
           std::logic_error, "The graph is globally indexed and "
           "k_rowPtrs_ has nonzero size " << k_rowPtrs_.extent(0)
           << ", but k_rowPtrs_(lclNumRows=" << lclNumRows << ")="
           << k_rowPtrs_(lclNumRows) << " != k_gblInds1D_.extent(0)="
           << k_gblInds1D_.extent(0) << "." << suffix);
      }
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (this->isLocallyIndexed () &&
         this->k_rowPtrs_.extent (0) != 0 &&
         (static_cast<size_t> (k_rowPtrs_.extent (0)) != static_cast<size_t> (lclNumRows + 1) ||
          this->k_rowPtrs_(lclNumRows) != static_cast<size_t> (this->k_lclInds1D_.extent (0))),
         std::logic_error, "If k_rowPtrs_ has nonzero size and "
         "the graph is locally indexed, then "
         "k_rowPtrs_ must have N+1 rows, and "
         "k_rowPtrs_(N) must equal k_lclInds1D_.extent(0)." << suffix);

      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (this->indicesAreAllocated () &&
         nodeAllocSize > 0 &&
         this->k_lclInds1D_.extent (0) == 0 &&
         this->k_gblInds1D_.extent (0) == 0,
         std::logic_error, "Graph is allocated nontrivially, but "
         "but 1-D allocations are not present." << suffix);

      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (! this->indicesAreAllocated () &&
         ((this->k_rowPtrs_.extent (0) != 0 ||
           this->k_numRowEntries_.extent (0) != 0) ||
          this->k_lclInds1D_.extent (0) != 0 ||
          this->k_gblInds1D_.extent (0) != 0),
         std::logic_error, "If indices are not allocated, "
         "then none of the buffers should be." << suffix);
      // indices may be local or global only if they are allocated
      // (numAllocated is redundant; could simply be indicesAreLocal_ ||
      // indicesAreGlobal_)
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        ((this->indicesAreLocal_ || this->indicesAreGlobal_) &&
         ! this->indicesAreAllocated_,
         std::logic_error, "Indices may be local or global only if they are "
         "allocated." << suffix);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (this->indicesAreLocal_ && this->indicesAreGlobal_,
         std::logic_error, "Indices may not be both local and global." << suffix);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (indicesAreLocal_ && k_gblInds1D_.extent (0) != 0,
         std::logic_error, "Indices are local, but "
         "k_gblInds1D_.extent(0) (= " << k_gblInds1D_.extent (0)
         << ") != 0.  In other words, if indices are local, then "
         "allocations of global indices should not be present."
         << suffix);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (indicesAreGlobal_ && k_lclInds1D_.extent (0) != 0,
         std::logic_error, "Indices are global, but "
         "k_lclInds1D_.extent(0) (= " << k_lclInds1D_.extent(0)
         << ") != 0.  In other words, if indices are global, "
         "then allocations for local indices should not be present."
         << suffix);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (indicesAreLocal_ && nodeAllocSize > 0 &&
         k_lclInds1D_.extent (0) == 0 && getNodeNumRows () > 0,
         std::logic_error, "Indices are local and "
         "getNodeAllocationSize() = " << nodeAllocSize << " > 0, but "
         "k_lclInds1D_.extent(0) = 0 and getNodeNumRows() = "
         << getNodeNumRows () << " > 0." << suffix);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (indicesAreGlobal_ && nodeAllocSize > 0 &&
         k_gblInds1D_.extent (0) == 0 && getNodeNumRows () > 0,
         std::logic_error, "Indices are global and "
         "getNodeAllocationSize() = " << nodeAllocSize << " > 0, but "
         "k_gblInds1D_.extent(0) = 0 and getNodeNumRows() = "
         << getNodeNumRows () << " > 0." << suffix);
      // check the actual allocations
      if (this->indicesAreAllocated () &&
          this->k_rowPtrs_.extent (0) != 0) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (static_cast<size_t> (this->k_rowPtrs_.extent (0)) !=
           this->getNodeNumRows () + 1,
           std::logic_error, "Indices are allocated and "
           "k_rowPtrs_ has nonzero length, but k_rowPtrs_.extent(0) = "
           << this->k_rowPtrs_.extent (0) << " != getNodeNumRows()+1 = "
           << (this->getNodeNumRows () + 1) << "." << suffix);
        const size_t actualNumAllocated =
          ::Tpetra::Details::getEntryOnHost (this->k_rowPtrs_, this->getNodeNumRows ());
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (this->isLocallyIndexed () &&
           static_cast<size_t> (this->k_lclInds1D_.extent (0)) != actualNumAllocated,
           std::logic_error, "Graph is locally indexed, indices are "
           "are allocated, and k_rowPtrs_ has nonzero length, but "
           "k_lclInds1D_.extent(0) = " << this->k_lclInds1D_.extent (0)
           << " != actualNumAllocated = " << actualNumAllocated << suffix);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (this->isGloballyIndexed () &&
           static_cast<size_t> (this->k_gblInds1D_.extent (0)) != actualNumAllocated,
           std::logic_error, "Graph is globally indexed, indices "
           "are allocated, and k_rowPtrs_ has nonzero length, but "
           "k_gblInds1D_.extent(0) = " << this->k_gblInds1D_.extent (0)
           << " != actualNumAllocated = " << actualNumAllocated << suffix);
      }

      if (verbose_) {
        std::ostringstream os;
        os << *prefix << "Done" << endl;
        std::cerr << os.str();
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getNumEntriesInGlobalRow (GlobalOrdinal globalRow) const
  {
    const RowInfo rowInfo = this->getRowInfoFromGlobalRowIndex (globalRow);
    if (rowInfo.localRow == Teuchos::OrdinalTraits<size_t>::invalid ()) {
      return Teuchos::OrdinalTraits<size_t>::invalid ();
    }
    else {
      return rowInfo.numEntries;
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getNumEntriesInLocalRow (LocalOrdinal localRow) const
  {
    const RowInfo rowInfo = this->getRowInfo (localRow);
    if (rowInfo.localRow == Teuchos::OrdinalTraits<size_t>::invalid ()) {
      return Teuchos::OrdinalTraits<size_t>::invalid ();
    }
    else {
      return rowInfo.numEntries;
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getNumAllocatedEntriesInGlobalRow (GlobalOrdinal globalRow) const
  {
    const RowInfo rowInfo = this->getRowInfoFromGlobalRowIndex (globalRow);
    if (rowInfo.localRow == Teuchos::OrdinalTraits<size_t>::invalid ()) {
      return Teuchos::OrdinalTraits<size_t>::invalid ();
    }
    else {
      return rowInfo.allocSize;
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getNumAllocatedEntriesInLocalRow (LocalOrdinal localRow) const
  {
    const RowInfo rowInfo = this->getRowInfo (localRow);
    if (rowInfo.localRow == Teuchos::OrdinalTraits<size_t>::invalid ()) {
      return Teuchos::OrdinalTraits<size_t>::invalid ();
    }
    else {
      return rowInfo.allocSize;
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<const size_t>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getNodeRowPtrs () const
  {
    using Kokkos::ViewAllocateWithoutInitializing;
    using Kokkos::create_mirror_view;
    using Teuchos::ArrayRCP;
    typedef typename local_graph_type::row_map_type row_map_type;
    typedef typename row_map_type::non_const_value_type row_offset_type;
    const char prefix[] = "Tpetra::CrsGraph::getNodeRowPtrs: ";
    const char suffix[] = "  Please report this bug to the Tpetra developers.";

    const size_t size = k_rowPtrs_.extent (0);
    constexpr bool same = std::is_same<size_t, row_offset_type>::value;

    if (size == 0) {
      return ArrayRCP<const size_t> ();
    }

    ArrayRCP<const row_offset_type> ptr_rot;
    ArrayRCP<const size_t> ptr_st;
    if (same) { // size_t == row_offset_type
      // NOTE (mfh 22 Mar 2015) In a debug build of Kokkos, the result
      // of create_mirror_view might actually be a new allocation.
      // This helps with debugging when there are two memory spaces.
      typename row_map_type::HostMirror ptr_h = create_mirror_view (k_rowPtrs_);
      Kokkos::deep_copy (ptr_h, k_rowPtrs_);
      if (debug_) {
        TEUCHOS_TEST_FOR_EXCEPTION
          (ptr_h.extent (0) != k_rowPtrs_.extent (0), std::logic_error,
           prefix << "size_t == row_offset_type, but ptr_h.extent(0) = "
           << ptr_h.extent (0) << " != k_rowPtrs_.extent(0) = "
           << k_rowPtrs_.extent (0) << ".");
        TEUCHOS_TEST_FOR_EXCEPTION
          (same && size != 0 && k_rowPtrs_.data () == nullptr, std::logic_error,
           prefix << "size_t == row_offset_type and k_rowPtrs_.extent(0) = "
           << size << " != 0, but k_rowPtrs_.data() == nullptr." << suffix);
        TEUCHOS_TEST_FOR_EXCEPTION
          (same && size != 0 && ptr_h.data () == nullptr, std::logic_error,
           prefix << "size_t == row_offset_type and k_rowPtrs_.extent(0) = "
           << size << " != 0, but create_mirror_view(k_rowPtrs_).data() "
           "== nullptr." << suffix);
      }
      ptr_rot = Kokkos::Compat::persistingView (ptr_h);
    }
    else { // size_t != row_offset_type
      typedef Kokkos::View<size_t*, device_type> ret_view_type;
      ret_view_type ptr_d (ViewAllocateWithoutInitializing ("ptr"), size);
      ::Tpetra::Details::copyOffsets (ptr_d, k_rowPtrs_);
      typename ret_view_type::HostMirror ptr_h = create_mirror_view (ptr_d);
      Kokkos::deep_copy (ptr_h, ptr_d);
      ptr_st = Kokkos::Compat::persistingView (ptr_h);
    }
    if (debug_) {
      TEUCHOS_TEST_FOR_EXCEPTION
        (same && size != 0 && ptr_rot.is_null (), std::logic_error,
         prefix << "size_t == row_offset_type and size = " << size
         << " != 0, but ptr_rot is null." << suffix);
      TEUCHOS_TEST_FOR_EXCEPTION
        (! same && size != 0 && ptr_st.is_null (), std::logic_error,
         prefix << "size_t != row_offset_type and size = " << size
         << " != 0, but ptr_st is null." << suffix);
    }

    // If size_t == row_offset_type, return a persisting host view of
    // k_rowPtrs_.  Otherwise, return a size_t host copy of k_rowPtrs_.
    ArrayRCP<const size_t> retval =
      Kokkos::Impl::if_c<same,
        ArrayRCP<const row_offset_type>,
        ArrayRCP<const size_t> >::select (ptr_rot, ptr_st);
    if (debug_) {
      TEUCHOS_TEST_FOR_EXCEPTION
        (size != 0 && retval.is_null (), std::logic_error,
         prefix << "size = " << size << " != 0, but retval is null." << suffix);
    }
    return retval;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayRCP<const LocalOrdinal>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getNodePackedIndices () const
  {
    return Kokkos::Compat::persistingView (k_lclInds1D_);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getLocalRowCopy (LocalOrdinal localRow,
                   const Teuchos::ArrayView<LocalOrdinal>&indices,
                   size_t& numEntries) const
  {
    using Teuchos::ArrayView;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    const char tfecfFuncName[] = "getLocalRowCopy: ";

    TEUCHOS_TEST_FOR_EXCEPTION(
      isGloballyIndexed () && ! hasColMap (), std::runtime_error,
      "Tpetra::CrsGraph::getLocalRowCopy: The graph is globally indexed and "
      "does not have a column Map yet.  That means we don't have local indices "
      "for columns yet, so it doesn't make sense to call this method.  If the "
      "graph doesn't have a column Map yet, you should call fillComplete on "
      "it first.");

    // This does the right thing (reports an empty row) if the input
    // row is invalid.
    const RowInfo rowinfo = this->getRowInfo (localRow);
    // No side effects on error.
    const size_t theNumEntries = rowinfo.numEntries;
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (static_cast<size_t> (indices.size ()) < theNumEntries, std::runtime_error,
       "Specified storage (size==" << indices.size () << ") does not suffice "
       "to hold all " << theNumEntries << " entry/ies for this row.");
    numEntries = theNumEntries;

    if (rowinfo.localRow != Teuchos::OrdinalTraits<size_t>::invalid ()) {
      if (isLocallyIndexed ()) {
        ArrayView<const LO> lview = getLocalView (rowinfo);
        for (size_t j = 0; j < theNumEntries; ++j) {
          indices[j] = lview[j];
        }
      }
      else if (isGloballyIndexed ()) {
        ArrayView<const GO> gview = getGlobalView (rowinfo);
        for (size_t j = 0; j < theNumEntries; ++j) {
          indices[j] = colMap_->getLocalElement (gview[j]);
        }
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getGlobalRowCopy (GlobalOrdinal globalRow,
                    const Teuchos::ArrayView<GlobalOrdinal>& indices,
                    size_t& numEntries) const
  {
    using Teuchos::ArrayView;
    const char tfecfFuncName[] = "getGlobalRowCopy: ";

    // This does the right thing (reports an empty row) if the input
    // row is invalid.
    const RowInfo rowinfo = getRowInfoFromGlobalRowIndex (globalRow);
    const size_t theNumEntries = rowinfo.numEntries;
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (indices.size ()) < theNumEntries, std::runtime_error,
      "Specified storage (size==" << indices.size () << ") does not suffice "
      "to hold all " << theNumEntries << " entry/ies for this row.");
    numEntries = theNumEntries; // first side effect

    if (rowinfo.localRow != Teuchos::OrdinalTraits<size_t>::invalid ()) {
      if (isLocallyIndexed ()) {
        ArrayView<const LocalOrdinal> lview = getLocalView (rowinfo);
        for (size_t j = 0; j < theNumEntries; ++j) {
          indices[j] = colMap_->getGlobalElement (lview[j]);
        }
      }
      else if (isGloballyIndexed ()) {
        ArrayView<const GlobalOrdinal> gview = getGlobalView (rowinfo);
        for (size_t j = 0; j < theNumEntries; ++j) {
          indices[j] = gview[j];
        }
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getLocalRowView (const LocalOrdinal localRow,
                   Teuchos::ArrayView<const LocalOrdinal>& indices) const
  {
    const char tfecfFuncName[] = "getLocalRowView: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (isGloballyIndexed (), std::runtime_error, "The graph's indices are "
       "currently stored as global indices, so we cannot return a view with "
       "local column indices, whether or not the graph has a column Map.  If "
       "the graph _does_ have a column Map, use getLocalRowCopy() instead.");

    // This does the right thing (reports an empty row) if the input
    // row is invalid.
    const RowInfo rowInfo = getRowInfo (localRow);
    indices = Teuchos::null;
    if (rowInfo.localRow != Teuchos::OrdinalTraits<size_t>::invalid () &&
        rowInfo.numEntries > 0) {
      indices = this->getLocalView (rowInfo);
      // getLocalView returns a view of the _entire_ row, including
      // any extra space at the end (which 1-D unpacked storage
      // might have, for example).  That's why we have to take a
      // subview of the returned view.
      indices = indices (0, rowInfo.numEntries);
    }

    if (debug_) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (static_cast<size_t> (indices.size ()) !=
         getNumEntriesInLocalRow (localRow), std::logic_error, "indices.size() "
         "= " << indices.size () << " != getNumEntriesInLocalRow(localRow=" <<
         localRow << ") = " << getNumEntriesInLocalRow (localRow) <<
         ".  Please report this bug to the Tpetra developers.");
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getGlobalRowView (const GlobalOrdinal globalRow,
                    Teuchos::ArrayView<const GlobalOrdinal>& indices) const
  {
    const char tfecfFuncName[] = "getGlobalRowView: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (isLocallyIndexed (), std::runtime_error, "The graph's indices are "
       "currently stored as local indices, so we cannot return a view with "
       "global column indices.  Use getGlobalRowCopy() instead.");

    // This does the right thing (reports an empty row) if the input
    // row is invalid.
    const RowInfo rowInfo = getRowInfoFromGlobalRowIndex (globalRow);
    indices = Teuchos::null;
    if (rowInfo.localRow != Teuchos::OrdinalTraits<size_t>::invalid () &&
        rowInfo.numEntries > 0) {
      indices = (this->getGlobalView (rowInfo)) (0, rowInfo.numEntries);
    }

    if (debug_) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (static_cast<size_t> (indices.size ()) !=
         getNumEntriesInGlobalRow (globalRow),
         std::logic_error, "indices.size() = " << indices.size ()
         << " != getNumEntriesInGlobalRow(globalRow=" << globalRow << ") = "
         << getNumEntriesInGlobalRow (globalRow)
         << ".  Please report this bug to the Tpetra developers.");
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  insertLocalIndices (const LocalOrdinal localRow,
                      const Teuchos::ArrayView<const LocalOrdinal>& indices)
  {
    const char tfecfFuncName[] = "insertLocalIndices: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! isFillActive (), std::runtime_error, "Fill must be active.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (isGloballyIndexed (), std::runtime_error,
       "Graph indices are global; use insertGlobalIndices().");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! hasColMap (), std::runtime_error,
       "Cannot insert local indices without a column Map.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! rowMap_->isNodeLocalElement (localRow), std::runtime_error,
       "Local row index " << localRow << " is not in the row Map "
       "on the calling process.");
    if (! indicesAreAllocated ()) {
      allocateIndices (LocalIndices, verbose_);
    }

    if (debug_) {
      // In debug mode, if the graph has a column Map, test whether any
      // of the given column indices are not in the column Map.  Keep
      // track of the invalid column indices so we can tell the user
      // about them.
      if (hasColMap ()) {
        using Teuchos::Array;
        using Teuchos::toString;
        using std::endl;
        typedef typename Teuchos::ArrayView<const LocalOrdinal>::size_type size_type;

        const map_type& colMap = *colMap_;
        Array<LocalOrdinal> badColInds;
        bool allInColMap = true;
        for (size_type k = 0; k < indices.size (); ++k) {
          if (! colMap.isNodeLocalElement (indices[k])) {
            allInColMap = false;
            badColInds.push_back (indices[k]);
          }
        }
        if (! allInColMap) {
          std::ostringstream os;
          os << "Tpetra::CrsGraph::insertLocalIndices: You attempted to insert "
            "entries in owned row " << localRow << ", at the following column "
            "indices: " << toString (indices) << "." << endl;
          os << "Of those, the following indices are not in the column Map on "
            "this process: " << toString (badColInds) << "." << endl << "Since "
            "the graph has a column Map already, it is invalid to insert entries "
            "at those locations.";
          TEUCHOS_TEST_FOR_EXCEPTION(! allInColMap, std::invalid_argument, os.str ());
        }
      }
    }

    insertLocalIndicesImpl (localRow, indices);

    if (debug_) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (! indicesAreAllocated () || ! isLocallyIndexed (), std::logic_error,
         "At the end of insertLocalIndices, ! indicesAreAllocated() || "
         "! isLocallyIndexed() is true.  Please report this bug to the "
         "Tpetra developers.");
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  insertLocalIndices (const LocalOrdinal localRow,
                      const LocalOrdinal numEnt,
                      const LocalOrdinal inds[])
  {
    Teuchos::ArrayView<const LocalOrdinal> indsT (inds, numEnt);
    this->insertLocalIndices (localRow, indsT);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  insertGlobalIndices (const GlobalOrdinal gblRow,
                       const LocalOrdinal numInputInds,
                       const GlobalOrdinal inputGblColInds[])
  {
    typedef LocalOrdinal LO;
    const char tfecfFuncName[] = "insertGlobalIndices: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (this->isLocallyIndexed (), std::runtime_error,
      "graph indices are local; use insertLocalIndices().");
    // This can't really be satisfied for now, because if we are
    // fillComplete(), then we are local.  In the future, this may
    // change.  However, the rule that modification require active
    // fill will not change.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! this->isFillActive (), std::runtime_error,
      "You are not allowed to call this method if fill is not active.  "
      "If fillComplete has been called, you must first call resumeFill "
      "before you may insert indices.");
    if (! indicesAreAllocated ()) {
      allocateIndices (GlobalIndices, verbose_);
    }
    const LO lclRow = this->rowMap_->getLocalElement (gblRow);
    if (lclRow != Tpetra::Details::OrdinalTraits<LO>::invalid ()) {
      if (debug_) {
        if (this->hasColMap ()) {
          using std::endl;
          const map_type& colMap = * (this->colMap_);
          // In a debug build, keep track of the nonowned ("bad") column
          // indices, so that we can display them in the exception
          // message.  In a release build, just ditch the loop early if
          // we encounter a nonowned column index.
          std::vector<GlobalOrdinal> badColInds;
          bool allInColMap = true;
          for (LO k = 0; k < numInputInds; ++k) {
            if (! colMap.isNodeGlobalElement (inputGblColInds[k])) {
              allInColMap = false;
              badColInds.push_back (inputGblColInds[k]);
            }
          }
          if (! allInColMap) {
            std::ostringstream os;
            os << "You attempted to insert entries in owned row " << gblRow
               << ", at the following column indices: [";
            for (LO k = 0; k < numInputInds; ++k) {
              os << inputGblColInds[k];
              if (k + static_cast<LO> (1) < numInputInds) {
                os << ",";
              }
            }
            os << "]." << endl << "Of those, the following indices are not in "
              "the column Map on this process: [";
            for (size_t k = 0; k < badColInds.size (); ++k) {
              os << badColInds[k];
              if (k + size_t (1) < badColInds.size ()) {
                os << ",";
              }
            }
            os << "]." << endl << "Since the matrix has a column Map already, "
              "it is invalid to insert entries at those locations.";
            TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
              (true, std::invalid_argument, os.str ());
          }
        }
      } // debug_
      this->insertGlobalIndicesImpl (lclRow, inputGblColInds, numInputInds);
    }
    else { // a nonlocal row
      this->insertGlobalIndicesIntoNonownedRows (gblRow, inputGblColInds,
                                                 numInputInds);
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  insertGlobalIndices (const GlobalOrdinal gblRow,
                       const Teuchos::ArrayView<const GlobalOrdinal>& inputGblColInds)
  {
    this->insertGlobalIndices (gblRow, inputGblColInds.size (),
                               inputGblColInds.getRawPtr ());
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  insertGlobalIndicesFiltered (const LocalOrdinal lclRow,
                               const GlobalOrdinal gblColInds[],
                               const LocalOrdinal numGblColInds)
  {
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    const char tfecfFuncName[] = "insertGlobalIndicesFiltered: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (this->isLocallyIndexed (), std::runtime_error,
       "Graph indices are local; use insertLocalIndices().");
    // This can't really be satisfied for now, because if we are
    // fillComplete(), then we are local.  In the future, this may
    // change.  However, the rule that modification require active
    // fill will not change.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! this->isFillActive (), std::runtime_error,
       "You are not allowed to call this method if fill is not active.  "
       "If fillComplete has been called, you must first call resumeFill "
       "before you may insert indices.");
    if (! indicesAreAllocated ()) {
      allocateIndices (GlobalIndices, verbose_);
    }

    Teuchos::ArrayView<const GO> gblColInds_av (gblColInds, numGblColInds);
    // If we have a column Map, use it to filter the entries.
    if (! colMap_.is_null ()) {
      const map_type& colMap = * (this->colMap_);

      LO curOffset = 0;
      while (curOffset < numGblColInds) {
        // Find a sequence of input indices that are in the column Map
        // on the calling process.  Doing a sequence at a time,
        // instead of one at a time, amortizes some overhead.
        LO endOffset = curOffset;
        for ( ; endOffset < numGblColInds; ++endOffset) {
          const LO lclCol = colMap.getLocalElement (gblColInds[endOffset]);
          if (lclCol == Tpetra::Details::OrdinalTraits<LO>::invalid ()) {
            break; // first entry, in current sequence, not in the column Map
          }
        }
        // curOffset, endOffset: half-exclusive range of indices in
        // the column Map on the calling process.  If endOffset ==
        // curOffset, the range is empty.
        const LO numIndInSeq = (endOffset - curOffset);
        if (numIndInSeq != 0) {
          this->insertGlobalIndicesImpl (lclRow, gblColInds + curOffset,
                                         numIndInSeq);
        }
        // Invariant before this line: Either endOffset ==
        // numGblColInds, or gblColInds[endOffset] is not in the
        // column Map on the calling process.
        curOffset = endOffset + 1;
      }
    }
    else {
      this->insertGlobalIndicesImpl (lclRow, gblColInds_av.getRawPtr (),
                                     gblColInds_av.size ());
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  insertGlobalIndicesIntoNonownedRows (const GlobalOrdinal gblRow,
                                       const GlobalOrdinal gblColInds[],
                                       const LocalOrdinal numGblColInds)
  {
    // This creates the std::vector if it doesn't exist yet.
    // std::map's operator[] does a lookup each time, so it's better
    // to pull nonlocals_[grow] out of the loop.
    std::vector<GlobalOrdinal>& nonlocalRow = this->nonlocals_[gblRow];
    for (LocalOrdinal k = 0; k < numGblColInds; ++k) {
      // FIXME (mfh 20 Jul 2017) Would be better to use a set, in
      // order to avoid duplicates.  globalAssemble() sorts these
      // anyway.
      nonlocalRow.push_back (gblColInds[k]);
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  removeLocalIndices (LocalOrdinal lrow)
  {
    const char tfecfFuncName[] = "removeLocalIndices: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isFillActive (), std::runtime_error, "requires that fill is active.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isStorageOptimized (), std::runtime_error,
      "cannot remove indices after optimizeStorage() has been called.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isGloballyIndexed (), std::runtime_error, "graph indices are global.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! rowMap_->isNodeLocalElement (lrow), std::runtime_error,
      "Local row " << lrow << " is not in the row Map on the calling process.");
    if (! indicesAreAllocated ()) {
      allocateIndices (LocalIndices, verbose_);
    }

    // FIXME (mfh 13 Aug 2014) What if they haven't been cleared on
    // all processes?
    clearGlobalConstants ();

    if (k_numRowEntries_.extent (0) != 0) {
      this->k_numRowEntries_(lrow) = 0;
    }

    if (debug_) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (getNumEntriesInLocalRow (lrow) != 0 ||
         ! indicesAreAllocated () ||
         ! isLocallyIndexed (), std::logic_error,
         "Violated stated post-conditions. Please contact Tpetra team.");
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  setAllIndices (const typename local_graph_type::row_map_type& rowPointers,
                 const typename local_graph_type::entries_type::non_const_type& columnIndices)
  {
    const char tfecfFuncName[] = "setAllIndices: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! hasColMap () || getColMap ().is_null (), std::runtime_error,
      "The graph must have a column Map before you may call this method.");
    LocalOrdinal numLocalRows = this->getNodeNumRows ();
    {
      LocalOrdinal rowPtrLen = rowPointers.size();
      if(numLocalRows == 0) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          rowPtrLen != 0 && rowPtrLen != 1,
          std::runtime_error, "Have 0 local rows, but rowPointers.size() is neither 0 nor 1.");
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          rowPtrLen != numLocalRows + 1,
          std::runtime_error, "rowPointers.size() = " << rowPtrLen <<
          " != this->getNodeNumRows()+1 = " << (numLocalRows + 1) << ".");
      }
    }

    if (debug_ && this->isSorted()) {
      // Verify that the local indices are actually sorted
      int notSorted = 0;
      using exec_space = typename local_graph_type::execution_space;
      using size_type = typename local_graph_type::size_type;
      Kokkos::parallel_reduce(Kokkos::RangePolicy<exec_space>(0, numLocalRows),
        KOKKOS_LAMBDA (const LocalOrdinal i, int& lNotSorted)
        {
          size_type rowBegin = rowPointers(i);
          size_type rowEnd = rowPointers(i + 1);
          for(size_type j = rowBegin + 1; j < rowEnd; j++)
          {
            if(columnIndices(j - 1) > columnIndices(j))
            {
              lNotSorted = 1;
            }
          }
        }, notSorted);
      //All-reduce notSorted to avoid rank divergence
      int globalNotSorted = 0;
      auto comm = this->getComm();
      Teuchos::reduceAll<int, int> (*comm, Teuchos::REDUCE_MAX, notSorted,
                           Teuchos::outArg (globalNotSorted));
      if (globalNotSorted)
      {
        std::string message;
        if (notSorted)
        {
          //Only print message from ranks with the problem
          message = std::string("ERROR, rank ") + std::to_string(comm->getRank()) + ", CrsGraph::setAllIndices(): provided columnIndices are not sorted!\n";
        }
        Details::gathervPrint(std::cout, message, *comm);
        throw std::invalid_argument("CrsGraph::setAllIndices(): provided columnIndices are not sorted within rows on at least one process.");
      }
    }

    // FIXME (mfh 07 Aug 2014) We need to relax this restriction,
    // since the future model will be allocation at construction, not
    // lazy allocation on first insert.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      ((this->k_lclInds1D_.extent (0) != 0 || this->k_gblInds1D_.extent (0) != 0),
       std::runtime_error, "You may not call this method if 1-D data "
       "structures are already allocated.");

    indicesAreAllocated_ = true;
    indicesAreLocal_     = true;
    indicesAreSorted_    = true;
    noRedundancies_      = true;
    k_lclInds1D_         = columnIndices;
    k_rowPtrs_           = rowPointers;
    // Storage MUST be packed, since the interface doesn't give any
    // way to indicate any extra space at the end of each row.
    storageStatus_       = Details::STORAGE_1D_PACKED;

    // Build the local graph.
    lclGraph_ = local_graph_type (k_lclInds1D_, k_rowPtrs_);

    // These normally get cleared out at the end of allocateIndices.
    // It makes sense to clear them out here, because at the end of
    // this method, the graph is allocated on the calling process.
    numAllocForAllRows_ = 0;
    k_numAllocPerRow_ = decltype (k_numAllocPerRow_) ();

    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  setAllIndices (const Teuchos::ArrayRCP<size_t>& rowPointers,
                 const Teuchos::ArrayRCP<LocalOrdinal>& columnIndices)
  {
    using Kokkos::View;
    typedef typename local_graph_type::row_map_type row_map_type;
    typedef typename row_map_type::array_layout layout_type;
    typedef typename row_map_type::non_const_value_type row_offset_type;
    typedef View<size_t*, layout_type , Kokkos::HostSpace,
      Kokkos::MemoryUnmanaged> input_view_type;
    typedef typename row_map_type::non_const_type nc_row_map_type;

    const size_t size = static_cast<size_t> (rowPointers.size ());
    constexpr bool same = std::is_same<size_t, row_offset_type>::value;
    input_view_type ptr_in (rowPointers.getRawPtr (), size);

    nc_row_map_type ptr_rot ("Tpetra::CrsGraph::ptr", size);

    if (same) { // size_t == row_offset_type
      // This compile-time logic ensures that the compiler never sees
      // an assignment of View<row_offset_type*, ...> to View<size_t*,
      // ...> unless size_t == row_offset_type.
      input_view_type ptr_decoy (rowPointers.getRawPtr (), size); // never used
      Kokkos::deep_copy (Kokkos::Impl::if_c<same,
                           nc_row_map_type,
                           input_view_type>::select (ptr_rot, ptr_decoy),
                         ptr_in);
    }
    else { // size_t != row_offset_type
      // CudaUvmSpace != HostSpace, so this will be false in that case.
      constexpr bool inHostMemory =
        std::is_same<typename row_map_type::memory_space,
          Kokkos::HostSpace>::value;
      if (inHostMemory) {
        // Copy (with cast from size_t to row_offset_type, with bounds
        // checking if necessary) to ptr_rot.
        ::Tpetra::Details::copyOffsets (ptr_rot, ptr_in);
      }
      else { // Copy input row offsets to device first.
        //
        // FIXME (mfh 24 Mar 2015) If CUDA UVM, running in the host's
        // execution space would avoid the double copy.
        //
        View<size_t*, layout_type ,execution_space > ptr_st ("Tpetra::CrsGraph::ptr", size);
        Kokkos::deep_copy (ptr_st, ptr_in);
        // Copy on device (casting from size_t to row_offset_type,
        // with bounds checking if necessary) to ptr_rot.  This
        // executes in the output View's execution space, which is the
        // same as execution_space.
        ::Tpetra::Details::copyOffsets (ptr_rot, ptr_st);
      }
    }

    Kokkos::View<LocalOrdinal*, layout_type , execution_space > k_ind =
      Kokkos::Compat::getKokkosViewDeepCopy<device_type> (columnIndices ());
    setAllIndices (ptr_rot, k_ind);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getNumEntriesPerLocalRowUpperBound (Teuchos::ArrayRCP<const size_t>& boundPerLocalRow,
                                      size_t& boundForAllLocalRows,
                                      bool& boundSameForAllLocalRows) const
  {
    const char tfecfFuncName[] = "getNumEntriesPerLocalRowUpperBound: ";
    const char suffix[] = "  Please report this bug to the Tpetra developers.";

    // The three output arguments.  We assign them to the actual
    // output arguments at the end, in order to implement
    // transactional semantics.
    Teuchos::ArrayRCP<const size_t> numEntriesPerRow;
    size_t numEntriesForAll = 0;
    bool allRowsSame = true;

    const ptrdiff_t numRows = static_cast<ptrdiff_t> (this->getNodeNumRows ());

    if (this->indicesAreAllocated ()) {
      if (this->isStorageOptimized ()) {
        // left with the case that we have optimized storage. in this
        // case, we have to construct a list of row sizes.
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (numRows != 0 && k_rowPtrs_.extent (0) == 0, std::logic_error,
           "The graph has " << numRows << " (> 0) row"
           << (numRows != 1 ? "s" : "") << " on the calling process, "
           "but the k_rowPtrs_ array has zero entries." << suffix);
        Teuchos::ArrayRCP<size_t> numEnt;
        if (numRows != 0) {
          numEnt = Teuchos::arcp<size_t> (numRows);
        }

        // We have to iterate through the row offsets anyway, so we
        // might as well check whether all rows' bounds are the same.
        bool allRowsReallySame = false;
        for (ptrdiff_t i = 0; i < numRows; ++i) {
          numEnt[i] = this->k_rowPtrs_(i+1) - this->k_rowPtrs_(i);
          if (i != 0 && numEnt[i] != numEnt[i-1]) {
            allRowsReallySame = false;
          }
        }
        if (allRowsReallySame) {
          if (numRows == 0) {
            numEntriesForAll = 0;
          } else {
            numEntriesForAll = numEnt[1] - numEnt[0];
          }
          allRowsSame = true;
        }
        else {
          numEntriesPerRow = numEnt; // Teuchos::arcp_const_cast<const size_t> (numEnt);
          allRowsSame = false; // conservatively; we don't check the array
        }
      }
      else if (k_numRowEntries_.extent (0) != 0) {
        // This is a shallow copy; the ArrayRCP wraps the View in a
        // custom destructor, which ensures correct deallocation if
        // that is the only reference to the View.  Furthermore, this
        // View is a host View, so this doesn't assume UVM.
        numEntriesPerRow = Kokkos::Compat::persistingView (k_numRowEntries_);
        allRowsSame = false; // conservatively; we don't check the array
      }
      else {
        numEntriesForAll = 0;
        allRowsSame = true;
      }
    }
    else { // indices not allocated
      if (k_numAllocPerRow_.extent (0) != 0) {
        // This is a shallow copy; the ArrayRCP wraps the View in a
        // custom destructor, which ensures correct deallocation if
        // that is the only reference to the View.  Furthermore, this
        // View is a host View, so this doesn't assume UVM.
        numEntriesPerRow = Kokkos::Compat::persistingView (k_numAllocPerRow_);
        allRowsSame = false; // conservatively; we don't check the array
      }
      else {
        numEntriesForAll = numAllocForAllRows_;
        allRowsSame = true;
      }
    }

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (numEntriesForAll != 0 && numEntriesPerRow.size () != 0, std::logic_error,
       "numEntriesForAll and numEntriesPerRow are not consistent.  The former "
       "is nonzero (" << numEntriesForAll << "), but the latter has nonzero "
       "size " << numEntriesPerRow.size () << "." << suffix);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (numEntriesForAll != 0 && ! allRowsSame, std::logic_error,
       "numEntriesForAll and allRowsSame are not consistent.  The former "
       "is nonzero (" << numEntriesForAll << "), but the latter is false."
       << suffix);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (numEntriesPerRow.size () != 0 && allRowsSame, std::logic_error,
       "numEntriesPerRow and allRowsSame are not consistent.  The former has "
       "nonzero length " << numEntriesForAll << ", but the latter is true."
       << suffix);

    boundPerLocalRow = numEntriesPerRow;
    boundForAllLocalRows = numEntriesForAll;
    boundSameForAllLocalRows = allRowsSame;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  globalAssemble ()
  {
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::REDUCE_MAX;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using std::endl;
    using crs_graph_type = CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
    using LO = local_ordinal_type;
    using GO = global_ordinal_type;
    using size_type = typename Teuchos::Array<GO>::size_type;
    const char tfecfFuncName[] = "globalAssemble: "; // for exception macro

    std::unique_ptr<std::string> prefix;
    if (verbose_) {
      prefix = this->createPrefix("CrsGraph", "globalAssemble");
      std::ostringstream os;
      os << *prefix << "Start" << endl;
      std::cerr << os.str();
    }
    RCP<const Comm<int> > comm = getComm ();

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! isFillActive (), std::runtime_error, "Fill must be active before "
       "you may call this method.");

    const size_t myNumNonlocalRows = this->nonlocals_.size ();

    // If no processes have nonlocal rows, then we don't have to do
    // anything.  Checking this is probably cheaper than constructing
    // the Map of nonlocal rows (see below) and noticing that it has
    // zero global entries.
    {
      const int iHaveNonlocalRows = (myNumNonlocalRows == 0) ? 0 : 1;
      int someoneHasNonlocalRows = 0;
      reduceAll<int, int> (*comm, REDUCE_MAX, iHaveNonlocalRows,
                           outArg (someoneHasNonlocalRows));
      if (someoneHasNonlocalRows == 0) {
        if (verbose_) {
          std::ostringstream os;
          os << *prefix << "Done: No nonlocal rows" << endl;
          std::cerr << os.str();
        }
        return;
      }
      else if (verbose_) {
        std::ostringstream os;
        os << *prefix << "At least 1 process has nonlocal rows"
           << endl;
        std::cerr << os.str();
      }
    }

    // 1. Create a list of the "nonlocal" rows on each process.  this
    //    requires iterating over nonlocals_, so while we do this,
    //    deduplicate the entries and get a count for each nonlocal
    //    row on this process.
    // 2. Construct a new row Map corresponding to those rows.  This
    //    Map is likely overlapping.  We know that the Map is not
    //    empty on all processes, because the above all-reduce and
    //    return exclude that case.

    RCP<const map_type> nonlocalRowMap;
    // Keep this for CrsGraph's constructor.
    Teuchos::Array<size_t> numEntPerNonlocalRow (myNumNonlocalRows);
    {
      Teuchos::Array<GO> myNonlocalGblRows (myNumNonlocalRows);
      size_type curPos = 0;
      for (auto mapIter = this->nonlocals_.begin ();
           mapIter != this->nonlocals_.end ();
           ++mapIter, ++curPos) {
        myNonlocalGblRows[curPos] = mapIter->first;
        std::vector<GO>& gblCols = mapIter->second; // by ref; change in place
        std::sort (gblCols.begin (), gblCols.end ());
        auto vecLast = std::unique (gblCols.begin (), gblCols.end ());
        gblCols.erase (vecLast, gblCols.end ());
        numEntPerNonlocalRow[curPos] = gblCols.size ();
      }

      // Currently, Map requires that its indexBase be the global min
      // of all its global indices.  Map won't compute this for us, so
      // we must do it.  If our process has no nonlocal rows, set the
      // "min" to the max possible GO value.  This ensures that if
      // some process has at least one nonlocal row, then it will pick
      // that up as the min.  We know that at least one process has a
      // nonlocal row, since the all-reduce and return at the top of
      // this method excluded that case.
      GO myMinNonlocalGblRow = std::numeric_limits<GO>::max ();
      {
        auto iter = std::min_element (myNonlocalGblRows.begin (),
                                      myNonlocalGblRows.end ());
        if (iter != myNonlocalGblRows.end ()) {
          myMinNonlocalGblRow = *iter;
        }
      }
      GO gblMinNonlocalGblRow = 0;
      reduceAll<int, GO> (*comm, REDUCE_MIN, myMinNonlocalGblRow,
                          outArg (gblMinNonlocalGblRow));
      const GO indexBase = gblMinNonlocalGblRow;
      const global_size_t INV = Teuchos::OrdinalTraits<global_size_t>::invalid ();
      nonlocalRowMap = rcp (new map_type (INV, myNonlocalGblRows (), indexBase, comm));
    }

    if (verbose_) {
      std::ostringstream os;
      os << *prefix << "nonlocalRowMap->getIndexBase()="
         << nonlocalRowMap->getIndexBase() << endl;
      std::cerr << os.str();
    }

    // 3. Use the column indices for each nonlocal row, as stored in
    //    nonlocals_, to construct a CrsGraph corresponding to
    //    nonlocal rows.  We need, but we have, exact counts of the
    //    number of entries in each nonlocal row.

    RCP<crs_graph_type> nonlocalGraph =
      rcp(new crs_graph_type(nonlocalRowMap, numEntPerNonlocalRow(),
                             StaticProfile));
    {
      size_type curPos = 0;
      for (auto mapIter = this->nonlocals_.begin ();
           mapIter != this->nonlocals_.end ();
           ++mapIter, ++curPos) {
        const GO gblRow = mapIter->first;
        std::vector<GO>& gblCols = mapIter->second; // by ref just to avoid copy
        const LO numEnt = static_cast<LO> (numEntPerNonlocalRow[curPos]);
        nonlocalGraph->insertGlobalIndices (gblRow, numEnt, gblCols.data ());
      }
    }
    if (verbose_) {
      std::ostringstream os;
      os << *prefix << "Built nonlocal graph" << endl;
      std::cerr << os.str();
    }
    // There's no need to fill-complete the nonlocals graph.
    // We just use it as a temporary container for the Export.

    // 4. If the original row Map is one to one, then we can Export
    //    directly from nonlocalGraph into this.  Otherwise, we have
    //    to create a temporary graph with a one-to-one row Map,
    //    Export into that, then Import from the temporary graph into
    //    *this.

    auto origRowMap = this->getRowMap ();
    const bool origRowMapIsOneToOne = origRowMap->isOneToOne ();

    if (origRowMapIsOneToOne) {
      if (verbose_) {
        std::ostringstream os;
        os << *prefix << "Original row Map is 1-to-1" << endl;
        std::cerr << os.str();
      }
      export_type exportToOrig (nonlocalRowMap, origRowMap);
      this->doExport (*nonlocalGraph, exportToOrig, Tpetra::INSERT);
      // We're done at this point!
    }
    else {
      if (verbose_) {
        std::ostringstream os;
        os << *prefix << "Original row Map is NOT 1-to-1" << endl;
        std::cerr << os.str();
      }
      // If you ask a Map whether it is one to one, it does some
      // communication and stashes intermediate results for later use
      // by createOneToOne.  Thus, calling createOneToOne doesn't cost
      // much more then the original cost of calling isOneToOne.
      auto oneToOneRowMap = Tpetra::createOneToOne (origRowMap);
      export_type exportToOneToOne (nonlocalRowMap, oneToOneRowMap);

      // Create a temporary graph with the one-to-one row Map.
      //
      // TODO (mfh 09 Sep 2016) Estimate the number of entries in each
      // row, to avoid reallocation during the Export operation.
      crs_graph_type oneToOneGraph (oneToOneRowMap, 0);

      // Export from graph of nonlocals into the temp one-to-one graph.
      if (verbose_) {
        std::ostringstream os;
        os << *prefix << "Export nonlocal graph" << endl;
        std::cerr << os.str();
      }
      oneToOneGraph.doExport (*nonlocalGraph, exportToOneToOne, Tpetra::INSERT);

      // We don't need the graph of nonlocals anymore, so get rid of
      // it, to keep the memory high-water mark down.
      nonlocalGraph = Teuchos::null;

      // Import from the one-to-one graph to the original graph.
      import_type importToOrig (oneToOneRowMap, origRowMap);
      if (verbose_) {
        std::ostringstream os;
        os << *prefix << "Import nonlocal graph" << endl;
        std::cerr << os.str();
      }
      this->doImport (oneToOneGraph, importToOrig, Tpetra::INSERT);
    }

    // It's safe now to clear out nonlocals_, since we've already
    // committed side effects to *this.  The standard idiom for
    // clearing a Container like std::map, is to swap it with an empty
    // Container and let the swapped Container fall out of scope.
    decltype (this->nonlocals_) newNonlocals;
    std::swap (this->nonlocals_, newNonlocals);

    checkInternalState ();
    if (verbose_) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str();
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  resumeFill (const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    clearGlobalConstants();
    if (params != Teuchos::null) this->setParameterList (params);
    // either still sorted/merged or initially sorted/merged
    indicesAreSorted_ = true;
    noRedundancies_ = true;
    fillComplete_ = false;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  fillComplete (const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    // If the graph already has domain and range Maps, don't clobber
    // them.  If it doesn't, use the current row Map for both the
    // domain and range Maps.
    //
    // NOTE (mfh 28 Sep 2014): If the graph was constructed without a
    // column Map, and column indices are inserted which are not in
    // the row Map on any process, this will cause troubles.  However,
    // that is not a common case for most applications that we
    // encounter, and checking for it might require more
    // communication.
    Teuchos::RCP<const map_type> domMap = this->getDomainMap ();
    if (domMap.is_null ()) {
      domMap = this->getRowMap ();
    }
    Teuchos::RCP<const map_type> ranMap = this->getRangeMap ();
    if (ranMap.is_null ()) {
      ranMap = this->getRowMap ();
    }
    this->fillComplete (domMap, ranMap, params);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  fillComplete (const Teuchos::RCP<const map_type>& domainMap,
                const Teuchos::RCP<const map_type>& rangeMap,
                const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    using std::endl;
    const char tfecfFuncName[] = "fillComplete: ";
    const bool verbose = verbose_;

    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsGraph", "fillComplete");
      std::ostringstream os;
      os << *prefix << "Start" << endl;
      std::cerr << os.str();
    }

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! isFillActive () || isFillComplete (), std::runtime_error,
       "Graph fill state must be active (isFillActive() "
       "must be true) before calling fillComplete().");

    const int numProcs = getComm ()->getSize ();

    //
    // Read and set parameters
    //

    // Does the caller want to sort remote GIDs (within those owned by
    // the same process) in makeColMap()?
    if (! params.is_null ()) {
      if (params->isParameter ("sort column map ghost gids")) {
        sortGhostsAssociatedWithEachProcessor_ =
          params->get<bool> ("sort column map ghost gids",
                             sortGhostsAssociatedWithEachProcessor_);
      }
      else if (params->isParameter ("Sort column Map ghost GIDs")) {
        sortGhostsAssociatedWithEachProcessor_ =
          params->get<bool> ("Sort column Map ghost GIDs",
                             sortGhostsAssociatedWithEachProcessor_);
      }
    }

    // If true, the caller promises that no process did nonlocal
    // changes since the last call to fillComplete.
    bool assertNoNonlocalInserts = false;
    if (! params.is_null ()) {
      assertNoNonlocalInserts =
        params->get<bool> ("No Nonlocal Changes", assertNoNonlocalInserts);
    }

    //
    // Allocate indices, if they haven't already been allocated
    //
    if (! indicesAreAllocated ()) {
      if (hasColMap ()) {
        // We have a column Map, so use local indices.
        allocateIndices (LocalIndices, verbose);
      } else {
        // We don't have a column Map, so use global indices.
        allocateIndices (GlobalIndices, verbose);
      }
    }

    //
    // Do global assembly, if requested and if the communicator
    // contains more than one process.
    //
    const bool mayNeedGlobalAssemble = ! assertNoNonlocalInserts && numProcs > 1;
    if (mayNeedGlobalAssemble) {
      // This first checks if we need to do global assembly.
      // The check costs a single all-reduce.
      globalAssemble ();
    }
    else {
      const size_t numNonlocals = nonlocals_.size();
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Do not need to call globalAssemble; "
          "assertNoNonlocalInserts="
           << (assertNoNonlocalInserts ? "true" : "false")
           << "numProcs=" << numProcs
           << ", nonlocals_.size()=" << numNonlocals << endl;
        std::cerr << os.str();
      }
      const int lclNeededGlobalAssemble =
        (numProcs > 1 && numNonlocals != 0) ? 1 : 0;
      if (lclNeededGlobalAssemble != 0 && verbose) {
        std::ostringstream os;
        os << *prefix;
        Details::Impl::verbosePrintMap(
          os, nonlocals_.begin(), nonlocals_.end(),
          nonlocals_.size(), "nonlocals_");
        std::cerr << os.str() << endl;
      }

      if (debug_) {
        auto map = this->getMap();
        auto comm = map.is_null() ? Teuchos::null : map->getComm();
        int gblNeededGlobalAssemble = lclNeededGlobalAssemble;
        if (! comm.is_null()) {
          using Teuchos::REDUCE_MAX;
          using Teuchos::reduceAll;
          reduceAll(*comm, REDUCE_MAX, lclNeededGlobalAssemble,
                    Teuchos::outArg(gblNeededGlobalAssemble));
        }
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (gblNeededGlobalAssemble != 0, std::runtime_error,
           "nonlocals_.size()=" << numNonlocals << " != 0 on at "
           "least one process in the CrsGraph's communicator.  This "
           "means either that you incorrectly set the "
           "\"No Nonlocal Changes\" fillComplete parameter to true, "
           "or that you inserted invalid entries.  "
           "Rerun with the environment variable TPETRA_VERBOSE="
           "CrsGraph set to see the entries of nonlocals_ on every "
           "MPI process (WARNING: lots of output).");
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (lclNeededGlobalAssemble != 0, std::runtime_error,
           "nonlocals_.size()=" << numNonlocals << " != 0 on the "
           "calling process.  This means either that you incorrectly "
           "set the \"No Nonlocal Changes\" fillComplete parameter "
           "to true, or that you inserted invalid entries.  "
           "Rerun with the environment "
           "variable TPETRA_VERBOSE=CrsGraph set to see the entries "
           "of nonlocals_ on every MPI process (WARNING: lots of "
           "output).");
      }
    }

    // Set domain and range Map.  This may clear the Import / Export
    // objects if the new Maps differ from any old ones.
    setDomainRangeMaps (domainMap, rangeMap);

    // If the graph does not already have a column Map (either from
    // the user constructor calling the version of the constructor
    // that takes a column Map, or from a previous fillComplete call),
    // then create it.
    Teuchos::Array<int> remotePIDs (0);
    const bool mustBuildColMap = ! this->hasColMap ();
    if (mustBuildColMap) {
      this->makeColMap (remotePIDs); // resized on output
    }

    // Make indices local, if they aren't already.
    // The method doesn't do any work if the indices are already local.
    const std::pair<size_t, std::string> makeIndicesLocalResult =
      this->makeIndicesLocal(verbose);
    if (debug_) {
      using Details::gathervPrint;
      using Teuchos::RCP;
      using Teuchos::REDUCE_MIN;
      using Teuchos::reduceAll;
      using Teuchos::outArg;

      RCP<const map_type> map = this->getMap ();
      RCP<const Teuchos::Comm<int> > comm;
      if (! map.is_null ()) {
        comm = map->getComm ();
      }
      if (comm.is_null ()) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (makeIndicesLocalResult.first != 0, std::runtime_error,
           makeIndicesLocalResult.second);
      }
      else {
        const int lclSuccess = (makeIndicesLocalResult.first == 0);
        int gblSuccess = 0; // output argument
        reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
        if (gblSuccess != 1) {
          std::ostringstream os;
          gathervPrint (os, makeIndicesLocalResult.second, *comm);
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (true, std::runtime_error, os.str ());
        }
      }
    }
    else {
      // TODO (mfh 20 Jul 2017) Instead of throwing here, pass along
      // the error state to makeImportExport or
      // computeGlobalConstants, which may do all-reduces and thus may
      // have the opportunity to communicate that error state.
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (makeIndicesLocalResult.first != 0, std::runtime_error,
         makeIndicesLocalResult.second);
    }

    // If this process has no indices, then CrsGraph considers it
    // already trivially sorted and merged.  Thus, this method need
    // not be called on all processes in the row Map's communicator.
    this->sortAndMergeAllIndices (this->isSorted (), this->isMerged ());

    // Make Import and Export objects, if they haven't been made
    // already.  If we made a column Map above, reuse information from
    // that process to avoid communiation in the Import setup.
    this->makeImportExport (remotePIDs, mustBuildColMap);

    // Create the Kokkos::StaticCrsGraph, if it doesn't already exist.
    this->fillLocalGraph (params);

    const bool callComputeGlobalConstants = params.get () == nullptr ||
      params->get ("compute global constants", true);
    if (callComputeGlobalConstants) {
      this->computeGlobalConstants ();
    }
    else {
      this->computeLocalConstants ();
    }
    this->fillComplete_ = true;
    this->checkInternalState ();

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str();
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  expertStaticFillComplete (const Teuchos::RCP<const map_type>& domainMap,
                            const Teuchos::RCP<const map_type>& rangeMap,
                            const Teuchos::RCP<const import_type>& importer,
                            const Teuchos::RCP<const export_type>& exporter,
                            const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    const char tfecfFuncName[] = "expertStaticFillComplete: ";
#ifdef HAVE_TPETRA_MMM_TIMINGS
    std::string label;
    if(!params.is_null())
      label = params->get("Timer Label",label);
    std::string prefix = std::string("Tpetra ")+ label + std::string(": ");
    using Teuchos::TimeMonitor;
    Teuchos::RCP<Teuchos::TimeMonitor> MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("ESFC-G-Setup"))));
#endif


    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      domainMap.is_null () || rangeMap.is_null (),
      std::runtime_error, "The input domain Map and range Map must be nonnull.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isFillComplete () || ! hasColMap (), std::runtime_error, "You may not "
      "call this method unless the graph has a column Map.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      getNodeNumRows () > 0 && k_rowPtrs_.extent (0) == 0,
      std::runtime_error, "The calling process has getNodeNumRows() = "
      << getNodeNumRows () << " > 0 rows, but the row offsets array has not "
      "been set.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (k_rowPtrs_.extent (0)) != getNodeNumRows () + 1,
      std::runtime_error, "The row offsets array has length " <<
      k_rowPtrs_.extent (0) << " != getNodeNumRows()+1 = " <<
      (getNodeNumRows () + 1) << ".");

    // Note: We don't need to do the following things which are normally done in fillComplete:
    // allocateIndices, globalAssemble, makeColMap, makeIndicesLocal, sortAndMergeAllIndices

    // Constants from allocateIndices
    //
    // mfh 08 Aug 2014: numAllocForAllRows_ and k_numAllocPerRow_ go
    // away once the graph is allocated.  expertStaticFillComplete
    // either presumes that the graph is allocated, or "allocates" it.
    //
    // FIXME (mfh 08 Aug 2014) The goal for the Kokkos refactor
    // version of CrsGraph is to allocate in the constructor, not
    // lazily on first insert.  That will make both
    // numAllocForAllRows_ and k_numAllocPerRow_ obsolete.
    numAllocForAllRows_  = 0;
    k_numAllocPerRow_    = decltype (k_numAllocPerRow_) ();
    indicesAreAllocated_ = true;

    // Constants from makeIndicesLocal
    //
    // The graph has a column Map, so its indices had better be local.
    indicesAreLocal_  = true;
    indicesAreGlobal_ = false;

    // set domain/range map: may clear the import/export objects
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("ESFC-G-Maps"))));
#endif
    setDomainRangeMaps (domainMap, rangeMap);

    // Presume the user sorted and merged the arrays first
    indicesAreSorted_ = true;
    noRedundancies_ = true;

    // makeImportExport won't create a new importer/exporter if I set one here first.
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("ESFC-G-mIXcheckI"))));
#endif

    importer_ = Teuchos::null;
    exporter_ = Teuchos::null;
    if (importer != Teuchos::null) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        ! importer->getSourceMap ()->isSameAs (*getDomainMap ()) ||
        ! importer->getTargetMap ()->isSameAs (*getColMap ()),
        std::invalid_argument,": importer does not match matrix maps.");
      importer_ = importer;

    }

#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("ESFC-G-mIXcheckE"))));
#endif

    if (exporter != Teuchos::null) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        ! exporter->getSourceMap ()->isSameAs (*getRowMap ()) ||
        ! exporter->getTargetMap ()->isSameAs (*getRangeMap ()),
        std::invalid_argument,": exporter does not match matrix maps.");
      exporter_ = exporter;
    }

#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("ESFC-G-mIXmake"))));
#endif
    Teuchos::Array<int> remotePIDs (0); // unused output argument
    this->makeImportExport (remotePIDs, false);

#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("ESFC-G-fLG"))));
#endif
    this->fillLocalGraph (params);

    const bool callComputeGlobalConstants = params.get () == nullptr ||
      params->get ("compute global constants", true);

    if (callComputeGlobalConstants) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("ESFC-G-cGC (const)"))));
#endif // HAVE_TPETRA_MMM_TIMINGS
      this->computeGlobalConstants ();
    }
    else {
#ifdef HAVE_TPETRA_MMM_TIMINGS
      MM = Teuchos::null;
      MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("ESFC-G-cGC (noconst)"))));
#endif // HAVE_TPETRA_MMM_TIMINGS
      this->computeLocalConstants ();
    }

    fillComplete_ = true;

#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("ESFC-G-cIS"))));
#endif
    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  fillLocalGraph (const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    using ::Tpetra::Details::computeOffsetsFromCounts;
    typedef decltype (k_numRowEntries_) row_entries_type;
    typedef typename local_graph_type::row_map_type row_map_type;
    typedef typename row_map_type::non_const_type non_const_row_map_type;
    typedef typename local_graph_type::entries_type::non_const_type lclinds_1d_type;
    const char tfecfFuncName[] = "fillLocalGraph (called from fillComplete or "
      "expertStaticFillComplete): ";
    const size_t lclNumRows = this->getNodeNumRows ();

    // This method's goal is to fill in the two arrays (compressed
    // sparse row format) that define the sparse graph's structure.
    //
    // Use the nonconst version of row_map_type for ptr_d, because
    // the latter is const and we need to modify ptr_d here.
    non_const_row_map_type ptr_d;
    row_map_type ptr_d_const;
    lclinds_1d_type ind_d;

    bool requestOptimizedStorage = true;
    if (! params.is_null () && ! params->get ("Optimize Storage", true)) {
      requestOptimizedStorage = false;
    }

    // The graph's column indices are currently stored in a 1-D
    // format, with row offsets in k_rowPtrs_ and local column indices
    // in k_lclInds1D_.

    if (debug_) {
      // The graph's array of row offsets must already be allocated.
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (k_rowPtrs_.extent (0) == 0, std::logic_error,
         "k_rowPtrs_ has size zero, but shouldn't");
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (k_rowPtrs_.extent (0) != lclNumRows + 1, std::logic_error,
         "k_rowPtrs_.extent(0) = "
         << k_rowPtrs_.extent (0) << " != (lclNumRows + 1) = "
         << (lclNumRows + 1) << ".");
      const size_t numOffsets = k_rowPtrs_.extent (0);
      const auto valToCheck =
        ::Tpetra::Details::getEntryOnHost (k_rowPtrs_, numOffsets - 1);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (numOffsets != 0 &&
         k_lclInds1D_.extent (0) != valToCheck,
         std::logic_error, "numOffsets=" << numOffsets << " != 0 "
         " and k_lclInds1D_.extent(0)=" << k_lclInds1D_.extent(0)
         << " != k_rowPtrs_(" << numOffsets << ")=" << valToCheck
         << ".");
    }

    size_t allocSize = 0;
    try {
      allocSize = this->getNodeAllocationSize ();
    }
    catch (std::logic_error& e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::logic_error, "getNodeAllocationSize threw "
         "std::logic_error: " << e.what ());
    }
    catch (std::runtime_error& e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::runtime_error, "getNodeAllocationSize threw "
         "std::runtime_error: " << e.what ());
    }
    catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::runtime_error, "getNodeAllocationSize threw "
         "std::exception: " << e.what ());
    }
    catch (...) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (true, std::runtime_error, "getNodeAllocationSize threw "
         "an exception not a subclass of std::exception.");
    }

    if (this->getNodeNumEntries () != allocSize) {
      // The graph's current 1-D storage is "unpacked."  This means
      // the row offsets may differ from what the final row offsets
      // should be.  This could happen, for example, if the user set
      // an upper bound on the number of entries in each row, but
      // didn't fill all those entries.

      if (debug_) {
        if (k_rowPtrs_.extent (0) != 0) {
          const size_t numOffsets =
            static_cast<size_t> (k_rowPtrs_.extent (0));
          const auto valToCheck =
            ::Tpetra::Details::getEntryOnHost (k_rowPtrs_, numOffsets - 1);
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (valToCheck != size_t(k_lclInds1D_.extent(0)),
             std::logic_error, "(Unpacked branch) Before allocating "
             "or packing, k_rowPtrs_(" << (numOffsets-1) << ")="
             << valToCheck << " != k_lclInds1D_.extent(0)="
             << k_lclInds1D_.extent (0) << ".");
        }
      }

      // Pack the row offsets into ptr_d, by doing a sum-scan of the
      // array of valid entry counts per row (k_numRowEntries_).

      // Total number of entries in the matrix on the calling
      // process.  We will compute this in the loop below.  It's
      // cheap to compute and useful as a sanity check.
      size_t lclTotalNumEntries = 0;
      {
        // Allocate the packed row offsets array.
        ptr_d = non_const_row_map_type ("Tpetra::CrsGraph::ptr", lclNumRows + 1);
        ptr_d_const = ptr_d;

        // It's ok that k_numRowEntries_ is a host View; the
        // function can handle this.
        typename row_entries_type::const_type numRowEnt_h = k_numRowEntries_;
        if (debug_) {
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (size_t(numRowEnt_h.extent (0)) != lclNumRows,
             std::logic_error, "(Unpacked branch) "
             "numRowEnt_h.extent(0)=" << numRowEnt_h.extent(0)
             << " != getNodeNumRows()=" << lclNumRows << "");
        }

        lclTotalNumEntries = computeOffsetsFromCounts (ptr_d, numRowEnt_h);

        if (debug_) {
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (static_cast<size_t> (ptr_d.extent (0)) != lclNumRows + 1,
             std::logic_error, "(Unpacked branch) After allocating "
             "ptr_d, ptr_d.extent(0) = " << ptr_d.extent(0)
             << " != lclNumRows+1 = " << (lclNumRows+1) << ".");
          const auto valToCheck =
            ::Tpetra::Details::getEntryOnHost (ptr_d, lclNumRows);
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (valToCheck != lclTotalNumEntries, std::logic_error,
             "Tpetra::CrsGraph::fillLocalGraph: In unpacked branch, "
             "after filling ptr_d, ptr_d(lclNumRows=" << lclNumRows
             << ") = " << valToCheck << " != total number of entries "
             "on the calling process = " << lclTotalNumEntries
             << ".");
        }
      }

      // Allocate the array of packed column indices.
      ind_d = lclinds_1d_type ("Tpetra::CrsGraph::ind", lclTotalNumEntries);

      // k_rowPtrs_ and k_lclInds1D_ are currently unpacked.  Pack
      // them, using the packed row offsets array ptr_d that we
      // created above.
      //
      // FIXME (mfh 08 Aug 2014) If "Optimize Storage" is false (in
      // CrsMatrix?), we need to keep around the unpacked row
      // offsets and column indices.

      // Pack the column indices from unpacked k_lclInds1D_ into
      // packed ind_d.  We will replace k_lclInds1D_ below.
      typedef pack_functor<
        typename local_graph_type::entries_type::non_const_type,
        row_map_type> inds_packer_type;
      inds_packer_type f (ind_d, k_lclInds1D_, ptr_d, k_rowPtrs_);
      {
        typedef typename decltype (ind_d)::execution_space exec_space;
        typedef Kokkos::RangePolicy<exec_space, LocalOrdinal> range_type;
        Kokkos::parallel_for (range_type (0, lclNumRows), f);
      }

      if (debug_) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (ptr_d.extent (0) == 0, std::logic_error,
           "(\"Optimize Storage\"=true branch) After packing, "
           "ptr_d.extent(0)=0.  This probably means k_rowPtrs_ was "
           "never allocated.");
        if (ptr_d.extent (0) != 0) {
          const size_t numOffsets = static_cast<size_t> (ptr_d.extent (0));
          const auto valToCheck =
            ::Tpetra::Details::getEntryOnHost (ptr_d, numOffsets - 1);
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (static_cast<size_t> (valToCheck) != ind_d.extent (0),
             std::logic_error, "(\"Optimize Storage\"=true branch) "
             "After packing, ptr_d(" << (numOffsets-1) << ")="
             << valToCheck << " != ind_d.extent(0)="
             << ind_d.extent(0) << ".");
        }
      }
    }
    else { // We don't have to pack, so just set the pointers.
      ptr_d_const = k_rowPtrs_;
      ind_d = k_lclInds1D_;

      if (debug_) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (ptr_d_const.extent (0) == 0, std::logic_error,
           "(\"Optimize Storage\"=false branch) "
           "ptr_d_const.extent(0) = 0.  This probably means that "
           "k_rowPtrs_ was never allocated.");
        if (ptr_d_const.extent (0) != 0) {
          const size_t numOffsets =
            static_cast<size_t> (ptr_d_const.extent (0));
          const size_t valToCheck =
            ::Tpetra::Details::getEntryOnHost (ptr_d_const, numOffsets - 1);
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (valToCheck != size_t(ind_d.extent (0)),
             std::logic_error, "(\"Optimize Storage\"=false branch) "
             "ptr_d_const(" << (numOffsets-1) << ")=" << valToCheck
             << " != ind_d.extent(0)=" << ind_d.extent (0) << ".");
        }
      }
    }

    if (debug_) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
        (static_cast<size_t> (ptr_d_const.extent (0)) != lclNumRows + 1,
         std::logic_error, "After packing, ptr_d_const.extent(0) = " <<
         ptr_d_const.extent (0) << " != lclNumRows+1 = " << (lclNumRows+1)
         << ".");
      if (ptr_d_const.extent (0) != 0) {
        const size_t numOffsets = static_cast<size_t> (ptr_d_const.extent (0));
        const auto valToCheck =
          ::Tpetra::Details::getEntryOnHost (ptr_d_const, numOffsets - 1);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (static_cast<size_t> (valToCheck) != ind_d.extent (0),
           std::logic_error, "After packing, ptr_d_const(" << (numOffsets-1)
           << ") = " << valToCheck << " != ind_d.extent(0) = "
           << ind_d.extent (0) << ".");
      }
    }

    if (requestOptimizedStorage) {
      // With optimized storage, we don't need to store
      // the array of row entry counts.

      // Free graph data structures that are only needed for
      // unpacked 1-D storage.
      k_numRowEntries_ = row_entries_type ();

      // Keep the new 1-D packed allocations.
      k_rowPtrs_   = ptr_d_const;
      k_lclInds1D_ = ind_d;

      storageStatus_ = Details::STORAGE_1D_PACKED;
    }

    // FIXME (mfh 28 Aug 2014) "Local Graph" sublist no longer used.

    // Build the local graph.
    lclGraph_ = local_graph_type (ind_d, ptr_d_const);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  replaceColMap (const Teuchos::RCP<const map_type>& newColMap)
  {
    // NOTE: This safety check matches the code, but not the documentation of Crsgraph
    //
    // FIXME (mfh 18 Aug 2014) This will break if the calling process
    // has no entries, because in that case, currently it is neither
    // locally nor globally indexed.  This will change once we get rid
    // of lazy allocation (so that the constructor allocates indices
    // and therefore commits to local vs. global).
    const char tfecfFuncName[] = "replaceColMap: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isLocallyIndexed () || isGloballyIndexed (), std::runtime_error,
      "Requires matching maps and non-static graph.");
    colMap_ = newColMap;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  reindexColumns (const Teuchos::RCP<const map_type>& newColMap,
                  const Teuchos::RCP<const import_type>& newImport,
                  const bool sortIndicesInEachRow)
  {
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    using Teuchos::RCP;
    typedef GlobalOrdinal GO;
    typedef LocalOrdinal LO;
    typedef typename local_graph_type::entries_type::non_const_type col_inds_type;
    const char tfecfFuncName[] = "reindexColumns: ";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isFillComplete (), std::runtime_error, "The graph is fill complete "
      "(isFillComplete() returns true).  You must call resumeFill() before "
      "you may call this method.");

    // mfh 19 Aug 2014: This method does NOT redistribute data; it
    // doesn't claim to do the work of an Import or Export.  This
    // means that for all processes, the calling process MUST own all
    // column indices, in both the old column Map (if it exists) and
    // the new column Map.  We check this via an all-reduce.
    //
    // Some processes may be globally indexed, others may be locally
    // indexed, and others (that have no graph entries) may be
    // neither.  This method will NOT change the graph's current
    // state.  If it's locally indexed, it will stay that way, and
    // vice versa.  It would easy to add an option to convert indices
    // from global to local, so as to save a global-to-local
    // conversion pass.  However, we don't do this here.  The intended
    // typical use case is that the graph already has a column Map and
    // is locally indexed, and this is the case for which we optimize.

    const LO lclNumRows = static_cast<LO> (this->getNodeNumRows ());

    // Attempt to convert indices to the new column Map's version of
    // local.  This will fail if on the calling process, the graph has
    // indices that are not on that process in the new column Map.
    // After the local conversion attempt, we will do an all-reduce to
    // see if any processes failed.

    // If this is false, then either the graph contains a column index
    // which is invalid in the CURRENT column Map, or the graph is
    // locally indexed but currently has no column Map.  In either
    // case, there is no way to convert the current local indices into
    // global indices, so that we can convert them into the new column
    // Map's local indices.  It's possible for this to be true on some
    // processes but not others, due to replaceColMap.
    bool allCurColIndsValid = true;
    // On the calling process, are all valid current column indices
    // also in the new column Map on the calling process?  In other
    // words, does local reindexing suffice, or should the user have
    // done an Import or Export instead?
    bool localSuffices = true;

    // Final arrays for the local indices.  We will allocate exactly
    // one of these ONLY if the graph is locally indexed on the
    // calling process, and ONLY if the graph has one or more entries
    // (is not empty) on the calling process.  In that case, we
    // allocate the first (1-D storage) if the graph has a static
    // profile, else we allocate the second (2-D storage).
    typename local_graph_type::entries_type::non_const_type newLclInds1D;
    Teuchos::ArrayRCP<Teuchos::Array<LO> > newLclInds2D;

    // If indices aren't allocated, that means the calling process
    // owns no entries in the graph.  Thus, there is nothing to
    // convert, and it trivially succeeds locally.
    if (indicesAreAllocated ()) {
      if (isLocallyIndexed ()) {
        if (hasColMap ()) { // locally indexed, and currently has a column Map
          const map_type& oldColMap = * (getColMap ());
          // Allocate storage for the new local indices.
          const size_t allocSize = this->getNodeAllocationSize ();
          newLclInds1D = col_inds_type ("Tpetra::CrsGraph::ind", allocSize);
          // Attempt to convert the new indices locally.
          for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
            const RowInfo rowInfo = this->getRowInfo (lclRow);
            const size_t beg = rowInfo.offset1D;
            const size_t end = beg + rowInfo.numEntries;
            for (size_t k = beg; k < end; ++k) {
              // FIXME (mfh 21 Aug 2014) This assumes UVM.  Should
              // use a DualView instead.
              const LO oldLclCol = k_lclInds1D_(k);
              if (oldLclCol == Teuchos::OrdinalTraits<LO>::invalid ()) {
                allCurColIndsValid = false;
                break; // Stop at the first invalid index
              }
              const GO gblCol = oldColMap.getGlobalElement (oldLclCol);

              // The above conversion MUST succeed.  Otherwise, the
              // current local index is invalid, which means that
              // the graph was constructed incorrectly.
              if (gblCol == Teuchos::OrdinalTraits<GO>::invalid ()) {
                allCurColIndsValid = false;
                break; // Stop at the first invalid index
              }
              else {
                const LO newLclCol = newColMap->getLocalElement (gblCol);
                if (newLclCol == Teuchos::OrdinalTraits<LO>::invalid ()) {
                  localSuffices = false;
                  break; // Stop at the first invalid index
                }
                // FIXME (mfh 21 Aug 2014) This assumes UVM.  Should
                // use a DualView instead.
                newLclInds1D(k) = newLclCol;
              }
            } // for each entry in the current row
          } // for each locally owned row
        }
        else { // locally indexed, but no column Map
          // This case is only possible if replaceColMap() was called
          // with a null argument on the calling process.  It's
          // possible, but it means that this method can't possibly
          // succeed, since we have no way of knowing how to convert
          // the current local indices to global indices.
          allCurColIndsValid = false;
        }
      }
      else { // globally indexed
        // If the graph is globally indexed, we don't need to save
        // local indices, but we _do_ need to know whether the current
        // global indices are valid in the new column Map.  We may
        // need to do a getRemoteIndexList call to find this out.
        //
        // In this case, it doesn't matter whether the graph currently
        // has a column Map.  We don't need the old column Map to
        // convert from global indices to the _new_ column Map's local
        // indices.  Furthermore, we can use the same code, whether
        // the graph is static or dynamic profile.

        // Test whether the current global indices are in the new
        // column Map on the calling process.
        for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
          const RowInfo rowInfo = this->getRowInfo (lclRow);
          Teuchos::ArrayView<const GO> oldGblRowView = getGlobalView (rowInfo);
          for (size_t k = 0; k < rowInfo.numEntries; ++k) {
            const GO gblCol = oldGblRowView[k];
            if (! newColMap->isNodeGlobalElement (gblCol)) {
              localSuffices = false;
              break; // Stop at the first invalid index
            }
          } // for each entry in the current row
        } // for each locally owned row
      } // locally or globally indexed
    } // whether indices are allocated

    // Do an all-reduce to check both possible error conditions.
    int lclSuccess[2];
    lclSuccess[0] = allCurColIndsValid ? 1 : 0;
    lclSuccess[1] = localSuffices ? 1 : 0;
    int gblSuccess[2];
    gblSuccess[0] = 0;
    gblSuccess[1] = 0;
    RCP<const Teuchos::Comm<int> > comm =
      getRowMap ().is_null () ? Teuchos::null : getRowMap ()->getComm ();
    if (! comm.is_null ()) {
      reduceAll<int, int> (*comm, REDUCE_MIN, 2, lclSuccess, gblSuccess);
    }

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      gblSuccess[0] == 0, std::runtime_error, "It is not possible to continue."
      "  The most likely reason is that the graph is locally indexed, but the "
      "column Map is missing (null) on some processes, due to a previous call "
      "to replaceColMap().");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      gblSuccess[1] == 0, std::runtime_error, "On some process, the graph "
      "contains column indices that are in the old column Map, but not in the "
      "new column Map (on that process).  This method does NOT redistribute "
      "data; it does not claim to do the work of an Import or Export operation."
      "  This means that for all processess, the calling process MUST own all "
      "column indices, in both the old column Map and the new column Map.  In "
      "this case, you will need to do an Import or Export operation to "
      "redistribute data.");

    // Commit the results.
    if (isLocallyIndexed ()) {
      k_lclInds1D_ = newLclInds1D;
      // We've reindexed, so we don't know if the indices are sorted.
      //
      // FIXME (mfh 17 Sep 2014) It could make sense to check this,
      // since we're already going through all the indices above.  We
      // could also sort each row in place; that way, we would only
      // have to make one pass over the rows.
      indicesAreSorted_ = false;
      if (sortIndicesInEachRow) {
        // NOTE (mfh 17 Sep 2014) The graph must be locally indexed in
        // order to call this method.
        //
        // FIXME (mfh 17 Sep 2014) This violates the strong exception
        // guarantee.  It would be better to sort the new index arrays
        // before committing them.
        const bool sorted = false; // need to resort
        const bool merged = true; // no need to merge, since no dups
        this->sortAndMergeAllIndices (sorted, merged);
      }
    }
    colMap_ = newColMap;

    if (newImport.is_null ()) {
      // FIXME (mfh 19 Aug 2014) Should use the above all-reduce to
      // check whether the input Import is null on any process.
      //
      // If the domain Map hasn't been set yet, we can't compute a new
      // Import object.  Leave it what it is; it should be null, but
      // it doesn't matter.  If the domain Map _has_ been set, then
      // compute a new Import object if necessary.
      if (! domainMap_.is_null ()) {
        if (! domainMap_->isSameAs (* newColMap)) {
          importer_ = Teuchos::rcp (new import_type (domainMap_, newColMap));
        } else {
          importer_ = Teuchos::null; // don't need an Import
        }
      }
    } else {
      // The caller gave us an Import object.  Assume that it's valid.
      importer_ = newImport;
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  replaceDomainMapAndImporter (const Teuchos::RCP<const map_type>& newDomainMap,
                               const Teuchos::RCP<const import_type>& newImporter)
  {
    const char prefix[] = "Tpetra::CrsGraph::replaceDomainMapAndImporter: ";
    TEUCHOS_TEST_FOR_EXCEPTION(
      colMap_.is_null (), std::invalid_argument, prefix << "You may not call "
      "this method unless the graph already has a column Map.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      newDomainMap.is_null (), std::invalid_argument,
      prefix << "The new domain Map must be nonnull.");

    if (debug_) {
      if (newImporter.is_null ()) {
        // It's not a good idea to put expensive operations in a macro
        // clause, even if they are side effect - free, because macros
        // don't promise that they won't evaluate their arguments more
        // than once.  It's polite for them to do so, but not required.
        const bool colSameAsDom = colMap_->isSameAs (*newDomainMap);
        TEUCHOS_TEST_FOR_EXCEPTION
          (colSameAsDom, std::invalid_argument, "If the new Import is null, "
           "then the new domain Map must be the same as the current column Map.");
      }
      else {
        const bool colSameAsTgt =
          colMap_->isSameAs (* (newImporter->getTargetMap ()));
        const bool newDomSameAsSrc =
          newDomainMap->isSameAs (* (newImporter->getSourceMap ()));
        TEUCHOS_TEST_FOR_EXCEPTION
          (! colSameAsTgt || ! newDomSameAsSrc, std::invalid_argument, "If the "
           "new Import is nonnull, then the current column Map must be the same "
           "as the new Import's target Map, and the new domain Map must be the "
           "same as the new Import's source Map.");
      }
    }

    domainMap_ = newDomainMap;
    importer_ = Teuchos::rcp_const_cast<import_type> (newImporter);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  typename CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::local_graph_type
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getLocalGraph () const
  {
    return lclGraph_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  computeGlobalConstants ()
  {
    using ::Tpetra::Details::ProfilingRegion;
    using Teuchos::ArrayView;
    using Teuchos::outArg;
    using Teuchos::reduceAll;
    typedef global_size_t GST;

    ProfilingRegion regionCGC ("Tpetra::CrsGraph::computeGlobalConstants");

    this->computeLocalConstants ();

    // Compute global constants from local constants.  Processes that
    // already have local constants still participate in the
    // all-reduces, using their previously computed values.
    if (! this->haveGlobalConstants_) {
      const Teuchos::Comm<int>& comm = * (this->getComm ());
      // Promote all the nodeNum* and nodeMaxNum* quantities from
      // size_t to global_size_t, when doing the all-reduces for
      // globalNum* / globalMaxNum* results.
      //
      // FIXME (mfh 07 May 2013) Unfortunately, we either have to do
      // this in two all-reduces (one for the sum and the other for
      // the max), or use a custom MPI_Op that combines the sum and
      // the max.  The latter might even be slower than two
      // all-reduces on modern network hardware.  It would also be a
      // good idea to use nonblocking all-reduces (MPI 3), so that we
      // don't have to wait around for the first one to finish before
      // starting the second one.
      GST lcl, gbl;
      lcl = static_cast<GST> (this->getNodeNumEntries ());

      reduceAll<int,GST> (comm, Teuchos::REDUCE_SUM, 1, &lcl, &gbl);
      this->globalNumEntries_ = gbl;

      const GST lclMaxNumRowEnt = static_cast<GST> (this->nodeMaxNumRowEntries_);
      reduceAll<int, GST> (comm, Teuchos::REDUCE_MAX, lclMaxNumRowEnt,
                           outArg (this->globalMaxNumRowEntries_));
      this->haveGlobalConstants_ = true;
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  computeLocalConstants ()
  {
    using ::Tpetra::Details::ProfilingRegion;

    ProfilingRegion regionCLC ("Tpetra::CrsGraph::computeLocalConstants");
    if (this->haveLocalConstants_) {
      return;
    }

    // Reset local properties
    this->nodeMaxNumRowEntries_ =
      Teuchos::OrdinalTraits<size_t>::invalid();

    using LO = local_ordinal_type;

    // KJ: This one is a bit different from the above. Conservatively thinking,
    // we also need the fence here as lclGraph_.row_map is on UVM and it can be
    // still updated. In practice, the local graph construction should be done
    // before this is called. This routine is computeLocalConstants. If we want
    // a better code, we need a flag stating that the local graph is completed
    // and safe to use it without fence.
    // For now, I recommend to put the fence. Defining the state of local
    // object can be improvements in the code.
    execution_space().fence ();

    auto ptr = this->lclGraph_.row_map;
    const LO lclNumRows = ptr.extent(0) == 0 ?
      static_cast<LO> (0) :
      (static_cast<LO> (ptr.extent(0)) - static_cast<LO> (1));

    const LO lclMaxNumRowEnt =
      ::Tpetra::Details::maxDifference ("Tpetra::CrsGraph: nodeMaxNumRowEntries",
                              ptr, lclNumRows);
    this->nodeMaxNumRowEntries_ = static_cast<size_t> (lclMaxNumRowEnt);
    this->haveLocalConstants_ = true;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  std::pair<size_t, std::string>
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  makeIndicesLocal (const bool verbose)
  {
    using Details::ProfilingRegion;
    using Teuchos::arcp;
    using Teuchos::Array;
    using std::endl;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef device_type DT;
    typedef typename local_graph_type::row_map_type::non_const_value_type offset_type;
    typedef decltype (k_numRowEntries_) row_entries_type;
    typedef typename row_entries_type::non_const_value_type num_ent_type;
    typedef typename local_graph_type::entries_type::non_const_type
      lcl_col_inds_type;
    typedef Kokkos::View<GO*, typename lcl_col_inds_type::array_layout,
      device_type> gbl_col_inds_type;
    const char tfecfFuncName[] = "makeIndicesLocal: ";
    ProfilingRegion regionMakeIndicesLocal ("Tpetra::CrsGraph::makeIndicesLocal");

    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsGraph", "makeIndicesLocal");
      std::ostringstream os;
      os << *prefix << "lclNumRows: " << getNodeNumRows() << endl;
      std::cerr << os.str();
    }

    // These are somewhat global properties, so it's safe to have
    // exception checks for them, rather than returning an error code.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! this->hasColMap (), std::logic_error, "The graph does not have a "
       "column Map yet.  This method should never be called in that case.  "
       "Please report this bug to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (this->getColMap ().is_null (), std::logic_error, "The graph claims "
       "that it has a column Map, because hasColMap() returns true.  However, "
       "the result of getColMap() is null.  This should never happen.  Please "
       "report this bug to the Tpetra developers.");

    // Return value 1: The number of column indices (counting
    // duplicates) that could not be converted to local indices,
    // because they were not in the column Map on the calling process.
    size_t lclNumErrs = 0;
    std::ostringstream errStrm; // for return value 2 (error string)

    const LO lclNumRows = static_cast<LO> (this->getNodeNumRows ());
    const map_type& colMap = * (this->getColMap ());

    if (this->isGloballyIndexed () && lclNumRows != 0) {
      // This is a host-accessible View.
      typename row_entries_type::const_type h_numRowEnt =
        this->k_numRowEntries_;

      // Allocate space for local indices.
      // If GO and LO are the same size, we can reuse the existing
      // array of 1-D index storage to convert column indices from
      // GO to LO.  Otherwise, we'll just allocate a new buffer.
      constexpr bool LO_GO_same = std::is_same<LO, GO>::value;
      if (LO_GO_same) {
        // This prevents a build error (illegal assignment) if
        // LO_GO_same is _not_ true.  Only the first branch
        // (returning k_gblInds1D_) should ever get taken.
        k_lclInds1D_ = Kokkos::Impl::if_c<LO_GO_same,
          t_GlobalOrdinal_1D,
          lcl_col_inds_type>::select (k_gblInds1D_, k_lclInds1D_);
      }
      else {
        if (k_rowPtrs_.extent (0) == 0) {
          errStrm << "k_rowPtrs_.extent(0) == 0.  This should never "
            "happen here.  Please report this bug to the Tpetra developers."
            << endl;
          // Need to return early.
          return std::make_pair (Tpetra::Details::OrdinalTraits<size_t>::invalid (),
                                 errStrm.str ());
        }
        const auto numEnt = ::Tpetra::Details::getEntryOnHost (k_rowPtrs_, lclNumRows);

        // mfh 17 Dec 2016: We don't need initial zero-fill of
        // k_lclInds1D_, because we will fill it below anyway.
        // AllowPadding would only help for aligned access (e.g.,
        // for vectorization) if we also were to pad each row to the
        // same alignment, so we'll skip AllowPadding for now.

        // using Kokkos::AllowPadding;
        using Kokkos::view_alloc;
        using Kokkos::WithoutInitializing;

        // When giving the label as an argument to
        // Kokkos::view_alloc, the label must be a string and not a
        // char*, else the code won't compile.  This is because
        // view_alloc also allows a raw pointer as its first
        // argument.  See
        // https://github.com/kokkos/kokkos/issues/434.  This is a
        // large allocation typically, so the overhead of creating
        // an std::string is minor.
        const std::string label ("Tpetra::CrsGraph::lclind");
        if (verbose) {
          std::ostringstream os;
          os << *prefix << "(Re)allocate k_lclInds1D_: old="
             << k_lclInds1D_.extent(0) << ", new=" << numEnt << endl;
          std::cerr << os.str();
        }
        k_lclInds1D_ =
          lcl_col_inds_type (view_alloc (label, WithoutInitializing), numEnt);
      }

      auto lclColMap = colMap.getLocalMap ();
      // This is a "device mirror" of the host View h_numRowEnt.
      //
      // NOTE (mfh 27 Sep 2016) Currently, the right way to get a
      // Device instance is to use its default constructor.  See the
      // following Kokkos issue:
      //
      // https://github.com/kokkos/kokkos/issues/442
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Allocate device mirror k_numRowEnt: "
           << h_numRowEnt.extent(0) << endl;
        std::cerr << os.str();
      }
      auto k_numRowEnt = Kokkos::create_mirror_view (device_type (), h_numRowEnt);

      using ::Tpetra::Details::convertColumnIndicesFromGlobalToLocal;
      lclNumErrs =
        convertColumnIndicesFromGlobalToLocal<LO, GO, DT, offset_type, num_ent_type> (k_lclInds1D_,
                                                                                      k_gblInds1D_,
                                                                                      k_rowPtrs_,
                                                                                      lclColMap,
                                                                                      k_numRowEnt);
      if (lclNumErrs != 0) {
        const int myRank = [this] () {
          auto map = this->getMap ();
          if (map.is_null ()) {
            return 0;
          }
          else {
            auto comm = map->getComm ();
            return comm.is_null () ? 0 : comm->getRank ();
          }
        } ();
        const bool pluralNumErrs = (lclNumErrs != static_cast<size_t> (1));
        errStrm << "(Process " << myRank << ") When converting column "
          "indices from global to local, we encountered " << lclNumErrs
          << " ind" << (pluralNumErrs ? "ices" : "ex")
          << " that do" << (pluralNumErrs ? "es" : "")
          << " not live in the column Map on this process." << endl;
      }

      // We've converted column indices from global to local, so we
      // can deallocate the global column indices (which we know are
      // in 1-D storage, because the graph has static profile).
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Free k_gblInds1D_: "
           << k_gblInds1D_.extent(0) << endl;
        std::cerr << os.str();
      }
      k_gblInds1D_ = gbl_col_inds_type ();
    } // globallyIndexed() && lclNumRows > 0

    this->lclGraph_ = local_graph_type (this->k_lclInds1D_, this->k_rowPtrs_);
    this->indicesAreLocal_  = true;
    this->indicesAreGlobal_ = false;
    this->checkInternalState ();

    return std::make_pair (lclNumErrs, errStrm.str ());
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  makeColMap (Teuchos::Array<int>& remotePIDs)
  {
    using Details::ProfilingRegion;
    using std::endl;
    const char tfecfFuncName[] = "makeColMap";

    ProfilingRegion regionSortAndMerge ("Tpetra::CrsGraph::makeColMap");
    std::unique_ptr<std::string> prefix;
    if (verbose_) {
      prefix = this->createPrefix("CrsGraph", tfecfFuncName);
      std::ostringstream os;
      os << *prefix << "Start" << endl;
      std::cerr << os.str();
    }

    // this->colMap_ should be null at this point, but we accept the
    // future possibility that it might not be (esp. if we decide
    // later to support graph structure changes after first
    // fillComplete, which CrsGraph does not currently (as of 12 Feb
    // 2017) support).
    Teuchos::RCP<const map_type> colMap = this->colMap_;
    const bool sortEachProcsGids =
      this->sortGhostsAssociatedWithEachProcessor_;

    // FIXME (mfh 12 Feb 2017) ::Tpetra::Details::makeColMap returns a
    // per-process error code.  If an error does occur on a process,
    // ::Tpetra::Details::makeColMap does NOT promise that all processes will
    // notice that error.  This is the caller's responsibility.  For
    // now, we only propagate (to all processes) and report the error
    // in debug mode.  In the future, we need to add the local/global
    // error handling scheme used in BlockCrsMatrix to this class.
    if (debug_) {
      using Teuchos::outArg;
      using Teuchos::REDUCE_MIN;
      using Teuchos::reduceAll;

      std::ostringstream errStrm;
      const int lclErrCode =
        Details::makeColMap (colMap, remotePIDs,
          getDomainMap (), *this, sortEachProcsGids, &errStrm);
      auto comm = this->getComm ();
      if (! comm.is_null ()) {
        const int lclSuccess = (lclErrCode == 0) ? 1 : 0;
        int gblSuccess = 0; // output argument
        reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess,
                             outArg (gblSuccess));
        if (gblSuccess != 1) {
          std::ostringstream os;
          Details::gathervPrint (os, errStrm.str (), *comm);
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
            (true, std::runtime_error, ": An error happened on at "
             "least one process in the CrsGraph's communicator.  "
             "Here are all processes' error messages:" << std::endl
             << os.str ());
        }
      }
    }
    else {
      (void) Details::makeColMap (colMap, remotePIDs,
        getDomainMap (), *this, sortEachProcsGids, nullptr);
    }
    // See above.  We want to admit the possibility of makeColMap
    // actually revising an existing column Map, even though that
    // doesn't currently (as of 10 May 2017) happen.
    this->colMap_ = colMap;

    checkInternalState ();
    if (verbose_) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str();
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  sortAndMergeAllIndices (const bool sorted, const bool merged)
  {
    using std::endl;
    using LO = LocalOrdinal;
    using host_execution_space =
      typename Kokkos::View<LO*, device_type>::HostMirror::
        execution_space;
    using range_type = Kokkos::RangePolicy<host_execution_space, LO>;
    const char tfecfFuncName[] = "sortAndMergeAllIndices";
    Details::ProfilingRegion regionSortAndMerge
      ("Tpetra::CrsGraph::sortAndMergeAllIndices");

    std::unique_ptr<std::string> prefix;
    if (verbose_) {
      prefix = this->createPrefix("CrsGraph", tfecfFuncName);
      std::ostringstream os;
      os << *prefix << "Start: "
         << "sorted=" << (sorted ? "true" : "false")
         << ", merged=" << (merged ? "true" : "false") << endl;
      std::cerr << os.str();
    }
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (this->isGloballyIndexed(), std::logic_error,
       "This method may only be called after makeIndicesLocal." );
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! merged && this->isStorageOptimized(), std::logic_error,
       "The graph is already storage optimized, so we shouldn't be "
       "merging any indices.  "
       "Please report this bug to the Tpetra developers.");

    if (! sorted || ! merged) {
      const LO lclNumRows(this->getNodeNumRows());
      auto range = range_type(0, lclNumRows);

      // FIXME (mfh 08 May 2017) Loops below assume CUDA UVM.
      if (verbose_) {
        size_t totalNumDups = 0;
        Kokkos::parallel_reduce(range,
          [this, sorted, merged] (const LO lclRow, size_t& numDups)
          {
            const RowInfo rowInfo = this->getRowInfo(lclRow);
            numDups += this->sortAndMergeRowIndices(rowInfo, sorted, merged);
          },
          totalNumDups);
        std::ostringstream os;
        os << *prefix << "totalNumDups=" << totalNumDups << endl;
        std::cerr << os.str();
      }
      else {
        // FIXME (mfh 08 May 2017) This may assume CUDA UVM.
        Kokkos::parallel_for(range,
          [this, sorted, merged] (const LO lclRow)
          {
            const RowInfo rowInfo = this->getRowInfo(lclRow);
            this->sortAndMergeRowIndices(rowInfo, sorted, merged);
          });
      }
      this->indicesAreSorted_ = true; // we just sorted every row
      this->noRedundancies_ = true; // we just merged every row
    }

    if (verbose_) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  makeImportExport (Teuchos::Array<int>& remotePIDs,
                    const bool useRemotePIDs)
  {
    using ::Tpetra::Details::ProfilingRegion;
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;
    const char tfecfFuncName[] = "makeImportExport: ";
    ProfilingRegion regionMIE ("Tpetra::CrsGraph::makeImportExport");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! this->hasColMap (), std::logic_error,
       "This method may not be called unless the graph has a column Map.");
    RCP<ParameterList> params = this->getNonconstParameterList (); // could be null

    // Don't do any checks to see if we need to create the Import, if
    // it exists already.
    //
    // FIXME (mfh 25 Mar 2013) This will become incorrect if we
    // change CrsGraph in the future to allow changing the column
    // Map after fillComplete.  For now, the column Map is fixed
    // after the first fillComplete call.
    if (importer_.is_null ()) {
      // Create the Import instance if necessary.
      if (domainMap_ != colMap_ && (! domainMap_->isSameAs (*colMap_))) {
        if (params.is_null () || ! params->isSublist ("Import")) {
          if (useRemotePIDs) {
            importer_ = rcp (new import_type (domainMap_, colMap_, remotePIDs));
          }
          else {
            importer_ = rcp (new import_type (domainMap_, colMap_));
          }
        }
        else {
          RCP<ParameterList> importSublist = sublist (params, "Import", true);
          if (useRemotePIDs) {
            RCP<import_type> newImp =
              rcp (new import_type (domainMap_, colMap_, remotePIDs,
                                    importSublist));
            importer_ = newImp;
          }
          else {
            importer_ = rcp (new import_type (domainMap_, colMap_, importSublist));
          }
        }
      }
    }

    // Don't do any checks to see if we need to create the Export, if
    // it exists already.
    if (exporter_.is_null ()) {
      // Create the Export instance if necessary.
      if (rangeMap_ != rowMap_ && ! rangeMap_->isSameAs (*rowMap_)) {
        if (params.is_null () || ! params->isSublist ("Export")) {
          exporter_ = rcp (new export_type (rowMap_, rangeMap_));
        }
        else {
          RCP<ParameterList> exportSublist = sublist (params, "Export", true);
          exporter_ = rcp (new export_type (rowMap_, rangeMap_, exportSublist));
        }
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  description () const
  {
    std::ostringstream oss;
    oss << dist_object_type::description ();
    if (isFillComplete ()) {
      oss << "{status = fill complete"
          << ", global rows = " << getGlobalNumRows()
          << ", global cols = " << getGlobalNumCols()
          << ", global num entries = " << getGlobalNumEntries()
          << "}";
    }
    else {
      oss << "{status = fill not complete"
          << ", global rows = " << getGlobalNumRows()
          << "}";
    }
    return oss.str();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    using Teuchos::ArrayView;
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    using std::endl;
    using std::setw;

    Teuchos::EVerbosityLevel vl = verbLevel;
    if (vl == VERB_DEFAULT) vl = VERB_LOW;
    RCP<const Comm<int> > comm = this->getComm();
    const int myImageID = comm->getRank(),
              numImages = comm->getSize();
    size_t width = 1;
    for (size_t dec=10; dec<getGlobalNumRows(); dec *= 10) {
      ++width;
    }
    width = std::max<size_t> (width, static_cast<size_t> (11)) + 2;
    Teuchos::OSTab tab (out);
    //    none: print nothing
    //     low: print O(1) info from node 0
    //  medium: print O(P) info, num entries per node
    //    high: print O(N) info, num entries per row
    // extreme: print O(NNZ) info: print graph indices
    //
    // for medium and higher, print constituent objects at specified verbLevel
    if (vl != VERB_NONE) {
      if (myImageID == 0) out << this->description() << std::endl;
      // O(1) globals, minus what was already printed by description()
      if (isFillComplete() && myImageID == 0) {
        out << "Global max number of row entries = " << globalMaxNumRowEntries_ << std::endl;
      }
      // constituent objects
      if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
        if (myImageID == 0) out << "\nRow map: " << std::endl;
        rowMap_->describe(out,vl);
        if (colMap_ != Teuchos::null) {
          if (myImageID == 0) out << "\nColumn map: " << std::endl;
          colMap_->describe(out,vl);
        }
        if (domainMap_ != Teuchos::null) {
          if (myImageID == 0) out << "\nDomain map: " << std::endl;
          domainMap_->describe(out,vl);
        }
        if (rangeMap_ != Teuchos::null) {
          if (myImageID == 0) out << "\nRange map: " << std::endl;
          rangeMap_->describe(out,vl);
        }
      }
      // O(P) data
      if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
        for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
          if (myImageID == imageCtr) {
            out << "Node ID = " << imageCtr << std::endl
                << "Node number of entries = " << this->getNodeNumEntries () << std::endl
                << "Node max number of entries = " << nodeMaxNumRowEntries_ << std::endl;
            if (! indicesAreAllocated ()) {
              out << "Indices are not allocated." << std::endl;
            }
          }
          comm->barrier();
          comm->barrier();
          comm->barrier();
        }
      }
      // O(N) and O(NNZ) data
      if (vl == VERB_HIGH || vl == VERB_EXTREME) {
        for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
          if (myImageID == imageCtr) {
            out << std::setw(width) << "Node ID"
                << std::setw(width) << "Global Row"
                << std::setw(width) << "Num Entries";
            if (vl == VERB_EXTREME) {
              out << " Entries";
            }
            out << std::endl;
            const LocalOrdinal lclNumRows =
              static_cast<LocalOrdinal> (this->getNodeNumRows ());
            for (LocalOrdinal r=0; r < lclNumRows; ++r) {
              const RowInfo rowinfo = this->getRowInfo (r);
              GlobalOrdinal gid = rowMap_->getGlobalElement(r);
              out << std::setw(width) << myImageID
                  << std::setw(width) << gid
                  << std::setw(width) << rowinfo.numEntries;
              if (vl == VERB_EXTREME) {
                out << " ";
                if (isGloballyIndexed()) {
                  ArrayView<const GlobalOrdinal> rowview = getGlobalView(rowinfo);
                  for (size_t j=0; j < rowinfo.numEntries; ++j) out << rowview[j] << " ";
                }
                else if (isLocallyIndexed()) {
                  ArrayView<const LocalOrdinal> rowview = getLocalView(rowinfo);
                  for (size_t j=0; j < rowinfo.numEntries; ++j) out << colMap_->getGlobalElement(rowview[j]) << " ";
                }
              }
              out << std::endl;
            }
          }
          comm->barrier();
          comm->barrier();
          comm->barrier();
        }
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  checkSizes (const SrcDistObject& /* source */)
  {
    // It's not clear what kind of compatibility checks on sizes can
    // be performed here.  Epetra_CrsGraph doesn't check any sizes for
    // compatibility.
    return true;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  copyAndPermute
  (const SrcDistObject& source,
   const size_t numSameIDs,
   const Kokkos::DualView<const local_ordinal_type*,
     buffer_device_type>& permuteToLIDs,
   const Kokkos::DualView<const local_ordinal_type*,
     buffer_device_type>& permuteFromLIDs)
  {
    using std::endl;
    using LO = local_ordinal_type;
    using GO = global_ordinal_type;
    using this_type = CrsGraph<LO, GO, node_type>;
    using row_graph_type = RowGraph<LO, GO, node_type>;
    const char tfecfFuncName[] = "copyAndPermute: ";
    const bool verbose = verbose_;

    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsGraph", "copyAndPermute");
      std::ostringstream os;
      os << *prefix << endl;
      std::cerr << os.str ();
    }

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (permuteToLIDs.extent (0) != permuteFromLIDs.extent (0),
       std::runtime_error, "permuteToLIDs.extent(0) = "
       << permuteToLIDs.extent (0) << " != permuteFromLIDs.extent(0) = "
       << permuteFromLIDs.extent (0) << ".");

    // We know from checkSizes that the source object is a
    // row_graph_type, so we don't need to check again.
    const row_graph_type& srcRowGraph =
      dynamic_cast<const row_graph_type&> (source);

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Compute padding" << endl;
      std::cerr << os.str ();
    }
    auto padding = computeCrsPadding(srcRowGraph, numSameIDs,
      permuteToLIDs, permuteFromLIDs, verbose);
    applyCrsPadding(*padding, verbose);

    // If the source object is actually a CrsGraph, we can use view
    // mode instead of copy mode to access the entries in each row,
    // if the graph is not fill complete.
    const this_type* srcCrsGraph =
      dynamic_cast<const this_type*> (&source);

    const map_type& srcRowMap = *(srcRowGraph.getRowMap());
    const map_type& tgtRowMap = *(getRowMap());
    const bool src_filled = srcRowGraph.isFillComplete();
    Teuchos::Array<GO> row_copy;
    LO myid = 0;

    //
    // "Copy" part of "copy and permute."
    //
    if (src_filled || srcCrsGraph == nullptr) {
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "src_filled || srcCrsGraph == nullptr" << endl;
        std::cerr << os.str ();
      }
      // If the source graph is fill complete, we can't use view mode,
      // because the data might be stored in a different format not
      // compatible with the expectations of view mode.  Also, if the
      // source graph is not a CrsGraph, we can't use view mode,
      // because RowGraph only provides copy mode access to the data.
      for (size_t i = 0; i < numSameIDs; ++i, ++myid) {
        const GO gid = srcRowMap.getGlobalElement (myid);
        size_t row_length = srcRowGraph.getNumEntriesInGlobalRow (gid);
        row_copy.resize (row_length);
        size_t check_row_length = 0;
        srcRowGraph.getGlobalRowCopy (gid, row_copy (), check_row_length);
        this->insertGlobalIndices (gid, row_copy ());
      }
    } else {
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "! src_filled && srcCrsGraph != nullptr" << endl;
        std::cerr << os.str ();
      }
      for (size_t i = 0; i < numSameIDs; ++i, ++myid) {
        const GO gid = srcRowMap.getGlobalElement (myid);
        Teuchos::ArrayView<const GO> row;
        srcCrsGraph->getGlobalRowView (gid, row);
        this->insertGlobalIndices (gid, row);
      }
    }

    //
    // "Permute" part of "copy and permute."
    //
    auto permuteToLIDs_h = permuteToLIDs.view_host ();
    auto permuteFromLIDs_h = permuteFromLIDs.view_host ();

    if (src_filled || srcCrsGraph == nullptr) {
      for (LO i = 0; i < static_cast<LO> (permuteToLIDs_h.extent (0)); ++i) {
        const GO mygid = tgtRowMap.getGlobalElement (permuteToLIDs_h[i]);
        const GO srcgid = srcRowMap.getGlobalElement (permuteFromLIDs_h[i]);
        size_t row_length = srcRowGraph.getNumEntriesInGlobalRow (srcgid);
        row_copy.resize (row_length);
        size_t check_row_length = 0;
        srcRowGraph.getGlobalRowCopy (srcgid, row_copy (), check_row_length);
        this->insertGlobalIndices (mygid, row_copy ());
      }
    } else {
      for (LO i = 0; i < static_cast<LO> (permuteToLIDs_h.extent (0)); ++i) {
        const GO mygid = tgtRowMap.getGlobalElement (permuteToLIDs_h[i]);
        const GO srcgid = srcRowMap.getGlobalElement (permuteFromLIDs_h[i]);
        Teuchos::ArrayView<const GO> row;
        srcCrsGraph->getGlobalRowView (srcgid, row);
        this->insertGlobalIndices (mygid, row);
      }
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str ();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  applyCrsPadding(const padding_type& padding,
                  const bool verbose)
  {
    using Details::ProfilingRegion;
    using Details::padCrsArrays;
    using std::endl;
    using LO = local_ordinal_type;
    using execution_space = typename device_type::execution_space;
    using row_ptrs_type =
      typename local_graph_type::row_map_type::non_const_type;
    using range_policy =
      Kokkos::RangePolicy<execution_space, Kokkos::IndexType<LO>>;
    const char tfecfFuncName[] = "applyCrsPadding";
    ProfilingRegion regionCAP("Tpetra::CrsGraph::applyCrsPadding");

    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsGraph", tfecfFuncName);
      std::ostringstream os;
      os << *prefix << "padding: ";
      padding.print(os);
      os << endl;
      std::cerr << os.str();
    }
    const int myRank = ! verbose ? -1 : [&] () {
      auto map = this->getMap();
      if (map.is_null()) {
        return -1;
      }
      auto comm = map->getComm();
      if (comm.is_null()) {
        return -1;
      }
      return comm->getRank();
    } ();

    // FIXME (mfh 10 Feb 2020) We shouldn't actually reallocate
    // row_ptrs_beg or allocate row_ptrs_end unless the allocation
    // size needs to increase.  That should be the job of
    // padCrsArrays.

    // Assume global indexing we don't have any indices yet
    if (! indicesAreAllocated()) {
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Call allocateIndices" << endl;
        std::cerr << os.str();
      }
      allocateIndices(GlobalIndices, verbose);
    }
    TEUCHOS_ASSERT( indicesAreAllocated() );

    // Making copies here because k_rowPtrs_ has a const type. Otherwise, we
    // would use it directly.

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Allocate row_ptrs_beg: "
         << k_rowPtrs_.extent(0) << endl;
      std::cerr << os.str();
    }
    using Kokkos::view_alloc;
    using Kokkos::WithoutInitializing;
    row_ptrs_type row_ptrs_beg(
      view_alloc("row_ptrs_beg", WithoutInitializing),
      k_rowPtrs_.extent(0));
    Kokkos::deep_copy(row_ptrs_beg, k_rowPtrs_);

    const size_t N = row_ptrs_beg.extent(0) == 0 ? size_t(0) :
      size_t(row_ptrs_beg.extent(0) - 1);
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Allocate row_ptrs_end: " << N << endl;
      std::cerr << os.str();
    }
    row_ptrs_type row_ptrs_end(
      view_alloc("row_ptrs_end", WithoutInitializing), N);

    const bool refill_num_row_entries = k_numRowEntries_.extent(0) != 0;
    if (refill_num_row_entries) { // Case 1: Unpacked storage
      // We can't assume correct *this capture until C++17, and it's
      // likely more efficient just to capture what we need anyway.
      auto num_row_entries = this->k_numRowEntries_;
      Kokkos::parallel_for
        ("Fill end row pointers", range_policy(0, N),
         KOKKOS_LAMBDA (const size_t i) {
          row_ptrs_end(i) = row_ptrs_beg(i) + num_row_entries(i);
        });
    }
    else {
      // FIXME (mfh 10 Feb 2020) Fix padCrsArrays so that if packed
      // storage, we don't need row_ptr_end to be separate allocation;
      // could just have it alias row_ptr_beg+1.
      Kokkos::parallel_for
        ("Fill end row pointers", range_policy(0, N),
         KOKKOS_LAMBDA (const size_t i) {
          row_ptrs_end(i) = row_ptrs_beg(i+1);
        });
    }

    if (isGloballyIndexed()) {
      padCrsArrays(row_ptrs_beg, row_ptrs_end, k_gblInds1D_,
                   padding, myRank, verbose);
    }
    else {
      padCrsArrays(row_ptrs_beg, row_ptrs_end, k_lclInds1D_,
                   padding, myRank, verbose);
    }

    if (refill_num_row_entries) {
      auto num_row_entries = this->k_numRowEntries_;
      Kokkos::parallel_for
        ("Fill num entries", range_policy(0, N),
         KOKKOS_LAMBDA (const size_t i) {
          num_row_entries(i) = row_ptrs_end(i) - row_ptrs_beg(i);
        });
    }
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Reassign k_rowPtrs_; old size: "
         << k_rowPtrs_.extent(0) << ", new size: "
         << row_ptrs_beg.extent(0) << endl;
      std::cerr << os.str();
      TEUCHOS_ASSERT( k_rowPtrs_.extent(0) == row_ptrs_beg.extent(0) );
    }
    this->k_rowPtrs_ = row_ptrs_beg;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  std::unique_ptr<
    typename CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::padding_type
    >
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  computeCrsPadding(
    const RowGraph<LocalOrdinal,GlobalOrdinal,Node>& source,
    const size_t numSameIDs,
    const Kokkos::DualView<const local_ordinal_type*,
      buffer_device_type>& permuteToLIDs,
    const Kokkos::DualView<const local_ordinal_type*,
      buffer_device_type>& permuteFromLIDs,
    const bool verbose) const
  {
    using LO = local_ordinal_type;
    using std::endl;

    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsGraph",
        "computeCrsPadding(same & permute)");
      std::ostringstream os;
      os << *prefix << "{numSameIDs: " << numSameIDs
         << ", numPermutes: " << permuteFromLIDs.extent(0) << "}"
         << endl;
      std::cerr << os.str();
    }

    const int myRank = [&] () {
      auto comm = rowMap_.is_null() ? Teuchos::null :
        rowMap_->getComm();
      return comm.is_null() ? -1 : comm->getRank();
    } ();
    std::unique_ptr<padding_type> padding(
      new padding_type(myRank, numSameIDs,
                       permuteFromLIDs.extent(0)));

    // We're accessing data on host, so make sure all device
    // computations on the graphs' data are done.
    //
    // NOTE (mfh 08 Feb 2020) If we ever get rid of this fence, look
    // carefully in computeCrsPaddingFor{Same,Permuted}IDs to see if
    // we need a fence there.
    Kokkos::fence();

    computeCrsPaddingForSameIDs(*padding, source,
                                static_cast<LO>(numSameIDs));
    computeCrsPaddingForPermutedIDs(*padding, source, permuteToLIDs,
                                    permuteFromLIDs);
    return padding;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  computeCrsPaddingForSameIDs(
    padding_type& padding,
    const RowGraph<local_ordinal_type, global_ordinal_type,
      node_type>& source,
    const local_ordinal_type numSameIDs) const
  {
    using LO = local_ordinal_type;
    using GO = global_ordinal_type;
    using Details::Impl::getRowGraphGlobalRow;
    using std::endl;
    const char tfecfFuncName[] = "computeCrsPaddingForSameIds";

    std::unique_ptr<std::string> prefix;
    const bool verbose = verbose_;
    if (verbose) {
      prefix = this->createPrefix("CrsGraph", tfecfFuncName);
      std::ostringstream os;
      os << *prefix << "numSameIDs: " << numSameIDs << endl;
      std::cerr << os.str();
    }

    if (numSameIDs == 0) {
      return;
    }

    const map_type& srcRowMap = *(source.getRowMap());
    const map_type& tgtRowMap = *rowMap_;
    using this_type = CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
    const this_type* srcCrs = dynamic_cast<const this_type*>(&source);
    const bool src_is_unique =
      srcCrs == nullptr ? false : srcCrs->isMerged();
    const bool tgt_is_unique = this->isMerged();

    std::vector<GO> srcGblColIndsScratch;
    std::vector<GO> tgtGblColIndsScratch;
    for (LO lclRowInd = 0; lclRowInd < numSameIDs; ++lclRowInd) {
      const GO srcGblRowInd = srcRowMap.getGlobalElement(lclRowInd);
      const GO tgtGblRowInd = tgtRowMap.getGlobalElement(lclRowInd);
      auto srcGblColInds = getRowGraphGlobalRow(
        srcGblColIndsScratch, source, srcGblRowInd);
      auto tgtGblColInds = getRowGraphGlobalRow(
        tgtGblColIndsScratch, *this, tgtGblRowInd);
      padding.update_same(lclRowInd, tgtGblColInds.getRawPtr(),
                          tgtGblColInds.size(), tgt_is_unique,
                          srcGblColInds.getRawPtr(),
                          srcGblColInds.size(), src_is_unique);
    }
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  computeCrsPaddingForPermutedIDs(
    padding_type& padding,
    const RowGraph<local_ordinal_type, global_ordinal_type,
      node_type>& source,
    const Kokkos::DualView<const local_ordinal_type*,
      buffer_device_type>& permuteToLIDs,
    const Kokkos::DualView<const local_ordinal_type*,
      buffer_device_type>& permuteFromLIDs) const
  {
    using LO = local_ordinal_type;
    using GO = global_ordinal_type;
    using Details::Impl::getRowGraphGlobalRow;
    using std::endl;
    const char tfecfFuncName[] = "computeCrsPaddingForPermutedIds";

    std::unique_ptr<std::string> prefix;
    const bool verbose = verbose_;
    if (verbose) {
      prefix = this->createPrefix("CrsGraph", tfecfFuncName);
      std::ostringstream os;
      os << *prefix << "permuteToLIDs.extent(0): "
         << permuteToLIDs.extent(0)
         << ", permuteFromLIDs.extent(0): "
         << permuteFromLIDs.extent(0) << endl;
      std::cerr << os.str();
    }

    if (permuteToLIDs.extent(0) == 0) {
      return;
    }

    const map_type& srcRowMap = *(source.getRowMap());
    const map_type& tgtRowMap = *rowMap_;
    using this_type = CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
    const this_type* srcCrs = dynamic_cast<const this_type*>(&source);
    const bool src_is_unique =
      srcCrs == nullptr ? false : srcCrs->isMerged();
    const bool tgt_is_unique = this->isMerged();

    TEUCHOS_ASSERT( ! permuteToLIDs.need_sync_host() );
    auto permuteToLIDs_h = permuteToLIDs.view_host();
    TEUCHOS_ASSERT( ! permuteFromLIDs.need_sync_host() );
    auto permuteFromLIDs_h = permuteFromLIDs.view_host();

    std::vector<GO> srcGblColIndsScratch;
    std::vector<GO> tgtGblColIndsScratch;
    const LO numPermutes = static_cast<LO>(permuteToLIDs_h.extent(0));
    for (LO whichPermute = 0; whichPermute < numPermutes; ++whichPermute) {
      const LO srcLclRowInd = permuteFromLIDs_h[whichPermute];
      const GO srcGblRowInd = srcRowMap.getGlobalElement(srcLclRowInd);
      auto srcGblColInds = getRowGraphGlobalRow(
        srcGblColIndsScratch, source, srcGblRowInd);
      const LO tgtLclRowInd = permuteToLIDs_h[whichPermute];
      const GO tgtGblRowInd = tgtRowMap.getGlobalElement(tgtLclRowInd);
      auto tgtGblColInds = getRowGraphGlobalRow(
        tgtGblColIndsScratch, *this, tgtGblRowInd);
      padding.update_permute(whichPermute, tgtLclRowInd,
                             tgtGblColInds.getRawPtr(),
                             tgtGblColInds.size(), tgt_is_unique,
                             srcGblColInds.getRawPtr(),
                             srcGblColInds.size(), src_is_unique);
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  std::unique_ptr<
    typename CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::padding_type
    >
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  computeCrsPaddingForImports(
    const Kokkos::DualView<const local_ordinal_type*,
      buffer_device_type>& importLIDs,
    Kokkos::DualView<packet_type*, buffer_device_type> imports,
    Kokkos::DualView<size_t*, buffer_device_type> numPacketsPerLID,
    const bool verbose) const
  {
    using Details::Impl::getRowGraphGlobalRow;
    using std::endl;
    using LO = local_ordinal_type;
    using GO = global_ordinal_type;
    const char tfecfFuncName[] = "computeCrsPaddingForImports";

    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsGraph", tfecfFuncName);
      std::ostringstream os;
      os << *prefix << "importLIDs.extent(0): "
         << importLIDs.extent(0)
         << ", imports.extent(0): "
         << imports.extent(0)
         << ", numPacketsPerLID.extent(0): "
         << numPacketsPerLID.extent(0) << endl;
      std::cerr << os.str();
    }

    const LO numImports = static_cast<LO>(importLIDs.extent(0));
    const int myRank = [&] () {
      auto comm = rowMap_.is_null() ? Teuchos::null :
        rowMap_->getComm();
      return comm.is_null() ? -1 : comm->getRank();
    } ();
    std::unique_ptr<padding_type> padding(
      new padding_type(myRank, numImports));
    Kokkos::fence(); // Make sure device sees changes made by host
    if (imports.need_sync_host()) {
      imports.sync_host();
    }
    auto imports_h = imports.view_host();
    if (numPacketsPerLID.need_sync_host ()) {
      numPacketsPerLID.sync_host();
    }
    auto numPacketsPerLID_h = numPacketsPerLID.view_host();

    TEUCHOS_ASSERT( ! importLIDs.need_sync_host() );
    auto importLIDs_h = importLIDs.view_host();

    const map_type& tgtRowMap = *rowMap_;
    // Always merge source column indices, since isMerged() is
    // per-process state, and we don't know its value on other
    // processes that sent us data.
    constexpr bool src_is_unique = false;
    const bool tgt_is_unique = isMerged();

    std::vector<GO> tgtGblColIndsScratch;
    size_t offset = 0;
    for (LO whichImport = 0; whichImport < numImports; ++whichImport) {
      // CrsGraph packs just global column indices, while CrsMatrix
      // packs bytes (first the number of entries in the row, then the
      // global column indices, then other stuff like the matrix
      // values in that row).
      const LO origSrcNumEnt =
        static_cast<LO>(numPacketsPerLID_h[whichImport]);
      GO* const srcGblColInds = imports_h.data() + offset;

      const LO tgtLclRowInd = importLIDs_h[whichImport];
      const GO tgtGblRowInd =
        tgtRowMap.getGlobalElement(tgtLclRowInd);
      auto tgtGblColInds = getRowGraphGlobalRow(
        tgtGblColIndsScratch, *this, tgtGblRowInd);
      const size_t origTgtNumEnt(tgtGblColInds.size());

      padding->update_import(whichImport, tgtLclRowInd,
                             tgtGblColInds.getRawPtr(),
                             origTgtNumEnt, tgt_is_unique,
                             srcGblColInds,
                             origSrcNumEnt, src_is_unique);
      offset += origSrcNumEnt;
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str();
    }
    return padding;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  std::unique_ptr<
    typename CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::padding_type
    >
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  computePaddingForCrsMatrixUnpack(
    const Kokkos::DualView<const local_ordinal_type*,
      buffer_device_type>& importLIDs,
    Kokkos::DualView<char*, buffer_device_type> imports,
    Kokkos::DualView<size_t*, buffer_device_type> numPacketsPerLID,
    const bool verbose) const
  {
    using Details::Impl::getRowGraphGlobalRow;
    using Details::PackTraits;
    using std::endl;
    using LO = local_ordinal_type;
    using GO = global_ordinal_type;
    const char tfecfFuncName[] = "computePaddingForCrsMatrixUnpack";

    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsGraph", tfecfFuncName);
      std::ostringstream os;
      os << *prefix << "importLIDs.extent(0): "
         << importLIDs.extent(0)
         << ", imports.extent(0): "
         << imports.extent(0)
         << ", numPacketsPerLID.extent(0): "
         << numPacketsPerLID.extent(0) << endl;
      std::cerr << os.str();
    }
    const bool extraVerbose =
      verbose && Details::Behavior::verbose("CrsPadding");

    const LO numImports = static_cast<LO>(importLIDs.extent(0));
    TEUCHOS_ASSERT( LO(numPacketsPerLID.extent(0)) >= numImports );
    const int myRank = [&] () {
      auto comm = rowMap_.is_null() ? Teuchos::null :
        rowMap_->getComm();
      return comm.is_null() ? -1 : comm->getRank();
    } ();
    std::unique_ptr<padding_type> padding(
      new padding_type(myRank, numImports));
    Kokkos::fence(); // Make sure host sees changes made by device
    if (imports.need_sync_host()) {
      imports.sync_host();
    }
    auto imports_h = imports.view_host();
    if (numPacketsPerLID.need_sync_host ()) {
      numPacketsPerLID.sync_host();
    }
    auto numPacketsPerLID_h = numPacketsPerLID.view_host();

    TEUCHOS_ASSERT( ! importLIDs.need_sync_host() );
    auto importLIDs_h = importLIDs.view_host();

    const map_type& tgtRowMap = *rowMap_;
    // Always merge source column indices, since isMerged() is
    // per-process state, and we don't know its value on other
    // processes that sent us data.
    constexpr bool src_is_unique = false;
    const bool tgt_is_unique = isMerged();

    std::vector<GO> srcGblColIndsScratch;
    std::vector<GO> tgtGblColIndsScratch;
    size_t offset = 0;
    for (LO whichImport = 0; whichImport < numImports; ++whichImport) {
      // CrsGraph packs just global column indices, while CrsMatrix
      // packs bytes (first the number of entries in the row, then the
      // global column indices, then other stuff like the matrix
      // values in that row).
      const size_t numBytes = numPacketsPerLID_h[whichImport];
      if (extraVerbose) {
        std::ostringstream os;
        os << *prefix << "whichImport=" << whichImport
           << ", numImports=" << numImports
           << ", numBytes=" << numBytes << endl;
        std::cerr << os.str();
      }
      if (numBytes == 0) {
        continue; // special case: no entries to unpack for this row
      }
      LO origSrcNumEnt = 0;
      const size_t numEntBeg = offset;
      const size_t numEntLen =
        PackTraits<LO>::packValueCount(origSrcNumEnt);
      TEUCHOS_ASSERT( numBytes >= numEntLen );
      TEUCHOS_ASSERT( imports_h.extent(0) >= numEntBeg + numEntLen );
      PackTraits<LO>::unpackValue(origSrcNumEnt,
                                  imports_h.data() + numEntBeg);
      if (extraVerbose) {
        std::ostringstream os;
        os << *prefix << "whichImport=" << whichImport
           << ", numImports=" << numImports
           << ", origSrcNumEnt=" << origSrcNumEnt << endl;
        std::cerr << os.str();
      }
      TEUCHOS_ASSERT( origSrcNumEnt >= LO(0) );
      TEUCHOS_ASSERT( numBytes >= size_t(numEntLen + origSrcNumEnt * sizeof(GO)) );
      const size_t gidsBeg = numEntBeg + numEntLen;
      if (srcGblColIndsScratch.size() < size_t(origSrcNumEnt)) {
        srcGblColIndsScratch.resize(origSrcNumEnt);
      }
      GO* const srcGblColInds = srcGblColIndsScratch.data();
      PackTraits<GO>::unpackArray(srcGblColInds,
                                  imports_h.data() + gidsBeg,
                                  origSrcNumEnt);
      const LO tgtLclRowInd = importLIDs_h[whichImport];
      const GO tgtGblRowInd =
        tgtRowMap.getGlobalElement(tgtLclRowInd);
      auto tgtGblColInds = getRowGraphGlobalRow(
        tgtGblColIndsScratch, *this, tgtGblRowInd);
      const size_t origNumTgtEnt(tgtGblColInds.size());

      if (extraVerbose) {
        std::ostringstream os;
        os << *prefix << "whichImport=" << whichImport
           << ", numImports=" << numImports
           << ": Call padding->update_import" << endl;
        std::cerr << os.str();
      }
      padding->update_import(whichImport, tgtLclRowInd,
                             tgtGblColInds.getRawPtr(),
                             origNumTgtEnt, tgt_is_unique,
                             srcGblColInds,
                             origSrcNumEnt, src_is_unique);
      offset += numBytes;
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str();
    }
    return padding;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  packAndPrepare
  (const SrcDistObject& source,
   const Kokkos::DualView<const local_ordinal_type*,
     buffer_device_type>& exportLIDs,
   Kokkos::DualView<packet_type*,
     buffer_device_type>& exports,
   Kokkos::DualView<size_t*,
     buffer_device_type> numPacketsPerLID,
   size_t& constantNumPackets,
   Distributor& distor)
  {
    using Tpetra::Details::ProfilingRegion;
    using GO = global_ordinal_type;
    using std::endl;
    using crs_graph_type =
      CrsGraph<local_ordinal_type, global_ordinal_type, node_type>;
    using row_graph_type =
      RowGraph<local_ordinal_type, global_ordinal_type, node_type>;
    const char tfecfFuncName[] = "packAndPrepare: ";
    ProfilingRegion region_papn ("Tpetra::CrsGraph::packAndPrepare");

    const bool verbose = verbose_;
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsGraph", "packAndPrepare");
      std::ostringstream os;
      os << *prefix << "Start" << endl;
      std::cerr << os.str();
    }

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (exportLIDs.extent (0) != numPacketsPerLID.extent (0),
       std::runtime_error,
       "exportLIDs.extent(0) = " << exportLIDs.extent (0)
       << " != numPacketsPerLID.extent(0) = " << numPacketsPerLID.extent (0)
       << ".");
    const row_graph_type* srcRowGraphPtr =
      dynamic_cast<const row_graph_type*> (&source);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (srcRowGraphPtr == nullptr, std::invalid_argument, "Source of an Export "
       "or Import operation to a CrsGraph must be a RowGraph with the same "
       "template parameters.");
    // We don't check whether src_graph has had fillComplete called,
    // because it doesn't matter whether the *source* graph has been
    // fillComplete'd. The target graph can not be fillComplete'd yet.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (this->isFillComplete (), std::runtime_error,
       "The target graph of an Import or Export must not be fill complete.");

    const crs_graph_type* srcCrsGraphPtr =
      dynamic_cast<const crs_graph_type*> (&source);

    if (srcCrsGraphPtr == nullptr) {
      using Teuchos::ArrayView;
      using LO = local_ordinal_type;

      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Source is a RowGraph but not a CrsGraph"
           << endl;
        std::cerr << os.str();
      }
      // RowGraph::pack serves the "old" DistObject interface.  It
      // takes Teuchos::ArrayView and Teuchos::Array&.  The latter
      // entails deep-copying the exports buffer on output.  RowGraph
      // is a convenience interface when not a CrsGraph, so we accept
      // the performance hit.
      TEUCHOS_ASSERT( ! exportLIDs.need_sync_host () );
      auto exportLIDs_h = exportLIDs.view_host ();
      ArrayView<const LO> exportLIDs_av (exportLIDs_h.data (),
                                         exportLIDs_h.extent (0));
      Teuchos::Array<GO> exports_a;

      numPacketsPerLID.clear_sync_state ();
      numPacketsPerLID.modify_host ();
      auto numPacketsPerLID_h = numPacketsPerLID.view_host ();
      ArrayView<size_t> numPacketsPerLID_av (numPacketsPerLID_h.data (),
                                             numPacketsPerLID_h.extent (0));
      srcRowGraphPtr->pack (exportLIDs_av, exports_a, numPacketsPerLID_av,
                            constantNumPackets, distor);
      const size_t newSize = static_cast<size_t> (exports_a.size ());
      if (static_cast<size_t> (exports.extent (0)) != newSize) {
        using exports_dv_type = Kokkos::DualView<packet_type*, buffer_device_type>;
        exports = exports_dv_type ("exports", newSize);
      }
      Kokkos::View<const packet_type*, Kokkos::HostSpace,
        Kokkos::MemoryUnmanaged> exports_a_h (exports_a.getRawPtr (), newSize);
      exports.clear_sync_state ();
      exports.modify_host ();
      Kokkos::deep_copy (exports.view_host (), exports_a_h);
    }
    // packCrsGraphNew requires a valid localGraph.
    else if (! getColMap ().is_null () &&
        (lclGraph_.row_map.extent (0) != 0 ||
         getRowMap ()->getNodeNumElements () == 0)) {
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "packCrsGraphNew path" << endl;
        std::cerr << os.str();
      }
      using export_pids_type =
        Kokkos::DualView<const int*, buffer_device_type>;
      export_pids_type exportPIDs; // not filling it; needed for syntax
      using LO = local_ordinal_type;
      using NT = node_type;
      using Tpetra::Details::packCrsGraphNew;
      packCrsGraphNew<LO,GO,NT> (*srcCrsGraphPtr, exportLIDs, exportPIDs,
                                 exports, numPacketsPerLID,
                                 constantNumPackets, false, distor);
    }
    else {
      srcCrsGraphPtr->packFillActiveNew (exportLIDs, exports, numPacketsPerLID,
                                         constantNumPackets, distor);
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  pack (const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
        Teuchos::Array<GlobalOrdinal>& exports,
        const Teuchos::ArrayView<size_t>& numPacketsPerLID,
        size_t& constantNumPackets,
        Distributor& distor) const
  {
    auto col_map = this->getColMap();
    // packCrsGraph requires a valid localGraph.
    if( !col_map.is_null() && (lclGraph_.row_map.extent(0) != 0  ||  getRowMap()->getNodeNumElements() ==0)) {
      using Tpetra::Details::packCrsGraph;
      packCrsGraph<LocalOrdinal,GlobalOrdinal,Node>(*this, exports, numPacketsPerLID,
                                                    exportLIDs, constantNumPackets, distor);
    }
    else {
      this->packFillActive(exportLIDs, exports, numPacketsPerLID,
                           constantNumPackets, distor);
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  packFillActive (const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
                  Teuchos::Array<GlobalOrdinal>& exports,
                  const Teuchos::ArrayView<size_t>& numPacketsPerLID,
                  size_t& constantNumPackets,
                  Distributor& /* distor */) const
  {
    using std::endl;
    using LO = LocalOrdinal;
    using GO = GlobalOrdinal;
    using host_execution_space =
      typename Kokkos::View<size_t*, device_type>::
        HostMirror::execution_space;
    using device_execution_space =
      typename device_type::execution_space;
    const char tfecfFuncName[] = "packFillActive: ";
    const bool verbose = verbose_;

    const auto numExportLIDs = exportLIDs.size ();
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsGraph", "allocateIndices");
      std::ostringstream os;
      os << *prefix << "numExportLIDs=" << numExportLIDs << endl;
      std::cerr << os.str();
    }
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (numExportLIDs != numPacketsPerLID.size (), std::runtime_error,
       "exportLIDs.size() = " << numExportLIDs << " != numPacketsPerLID.size()"
       " = " << numPacketsPerLID.size () << ".");

    // We may be accessing UVM data on host below, so ensure that the
    // device is done accessing it.
    device_execution_space().fence ();

    const map_type& rowMap = * (this->getRowMap ());
    const map_type* const colMapPtr = this->colMap_.getRawPtr ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (this->isLocallyIndexed () && colMapPtr == nullptr, std::logic_error,
       "This graph claims to be locally indexed, but its column Map is nullptr.  "
       "This should never happen.  Please report this bug to the Tpetra "
       "developers.");

    // We may pack different amounts of data for different rows.
    constantNumPackets = 0;

    // mfh 20 Sep 2017: Teuchos::ArrayView isn't thread safe (well,
    // it might be now, but we might as well be safe).
    size_t* const numPacketsPerLID_raw = numPacketsPerLID.getRawPtr ();
    const LO* const exportLIDs_raw = exportLIDs.getRawPtr ();

    // Count the total number of packets (column indices, in the case
    // of a CrsGraph) to pack.  While doing so, set
    // numPacketsPerLID[i] to the number of entries owned by the
    // calling process in (local) row exportLIDs[i] of the graph, that
    // the caller wants us to send out.
    Kokkos::RangePolicy<host_execution_space, LO> inputRange (0, numExportLIDs);
    size_t totalNumPackets = 0;
    size_t errCount = 0;
    // lambdas turn what they capture const, so we can't
    // atomic_add(&errCount,1).  Instead, we need a View to modify.
    typedef Kokkos::Device<host_execution_space, Kokkos::HostSpace>
      host_device_type;
    Kokkos::View<size_t, host_device_type> errCountView (&errCount);
    constexpr size_t ONE = 1;

    Kokkos::parallel_reduce ("Tpetra::CrsGraph::pack: totalNumPackets",
      inputRange,
      [=] (const LO& i, size_t& curTotalNumPackets) {
        const GO gblRow = rowMap.getGlobalElement (exportLIDs_raw[i]);
        if (gblRow == Tpetra::Details::OrdinalTraits<GO>::invalid ()) {
          Kokkos::atomic_add (&errCountView(), ONE);
          numPacketsPerLID_raw[i] = 0;
        }
        else {
          const size_t numEnt = this->getNumEntriesInGlobalRow (gblRow);
          numPacketsPerLID_raw[i] = numEnt;
          curTotalNumPackets += numEnt;
        }
      },
      totalNumPackets);

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "totalNumPackets=" << totalNumPackets << endl;
      std::cerr << os.str();
    }
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (errCount != 0, std::logic_error, "totalNumPackets count encountered "
       "one or more errors!  errCount = " << errCount
       << ", totalNumPackets = " << totalNumPackets << ".");
    errCount = 0;

    // Allocate space for all the column indices to pack.
    exports.resize (totalNumPackets);

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! this->supportsRowViews (), std::logic_error,
       "this->supportsRowViews() returns false; this should never happen.  "
       "Please report this bug to the Tpetra developers.");

    // Loop again over the rows to export, and pack rows of indices
    // into the output buffer.

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Pack into exports" << endl;
      std::cerr << os.str();
    }

    // Teuchos::ArrayView may not be thread safe, or may not be
    // efficiently thread safe.  Better to use the raw pointer.
    GO* const exports_raw = exports.getRawPtr ();
    errCount = 0;
    Kokkos::parallel_scan ("Tpetra::CrsGraph::pack: pack from views",
      inputRange, [=, &prefix]
      (const LO i, size_t& exportsOffset, const bool final) {
        const size_t curOffset = exportsOffset;
        const GO gblRow = rowMap.getGlobalElement (exportLIDs_raw[i]);
        const RowInfo rowInfo =
          this->getRowInfoFromGlobalRowIndex (gblRow);

        using TDO = Tpetra::Details::OrdinalTraits<size_t>;
        if (rowInfo.localRow == TDO::invalid ()) {
          if (verbose) {
            std::ostringstream os;
            os << *prefix << ": INVALID rowInfo: i=" << i
               << ", lclRow=" << exportLIDs_raw[i] << endl;
            std::cerr << os.str();
          }
          Kokkos::atomic_add (&errCountView(), ONE);
        }
        else if (curOffset + rowInfo.numEntries > totalNumPackets) {
          if (verbose) {
            std::ostringstream os;
            os << *prefix << ": UH OH!  For i=" << i << ", lclRow="
               << exportLIDs_raw[i] << ", gblRow=" << gblRow << ", curOffset "
              "(= " << curOffset << ") + numEnt (= " << rowInfo.numEntries
               << ") > totalNumPackets (= " << totalNumPackets << ")."
               << endl;
            std::cerr << os.str();
          }
          Kokkos::atomic_add (&errCountView(), ONE);
        }
        else {
          const LO numEnt = static_cast<LO> (rowInfo.numEntries);
          if (this->isLocallyIndexed ()) {
            const LO* lclColInds = nullptr;
            LO capacity = 0;
            const LO errCode =
              this->getLocalViewRawConst (lclColInds, capacity, rowInfo);
            if (errCode == 0) {
              if (final) {
                for (LO k = 0; k < numEnt; ++k) {
                  const LO lclColInd = lclColInds[k];
                  const GO gblColInd = colMapPtr->getGlobalElement (lclColInd);
                  // Pack it, even if it's wrong.  Let the receiving
                  // process deal with it.  Otherwise, we'll miss out
                  // on any correct data.
                  exports_raw[curOffset + k] = gblColInd;
                } // for each entry in the row
              } // final pass?
              exportsOffset = curOffset + numEnt;
            }
            else { // error in getting local row view
              Kokkos::atomic_add (&errCountView(), ONE);
            }
          }
          else if (this->isGloballyIndexed ()) {
            const GO* gblColInds = nullptr;
            LO capacity = 0;
            const LO errCode =
              this->getGlobalViewRawConst (gblColInds, capacity, rowInfo);
            if (errCode == 0) {
              if (final) {
                for (LO k = 0; k < numEnt; ++k) {
                  const GO gblColInd = gblColInds[k];
                  // Pack it, even if it's wrong.  Let the receiving
                  // process deal with it.  Otherwise, we'll miss out
                  // on any correct data.
                  exports_raw[curOffset + k] = gblColInd;
                } // for each entry in the row
              } // final pass?
              exportsOffset = curOffset + numEnt;
            }
            else { // error in getting global row view
              Kokkos::atomic_add (&errCountView(), ONE);
            }
          }
          // If neither globally nor locally indexed, then the graph
          // has no entries in this row (or indeed, in any row on this
          // process) to pack.
        }
      });

    // We may have accessed UVM data on host above, so ensure that the
    // device sees these changes.
    device_execution_space().fence ();

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (errCount != 0, std::logic_error, "Packing encountered "
       "one or more errors!  errCount = " << errCount
       << ", totalNumPackets = " << totalNumPackets << ".");

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  packFillActiveNew (const Kokkos::DualView<const local_ordinal_type*,
                       buffer_device_type>& exportLIDs,
                     Kokkos::DualView<packet_type*,
                       buffer_device_type>& exports,
                     Kokkos::DualView<size_t*,
                       buffer_device_type> numPacketsPerLID,
                     size_t& constantNumPackets,
                     Distributor& distor) const
  {
    using std::endl;
    using LO = local_ordinal_type;
    using GO = global_ordinal_type;
    using host_execution_space = typename Kokkos::View<size_t*,
      device_type>::HostMirror::execution_space;
    using host_device_type =
      Kokkos::Device<host_execution_space, Kokkos::HostSpace>;
    using device_execution_space = typename device_type::execution_space;
    using exports_dv_type =
      Kokkos::DualView<packet_type*, buffer_device_type>;
    const char tfecfFuncName[] = "packFillActiveNew: ";
    const bool verbose = verbose_;

    const auto numExportLIDs = exportLIDs.extent (0);
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsGraph", "packFillActiveNew");
      std::ostringstream os;
      os << *prefix << "numExportLIDs: " << numExportLIDs
         << ", numPacketsPerLID.extent(0): "
         << numPacketsPerLID.extent(0) << endl;
      std::cerr << os.str();
    }
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (numExportLIDs != numPacketsPerLID.extent (0), std::runtime_error,
       "exportLIDs.extent(0) = " << numExportLIDs
       << " != numPacketsPerLID.extent(0) = "
       << numPacketsPerLID.extent (0) << ".");
    TEUCHOS_ASSERT( ! exportLIDs.need_sync_host () );
    auto exportLIDs_h = exportLIDs.view_host ();

    // We may be accessing UVM data on host below, so ensure that the
    // device is done accessing it.
    device_execution_space().fence ();

    const map_type& rowMap = * (this->getRowMap ());
    const map_type* const colMapPtr = this->colMap_.getRawPtr ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (this->isLocallyIndexed () && colMapPtr == nullptr, std::logic_error,
       "This graph claims to be locally indexed, but its column Map is nullptr.  "
       "This should never happen.  Please report this bug to the Tpetra "
       "developers.");

    // We may pack different amounts of data for different rows.
    constantNumPackets = 0;

    numPacketsPerLID.clear_sync_state ();
    numPacketsPerLID.modify_host ();
    auto numPacketsPerLID_h = numPacketsPerLID.view_host ();

    // Count the total number of packets (column indices, in the case
    // of a CrsGraph) to pack.  While doing so, set
    // numPacketsPerLID[i] to the number of entries owned by the
    // calling process in (local) row exportLIDs[i] of the graph, that
    // the caller wants us to send out.
    using range_type = Kokkos::RangePolicy<host_execution_space, LO>;
    range_type inputRange (0, numExportLIDs);
    size_t totalNumPackets = 0;
    size_t errCount = 0;
    // lambdas turn what they capture const, so we can't
    // atomic_add(&errCount,1).  Instead, we need a View to modify.
    Kokkos::View<size_t, host_device_type> errCountView (&errCount);
    constexpr size_t ONE = 1;

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Compute totalNumPackets" << endl;
      std::cerr << os.str ();
    }

    Kokkos::parallel_reduce
      ("Tpetra::CrsGraph::pack: totalNumPackets",
       inputRange,
       [=, &prefix] (const LO i, size_t& curTotalNumPackets) {
         const LO lclRow = exportLIDs_h[i];
         const GO gblRow = rowMap.getGlobalElement (lclRow);
         if (gblRow == Tpetra::Details::OrdinalTraits<GO>::invalid ()) {
           if (verbose) {
             std::ostringstream os;
             os << *prefix << "For i=" << i << ", lclRow=" << lclRow
                << " not in row Map on this process" << endl;
             std::cerr << os.str();
           }
           Kokkos::atomic_add (&errCountView(), ONE);
           numPacketsPerLID_h(i) = 0;
         }
         else {
           const size_t numEnt = this->getNumEntriesInGlobalRow (gblRow);
           numPacketsPerLID_h(i) = numEnt;
           curTotalNumPackets += numEnt;
         }
      },
      totalNumPackets);

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "totalNumPackets: " << totalNumPackets
         << ", errCount: " << errCount << endl;
      std::cerr << os.str ();
    }
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (errCount != 0, std::logic_error, "totalNumPackets count encountered "
       "one or more errors! totalNumPackets: " << totalNumPackets
       << ", errCount: " << errCount << ".");

    // Allocate space for all the column indices to pack.
    if (size_t(exports.extent (0)) < totalNumPackets) {
      // FIXME (mfh 09 Apr 2019) Create without initializing.
      exports = exports_dv_type ("exports", totalNumPackets);
    }

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! this->supportsRowViews (), std::logic_error,
       "this->supportsRowViews() returns false; this should never happen.  "
       "Please report this bug to the Tpetra developers.");

    // Loop again over the rows to export, and pack rows of indices
    // into the output buffer.

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Pack into exports buffer" << endl;
      std::cerr << os.str();
    }

    exports.clear_sync_state ();
    exports.modify_host ();
    auto exports_h = exports.view_host ();

    // The graph may store its data in UVM memory, so make sure that
    // any device kernels are done modifying the graph's data before
    // reading the data.
    device_execution_space().fence ();

    errCount = 0;
    Kokkos::parallel_scan
      ("Tpetra::CrsGraph::packFillActiveNew: Pack exports",
       inputRange, [=, &prefix]
       (const LO i, size_t& exportsOffset, const bool final) {
         const size_t curOffset = exportsOffset;
         const LO lclRow = exportLIDs_h(i);
         const GO gblRow = rowMap.getGlobalElement (lclRow);
         if (gblRow == Details::OrdinalTraits<GO>::invalid ()) {
           if (verbose) {
             std::ostringstream os;
             os << *prefix << "For i=" << i << ", lclRow=" << lclRow
                << " not in row Map on this process" << endl;
             std::cerr << os.str();
           }
           Kokkos::atomic_add (&errCountView(), ONE);
           return;
         }

         const RowInfo rowInfo = this->getRowInfoFromGlobalRowIndex (gblRow);
         if (rowInfo.localRow == Details::OrdinalTraits<size_t>::invalid ()) {
           if (verbose) {
             std::ostringstream os;
             os << *prefix << "For i=" << i << ", lclRow=" << lclRow
                << ", gblRow=" << gblRow << ": invalid rowInfo"
                << endl;
             std::cerr << os.str();
           }
           Kokkos::atomic_add (&errCountView(), ONE);
           return;
         }

         if (curOffset + rowInfo.numEntries > totalNumPackets) {
           if (verbose) {
             std::ostringstream os;
             os << *prefix << "For i=" << i << ", lclRow=" << lclRow
                << ", gblRow=" << gblRow << ", curOffset (= "
                << curOffset << ") + numEnt (= " << rowInfo.numEntries
                << ") > totalNumPackets (= " << totalNumPackets
                << ")." << endl;
             std::cerr << os.str();
           }
           Kokkos::atomic_add (&errCountView(), ONE);
           return;
         }

         const LO numEnt = static_cast<LO> (rowInfo.numEntries);
         if (this->isLocallyIndexed ()) {
           const LO* lclColInds = nullptr;
           LO capacity = 0;
           const LO errCode =
             this->getLocalViewRawConst (lclColInds, capacity, rowInfo);
           if (errCode == 0) {
             if (final) {
               for (LO k = 0; k < numEnt; ++k) {
                 const LO lclColInd = lclColInds[k];
                 const GO gblColInd = colMapPtr->getGlobalElement (lclColInd);
                 // Pack it, even if it's wrong.  Let the receiving
                 // process deal with it.  Otherwise, we'll miss out
                 // on any correct data.
                 exports_h(curOffset + k) = gblColInd;
               } // for each entry in the row
             } // final pass?
             exportsOffset = curOffset + numEnt;
           }
           else { // error in getting local row view
             if (verbose) {
               std::ostringstream os;
               os << *prefix << "For i=" << i << ", lclRow=" << lclRow
                  << ", gblRow=" << gblRow << ": getLocalViewRawConst"
                 "returned nonzero error code " << errCode << endl;
               std::cerr << os.str();
             }
             Kokkos::atomic_add (&errCountView(), ONE);
           }
         }
         else if (this->isGloballyIndexed ()) {
           const GO* gblColInds = nullptr;
           LO capacity = 0;
           const LO errCode =
             this->getGlobalViewRawConst (gblColInds, capacity, rowInfo);
           if (errCode == 0) {
             if (final) {
               for (LO k = 0; k < numEnt; ++k) {
                 const GO gblColInd = gblColInds[k];
                 // Pack it, even if it's wrong.  Let the receiving
                 // process deal with it.  Otherwise, we'll miss out
                 // on any correct data.
                 exports_h(curOffset + k) = gblColInd;
               } // for each entry in the row
             } // final pass?
             exportsOffset = curOffset + numEnt;
           }
           else { // error in getting global row view
             if (verbose) {
               std::ostringstream os;
               os << *prefix << "For i=" << i << ", lclRow=" << lclRow
                  << ", gblRow=" << gblRow << ": "
                 "getGlobalViewRawConst returned nonzero error code "
                  << errCode << endl;
               std::cerr << os.str();
             }
             Kokkos::atomic_add (&errCountView(), ONE);
           }
         }
         // If neither globally nor locally indexed, then the graph
         // has no entries in this row (or indeed, in any row on this
         // process) to pack.
      });

    // TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
    //   (errCount != 0, std::logic_error, "Packing encountered "
    //    "one or more errors!  errCount = " << errCount
    //    << ", totalNumPackets = " << totalNumPackets << ".");

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "errCount=" << errCount << "; Done" << endl;
      std::cerr << os.str();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  unpackAndCombine
  (const Kokkos::DualView<const local_ordinal_type*,
     buffer_device_type>& importLIDs,
   Kokkos::DualView<packet_type*,
     buffer_device_type> imports,
   Kokkos::DualView<size_t*,
     buffer_device_type> numPacketsPerLID,
   const size_t /* constantNumPackets */,
   Distributor& /* distor */,
   const CombineMode /* combineMode */ )
  {
    using Details::ProfilingRegion;
    using std::endl;
    using LO = local_ordinal_type;
    using GO = global_ordinal_type;
    const char tfecfFuncName[] = "unpackAndCombine";

    ProfilingRegion regionCGC("Tpetra::CrsGraph::unpackAndCombine");
    const bool verbose = verbose_;

    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsGraph", tfecfFuncName);
      std::ostringstream os;
      os << *prefix << "Start" << endl;
      std::cerr << os.str ();
    }
    {
      auto padding = computeCrsPaddingForImports(
        importLIDs, imports, numPacketsPerLID, verbose);
      applyCrsPadding(*padding, verbose);
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Done computing & applying padding" << endl;
        std::cerr << os.str ();
      }
    }

    // FIXME (mfh 02 Apr 2012) REPLACE combine mode has a perfectly
    // reasonable meaning, whether or not the matrix is fill complete.
    // It's just more work to implement.

    // We are not checking the value of the CombineMode input
    // argument.  For CrsGraph, we only support import/export
    // operations if fillComplete has not yet been called.  Any
    // incoming column-indices are inserted into the target graph. In
    // this context, CombineMode values of ADD vs INSERT are
    // equivalent. What is the meaning of REPLACE for CrsGraph? If a
    // duplicate column-index is inserted, it will be compressed out
    // when fillComplete is called.
    //
    // Note: I think REPLACE means that an existing row is replaced by
    // the imported row, i.e., the existing indices are cleared. CGB,
    // 6/17/2010

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (importLIDs.extent (0) != numPacketsPerLID.extent (0),
       std::runtime_error, ": importLIDs.extent(0) = "
       << importLIDs.extent (0) << " != numPacketsPerLID.extent(0) = "
       << numPacketsPerLID.extent (0) << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (isFillComplete (), std::runtime_error,
       ": Import or Export operations are not allowed on a target "
       "CrsGraph that is fillComplete.");

    const size_t numImportLIDs(importLIDs.extent(0));
    if (numPacketsPerLID.need_sync_host()) {
      numPacketsPerLID.sync_host();
    }
    auto numPacketsPerLID_h = numPacketsPerLID.view_host();
    if (imports.need_sync_host()) {
      imports.sync_host();
    }
    auto imports_h = imports.view_host();
    TEUCHOS_ASSERT( ! importLIDs.need_sync_host() );
    auto importLIDs_h = importLIDs.view_host();

    // If we're inserting in local indices, let's pre-allocate
    Teuchos::Array<LO> lclColInds;
    if (isLocallyIndexed()) {
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Preallocate local indices scratch" << endl;
        std::cerr << os.str();
      }
      size_t maxNumInserts = 0;
      for (size_t i = 0; i < numImportLIDs; ++i) {
        maxNumInserts = std::max (maxNumInserts, numPacketsPerLID_h[i]);
      }
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Local indices scratch size: "
           << maxNumInserts << endl;
        std::cerr << os.str();
      }
      lclColInds.resize (maxNumInserts);
    }
    else {
      if (verbose) {
        std::ostringstream os;
        os << *prefix;
        if (isGloballyIndexed()) {
          os << "Graph is globally indexed";
        }
        else {
          os << "Graph is neither locally nor globally indexed";
        }
        os << endl;
        std::cerr << os.str();
      }
    }

    TEUCHOS_ASSERT( ! rowMap_.is_null() );
    const map_type& rowMap = *rowMap_;

    try {
      size_t importsOffset = 0;
      for (size_t i = 0; i < numImportLIDs; ++i) {
        if (verbose) {
          std::ostringstream os;
          os << *prefix << "i=" << i << ", numImportLIDs="
             << numImportLIDs << endl;
          std::cerr << os.str();
        }
        // We can only unpack into owned rows, since we only have
        // local row indices.
        const LO lclRow = importLIDs_h[i];
        const GO gblRow = rowMap.getGlobalElement(lclRow);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
          (gblRow == Teuchos::OrdinalTraits<GO>::invalid(),
           std::logic_error, "importLIDs[i=" << i << "]="
           << lclRow << " is not in the row Map on the calling "
           "process.");
        const LO numEnt = numPacketsPerLID_h[i];
        const GO* const gblColInds = (numEnt == 0) ? nullptr :
          imports_h.data() + importsOffset;
        if (! isLocallyIndexed()) {
          insertGlobalIndicesFiltered(lclRow, gblColInds, numEnt);
        }
        else {
          // FIXME (mfh 09 Feb 2020) Now would be a good time to do
          // column Map filtering.
          for (LO j = 0; j < numEnt; j++)  {
            lclColInds[j] = colMap_->getLocalElement(gblColInds[j]);
          }
          insertLocalIndices(lclRow, numEnt, lclColInds.data());
        }
        importsOffset += numEnt;
      }
    }
    catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::runtime_error,
         "Tpetra::CrsGraph::unpackAndCombine: Insert loop threw an "
         "exception: " << endl << e.what());
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str();
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  removeEmptyProcessesInPlace (const Teuchos::RCP<const map_type>& newMap)
  {
    using Teuchos::Comm;
    using Teuchos::null;
    using Teuchos::ParameterList;
    using Teuchos::RCP;

    // We'll set all the state "transactionally," so that this method
    // satisfies the strong exception guarantee.  This object's state
    // won't be modified until the end of this method.
    RCP<const map_type> rowMap, domainMap, rangeMap, colMap;
    RCP<import_type> importer;
    RCP<export_type> exporter;

    rowMap = newMap;
    RCP<const Comm<int> > newComm =
      (newMap.is_null ()) ? null : newMap->getComm ();

    if (! domainMap_.is_null ()) {
      if (domainMap_.getRawPtr () == rowMap_.getRawPtr ()) {
        // Common case: original domain and row Maps are identical.
        // In that case, we need only replace the original domain Map
        // with the new Map.  This ensures that the new domain and row
        // Maps _stay_ identical.
        domainMap = newMap;
      } else {
        domainMap = domainMap_->replaceCommWithSubset (newComm);
      }
    }
    if (! rangeMap_.is_null ()) {
      if (rangeMap_.getRawPtr () == rowMap_.getRawPtr ()) {
        // Common case: original range and row Maps are identical.  In
        // that case, we need only replace the original range Map with
        // the new Map.  This ensures that the new range and row Maps
        // _stay_ identical.
        rangeMap = newMap;
      } else {
        rangeMap = rangeMap_->replaceCommWithSubset (newComm);
      }
    }
    if (! colMap.is_null ()) {
      colMap = colMap_->replaceCommWithSubset (newComm);
    }

    // (Re)create the Export and / or Import if necessary.
    if (! newComm.is_null ()) {
      RCP<ParameterList> params = this->getNonconstParameterList (); // could be null
      //
      // The operations below are collective on the new communicator.
      //
      // (Re)create the Export object if necessary.  If I haven't
      // called fillComplete yet, I don't have a rangeMap, so I must
      // first check if the _original_ rangeMap is not null.  Ditto
      // for the Import object and the domain Map.
      if (! rangeMap_.is_null () &&
          rangeMap != rowMap &&
          ! rangeMap->isSameAs (*rowMap)) {
        if (params.is_null () || ! params->isSublist ("Export")) {
          exporter = rcp (new export_type (rowMap, rangeMap));
        }
        else {
          RCP<ParameterList> exportSublist = sublist (params, "Export", true);
          exporter = rcp (new export_type (rowMap, rangeMap, exportSublist));
        }
      }
      // (Re)create the Import object if necessary.
      if (! domainMap_.is_null () &&
          domainMap != colMap &&
          ! domainMap->isSameAs (*colMap)) {
        if (params.is_null () || ! params->isSublist ("Import")) {
          importer = rcp (new import_type (domainMap, colMap));
        } else {
          RCP<ParameterList> importSublist = sublist (params, "Import", true);
          importer = rcp (new import_type (domainMap, colMap, importSublist));
        }
      }
    } // if newComm is not null

    // Defer side effects until the end.  If no destructors throw
    // exceptions (they shouldn't anyway), then this method satisfies
    // the strong exception guarantee.
    exporter_ = exporter;
    importer_ = importer;
    rowMap_ = rowMap;
    // mfh 31 Mar 2013: DistObject's map_ is the row Map of a CrsGraph
    // or CrsMatrix.  CrsGraph keeps a redundant pointer (rowMap_) to
    // the same object.  We might want to get rid of this redundant
    // pointer sometime, but for now, we'll leave it alone and just
    // set map_ to the same object.
    this->map_ = rowMap;
    domainMap_ = domainMap;
    rangeMap_ = rangeMap;
    colMap_ = colMap;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getLocalDiagOffsets (const Kokkos::View<size_t*, device_type, Kokkos::MemoryUnmanaged>& offsets) const
  {
    using std::endl;
    using LO = LocalOrdinal;
    using GO = GlobalOrdinal;
    const char tfecfFuncName[] = "getLocalDiagOffsets: ";
    const bool verbose = verbose_;

    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("CrsGraph", "getLocalDiagOffsets");
      std::ostringstream os;
      os << *prefix << "offsets.extent(0)=" << offsets.extent(0)
         << endl;
      std::cerr << os.str();
    }

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! hasColMap (), std::runtime_error, "The graph must have a column Map.");
    const LO lclNumRows = static_cast<LO> (this->getNodeNumRows ());
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (static_cast<LO> (offsets.extent (0)) < lclNumRows,
       std::invalid_argument, "offsets.extent(0) = " <<
       offsets.extent (0) << " < getNodeNumRows() = " << lclNumRows << ".");

    const map_type& rowMap = * (this->getRowMap ());
    const map_type& colMap = * (this->getColMap ());

    // We only use these in debug mode, but since debug mode is a
    // run-time option, they need to exist here.  That's why we create
    // the vector with explicit size zero, to avoid overhead if debug
    // mode is off.
    bool allRowMapDiagEntriesInColMap = true;
    bool allDiagEntriesFound = true;
    bool allOffsetsCorrect = true;
    bool noOtherWeirdness = true;
    using wrong_offsets_type = std::vector<std::pair<LO, size_t> >;
    wrong_offsets_type wrongOffsets(0);

    // mfh 12 Mar 2016: LocalMap works on (CUDA) device.  It has just
    // the subset of Map functionality that we need below.
    auto lclRowMap = rowMap.getLocalMap ();
    auto lclColMap = colMap.getLocalMap ();

    // FIXME (mfh 16 Dec 2015) It's easy to thread-parallelize this
    // setup, at least on the host.  For CUDA, we have to use LocalMap
    // (that comes from each of the two Maps).

    const bool sorted = this->isSorted ();
    if (isFillComplete ()) {
      auto lclGraph = this->getLocalGraph ();
      ::Tpetra::Details::getGraphDiagOffsets (offsets, lclRowMap, lclColMap,
                                              lclGraph.row_map,
                                              lclGraph.entries, sorted);
    }
    else {
      // NOTE (mfh 22 Feb 2017): We have to run this code on host,
      // since the graph is not fill complete.  The previous version
      // of this code assumed UVM; this version does not.
      auto offsets_h = Kokkos::create_mirror_view (offsets);

      for (LO lclRowInd = 0; lclRowInd < lclNumRows; ++lclRowInd) {
        // Find the diagonal entry.  Since the row Map and column Map
        // may differ, we have to compare global row and column
        // indices, not local.
        const GO gblRowInd = lclRowMap.getGlobalElement (lclRowInd);
        const GO gblColInd = gblRowInd;
        const LO lclColInd = lclColMap.getLocalElement (gblColInd);

        if (lclColInd == Tpetra::Details::OrdinalTraits<LO>::invalid ()) {
          allRowMapDiagEntriesInColMap = false;
          offsets_h(lclRowInd) = Tpetra::Details::OrdinalTraits<size_t>::invalid ();
        }
        else {
          const RowInfo rowInfo = this->getRowInfo (lclRowInd);
          if (static_cast<LO> (rowInfo.localRow) == lclRowInd &&
              rowInfo.numEntries > 0) {

            auto colInds = this->getLocalKokkosRowView (rowInfo);
            const size_t hint = 0; // not needed for this algorithm
            const size_t offset =
              KokkosSparse::findRelOffset (colInds, rowInfo.numEntries,
                                           lclColInd, hint, sorted);
            offsets_h(lclRowInd) = offset;

            if (debug_) {
              // Now that we have what we think is an offset, make sure
              // that it really does point to the diagonal entry.  Offsets
              // are _relative_ to each row, not absolute (for the whole
              // (local) graph).
              Teuchos::ArrayView<const LO> lclColInds;
              try {
                this->getLocalRowView (lclRowInd, lclColInds);
              }
              catch (...) {
                noOtherWeirdness = false;
              }
              // Don't continue with error checking if the above failed.
              if (noOtherWeirdness) {
                const size_t numEnt = lclColInds.size ();
                if (offset >= numEnt) {
                  // Offsets are relative to each row, so this means that
                  // the offset is out of bounds.
                  allOffsetsCorrect = false;
                  wrongOffsets.push_back (std::make_pair (lclRowInd, offset));
                } else {
                  const LO actualLclColInd = lclColInds[offset];
                  const GO actualGblColInd = lclColMap.getGlobalElement (actualLclColInd);
                  if (actualGblColInd != gblColInd) {
                    allOffsetsCorrect = false;
                    wrongOffsets.push_back (std::make_pair (lclRowInd, offset));
                  }
                }
              }
            } // debug_
          }
          else { // either row is empty, or something went wrong w/ getRowInfo()
            offsets_h(lclRowInd) = Tpetra::Details::OrdinalTraits<size_t>::invalid ();
            allDiagEntriesFound = false;
          }
        } // whether lclColInd is a valid local column index
      } // for each local row

      Kokkos::deep_copy (offsets, offsets_h);
    } // whether the graph is fill complete

    if (verbose && wrongOffsets.size () != 0) {
      std::ostringstream os;
      os << *prefix << "Wrong offsets: [";
      for (size_t k = 0; k < wrongOffsets.size (); ++k) {
        os << "(" << wrongOffsets[k].first << ","
           << wrongOffsets[k].second << ")";
        if (k + 1 < wrongOffsets.size ()) {
          os << ", ";
        }
      }
      os << "]" << endl;
      std::cerr << os.str();
    }

    if (debug_) {
      using Teuchos::reduceAll;
      using std::endl;
      Teuchos::RCP<const Teuchos::Comm<int> > comm = this->getComm ();
      const bool localSuccess =
        allRowMapDiagEntriesInColMap && allDiagEntriesFound && allOffsetsCorrect;
      const int numResults = 5;
      int lclResults[5];
      lclResults[0] = allRowMapDiagEntriesInColMap ? 1 : 0;
      lclResults[1] = allDiagEntriesFound ? 1 : 0;
      lclResults[2] = allOffsetsCorrect ? 1 : 0;
      lclResults[3] = noOtherWeirdness ? 1 : 0;
      // min-all-reduce will compute least rank of all the processes
      // that didn't succeed.
      lclResults[4] = ! localSuccess ? comm->getRank () : comm->getSize ();

      int gblResults[5];
      gblResults[0] = 0;
      gblResults[1] = 0;
      gblResults[2] = 0;
      gblResults[3] = 0;
      gblResults[4] = 0;
      reduceAll<int, int> (*comm, Teuchos::REDUCE_MIN,
                           numResults, lclResults, gblResults);

      if (gblResults[0] != 1 || gblResults[1] != 1 || gblResults[2] != 1
          || gblResults[3] != 1) {
        std::ostringstream os; // build error message
        os << "Issue(s) that we noticed (on Process " << gblResults[4] << ", "
          "possibly among others): " << endl;
        if (gblResults[0] == 0) {
          os << "  - The column Map does not contain at least one diagonal entry "
            "of the graph." << endl;
        }
        if (gblResults[1] == 0) {
          os << "  - On one or more processes, some row does not contain a "
            "diagonal entry." << endl;
        }
        if (gblResults[2] == 0) {
          os << "  - On one or more processes, some offsets are incorrect."
             << endl;
        }
        if (gblResults[3] == 0) {
          os << "  - One or more processes had some other error."
             << endl;
        }
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(true, std::runtime_error, os.str());
      }
    } // debug_
  }

  namespace { // (anonymous)

    // mfh 21 Jan 2016: This is useful for getLocalDiagOffsets (see
    // below).  The point is to avoid the deep copy between the input
    // Teuchos::ArrayRCP and the internally used Kokkos::View.  We
    // can't use UVM to avoid the deep copy with CUDA, because the
    // ArrayRCP is a host pointer, while the input to the graph's
    // getLocalDiagOffsets method is a device pointer.  Assigning a
    // host pointer to a device pointer is incorrect unless the host
    // pointer points to host pinned memory.  The goal is to get rid
    // of the Teuchos::ArrayRCP overload anyway, so we accept the deep
    // copy for backwards compatibility.
    //
    // We have to use template magic because
    // "staticGraph_->getLocalDiagOffsets(offsetsHosts)" won't compile
    // if device_type::memory_space is not Kokkos::HostSpace (as is
    // the case with CUDA).

    template<class DeviceType,
             const bool memSpaceIsHostSpace =
               std::is_same<typename DeviceType::memory_space,
                            Kokkos::HostSpace>::value>
    struct HelpGetLocalDiagOffsets {};

    template<class DeviceType>
    struct HelpGetLocalDiagOffsets<DeviceType, true> {
      typedef DeviceType device_type;
      typedef Kokkos::View<size_t*, Kokkos::HostSpace,
                           Kokkos::MemoryUnmanaged> device_offsets_type;
      typedef Kokkos::View<size_t*, Kokkos::HostSpace,
                           Kokkos::MemoryUnmanaged> host_offsets_type;

      static device_offsets_type
      getDeviceOffsets (const host_offsets_type& hostOffsets)
      {
        // Host and device are the same; no need to allocate a
        // temporary device View.
        return hostOffsets;
      }

      static void
      copyBackIfNeeded (const host_offsets_type& /* hostOffsets */,
                        const device_offsets_type& /* deviceOffsets */)
      { /* copy back not needed; host and device are the same */ }
    };

    template<class DeviceType>
    struct HelpGetLocalDiagOffsets<DeviceType, false> {
      typedef DeviceType device_type;
      // We have to do a deep copy, since host memory space != device
      // memory space.  Thus, the device View is managed (we need to
      // allocate a temporary device View).
      typedef Kokkos::View<size_t*, device_type> device_offsets_type;
      typedef Kokkos::View<size_t*, Kokkos::HostSpace,
                           Kokkos::MemoryUnmanaged> host_offsets_type;

      static device_offsets_type
      getDeviceOffsets (const host_offsets_type& hostOffsets)
      {
        // Host memory space != device memory space, so we must
        // allocate a temporary device View for the graph.
        return device_offsets_type ("offsets", hostOffsets.extent (0));
      }

      static void
      copyBackIfNeeded (const host_offsets_type& hostOffsets,
                        const device_offsets_type& deviceOffsets)
      {
        Kokkos::deep_copy (hostOffsets, deviceOffsets);
      }
    };
  } // namespace (anonymous)


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  getLocalDiagOffsets (Teuchos::ArrayRCP<size_t>& offsets) const
  {
    typedef LocalOrdinal LO;
    const char tfecfFuncName[] = "getLocalDiagOffsets: ";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
      (! this->hasColMap (), std::runtime_error,
       "The graph does not yet have a column Map.");
    const LO myNumRows = static_cast<LO> (this->getNodeNumRows ());
    if (static_cast<LO> (offsets.size ()) != myNumRows) {
      // NOTE (mfh 21 Jan 2016) This means that the method does not
      // satisfy the strong exception guarantee (no side effects
      // unless successful).
      offsets.resize (myNumRows);
    }

    // mfh 21 Jan 2016: This method unfortunately takes a
    // Teuchos::ArrayRCP, which is host memory.  The graph wants a
    // device pointer.  We can't access host memory from the device;
    // that's the wrong direction for UVM.  (It's the right direction
    // for inefficient host pinned memory, but we don't want to use
    // that here.)  Thus, if device memory space != host memory space,
    // we allocate and use a temporary device View to get the offsets.
    // If the two spaces are equal, the template magic makes the deep
    // copy go away.
    typedef HelpGetLocalDiagOffsets<device_type> helper_type;
    typedef typename helper_type::host_offsets_type host_offsets_type;
    // Unmanaged host View that views the output array.
    host_offsets_type hostOffsets (offsets.getRawPtr (), myNumRows);
    // Allocate temp device View if host != device, else reuse host array.
    auto deviceOffsets = helper_type::getDeviceOffsets (hostOffsets);
    // NOT recursion; this calls the overload that takes a device View.
    this->getLocalDiagOffsets (deviceOffsets);
    helper_type::copyBackIfNeeded (hostOffsets, deviceOffsets);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  supportsRowViews () const {
    return true;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  transferAndFillComplete (Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >& destGraph,
                           const ::Tpetra::Details::Transfer<LocalOrdinal, GlobalOrdinal, Node>& rowTransfer,
                           const Teuchos::RCP<const ::Tpetra::Details::Transfer<LocalOrdinal, GlobalOrdinal, Node> > & domainTransfer,
                           const Teuchos::RCP<const map_type>& domainMap,
                           const Teuchos::RCP<const map_type>& rangeMap,
                           const Teuchos::RCP<Teuchos::ParameterList>& params) const
  {
    using Tpetra::Details::getArrayViewFromDualView;
    using Tpetra::Details::packCrsGraphWithOwningPIDs;
    using Tpetra::Details::unpackAndCombineWithOwningPIDsCount;
    using Tpetra::Details::unpackAndCombineIntoCrsArrays;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::Comm;
    using Teuchos::ParameterList;
    using Teuchos::rcp;
    using Teuchos::RCP;
#ifdef HAVE_TPETRA_MMM_TIMINGS
    using std::string;
    using Teuchos::TimeMonitor;
#endif

    using LO = LocalOrdinal;
    using GO = GlobalOrdinal;
    using NT = node_type;
    using this_type = CrsGraph<LO, GO, NT>;
    using ivector_type = Vector<int, LO, GO, NT>;
    using packet_type = typename this_type::packet_type;

    const char* prefix = "Tpetra::CrsGraph::transferAndFillComplete: ";

#ifdef HAVE_TPETRA_MMM_TIMINGS
    string label;
    if(!params.is_null()) label = params->get("Timer Label", label);
    string prefix2 = string("Tpetra ")+ label + std::string(": CrsGraph TAFC ");
    RCP<TimeMonitor> MM =
      rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix2+string("Pack-1"))));
#endif

    // Make sure that the input argument rowTransfer is either an
    // Import or an Export.  Import and Export are the only two
    // subclasses of Transfer that we defined, but users might
    // (unwisely, for now at least) decide to implement their own
    // subclasses.  Exclude this possibility.
    const import_type* xferAsImport = dynamic_cast<const import_type*>(&rowTransfer);
    const export_type* xferAsExport = dynamic_cast<const export_type*>(&rowTransfer);
    TEUCHOS_TEST_FOR_EXCEPTION(
      xferAsImport == nullptr && xferAsExport == nullptr, std::invalid_argument,
      prefix << "The 'rowTransfer' input argument must be either an Import or "
      "an Export, and its template parameters must match the corresponding "
      "template parameters of the CrsGraph.");

    // Make sure that the input argument domainTransfer is either an
    // Import or an Export.  Import and Export are the only two
    // subclasses of Transfer that we defined, but users might
    // (unwisely, for now at least) decide to implement their own
    // subclasses.  Exclude this possibility.
    Teuchos::RCP<const import_type> xferDomainAsImport =
      Teuchos::rcp_dynamic_cast<const import_type>(domainTransfer);
    Teuchos::RCP<const export_type> xferDomainAsExport =
      Teuchos::rcp_dynamic_cast<const export_type>(domainTransfer);

    if(! domainTransfer.is_null()) {

      TEUCHOS_TEST_FOR_EXCEPTION(
         (xferDomainAsImport.is_null() && xferDomainAsExport.is_null()), std::invalid_argument,
         prefix << "The 'domainTransfer' input argument must be either an "
         "Import or an Export, and its template parameters must match the "
         "corresponding template parameters of the CrsGraph.");

      TEUCHOS_TEST_FOR_EXCEPTION(
         ( xferAsImport != nullptr || ! xferDomainAsImport.is_null() ) &&
         (( xferAsImport != nullptr &&   xferDomainAsImport.is_null() ) ||
          ( xferAsImport == nullptr && ! xferDomainAsImport.is_null() )), std::invalid_argument,
         prefix << "The 'rowTransfer' and 'domainTransfer' input arguments "
         "must be of the same type (either Import or Export).");

      TEUCHOS_TEST_FOR_EXCEPTION(
         ( xferAsExport != nullptr || ! xferDomainAsExport.is_null() ) &&
         (( xferAsExport != nullptr &&   xferDomainAsExport.is_null() ) ||
          ( xferAsExport == nullptr && ! xferDomainAsExport.is_null() )), std::invalid_argument,
         prefix << "The 'rowTransfer' and 'domainTransfer' input arguments "
         "must be of the same type (either Import or Export).");

    } // domainTransfer != null


    // FIXME (mfh 15 May 2014) Wouldn't communication still be needed,
    // if the source Map is not distributed but the target Map is?
    const bool communication_needed = rowTransfer.getSourceMap()->isDistributed();

    //
    // Get the caller's parameters
    //

    bool reverseMode = false; // Are we in reverse mode?
    bool restrictComm = false; // Do we need to restrict the communicator?
    RCP<ParameterList> graphparams; // parameters for the destination graph
    if (! params.is_null()) {
      reverseMode = params->get("Reverse Mode", reverseMode);
      restrictComm = params->get("Restrict Communicator", restrictComm);
      graphparams = sublist(params, "CrsGraph");
    }

    // Get the new domain and range Maps.  We need some of them for error
    // checking, now that we have the reverseMode parameter.
    RCP<const map_type> MyRowMap = reverseMode ?
      rowTransfer.getSourceMap() : rowTransfer.getTargetMap();
    RCP<const map_type> MyColMap; // create this below
    RCP<const map_type> MyDomainMap = ! domainMap.is_null() ? domainMap : getDomainMap();
    RCP<const map_type> MyRangeMap = ! rangeMap.is_null() ? rangeMap : getRangeMap();
    RCP<const map_type> BaseRowMap = MyRowMap;
    RCP<const map_type> BaseDomainMap = MyDomainMap;

    // If the user gave us a nonnull destGraph, then check whether it's
    // "pristine."  That means that it has no entries.
    //
    // FIXME (mfh 15 May 2014) If this is not true on all processes,
    // then this exception test may hang.  It would be better to
    // forward an error flag to the next communication phase.
    if (! destGraph.is_null()) {
      // FIXME (mfh 15 May 2014): The Epetra idiom for checking
      // whether a graph or matrix has no entries on the calling
      // process, is that it is neither locally nor globally indexed.
      // This may change eventually with the Kokkos refactor version
      // of Tpetra, so it would be better just to check the quantity
      // of interest directly.  Note that with the Kokkos refactor
      // version of Tpetra, asking for the total number of entries in
      // a graph or matrix that is not fill complete might require
      // computation (kernel launch), since it is not thread scalable
      // to update a count every time an entry is inserted.
      const bool NewFlag =
        ! destGraph->isLocallyIndexed() && ! destGraph->isGloballyIndexed();
      TEUCHOS_TEST_FOR_EXCEPTION(! NewFlag, std::invalid_argument,
        prefix << "The input argument 'destGraph' is only allowed to be nonnull, "
        "if its graph is empty (neither locally nor globally indexed).");

      // FIXME (mfh 15 May 2014) At some point, we want to change
      // graphs and matrices so that their DistObject Map
      // (this->getMap()) may differ from their row Map.  This will
      // make redistribution for 2-D distributions more efficient.  I
      // hesitate to change this check, because I'm not sure how much
      // the code here depends on getMap() and getRowMap() being the
      // same.
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! destGraph->getRowMap()->isSameAs(*MyRowMap), std::invalid_argument,
        prefix << "The (row) Map of the input argument 'destGraph' is not the "
        "same as the (row) Map specified by the input argument 'rowTransfer'.");

      TEUCHOS_TEST_FOR_EXCEPTION(
        ! destGraph->checkSizes(*this), std::invalid_argument,
        prefix << "You provided a nonnull destination graph, but checkSizes() "
        "indicates that it is not a legal legal target for redistribution from "
        "the source graph (*this).  This may mean that they do not have the "
        "same dimensions.");
    }

    // If forward mode (the default), then *this's (row) Map must be
    // the same as the source Map of the Transfer.  If reverse mode,
    // then *this's (row) Map must be the same as the target Map of
    // the Transfer.
    //
    // FIXME (mfh 15 May 2014) At some point, we want to change graphs
    // and matrices so that their DistObject Map (this->getMap()) may
    // differ from their row Map.  This will make redistribution for
    // 2-D distributions more efficient.  I hesitate to change this
    // check, because I'm not sure how much the code here depends on
    // getMap() and getRowMap() being the same.
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! (reverseMode || getRowMap()->isSameAs(*rowTransfer.getSourceMap())),
      std::invalid_argument, prefix <<
      "rowTransfer->getSourceMap() must match this->getRowMap() in forward mode.");

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! (! reverseMode || getRowMap()->isSameAs(*rowTransfer.getTargetMap())),
      std::invalid_argument, prefix <<
      "rowTransfer->getTargetMap() must match this->getRowMap() in reverse mode.");

    // checks for domainTransfer
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! xferDomainAsImport.is_null() && ! xferDomainAsImport->getTargetMap()->isSameAs(*domainMap),
      std::invalid_argument,
      prefix << "The target map of the 'domainTransfer' input argument must be "
      "the same as the rebalanced domain map 'domainMap'");

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! xferDomainAsExport.is_null() && ! xferDomainAsExport->getSourceMap()->isSameAs(*domainMap),
      std::invalid_argument,
      prefix << "The source map of the 'domainTransfer' input argument must be "
      "the same as the rebalanced domain map 'domainMap'");

    // The basic algorithm here is:
    //
    // 1. Call the moral equivalent of "distor.do" to handle the import.
    // 2. Copy all the Imported and Copy/Permuted data into the raw
    //    CrsGraph pointers, still using GIDs.
    // 3. Call an optimized version of MakeColMap that avoids the
    //    Directory lookups (since the importer knows who owns all the
    //    GIDs) AND reindexes to LIDs.
    // 4. Call expertStaticFillComplete()

    // Get information from the Importer
    const size_t NumSameIDs = rowTransfer.getNumSameIDs();
    ArrayView<const LO> ExportLIDs = reverseMode ?
      rowTransfer.getRemoteLIDs() : rowTransfer.getExportLIDs();
    ArrayView<const LO> RemoteLIDs = reverseMode ?
      rowTransfer.getExportLIDs() : rowTransfer.getRemoteLIDs();
    ArrayView<const LO> PermuteToLIDs = reverseMode ?
      rowTransfer.getPermuteFromLIDs() : rowTransfer.getPermuteToLIDs();
    ArrayView<const LO> PermuteFromLIDs = reverseMode ?
      rowTransfer.getPermuteToLIDs() : rowTransfer.getPermuteFromLIDs();
    Distributor& Distor = rowTransfer.getDistributor();

    // Owning PIDs
    Teuchos::Array<int> SourcePids;
    Teuchos::Array<int> TargetPids;
    int MyPID = getComm()->getRank();

    // Temp variables for sub-communicators
    RCP<const map_type> ReducedRowMap, ReducedColMap,
      ReducedDomainMap, ReducedRangeMap;
    RCP<const Comm<int> > ReducedComm;

    // If the user gave us a null destGraph, then construct the new
    // destination graph.  We will replace its column Map later.
    if (destGraph.is_null()) {
      destGraph = rcp(new this_type(MyRowMap, 0, StaticProfile, graphparams));
    }

    /***************************************************/
    /***** 1) First communicator restriction phase ****/
    /***************************************************/
    if (restrictComm) {
      ReducedRowMap = MyRowMap->removeEmptyProcesses();
      ReducedComm = ReducedRowMap.is_null() ?
        Teuchos::null :
        ReducedRowMap->getComm();
      destGraph->removeEmptyProcessesInPlace(ReducedRowMap);

      ReducedDomainMap = MyRowMap.getRawPtr() == MyDomainMap.getRawPtr() ?
        ReducedRowMap :
        MyDomainMap->replaceCommWithSubset(ReducedComm);
      ReducedRangeMap = MyRowMap.getRawPtr() == MyRangeMap.getRawPtr() ?
        ReducedRowMap :
        MyRangeMap->replaceCommWithSubset(ReducedComm);

      // Reset the "my" maps
      MyRowMap    = ReducedRowMap;
      MyDomainMap = ReducedDomainMap;
      MyRangeMap  = ReducedRangeMap;

      // Update my PID, if we've restricted the communicator
      if (! ReducedComm.is_null()) {
        MyPID = ReducedComm->getRank();
      }
      else {
        MyPID = -2; // For debugging
      }
    }
    else {
      ReducedComm = MyRowMap->getComm();
    }

    /***************************************************/
    /***** 2) From Tpera::DistObject::doTransfer() ****/
    /***************************************************/
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix2+string("ImportSetup"))));
#endif
    // Get the owning PIDs
    RCP<const import_type> MyImporter = getImporter();

    // check whether domain maps of source graph and base domain map is the same
    bool bSameDomainMap = BaseDomainMap->isSameAs(*getDomainMap());

    if (! restrictComm && ! MyImporter.is_null() && bSameDomainMap ) {
      // Same domain map as source graph
      //
      // NOTE: This won't work for restrictComm (because the Import
      // doesn't know the restricted PIDs), though writing an
      // optimized version for that case would be easy (Import an
      // IntVector of the new PIDs).  Might want to add this later.
      Import_Util::getPids(*MyImporter, SourcePids, false);
    }
    else if (restrictComm && ! MyImporter.is_null() && bSameDomainMap) {
      // Same domain map as source graph (restricted communicator)
      // We need one import from the domain to the column map
      ivector_type SourceDomain_pids(getDomainMap(),true);
      ivector_type SourceCol_pids(getColMap());
      // SourceDomain_pids contains the restricted pids
      SourceDomain_pids.putScalar(MyPID);

      SourceCol_pids.doImport(SourceDomain_pids, *MyImporter, INSERT);
      SourcePids.resize(getColMap()->getNodeNumElements());
      SourceCol_pids.get1dCopy(SourcePids());
    }
    else if (MyImporter.is_null() && bSameDomainMap) {
      // Graph has no off-process entries
      SourcePids.resize(getColMap()->getNodeNumElements());
      SourcePids.assign(getColMap()->getNodeNumElements(), MyPID);
    }
    else if ( ! MyImporter.is_null() &&
              ! domainTransfer.is_null() ) {
      // general implementation for rectangular matrices with
      // domain map different than SourceGraph domain map.
      // User has to provide a DomainTransfer object. We need
      // to communications (import/export)

      // TargetDomain_pids lives on the rebalanced new domain map
      ivector_type TargetDomain_pids(domainMap);
      TargetDomain_pids.putScalar(MyPID);

      // SourceDomain_pids lives on the non-rebalanced old domain map
      ivector_type SourceDomain_pids(getDomainMap());

      // SourceCol_pids lives on the non-rebalanced old column map
      ivector_type SourceCol_pids(getColMap());

      if (! reverseMode && ! xferDomainAsImport.is_null() ) {
        SourceDomain_pids.doExport(TargetDomain_pids, *xferDomainAsImport, INSERT);
      }
      else if (reverseMode && ! xferDomainAsExport.is_null() ) {
        SourceDomain_pids.doExport(TargetDomain_pids, *xferDomainAsExport, INSERT);
      }
      else if (! reverseMode && ! xferDomainAsExport.is_null() ) {
        SourceDomain_pids.doImport(TargetDomain_pids, *xferDomainAsExport, INSERT);
      }
      else if (reverseMode && ! xferDomainAsImport.is_null() ) {
        SourceDomain_pids.doImport(TargetDomain_pids, *xferDomainAsImport, INSERT);
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error,
          prefix << "Should never get here!  Please report this bug to a Tpetra developer.");
      }
      SourceCol_pids.doImport(SourceDomain_pids, *MyImporter, INSERT);
      SourcePids.resize(getColMap()->getNodeNumElements());
      SourceCol_pids.get1dCopy(SourcePids());
    }
    else if (BaseDomainMap->isSameAs(*BaseRowMap) &&
             getDomainMap()->isSameAs(*getRowMap())) {
      // We can use the rowTransfer + SourceGraph's Import to find out who owns what.
      ivector_type TargetRow_pids(domainMap);
      ivector_type SourceRow_pids(getRowMap());
      ivector_type SourceCol_pids(getColMap());

      TargetRow_pids.putScalar(MyPID);
      if (! reverseMode && xferAsImport != nullptr) {
        SourceRow_pids.doExport(TargetRow_pids, *xferAsImport, INSERT);
      }
      else if (reverseMode && xferAsExport != nullptr) {
        SourceRow_pids.doExport(TargetRow_pids, *xferAsExport, INSERT);
      }
      else if (! reverseMode && xferAsExport != nullptr) {
        SourceRow_pids.doImport(TargetRow_pids, *xferAsExport, INSERT);
      }
      else if (reverseMode && xferAsImport != nullptr) {
        SourceRow_pids.doImport(TargetRow_pids, *xferAsImport, INSERT);
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error,
          prefix << "Should never get here!  Please report this bug to a Tpetra developer.");
      }
      SourceCol_pids.doImport(SourceRow_pids, *MyImporter, INSERT);
      SourcePids.resize(getColMap()->getNodeNumElements());
      SourceCol_pids.get1dCopy(SourcePids());
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument,
        prefix << "This method only allows either domainMap == getDomainMap(), "
        "or (domainMap == rowTransfer.getTargetMap() and getDomainMap() == getRowMap()).");
    }

    // Tpetra-specific stuff
    size_t constantNumPackets = destGraph->constantNumberOfPackets();
    if (constantNumPackets == 0) {
      destGraph->reallocArraysForNumPacketsPerLid(ExportLIDs.size(),
                                                 RemoteLIDs.size());
    }
    else {
      // There are a constant number of packets per element.  We
      // already know (from the number of "remote" (incoming)
      // elements) how many incoming elements we expect, so we can
      // resize the buffer accordingly.
      const size_t rbufLen = RemoteLIDs.size() * constantNumPackets;
      destGraph->reallocImportsIfNeeded(rbufLen, false, nullptr);
    }

    {
      // packAndPrepare* methods modify numExportPacketsPerLID_.
      destGraph->numExportPacketsPerLID_.modify_host();
      Teuchos::ArrayView<size_t> numExportPacketsPerLID =
        getArrayViewFromDualView(destGraph->numExportPacketsPerLID_);

      // Pack & Prepare w/ owning PIDs
      packCrsGraphWithOwningPIDs(*this, destGraph->exports_,
                                 numExportPacketsPerLID, ExportLIDs,
                                 SourcePids, constantNumPackets, Distor);
    }

    // Do the exchange of remote data.
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix2+string("Transfer"))));
#endif

    if (communication_needed) {
      if (reverseMode) {
        if (constantNumPackets == 0) { // variable number of packets per LID
          // Make sure that host has the latest version, since we're
          // using the version on host.  If host has the latest
          // version, syncing to host does nothing.
          destGraph->numExportPacketsPerLID_.sync_host();
          Teuchos::ArrayView<const size_t> numExportPacketsPerLID =
            getArrayViewFromDualView(destGraph->numExportPacketsPerLID_);
          destGraph->numImportPacketsPerLID_.sync_host();
          Teuchos::ArrayView<size_t> numImportPacketsPerLID =
            getArrayViewFromDualView(destGraph->numImportPacketsPerLID_);
          Distor.doReversePostsAndWaits(numExportPacketsPerLID, 1,
                                         numImportPacketsPerLID);
          size_t totalImportPackets = 0;
          for (Array_size_type i = 0; i < numImportPacketsPerLID.size(); ++i) {
            totalImportPackets += numImportPacketsPerLID[i];
          }

          // Reallocation MUST go before setting the modified flag,
          // because it may clear out the flags.
          destGraph->reallocImportsIfNeeded(totalImportPackets, false, nullptr);
          destGraph->imports_.modify_host();
          Teuchos::ArrayView<packet_type> hostImports =
            getArrayViewFromDualView(destGraph->imports_);
          // This is a legacy host pack/unpack path, so use the host
          // version of exports_.
          destGraph->exports_.sync_host();
          Teuchos::ArrayView<const packet_type> hostExports =
            getArrayViewFromDualView(destGraph->exports_);
          Distor.doReversePostsAndWaits(hostExports,
                                         numExportPacketsPerLID,
                                         hostImports,
                                         numImportPacketsPerLID);
        }
        else { // constant number of packets per LI
          destGraph->imports_.modify_host();
          Teuchos::ArrayView<packet_type> hostImports =
            getArrayViewFromDualView(destGraph->imports_);
          // This is a legacy host pack/unpack path, so use the host
          // version of exports_.
          destGraph->exports_.sync_host();
          Teuchos::ArrayView<const packet_type> hostExports =
            getArrayViewFromDualView(destGraph->exports_);
          Distor.doReversePostsAndWaits(hostExports,
                                         constantNumPackets,
                                         hostImports);
        }
      }
      else { // forward mode (the default)
        if (constantNumPackets == 0) { // variable number of packets per LID
          // Make sure that host has the latest version, since we're
          // using the version on host.  If host has the latest
          // version, syncing to host does nothing.
          destGraph->numExportPacketsPerLID_.sync_host();
          Teuchos::ArrayView<const size_t> numExportPacketsPerLID =
            getArrayViewFromDualView(destGraph->numExportPacketsPerLID_);
          destGraph->numImportPacketsPerLID_.sync_host();
          Teuchos::ArrayView<size_t> numImportPacketsPerLID =
            getArrayViewFromDualView(destGraph->numImportPacketsPerLID_);
          Distor.doPostsAndWaits(numExportPacketsPerLID, 1,
                                  numImportPacketsPerLID);
          size_t totalImportPackets = 0;
          for (Array_size_type i = 0; i < numImportPacketsPerLID.size(); ++i) {
            totalImportPackets += numImportPacketsPerLID[i];
          }

          // Reallocation MUST go before setting the modified flag,
          // because it may clear out the flags.
          destGraph->reallocImportsIfNeeded(totalImportPackets, false, nullptr);
          destGraph->imports_.modify_host();
          Teuchos::ArrayView<packet_type> hostImports =
            getArrayViewFromDualView(destGraph->imports_);
          // This is a legacy host pack/unpack path, so use the host
          // version of exports_.
          destGraph->exports_.sync_host();
          Teuchos::ArrayView<const packet_type> hostExports =
            getArrayViewFromDualView(destGraph->exports_);
          Distor.doPostsAndWaits(hostExports,
                                  numExportPacketsPerLID,
                                  hostImports,
                                  numImportPacketsPerLID);
        }
        else { // constant number of packets per LID
          destGraph->imports_.modify_host();
          Teuchos::ArrayView<packet_type> hostImports =
            getArrayViewFromDualView(destGraph->imports_);
          // This is a legacy host pack/unpack path, so use the host
          // version of exports_.
          destGraph->exports_.sync_host();
          Teuchos::ArrayView<const packet_type> hostExports =
            getArrayViewFromDualView(destGraph->exports_);
          Distor.doPostsAndWaits(hostExports,
                                  constantNumPackets,
                                  hostImports);
        }
      }
    }

    /*********************************************************************/
    /**** 3) Copy all of the Same/Permute/Remote data into CSR_arrays ****/
    /*********************************************************************/

#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix2+string("Unpack-1"))));
#endif

    // Backwards compatibility measure.  We'll use this again below.
    destGraph->numImportPacketsPerLID_.sync_host();
    Teuchos::ArrayView<const size_t> numImportPacketsPerLID =
      getArrayViewFromDualView(destGraph->numImportPacketsPerLID_);
    destGraph->imports_.sync_host();
    Teuchos::ArrayView<const packet_type> hostImports =
      getArrayViewFromDualView(destGraph->imports_);
    size_t mynnz =
      unpackAndCombineWithOwningPIDsCount(*this, RemoteLIDs, hostImports,
                                           numImportPacketsPerLID,
                                           constantNumPackets, Distor, INSERT,
                                           NumSameIDs, PermuteToLIDs, PermuteFromLIDs);
    size_t N = BaseRowMap->getNodeNumElements();

    // Allocations
    ArrayRCP<size_t> CSR_rowptr(N+1);
    ArrayRCP<GO> CSR_colind_GID;
    ArrayRCP<LO> CSR_colind_LID;
    CSR_colind_GID.resize(mynnz);

    // If LO and GO are the same, we can reuse memory when
    // converting the column indices from global to local indices.
    if (typeid(LO) == typeid(GO)) {
      CSR_colind_LID = Teuchos::arcp_reinterpret_cast<LO>(CSR_colind_GID);
    }
    else {
      CSR_colind_LID.resize(mynnz);
    }

    // FIXME (mfh 15 May 2014) Why can't we abstract this out as an
    // unpackAndCombine method on a "CrsArrays" object?  This passing
    // in a huge list of arrays is icky.  Can't we have a bit of an
    // abstraction?  Implementing a concrete DistObject subclass only
    // takes five methods.
    unpackAndCombineIntoCrsArrays(*this, RemoteLIDs, hostImports,
                                  numImportPacketsPerLID, constantNumPackets,
                                  Distor, INSERT, NumSameIDs, PermuteToLIDs,
                                  PermuteFromLIDs, N, mynnz, MyPID,
                                  CSR_rowptr(), CSR_colind_GID(),
                                  SourcePids(), TargetPids);

    /**************************************************************/
    /**** 4) Call Optimized MakeColMap w/ no Directory Lookups ****/
    /**************************************************************/
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix2+string("Unpack-2"))));
#endif
    // Call an optimized version of makeColMap that avoids the
    // Directory lookups (since the Import object knows who owns all
    // the GIDs).
    Teuchos::Array<int> RemotePids;
    Import_Util::lowCommunicationMakeColMapAndReindex(CSR_rowptr(),
                                                       CSR_colind_LID(),
                                                       CSR_colind_GID(),
                                                       BaseDomainMap,
                                                       TargetPids, RemotePids,
                                                       MyColMap);

    /*******************************************************/
    /**** 4) Second communicator restriction phase      ****/
    /*******************************************************/
    if (restrictComm) {
      ReducedColMap = (MyRowMap.getRawPtr() == MyColMap.getRawPtr()) ?
        ReducedRowMap :
        MyColMap->replaceCommWithSubset(ReducedComm);
      MyColMap = ReducedColMap; // Reset the "my" maps
    }

    // Replace the col map
    destGraph->replaceColMap(MyColMap);

    // Short circuit if the processor is no longer in the communicator
    //
    // NOTE: Epetra replaces modifies all "removed" processes so they
    // have a dummy (serial) Map that doesn't touch the original
    // communicator.  Duplicating that here might be a good idea.
    if (ReducedComm.is_null()) {
      return;
    }

    /***************************************************/
    /**** 5) Sort                                   ****/
    /***************************************************/
    if ((! reverseMode && xferAsImport != nullptr) ||
        (reverseMode && xferAsExport != nullptr)) {
      Import_Util::sortCrsEntries(CSR_rowptr(),
                                   CSR_colind_LID());
    }
    else if ((! reverseMode && xferAsExport != nullptr) ||
             (reverseMode && xferAsImport != nullptr)) {
      Import_Util::sortAndMergeCrsEntries(CSR_rowptr(),
                                           CSR_colind_LID());
      if (CSR_rowptr[N] != mynnz) {
        CSR_colind_LID.resize(CSR_rowptr[N]);
      }
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        prefix << "Should never get here!  Please report this bug to a Tpetra developer.");
    }
    /***************************************************/
    /**** 6) Reset the colmap and the arrays        ****/
    /***************************************************/

    // Call constructor for the new graph (restricted as needed)
    //
    destGraph->setAllIndices(CSR_rowptr, CSR_colind_LID);

    /***************************************************/
    /**** 7) Build Importer & Call ESFC             ****/
    /***************************************************/
    // Pre-build the importer using the existing PIDs
    Teuchos::ParameterList esfc_params;
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix2+string("CreateImporter"))));
#endif
    RCP<import_type> MyImport = rcp(new import_type(MyDomainMap, MyColMap, RemotePids));
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix2+string("ESFC"))));

    esfc_params.set("Timer Label",prefix + std::string("TAFC"));
#endif
    if(!params.is_null())
      esfc_params.set("compute global constants",params->get("compute global constants",true));

    destGraph->expertStaticFillComplete(MyDomainMap, MyRangeMap,
        MyImport, Teuchos::null, rcp(&esfc_params,false));

  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  importAndFillComplete(Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >& destGraph,
                         const import_type& importer,
                         const Teuchos::RCP<const map_type>& domainMap,
                         const Teuchos::RCP<const map_type>& rangeMap,
                         const Teuchos::RCP<Teuchos::ParameterList>& params) const
  {
    transferAndFillComplete(destGraph, importer, Teuchos::null, domainMap, rangeMap, params);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  importAndFillComplete(Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >& destGraph,
                         const import_type& rowImporter,
                         const import_type& domainImporter,
                         const Teuchos::RCP<const map_type>& domainMap,
                         const Teuchos::RCP<const map_type>& rangeMap,
                         const Teuchos::RCP<Teuchos::ParameterList>& params) const
  {
    transferAndFillComplete(destGraph, rowImporter, Teuchos::rcpFromRef(domainImporter), domainMap, rangeMap, params);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  exportAndFillComplete(Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >& destGraph,
                         const export_type& exporter,
                         const Teuchos::RCP<const map_type>& domainMap,
                         const Teuchos::RCP<const map_type>& rangeMap,
                         const Teuchos::RCP<Teuchos::ParameterList>& params) const
  {
    transferAndFillComplete(destGraph, exporter, Teuchos::null, domainMap, rangeMap, params);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  exportAndFillComplete(Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >& destGraph,
                         const export_type& rowExporter,
                         const export_type& domainExporter,
                         const Teuchos::RCP<const map_type>& domainMap,
                         const Teuchos::RCP<const map_type>& rangeMap,
                         const Teuchos::RCP<Teuchos::ParameterList>& params) const
  {
    transferAndFillComplete(destGraph, rowExporter, Teuchos::rcpFromRef(domainExporter), domainMap, rangeMap, params);
  }


  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  swap(CrsGraph<LocalOrdinal, GlobalOrdinal, Node>& graph)
  {
    std::swap(graph.rowMap_, this->rowMap_);
    std::swap(graph.colMap_, this->colMap_);
    std::swap(graph.rangeMap_, this->rangeMap_);
    std::swap(graph.domainMap_, this->domainMap_);

    std::swap(graph.importer_, this->importer_);
    std::swap(graph.exporter_, this->exporter_);

    std::swap(graph.lclGraph_, this->lclGraph_);

    std::swap(graph.nodeMaxNumRowEntries_, this->nodeMaxNumRowEntries_);

    std::swap(graph.globalNumEntries_, this->globalNumEntries_);
    std::swap(graph.globalMaxNumRowEntries_, this->globalMaxNumRowEntries_);

    std::swap(graph.numAllocForAllRows_, this->numAllocForAllRows_);

    std::swap(graph.k_rowPtrs_, this->k_rowPtrs_);

    std::swap(graph.k_lclInds1D_, this->k_lclInds1D_);
    std::swap(graph.k_gblInds1D_, this->k_gblInds1D_);

    std::swap(graph.storageStatus_, this->storageStatus_);

    std::swap(graph.indicesAreAllocated_, this->indicesAreAllocated_);
    std::swap(graph.indicesAreLocal_, this->indicesAreLocal_);
    std::swap(graph.indicesAreGlobal_, this->indicesAreGlobal_);
    std::swap(graph.fillComplete_, this->fillComplete_);
    std::swap(graph.indicesAreSorted_, this->indicesAreSorted_);
    std::swap(graph.noRedundancies_, this->noRedundancies_);
    std::swap(graph.haveLocalConstants_, this->haveLocalConstants_);
    std::swap(graph.haveGlobalConstants_, this->haveGlobalConstants_);

    std::swap(graph.sortGhostsAssociatedWithEachProcessor_, this->sortGhostsAssociatedWithEachProcessor_);

    std::swap(graph.k_numAllocPerRow_, this->k_numAllocPerRow_);  // View
    std::swap(graph.k_numRowEntries_, this->k_numRowEntries_);    // View
    std::swap(graph.nonlocals_, this->nonlocals_);                // std::map
  }


  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
  isIdenticalTo(const CrsGraph<LocalOrdinal, GlobalOrdinal, Node> & graph) const
  {
    auto compare_nonlocals = [&] (const nonlocals_type & m1, const nonlocals_type & m2) {
      bool output = true;
      output = m1.size() == m2.size() ? output : false;
      for(auto & it_m: m1)
      {
        size_t key = it_m.first;
        output = m2.find(key) != m2.end() ? output : false;
        if(output)
        {
          auto v1 = m1.find(key)->second;
          auto v2 = m2.find(key)->second;
          std::sort(v1.begin(), v1.end());
          std::sort(v2.begin(), v2.end());

          output = v1.size() == v2.size() ? output : false;
          for(size_t i=0; output && i<v1.size(); i++)
          {
            output = v1[i]==v2[i] ? output : false;
          }
        }
      }
      return output;
    };

    bool output = true;

    output = this->rowMap_->isSameAs( *(graph.rowMap_) ) ? output : false;
    output = this->colMap_->isSameAs( *(graph.colMap_) ) ? output : false;
    output = this->rangeMap_->isSameAs( *(graph.rangeMap_) ) ? output : false;
    output = this->domainMap_->isSameAs( *(graph.domainMap_) ) ? output : false;

    output = this->nodeMaxNumRowEntries_ == graph.nodeMaxNumRowEntries_ ? output : false;

    output = this->globalNumEntries_ == graph.globalNumEntries_ ? output : false;
    output = this->globalMaxNumRowEntries_ == graph.globalMaxNumRowEntries_ ? output : false;

    output = this->numAllocForAllRows_ == graph.numAllocForAllRows_ ? output : false;

    output = this->storageStatus_ == graph.storageStatus_ ? output : false;  // EStorageStatus is an enum

    output = this->indicesAreAllocated_ == graph.indicesAreAllocated_ ? output : false;
    output = this->indicesAreLocal_ == graph.indicesAreLocal_ ? output : false;
    output = this->indicesAreGlobal_ == graph.indicesAreGlobal_ ? output : false;
    output = this->fillComplete_ == graph.fillComplete_ ? output : false;
    output = this->indicesAreSorted_ == graph.indicesAreSorted_ ? output : false;
    output = this->noRedundancies_ == graph.noRedundancies_ ? output : false;
    output = this->haveLocalConstants_ == graph.haveLocalConstants_ ? output : false;
    output = this->haveGlobalConstants_ == graph.haveGlobalConstants_ ? output : false;
    output = this->sortGhostsAssociatedWithEachProcessor_ == this->sortGhostsAssociatedWithEachProcessor_ ? output : false;

    // Compare nonlocals_ -- std::map<GlobalOrdinal, std::vector<GlobalOrdinal> >
    // nonlocals_ isa std::map<GO, std::vector<GO> >
    output = compare_nonlocals(this->nonlocals_, graph.nonlocals_) ? output : false;

    // Compare k_numAllocPerRow_ isa Kokkos::View::HostMirror
    // - since this is a HostMirror type, it should be in host memory already
    output = this->k_numAllocPerRow_.extent(0) == graph.k_numAllocPerRow_.extent(0) ? output : false;
    if(output && this->k_numAllocPerRow_.extent(0) > 0)
    {
      for(size_t i=0; output && i<this->k_numAllocPerRow_.extent(0); i++)
        output = this->k_numAllocPerRow_(i) == graph.k_numAllocPerRow_(i) ? output : false;
    }

    // Compare k_numRowEntries_ isa Kokkos::View::HostMirror
    // - since this is a HostMirror type, it should be in host memory already
    output = this->k_numRowEntries_.extent(0) == graph.k_numRowEntries_.extent(0) ? output : false;
    if(output && this->k_numRowEntries_.extent(0) > 0)
    {
      for(size_t i = 0; output && i < this->k_numRowEntries_.extent(0); i++)
        output = this->k_numRowEntries_(i) == graph.k_numRowEntries_(i) ? output : false;
    }

    // Compare this->k_rowPtrs_ isa Kokkos::View<LocalOrdinal*, ...>
    output = this->k_rowPtrs_.extent(0) == graph.k_rowPtrs_.extent(0) ? output : false;
    if(output && this->k_rowPtrs_.extent(0) > 0)
    {
      typename local_graph_type::row_map_type::const_type::HostMirror k_rowPtrs_host_this = Kokkos::create_mirror_view(this->k_rowPtrs_);
      typename local_graph_type::row_map_type::const_type::HostMirror k_rowPtrs_host_graph= Kokkos::create_mirror_view(graph.k_rowPtrs_);
      Kokkos::deep_copy(k_rowPtrs_host_this, this->k_rowPtrs_);
      Kokkos::deep_copy(k_rowPtrs_host_graph, graph.k_rowPtrs_);
      for(size_t i=0; output && i<k_rowPtrs_host_this.extent(0); i++)
        output = k_rowPtrs_host_this(i) == k_rowPtrs_host_graph(i) ? output : false;
    }

    // Compare k_lclInds1D_ isa Kokkos::View<LocalOrdinal*, ...>
    output = this->k_lclInds1D_.extent(0) == graph.k_lclInds1D_.extent(0) ? output : false;
    if(output && this->k_lclInds1D_.extent(0) > 0)
    {
      typename local_graph_type::entries_type::non_const_type::HostMirror k_lclInds1D_host_this = Kokkos::create_mirror_view(this->k_lclInds1D_);
      typename local_graph_type::entries_type::non_const_type::HostMirror k_lclInds1D_host_graph= Kokkos::create_mirror_view(graph.k_lclInds1D_);
      Kokkos::deep_copy(k_lclInds1D_host_this, this->k_lclInds1D_);
      Kokkos::deep_copy(k_lclInds1D_host_graph, graph.k_lclInds1D_);
      for(size_t i=0; output && i < k_lclInds1D_host_this.extent(0); i++)
        output = k_lclInds1D_host_this(i) == k_lclInds1D_host_graph(i) ? output : false;
    }

    // Compare k_gblInds1D_ isa Kokkos::View<GlobalOrdinal*, ...>
    output = this->k_gblInds1D_.extent(0) == graph.k_gblInds1D_.extent(0) ? output : false;
    if(output && this->k_gblInds1D_.extent(0) > 0)
    {
      typename t_GlobalOrdinal_1D::HostMirror k_gblInds1D_host_this  = Kokkos::create_mirror_view(this->k_gblInds1D_);
      typename t_GlobalOrdinal_1D::HostMirror k_gblInds1D_host_graph = Kokkos::create_mirror_view(graph.k_gblInds1D_);
      Kokkos::deep_copy(k_gblInds1D_host_this, this->k_gblInds1D_);
      Kokkos::deep_copy(k_gblInds1D_host_graph, graph.k_gblInds1D_);
      for(size_t i=0; output && i<k_gblInds1D_host_this.extent(0); i++)
        output = k_gblInds1D_host_this(i) == k_gblInds1D_host_graph(i) ? output : false;
    }

    // Check lclGraph_      // isa Kokkos::StaticCrsGraph<LocalOrdinal, Kokkos::LayoutLeft, execution_space>
    // Kokkos::StaticCrsGraph has 3 data members in it:
    //   Kokkos::View<size_type*, ...> row_map            (local_graph_type::row_map_type)
    //   Kokkos::View<data_type*, ...> entries            (local_graph_type::entries_type)
    //   Kokkos::View<size_type*, ...> row_block_offsets  (local_graph_type::row_block_type)
    // There is currently no Kokkos::StaticCrsGraph comparison function that's built-in, so we will just compare
    // the three data items here. This can be replaced if Kokkos ever puts in its own comparison routine.
    output = this->lclGraph_.row_map.extent(0) == graph.lclGraph_.row_map.extent(0) ? output : false;
    if(output && this->lclGraph_.row_map.extent(0) > 0)
    {
      typename local_graph_type::row_map_type::HostMirror lclGraph_rowmap_host_this  = Kokkos::create_mirror_view(this->lclGraph_.row_map);
      typename local_graph_type::row_map_type::HostMirror lclGraph_rowmap_host_graph = Kokkos::create_mirror_view(graph.lclGraph_.row_map);
      Kokkos::deep_copy(lclGraph_rowmap_host_this, this->lclGraph_.row_map);
      Kokkos::deep_copy(lclGraph_rowmap_host_graph, graph.lclGraph_.row_map);
      for(size_t i=0; output && i<lclGraph_rowmap_host_this.extent(0); i++)
        output = lclGraph_rowmap_host_this(i) == lclGraph_rowmap_host_graph(i) ? output : false;
    }

    output = this->lclGraph_.entries.extent(0) == graph.lclGraph_.entries.extent(0) ? output : false;
    if(output && this->lclGraph_.entries.extent(0) > 0)
    {
      typename local_graph_type::entries_type::HostMirror lclGraph_entries_host_this = Kokkos::create_mirror_view(this->lclGraph_.entries);
      typename local_graph_type::entries_type::HostMirror lclGraph_entries_host_graph = Kokkos::create_mirror_view(graph.lclGraph_.entries);
      Kokkos::deep_copy(lclGraph_entries_host_this, this->lclGraph_.entries);
      Kokkos::deep_copy(lclGraph_entries_host_graph, graph.lclGraph_.entries);
      for(size_t i=0; output && i<lclGraph_entries_host_this.extent(0); i++)
        output = lclGraph_entries_host_this(i) == lclGraph_entries_host_graph(i) ? output : false;
    }

    output = this->lclGraph_.row_block_offsets.extent(0) == graph.lclGraph_.row_block_offsets.extent(0) ? output : false;
    if(output && this->lclGraph_.row_block_offsets.extent(0) > 0)
    {
      typename local_graph_type::row_block_type::HostMirror lclGraph_rbo_host_this = Kokkos::create_mirror_view(this->lclGraph_.row_block_offsets);
      typename local_graph_type::row_block_type::HostMirror lclGraph_rbo_host_graph = Kokkos::create_mirror_view(graph.lclGraph_.row_block_offsets);
      Kokkos::deep_copy(lclGraph_rbo_host_this, this->lclGraph_.row_block_offsets);
      Kokkos::deep_copy(lclGraph_rbo_host_graph, graph.lclGraph_.row_block_offsets);
      for(size_t i=0; output && i < lclGraph_rbo_host_this.extent(0); i++)
        output = lclGraph_rbo_host_this(i) == lclGraph_rbo_host_graph(i) ? output : false;
    }

    // For the Importer and Exporter, we shouldn't need to explicitly check them since
    // they will be consistent with the maps.
    // Note: importer_  isa Teuchos::RCP<const import_type>
    //       exporter_  isa Teuchos::RCP<const export_type>

    return output;
  }



} // namespace Tpetra

//
// Explicit instantiation macros
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_CRSGRAPH_IMPORT_AND_FILL_COMPLETE_INSTANT(LO,GO,NODE) \
  template<>                                                                        \
  Teuchos::RCP<CrsGraph<LO,GO,NODE> >                        \
  importAndFillCompleteCrsGraph(const Teuchos::RCP<const CrsGraph<LO,GO,NODE> >& sourceGraph, \
                                  const Import<CrsGraph<LO,GO,NODE>::local_ordinal_type,  \
                                               CrsGraph<LO,GO,NODE>::global_ordinal_type,  \
                                               CrsGraph<LO,GO,NODE>::node_type>& importer, \
                                  const Teuchos::RCP<const Map<CrsGraph<LO,GO,NODE>::local_ordinal_type,      \
                                                               CrsGraph<LO,GO,NODE>::global_ordinal_type,     \
                                                               CrsGraph<LO,GO,NODE>::node_type> >& domainMap, \
                                  const Teuchos::RCP<const Map<CrsGraph<LO,GO,NODE>::local_ordinal_type,      \
                                                               CrsGraph<LO,GO,NODE>::global_ordinal_type,     \
                                                               CrsGraph<LO,GO,NODE>::node_type> >& rangeMap,  \
                                                               const Teuchos::RCP<Teuchos::ParameterList>& params);

#define TPETRA_CRSGRAPH_IMPORT_AND_FILL_COMPLETE_INSTANT_TWO(LO,GO,NODE) \
  template<>                                                                        \
  Teuchos::RCP<CrsGraph<LO,GO,NODE> >                        \
  importAndFillCompleteCrsGraph(const Teuchos::RCP<const CrsGraph<LO,GO,NODE> >& sourceGraph, \
                                  const Import<CrsGraph<LO,GO,NODE>::local_ordinal_type,  \
                                               CrsGraph<LO,GO,NODE>::global_ordinal_type,  \
                                               CrsGraph<LO,GO,NODE>::node_type>& rowImporter, \
                                  const Import<CrsGraph<LO,GO,NODE>::local_ordinal_type,  \
                                               CrsGraph<LO,GO,NODE>::global_ordinal_type,  \
                                               CrsGraph<LO,GO,NODE>::node_type>& domainImporter, \
                                  const Teuchos::RCP<const Map<CrsGraph<LO,GO,NODE>::local_ordinal_type,      \
                                                               CrsGraph<LO,GO,NODE>::global_ordinal_type,     \
                                                               CrsGraph<LO,GO,NODE>::node_type> >& domainMap, \
                                  const Teuchos::RCP<const Map<CrsGraph<LO,GO,NODE>::local_ordinal_type,      \
                                                               CrsGraph<LO,GO,NODE>::global_ordinal_type,     \
                                                               CrsGraph<LO,GO,NODE>::node_type> >& rangeMap,  \
                                                               const Teuchos::RCP<Teuchos::ParameterList>& params);


#define TPETRA_CRSGRAPH_EXPORT_AND_FILL_COMPLETE_INSTANT(LO,GO,NODE) \
  template<>                                                                        \
  Teuchos::RCP<CrsGraph<LO,GO,NODE> >                        \
  exportAndFillCompleteCrsGraph(const Teuchos::RCP<const CrsGraph<LO,GO,NODE> >& sourceGraph, \
                                  const Export<CrsGraph<LO,GO,NODE>::local_ordinal_type,  \
                                               CrsGraph<LO,GO,NODE>::global_ordinal_type,  \
                                               CrsGraph<LO,GO,NODE>::node_type>& exporter, \
                                  const Teuchos::RCP<const Map<CrsGraph<LO,GO,NODE>::local_ordinal_type,      \
                                                               CrsGraph<LO,GO,NODE>::global_ordinal_type,     \
                                                               CrsGraph<LO,GO,NODE>::node_type> >& domainMap, \
                                  const Teuchos::RCP<const Map<CrsGraph<LO,GO,NODE>::local_ordinal_type,      \
                                                               CrsGraph<LO,GO,NODE>::global_ordinal_type,     \
                                                               CrsGraph<LO,GO,NODE>::node_type> >& rangeMap,  \
                                                               const Teuchos::RCP<Teuchos::ParameterList>& params);

#define TPETRA_CRSGRAPH_EXPORT_AND_FILL_COMPLETE_INSTANT_TWO(LO,GO,NODE) \
  template<>                                                                        \
  Teuchos::RCP<CrsGraph<LO,GO,NODE> >                        \
  exportAndFillCompleteCrsGraph(const Teuchos::RCP<const CrsGraph<LO,GO,NODE> >& sourceGraph, \
                                  const Export<CrsGraph<LO,GO,NODE>::local_ordinal_type,  \
                                               CrsGraph<LO,GO,NODE>::global_ordinal_type,  \
                                               CrsGraph<LO,GO,NODE>::node_type>& rowExporter, \
                                  const Export<CrsGraph<LO,GO,NODE>::local_ordinal_type,  \
                                               CrsGraph<LO,GO,NODE>::global_ordinal_type,  \
                                               CrsGraph<LO,GO,NODE>::node_type>& domainExporter, \
                                  const Teuchos::RCP<const Map<CrsGraph<LO,GO,NODE>::local_ordinal_type,      \
                                                               CrsGraph<LO,GO,NODE>::global_ordinal_type,     \
                                                               CrsGraph<LO,GO,NODE>::node_type> >& domainMap, \
                                  const Teuchos::RCP<const Map<CrsGraph<LO,GO,NODE>::local_ordinal_type,      \
                                                               CrsGraph<LO,GO,NODE>::global_ordinal_type,     \
                                                               CrsGraph<LO,GO,NODE>::node_type> >& rangeMap,  \
                                                               const Teuchos::RCP<Teuchos::ParameterList>& params);


#define TPETRA_CRSGRAPH_INSTANT( LO, GO, NODE ) \
  template class CrsGraph<LO, GO, NODE>; \
  TPETRA_CRSGRAPH_IMPORT_AND_FILL_COMPLETE_INSTANT(LO,GO,NODE) \
  TPETRA_CRSGRAPH_EXPORT_AND_FILL_COMPLETE_INSTANT(LO,GO,NODE) \
  TPETRA_CRSGRAPH_IMPORT_AND_FILL_COMPLETE_INSTANT_TWO(LO,GO,NODE) \
  TPETRA_CRSGRAPH_EXPORT_AND_FILL_COMPLETE_INSTANT_TWO(LO,GO,NODE)


#endif // TPETRA_CRSGRAPH_DEF_HPP
