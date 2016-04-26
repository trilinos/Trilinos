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

#ifndef TPETRA_DISTOBJECT_DEF_HPP
#define TPETRA_DISTOBJECT_DEF_HPP

/// \file Tpetra_DistObject_def.hpp
/// \brief Definition of the Tpetra::DistObject class
///
/// If you want to use Tpetra::DistObject, include
/// "Tpetra_DistObject.hpp" (a file which CMake generates and installs
/// for you).  If you only want the declaration of Tpetra::DistObject,
/// include "Tpetra_DistObject_decl.hpp".

#include "Tpetra_Distributor.hpp"

namespace { // (anonymous)

  // Get a Teuchos::ArrayView which views the host Kokkos::View of the
  // input 1-D Kokkos::DualView.
  template<class DualViewType>
  Teuchos::ArrayView<typename DualViewType::t_dev::value_type>
  getArrayViewFromDualView (const DualViewType& x)
  {
    static_assert (static_cast<int> (DualViewType::t_dev::rank) == 1,
                   "The input DualView must have rank 1.");
    TEUCHOS_TEST_FOR_EXCEPTION
      (x.template need_sync<Kokkos::HostSpace> (), std::logic_error, "The "
       "input Kokkos::DualView was most recently modified on device, but this "
       "function needs the host view of the data to be the most recently "
       "modified.");

    auto x_host = x.template view<Kokkos::HostSpace> ();
    typedef typename DualViewType::t_dev::value_type value_type;
    return Teuchos::ArrayView<value_type> (x_host.ptr_on_device (),
                                           x_host.dimension_0 ());
  }

} // namespace (anonymous)

namespace Tpetra {

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  DistObject (const Teuchos::RCP<const map_type>& map) :
    map_ (map)
  {
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
    using Teuchos::RCP;
    using Teuchos::Time;
    using Teuchos::TimeMonitor;

    RCP<Time> doXferTimer =
      TimeMonitor::lookupCounter ("Tpetra::DistObject::doTransfer");
    if (doXferTimer.is_null ()) {
      doXferTimer =
        TimeMonitor::getNewCounter ("Tpetra::DistObject::doTransfer");
    }
    doXferTimer_ = doXferTimer;

    RCP<Time> copyAndPermuteTimer =
      TimeMonitor::lookupCounter ("Tpetra::DistObject::copyAndPermute");
    if (copyAndPermuteTimer.is_null ()) {
      copyAndPermuteTimer =
        TimeMonitor::getNewCounter ("Tpetra::DistObject::copyAndPermute");
    }
    copyAndPermuteTimer_ = copyAndPermuteTimer;

    RCP<Time> packAndPrepareTimer =
      TimeMonitor::lookupCounter ("Tpetra::DistObject::packAndPrepare");
    if (packAndPrepareTimer.is_null ()) {
      packAndPrepareTimer =
        TimeMonitor::getNewCounter ("Tpetra::DistObject::packAndPrepare");
    }
    packAndPrepareTimer_ = packAndPrepareTimer;

    RCP<Time> doPostsAndWaitsTimer =
      TimeMonitor::lookupCounter ("Tpetra::DistObject::doPostsAndWaits");
    if (doPostsAndWaitsTimer.is_null ()) {
      doPostsAndWaitsTimer =
        TimeMonitor::getNewCounter ("Tpetra::DistObject::doPostsAndWaits");
    }
    doPostsAndWaitsTimer_ = doPostsAndWaitsTimer;

    RCP<Time> unpackAndCombineTimer =
      TimeMonitor::lookupCounter ("Tpetra::DistObject::unpackAndCombine");
    if (unpackAndCombineTimer.is_null ()) {
      unpackAndCombineTimer =
        TimeMonitor::getNewCounter ("Tpetra::DistObject::unpackAndCombine");
    }
    unpackAndCombineTimer_ = unpackAndCombineTimer;
#endif // HAVE_TPETRA_TRANSFER_TIMERS
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  DistObject (const DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>& rhs)
    : map_ (rhs.map_)
  {}

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  ~DistObject ()
  {}

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  std::string
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  description () const
  {
    using Teuchos::TypeNameTraits;

    std::ostringstream os;
    os << "\"Tpetra::DistObject\": {"
       << "Packet: " << TypeNameTraits<packet_type>::name ()
       << ", LocalOrdinal: " << TypeNameTraits<local_ordinal_type>::name ()
       << ", GlobalOrdinal: " << TypeNameTraits<global_ordinal_type>::name ()
       << ", Node: " << TypeNameTraits<Node>::name ();
    if (this->getObjectLabel () != "") {
      os << "Label: \"" << this->getObjectLabel () << "\"";
    }
    os << "}";
    return os.str ();
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    using Teuchos::rcpFromRef;
    using Teuchos::TypeNameTraits;
    using std::endl;
    const Teuchos::EVerbosityLevel vl = (verbLevel == Teuchos::VERB_DEFAULT) ?
      Teuchos::VERB_LOW : verbLevel;
    Teuchos::RCP<const Teuchos::Comm<int> > comm = this->getMap ()->getComm ();
    const int myRank = comm.is_null () ? 0 : comm->getRank ();
    const int numProcs = comm.is_null () ? 1 : comm->getSize ();

    if (vl != Teuchos::VERB_NONE) {
      Teuchos::OSTab tab0 (out);
      if (myRank == 0) {
        out << "\"Tpetra::DistObject\":" << endl;
      }
      Teuchos::OSTab tab1 (out);
      if (myRank == 0) {
        out << "Template parameters:" << endl;
        {
          Teuchos::OSTab tab2 (out);
          out << "Packet: " << TypeNameTraits<packet_type>::name () << endl
              << "LocalOrdinal: " << TypeNameTraits<local_ordinal_type>::name () << endl
              << "GlobalOrdinal: " << TypeNameTraits<global_ordinal_type>::name () << endl
              << "Node: " << TypeNameTraits<node_type>::name () << endl;
        }
        if (this->getObjectLabel () != "") {
          out << "Label: \"" << this->getObjectLabel () << "\"" << endl;
        }
      } // if myRank == 0

      // Describe the Map.
      {
        if (myRank == 0) {
          out << "Map:" << endl;
        }
        Teuchos::OSTab tab2 (out);
        map_->describe (out, vl);
      }

      // At verbosity > VERB_LOW, each process prints something.
      if (vl > Teuchos::VERB_LOW) {
        for (int p = 0; p < numProcs; ++p) {
          if (myRank == p) {
            out << "Process " << myRank << ":" << endl;
            Teuchos::OSTab tab2 (out);
            out << "Export buffer size (in packets): "
                << exports_.dimension_0 ()
                << endl
                << "Import buffer size (in packets): "
                << imports_.dimension_0 ()
                << endl;
          }
          if (! comm.is_null ()) {
            comm->barrier (); // give output time to finish
            comm->barrier ();
            comm->barrier ();
          }
        } // for each process rank p
      } // if vl > VERB_LOW
    } // if vl != VERB_NONE
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  removeEmptyProcessesInPlace (const Teuchos::RCP<const map_type>& newMap)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
      "Tpetra::DistObject::removeEmptyProcessesInPlace: Not implemented");
  }

  /* These are provided in base DistObject template
  template<class DistObjectType>
  void
  removeEmptyProcessesInPlace (Teuchos::RCP<DistObjectType>& input,
                               const Teuchos::RCP<const Map<typename DistObjectType::local_ordinal_type,
                                                            typename DistObjectType::global_ordinal_type,
                                                            typename DistObjectType::node_type> >& newMap)
  {
    input->removeEmptyProcessesInPlace (newMap);
    if (newMap.is_null ()) { // my process is excluded
      input = Teuchos::null;
    }
  }

  template<class DistObjectType>
  void
  removeEmptyProcessesInPlace (Teuchos::RCP<DistObjectType>& input)
  {
    using Teuchos::RCP;
    typedef typename DistObjectType::local_ordinal_type LO;
    typedef typename DistObjectType::global_ordinal_type GO;
    typedef typename DistObjectType::node_type NT;
    typedef Map<LO, GO, NT> map_type;

    RCP<const map_type> newMap = input->getMap ()->removeEmptyProcesses ();
    removeEmptyProcessesInPlace<DistObjectType> (input, newMap);
  }
  */

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  doImport (const SrcDistObject& source,
            const Import<LocalOrdinal, GlobalOrdinal, Node>& importer,
            CombineMode CM)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(*getMap() != *importer.getTargetMap(),
      std::invalid_argument, "doImport: The target DistObject's Map is not "
      "identical to the Import's target Map.");
#ifdef HAVE_TPETRA_DEBUG
    {
      const this_type* srcDistObj = dynamic_cast<const this_type*> (&source);
      TEUCHOS_TEST_FOR_EXCEPTION(
        srcDistObj != NULL && * (srcDistObj->getMap ()) != *importer.getSourceMap(),
        std::invalid_argument, "doImport: The source is a DistObject, yet its "
        "Map is not identical to the Import's source Map.");
    }
#endif // HAVE_TPETRA_DEBUG
    size_t numSameIDs = importer.getNumSameIDs ();

    typedef Teuchos::ArrayView<const LocalOrdinal> view_type;
    const view_type exportLIDs      = importer.getExportLIDs();
    const view_type remoteLIDs      = importer.getRemoteLIDs();
    const view_type permuteToLIDs   = importer.getPermuteToLIDs();
    const view_type permuteFromLIDs = importer.getPermuteFromLIDs();
    this->doTransfer (source, CM, numSameIDs, permuteToLIDs, permuteFromLIDs,
                      remoteLIDs, exportLIDs, importer.getDistributor (),
                      DoForward);
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  doExport (const SrcDistObject& source,
            const Export<LocalOrdinal, GlobalOrdinal, Node>& exporter,
            CombineMode CM)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      *getMap() != *exporter.getTargetMap(), std::invalid_argument,
      "doExport: The target DistObject's Map is not identical to the Export's "
      "target Map.");
#ifdef HAVE_TPETRA_DEBUG
    {
      const this_type* srcDistObj = dynamic_cast<const this_type*> (&source);
      TEUCHOS_TEST_FOR_EXCEPTION(
        srcDistObj != NULL && * (srcDistObj->getMap ()) != *exporter.getSourceMap(),
        std::invalid_argument, "doExport: The source is a DistObject, yet its "
        "Map is not identical to the Export's source Map.");
    }
#endif // HAVE_TPETRA_DEBUG
    size_t numSameIDs = exporter.getNumSameIDs();

    typedef ArrayView<const LocalOrdinal> view_type;
    view_type exportLIDs      = exporter.getExportLIDs();
    view_type remoteLIDs      = exporter.getRemoteLIDs();
    view_type permuteToLIDs   = exporter.getPermuteToLIDs();
    view_type permuteFromLIDs = exporter.getPermuteFromLIDs();
    doTransfer (source, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs,
                exportLIDs, exporter.getDistributor (), DoForward);
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  doImport (const SrcDistObject& source,
            const Export<LocalOrdinal, GlobalOrdinal, Node> & exporter,
            CombineMode CM)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      *getMap() != *exporter.getSourceMap(), std::invalid_argument,
      "doImport (reverse mode): The target DistObject's Map is not identical "
      "to the Export's source Map.");
#ifdef HAVE_TPETRA_DEBUG
    {
      const this_type* srcDistObj = dynamic_cast<const this_type*> (&source);
      TEUCHOS_TEST_FOR_EXCEPTION(
        srcDistObj != NULL && * (srcDistObj->getMap ()) != *exporter.getTargetMap(),
        std::invalid_argument,
        "doImport (reverse mode): The source is a DistObject, yet its "
        "Map is not identical to the Export's target Map.");
    }
#endif // HAVE_TPETRA_DEBUG
    size_t numSameIDs = exporter.getNumSameIDs();

    typedef ArrayView<const LocalOrdinal> view_type;
    view_type exportLIDs      = exporter.getRemoteLIDs();
    view_type remoteLIDs      = exporter.getExportLIDs();
    view_type permuteToLIDs   = exporter.getPermuteFromLIDs();
    view_type permuteFromLIDs = exporter.getPermuteToLIDs();
    doTransfer (source, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs,
                exportLIDs, exporter.getDistributor (), DoReverse);
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  doExport (const SrcDistObject& source,
            const Import<LocalOrdinal, GlobalOrdinal, Node> & importer,
            CombineMode CM)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      *getMap() != *importer.getSourceMap(), std::invalid_argument,
      "doExport (reverse mode): The target object's Map "
      "is not identical to the Import's source Map.");
#ifdef HAVE_TPETRA_DEBUG
    {
      const this_type* srcDistObj = dynamic_cast<const this_type*> (&source);
      TEUCHOS_TEST_FOR_EXCEPTION(
        srcDistObj != NULL && * (srcDistObj->getMap ()) != *importer.getTargetMap(),
        std::invalid_argument,
        "doExport (reverse mode): The source is a DistObject, yet its "
        "Map is not identical to the Import's target Map.");
    }
#endif // HAVE_TPETRA_DEBUG
    size_t numSameIDs = importer.getNumSameIDs();

    typedef ArrayView<const LocalOrdinal> view_type;
    view_type exportLIDs      = importer.getRemoteLIDs();
    view_type remoteLIDs      = importer.getExportLIDs();
    view_type permuteToLIDs   = importer.getPermuteFromLIDs();
    view_type permuteFromLIDs = importer.getPermuteToLIDs();
    doTransfer (source, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs,
                exportLIDs, importer.getDistributor (), DoReverse);
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  bool
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  isDistributed () const {
    return map_->isDistributed ();
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  size_t
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  constantNumberOfPackets () const {
    return 0; // default implementation; subclasses may override
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  doTransfer (const SrcDistObject& src,
              CombineMode CM,
              size_t numSameIDs,
              const Teuchos::ArrayView<const LocalOrdinal>& permuteToLIDs_,
              const Teuchos::ArrayView<const LocalOrdinal>& permuteFromLIDs_,
              const Teuchos::ArrayView<const LocalOrdinal>& remoteLIDs_,
              const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs_,
              Distributor &distor,
              ReverseOption revOp)
  {
    if (this->useNewInterface()) {
      doTransferNew(src, CM, numSameIDs, permuteToLIDs_, permuteFromLIDs_, remoteLIDs_, exportLIDs_, distor, revOp);
    }
    else {
      doTransferOld(src, CM, numSameIDs, permuteToLIDs_, permuteFromLIDs_, remoteLIDs_, exportLIDs_, distor, revOp);
    }
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  reallocImportsIfNeeded (const size_t newSize, const bool debug)
  {
    if (static_cast<size_t> (imports_.dimension_0 ()) != newSize) {
      if (debug) {
        std::ostringstream os;
        os << "*** Realloc imports_ from " << imports_.dimension_0 () << " to "
           << newSize << std::endl;
        std::cerr << os.str ();
      }
      // FIXME (mfh 28 Mar 2016, 25 Apr 2016) Fences around (UVM)
      // allocations are for #227 debugging, but shouldn't be needed
      // once #227 is fixed.
      execution_space::fence ();
      imports_ = decltype (imports_) ("imports", newSize);
      execution_space::fence ();
      TEUCHOS_TEST_FOR_EXCEPTION
        (static_cast<size_t> (imports_.dimension_0 ()) != newSize,
         std::logic_error, "DualView reallocation failed: "
         "imports_.dimension_0() = " << imports_.dimension_0 ()
         << " != " << newSize << ".");
    }
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  doTransferOld (const SrcDistObject& src,
              CombineMode CM,
              size_t numSameIDs,
              const Teuchos::ArrayView<const LocalOrdinal>& permuteToLIDs,
              const Teuchos::ArrayView<const LocalOrdinal>& permuteFromLIDs,
              const Teuchos::ArrayView<const LocalOrdinal>& remoteLIDs,
              const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
              Distributor &distor,
              ReverseOption revOp)
  {
    const bool debug = false;

#ifdef HAVE_TPETRA_TRANSFER_TIMERS
    Teuchos::TimeMonitor doXferMon (*doXferTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! checkSizes (src), std::invalid_argument,
      "Tpetra::DistObject::doTransfer(): checkSizes() indicates that the "
      "destination object is not a legal target for redistribution from the "
      "source object.  This probably means that they do not have the same "
      "dimensions.  For example, MultiVectors must have the same number of "
      "rows and columns.");
    KokkosClassic::ReadWriteOption rwo = KokkosClassic::ReadWrite;
    if (CM == INSERT || CM == REPLACE) {
      const size_t numIDsToWrite = numSameIDs +
        static_cast<size_t> (permuteToLIDs.size ()) +
        static_cast<size_t> (remoteLIDs.size ());
      if (numIDsToWrite == this->getMap ()->getNodeNumElements ()) {
        // We're overwriting all of our local data in the destination
        // object, so a write-only view suffices.
        //
        // FIXME (mfh 10 Apr 2012) This doesn't make sense for a
        // CrsMatrix with a dynamic graph.  INSERT mode could mean
        // that we're adding new entries to the object, but we don't
        // want to get rid of the old ones.
        rwo = KokkosClassic::WriteOnly;
      }
    }
    // Tell the source to create a read-only view of its data.  On a
    // discrete accelerator such as a GPU, this brings EVERYTHING from
    // device memory to host memory.
    //
    // FIXME (mfh 23 Mar 2012) By passing in the list of GIDs (or
    // rather, local LIDs to send) and packet counts, createViews()
    // could create a "sparse view" that only brings in the necessary
    // data from device to host memory.
    const this_type* srcDistObj = dynamic_cast<const this_type*> (&src);
    if (srcDistObj != NULL) {
      srcDistObj->createViews ();
    }

    // Tell the target to create a view of its data.  Depending on
    // rwo, this could be a write-only view or a read-and-write view.
    // On a discrete accelerator such as a GPU, a write-only view only
    // requires a transfer from host to device memory.  A
    // read-and-write view requires a two-way transfer.  This has the
    // same problem as createViews(): it transfers EVERYTHING, not
    // just the necessary data.
    //
    // FIXME (mfh 23 Mar 2012) By passing in the list of GIDs (or
    // rather, local LIDs into which to receive) and packet counts,
    // createViewsNonConst() could create a "sparse view" that only
    // transfers the necessary data.
    this->createViewsNonConst (rwo);

    if (numSameIDs + permuteToLIDs.size()) {
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
      Teuchos::TimeMonitor copyAndPermuteMon (*copyAndPermuteTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
      // There is at least one GID to copy or permute.
      copyAndPermute (src, numSameIDs, permuteToLIDs, permuteFromLIDs);
    }

    // The method may return zero even if the implementation actually
    // does have a constant number of packets per LID.  However, if it
    // returns nonzero, we may use this information to avoid
    // (re)allocating num{Ex,Im}portPacketsPerLID_.  packAndPrepare()
    // will set this to its final value.
    //
    // We only need this if CM != ZERO, but it has to be lifted out of
    // that scope because there are multiple tests for CM != ZERO.
    size_t constantNumPackets = this->constantNumberOfPackets ();

    // We only need to pack communication buffers if the combine mode
    // is not ZERO. A "ZERO combine mode" means that the results are
    // the same as if we had received all zeros, and added them to the
    // existing values. That means we don't need to communicate.
    if (CM != ZERO) {
      if (constantNumPackets == 0) {
        // FIXME (mfh 25 Apr 2016) Fences around (UVM) allocations
        // facilitate #227 debugging, but we shouldn't need them.
        execution_space::fence ();
        numExportPacketsPerLID_ =
          decltype (numExportPacketsPerLID_) ("numExportPacketsPerLID",
                                              exportLIDs.size ());
        execution_space::fence ();
        numImportPacketsPerLID_ =
          decltype (numImportPacketsPerLID_) ("numImportPacketsPerLID",
                                              remoteLIDs.size ());
        execution_space::fence ();
      }

      {
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
        Teuchos::TimeMonitor packAndPrepareMon (*packAndPrepareTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
        // Ask the source to pack data.  Also ask it whether there are a
        // constant number of packets per element (constantNumPackets is
        // an output argument).  If there are, constantNumPackets will
        // come back nonzero.  Otherwise, the source will fill the
        // numExportPacketsPerLID_ array.
        numExportPacketsPerLID_.template modify<Kokkos::HostSpace> ();
        Teuchos::ArrayView<size_t> numExportPacketsPerLID =
          getArrayViewFromDualView (numExportPacketsPerLID_);

        // FIXME (mfh 26 Apr 2016) For backwards compatibility, use
        // the old packAndPrepare interface that takes and resizes the
        // exports buffer as a Teuchos::Array<packet_type>.  Then,
        // copy out that buffer into the host version of exports_.

        Teuchos::Array<packet_type> exportsOld;
        packAndPrepare (src, exportLIDs, exportsOld, numExportPacketsPerLID,
                        constantNumPackets, distor);
        const size_t exportsLen = static_cast<size_t> (exportsOld.size ());
        if (static_cast<size_t> (exports_.dimension_0 ()) != exportsLen) {
          // FIXME (mfh 26 Apr 2016) Fences around (UVM) allocations
          // facilitate #227 debugging, but we shouldn't need them.
          execution_space::fence ();
          exports_ = decltype (exports_) ("exports", exportsLen);
          execution_space::fence ();
        }
        Kokkos::View<const packet_type*, Kokkos::HostSpace,
          Kokkos::MemoryUnmanaged> exportsOldK (exportsOld.getRawPtr (),
                                                exportsLen);
        exports_.template modify<Kokkos::HostSpace> ();
        Kokkos::deep_copy (exports_.template view<Kokkos::HostSpace> (),
                           exportsOldK);
      }
    }

    // We don't need the source's data anymore, so it can let go of
    // its views.  On an accelerator device with a separate memory
    // space (like a GPU), this frees host memory, since device memory
    // has the "master" version of the data.
    if (srcDistObj != NULL) {
      srcDistObj->releaseViews ();
    }

    // We only need to send data if the combine mode is not ZERO.
    if (CM != ZERO) {
      if (constantNumPackets != 0) {
        // There are a constant number of packets per element.  We
        // already know (from the number of "remote" (incoming)
        // elements) how many incoming elements we expect, so we can
        // resize the buffer accordingly.
        const size_t rbufLen = remoteLIDs.size() * constantNumPackets;
        if (debug) {
          std::ostringstream os;
          os << "*** doTransferOld: Const # packets: imports_.dimension_0() = "
             << imports_.dimension_0 () << ", rbufLen = " << rbufLen
             << std::endl;
          std::cerr << os.str ();
        }
        reallocImportsIfNeeded (rbufLen, debug);
      }

      // Do we need to do communication (via doPostsAndWaits)?
      bool needCommunication = true;
      if (revOp == DoReverse && ! isDistributed ()) {
        needCommunication = false;
      }
      // FIXME (mfh 30 Jun 2013): Checking whether the source object
      // is distributed requires a cast to DistObject.  If it's not a
      // DistObject, then I'm not quite sure what to do.  Perhaps it
      // would be more appropriate for SrcDistObject to have an
      // isDistributed() method.  For now, I'll just assume that we
      // need to do communication unless the cast succeeds and the
      // source is not distributed.
      else if (revOp == DoForward && srcDistObj != NULL &&
               ! srcDistObj->isDistributed ()) {
        needCommunication = false;
      }

      if (needCommunication) {
        if (revOp == DoReverse) {
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
          Teuchos::TimeMonitor doPostsAndWaitsMon (*doPostsAndWaitsTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
          if (constantNumPackets == 0) { //variable num-packets-per-LID:
            // First communicate the number of packets per LID to receive.

            // Make sure that host has the latest version, since we're
            // using the version on host.  If host has the latest
            // version already, syncing to host does nothing.
            numExportPacketsPerLID_.template sync<Kokkos::HostSpace> ();
            Teuchos::ArrayView<const size_t> numExportPacketsPerLID =
              getArrayViewFromDualView (numExportPacketsPerLID_);

            // numImportPacketsPerLID_ is the output array here, so
            // mark it as modified.  It's strictly output, so we don't
            // have to sync from device.
            //numImportPacketsPerLID_.template sync<Kokkos::HostSpace> ();
            numImportPacketsPerLID_.template modify<Kokkos::HostSpace> ();
            Teuchos::ArrayView<size_t> numImportPacketsPerLID =
              getArrayViewFromDualView (numImportPacketsPerLID_);
            distor.doReversePostsAndWaits (numExportPacketsPerLID, 1,
                                           numImportPacketsPerLID);
            size_t totalImportPackets = 0;
            for (Array_size_type i = 0; i < numImportPacketsPerLID.size (); ++i) {
              totalImportPackets += numImportPacketsPerLID[i];
            }

            reallocImportsIfNeeded (totalImportPackets, debug);

            // We don't need to sync imports_, because it is only for
            // output here.  Similarly, we don't need to mark exports_
            // as modified, since it is read only here. This legacy
            // version of doTransfer only uses host arrays.
            imports_.template modify<Kokkos::HostSpace> ();
            Teuchos::ArrayView<packet_type> hostImports =
              getArrayViewFromDualView (imports_);
            exports_.template sync<Kokkos::HostSpace> ();
            Teuchos::ArrayView<const packet_type> hostExports =
              getArrayViewFromDualView (exports_);
            distor.doReversePostsAndWaits (hostExports,
                                           numExportPacketsPerLID,
                                           hostImports,
                                           numImportPacketsPerLID);
          }
          else {
            // We don't need to sync imports_, because it is only for
            // output here.  Similarly, we don't need to mark exports_
            // as modified, since it is read only here. This legacy
            // version of doTransfer only uses host arrays.
            imports_.template modify<Kokkos::HostSpace> ();
            Teuchos::ArrayView<packet_type> hostImports =
              getArrayViewFromDualView (imports_);
            exports_.template sync<Kokkos::HostSpace> ();
            Teuchos::ArrayView<const packet_type> hostExports =
              getArrayViewFromDualView (exports_);
            distor.doReversePostsAndWaits (hostExports,
                                           constantNumPackets,
                                           hostImports);
          }
        }
        else { // revOp == DoForward
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
          Teuchos::TimeMonitor doPostsAndWaitsMon (*doPostsAndWaitsTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
          if (constantNumPackets == 0) { //variable num-packets-per-LID:
            // First communicate the number of packets per LID to receive.

            // Make sure that host has the latest version, since we're
            // using the version on host.  If host has the latest
            // version already, syncing to host does nothing.
            numExportPacketsPerLID_.template sync<Kokkos::HostSpace> ();
            Teuchos::ArrayView<const size_t> numExportPacketsPerLID =
              getArrayViewFromDualView (numExportPacketsPerLID_);

            // numImportPacketsPerLID_ is the output array here, so
            // mark it as modified.  It's strictly output, so we don't
            // have to sync from device.
            //numImportPacketsPerLID_.template sync<Kokkos::HostSpace> ();
            numImportPacketsPerLID_.template modify<Kokkos::HostSpace> ();
            Teuchos::ArrayView<size_t> numImportPacketsPerLID =
              getArrayViewFromDualView (numImportPacketsPerLID_);
            distor.doPostsAndWaits (numExportPacketsPerLID, 1,
                                    numImportPacketsPerLID);
            size_t totalImportPackets = 0;
            for (Array_size_type i = 0; i < numImportPacketsPerLID.size (); ++i) {
              totalImportPackets += numImportPacketsPerLID[i];
            }

            reallocImportsIfNeeded (totalImportPackets, debug);

            // We don't need to sync imports_, because it is only for
            // output here.  Similarly, we don't need to mark exports_
            // as modified, since it is read only here. This legacy
            // version of doTransfer only uses host arrays.
            imports_.template modify<Kokkos::HostSpace> ();
            Teuchos::ArrayView<packet_type> hostImports =
              getArrayViewFromDualView (imports_);
            exports_.template sync<Kokkos::HostSpace> ();
            Teuchos::ArrayView<const packet_type> hostExports =
              getArrayViewFromDualView (exports_);
            distor.doPostsAndWaits (hostExports,
                                    numExportPacketsPerLID,
                                    hostImports,
                                    numImportPacketsPerLID);
          }
          else {
            // We don't need to sync imports_, because it is only for
            // output here.  Similarly, we don't need to mark exports_
            // as modified, since it is read only here. This legacy
            // version of doTransfer only uses host arrays.
            imports_.template modify<Kokkos::HostSpace> ();
            Teuchos::ArrayView<packet_type> hostImports =
              getArrayViewFromDualView (imports_);
            exports_.template sync<Kokkos::HostSpace> ();
            Teuchos::ArrayView<const packet_type> hostExports =
              getArrayViewFromDualView (exports_);
            distor.doPostsAndWaits (hostExports,
                                    constantNumPackets,
                                    hostImports);
          }
        }
        {
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
          Teuchos::TimeMonitor unpackAndCombineMon (*unpackAndCombineTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS

          // We don't need to sync imports_, because it is only for
          // output here.  This legacy version of doTransfer only uses
          // host arrays.
          imports_.template modify<Kokkos::HostSpace> ();
          Teuchos::ArrayView<packet_type> hostImports =
            getArrayViewFromDualView (imports_);
          // NOTE (mfh 25 Apr 2016) unpackAndCombine doesn't actually
          // change its numImportPacketsPerLID argument, so we don't
          // have to mark it modified here.
          numImportPacketsPerLID_.template sync<Kokkos::HostSpace> ();
          // FIXME (mfh 25 Apr 2016) unpackAndCombine doesn't actually
          // change its numImportPacketsPerLID argument, so we should
          // be able to use a const Teuchos::ArrayView here.
          Teuchos::ArrayView<size_t> numImportPacketsPerLID =
            getArrayViewFromDualView (numImportPacketsPerLID_);
          unpackAndCombine (remoteLIDs, hostImports, numImportPacketsPerLID,
                            constantNumPackets, distor, CM);
        }
      }
    } // if (CM != ZERO)

    this->releaseViews ();
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  doTransferNew (const SrcDistObject& src,
                 CombineMode CM,
                 size_t numSameIDs,
                 const Teuchos::ArrayView<const LocalOrdinal>& permuteToLIDs_,
                 const Teuchos::ArrayView<const LocalOrdinal>& permuteFromLIDs_,
                 const Teuchos::ArrayView<const LocalOrdinal>& remoteLIDs_,
                 const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs_,
                 Distributor &distor,
                 ReverseOption revOp)
  {
    using Kokkos::Compat::getArrayView;
    using Kokkos::Compat::getConstArrayView;
    using Kokkos::Compat::getKokkosViewDeepCopy;
    using Kokkos::Compat::create_const_view;
    typedef LocalOrdinal LO;
    typedef typename Kokkos::DualView<LO*, device_type>::t_dev::memory_space dev_memory_space;
    typedef typename Kokkos::DualView<LO*, device_type>::t_host::memory_space host_memory_space;

    const bool packOnHost = false;
    const bool debug = false;

    if (debug) {
      std::ostringstream os;
      os << ">>> DistObject::doTransferNew: remoteLIDs_.size() = "
         << remoteLIDs_.size () << std::endl;
      std::cerr << os.str ();
    }

#ifdef HAVE_TPETRA_TRANSFER_TIMERS
    Teuchos::TimeMonitor doXferMon (*doXferTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS

    // Convert arguments to Kokkos::View's (involves deep copy to device)
    //
    // FIXME (mfh 17 Feb 2014) getKokkosViewDeepCopy _always_ does a
    // deep copy.  It has to, since it returns a managed Kokkos::View,
    // but the input Teuchos::ArrayView is unmanaged.  One way to fix
    // this would be to change the interface by replacing the
    // Teuchos::ArrayView inputs with Kokkos::View inputs.  This is
    // the approach taken by Kokkos::DistObjectKA.  Some points:
    //
    //   1. It's perfectly OK to change the interface of doTransfer.
    //      It should take Kokkos::View by default, and convert to
    //      Teuchos::Array{RCP,View} if needed internally.
    //   2. Recall that Teuchos::ArrayView is an unmanaged view.
    //   3. If DistObject ever gets a nonblocking interface, that
    //      interface should take managed views.
    typedef Kokkos::View<const LO*, execution_space> lo_const_view_type;
    lo_const_view_type permuteToLIDs =
      getKokkosViewDeepCopy<execution_space> (permuteToLIDs_);
    lo_const_view_type permuteFromLIDs =
      getKokkosViewDeepCopy<execution_space> (permuteFromLIDs_);

    // No need to sync this.  packAndPrepareNew will use it to
    // determine where to pack (in host or device memory).
    Kokkos::DualView<LO*, device_type> remoteLIDs ("remoteLIDs",
                                                   remoteLIDs_.size ());
    {
      Kokkos::View<const LO*, host_memory_space,
        Kokkos::MemoryUnmanaged> remoteLIDs_host_in (remoteLIDs_.getRawPtr (),
                                                     remoteLIDs_.size ());
      if (packOnHost) {
        remoteLIDs.template modify<host_memory_space> ();
        Kokkos::deep_copy (remoteLIDs.template view<host_memory_space> (),
                           remoteLIDs_host_in);
      }
      else { // pack on device
        remoteLIDs.template modify<dev_memory_space> ();
        Kokkos::deep_copy (remoteLIDs.template view<dev_memory_space> (),
                           remoteLIDs_host_in);
      }
    }

    if (debug) {
      TEUCHOS_TEST_FOR_EXCEPTION
        (static_cast<size_t> (permuteFromLIDs.size ()) !=
         static_cast<size_t> (permuteFromLIDs_.size ()), std::logic_error,
         "permuteFromLIDs.size() = " << permuteFromLIDs.size () <<
         " != permuteFromLIDs_.size() = " << permuteFromLIDs_.size () << ".");
    }

    Kokkos::DualView<LO*, device_type> exportLIDs ("exportLIDs",
                                                   exportLIDs_.size ());
    {
      Kokkos::View<const LO*, host_memory_space,
        Kokkos::MemoryUnmanaged> exportLIDs_host_in (exportLIDs_.getRawPtr (),
                                                     exportLIDs_.size ());

      // mfh 12 Apr 2016: packAndPrepareNew decides where to pack
      // based on the memory space in which exportLIDs was last
      // modified.
      if (packOnHost) {
        exportLIDs.template modify<host_memory_space> ();
        Kokkos::deep_copy (exportLIDs.template view<host_memory_space> (),
                           exportLIDs_host_in);
      }
      else {
        exportLIDs.template modify<dev_memory_space> ();
        Kokkos::deep_copy (exportLIDs.template view<dev_memory_space> (),
                           exportLIDs_host_in);
      }
    }

    if (debug) {
      std::cerr << ">>> 1. checkSizes" << std::endl;
    }

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! checkSizes (src), std::invalid_argument,
      "Tpetra::DistObject::doTransfer(): checkSizes() indicates that the "
      "destination object is not a legal target for redistribution from the "
      "source object.  This probably means that they do not have the same "
      "dimensions.  For example, MultiVectors must have the same number of "
      "rows and columns.");
    KokkosClassic::ReadWriteOption rwo = KokkosClassic::ReadWrite;
    if (CM == INSERT || CM == REPLACE) {
      const size_t numIDsToWrite = numSameIDs +
        static_cast<size_t> (permuteToLIDs.size ()) +
        static_cast<size_t> (remoteLIDs.dimension_0 ());
      if (numIDsToWrite == this->getMap ()->getNodeNumElements ()) {
        // We're overwriting all of our local data in the destination
        // object, so a write-only view suffices.
        //
        // FIXME (mfh 10 Apr 2012) This doesn't make sense for a
        // CrsMatrix with a dynamic graph.  INSERT mode could mean
        // that we're adding new entries to the object, but we don't
        // want to get rid of the old ones.
        rwo = KokkosClassic::WriteOnly;
      }
    }

    if (debug) {
      std::cerr << ">>> 2. createViews" << std::endl;
    }

    // FIXME (mfh 17 Feb 2014) We're assuming that MPI can read device
    // memory (that's even pre-UVM), so there is no need to create
    // host views of the source object's data.

    // Tell the source to create a read-only view of its data.  On a
    // discrete accelerator such as a GPU, this brings EVERYTHING from
    // device memory to host memory.
    //
    // FIXME (mfh 23 Mar 2012) By passing in the list of GIDs (or
    // rather, local LIDs to send) and packet counts, createViews()
    // could create a "sparse view" that only brings in the necessary
    // data from device to host memory.
    const this_type* srcDistObj = dynamic_cast<const this_type*> (&src);
    if (srcDistObj != NULL) {
      srcDistObj->createViews ();
    }

    // FIXME (mfh 17 Feb 2014) We're assuming that MPI can read device
    // memory (that's even pre-UVM), so there is no need to create
    // host views of the target object's data.

    // Tell the target to create a view of its data.  Depending on
    // rwo, this could be a write-only view or a read-and-write view.
    // On a discrete accelerator such as a GPU, a write-only view only
    // requires a transfer from host to device memory.  A
    // read-and-write view requires a two-way transfer.  This has the
    // same problem as createViews(): it transfers EVERYTHING, not
    // just the necessary data.
    //
    // FIXME (mfh 23 Mar 2012) By passing in the list of GIDs (or
    // rather, local LIDs into which to receive) and packet counts,
    // createViewsNonConst() could create a "sparse view" that only
    // transfers the necessary data.
    this->createViewsNonConst (rwo);

    if (debug) {
      std::cerr << ">>> 3. copyAndPermuteNew" << std::endl;
    }

    if (numSameIDs + permuteToLIDs.size()) {
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
      Teuchos::TimeMonitor copyAndPermuteMon (*copyAndPermuteTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS

      // FIXME (mfh 17 Feb 2014) Nobody implements DistObject
      // subclasses but Tpetra developers anyway, so don't bother with
      // this "new" business.  Just write the interface you want.

      // There is at least one GID to copy or permute.
      copyAndPermuteNew (src, numSameIDs, permuteToLIDs, permuteFromLIDs);
    }

    // The method may return zero even if the implementation actually
    // does have a constant number of packets per LID.  However, if it
    // returns nonzero, we may use this information to avoid
    // (re)allocating num{Ex,Im}portPacketsPerLID_.  packAndPrepare()
    // will set this to its final value.
    //
    // We only need this if CM != ZERO, but it has to be lifted out of
    // that scope because there are multiple tests for CM != ZERO.
    size_t constantNumPackets = this->constantNumberOfPackets ();

    // We only need to pack communication buffers if the combine mode
    // is not ZERO. A "ZERO combine mode" means that the results are
    // the same as if we had received all zeros, and added them to the
    // existing values. That means we don't need to communicate.
    if (CM != ZERO) {
      if (constantNumPackets == 0) {
        if (debug) {
          std::cerr << ">>> 4. Allocate num{Ex,Im}portPacketsPerLID" << std::endl;
        }

        // Don't "realloc" unless you really need to.  That will avoid
        // a bit of time for reinitializing the Views' data.

        //Kokkos::realloc (numExportPacketsPerLID_, exportLIDs.dimension_0 ());
        if (numExportPacketsPerLID_.dimension_0 () != exportLIDs.dimension_0 ()) {
          execution_space::fence ();
          numExportPacketsPerLID_ =
            decltype (numExportPacketsPerLID_) ("numExportPacketsPerLID",
                                                exportLIDs.dimension_0 ());
          execution_space::fence ();
        }
        if (numImportPacketsPerLID_.dimension_0 () != remoteLIDs.dimension_0 ()) {
          // FIXME (mfh 25 Apr 2016) Fences around (UVM) allocations
          // facilitate #227 debugging, but shouldn't be needed.
          execution_space::fence ();
          numImportPacketsPerLID_ =
            decltype (numImportPacketsPerLID_) ("numImportPacketsPerLID",
                                                remoteLIDs.dimension_0 ());
          execution_space::fence ();
        }
      }

      if (debug) {
        std::cerr << ">>> 5. packAndPrepareNew" << std::endl;
      }

      {
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
        Teuchos::TimeMonitor packAndPrepareMon (*packAndPrepareTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS

        if (debug) {
          std::ostringstream os;
          const int myRank = this->getMap ()->getComm ()->getRank ();
          os << ">>> (Proc " << myRank << ") 5.0. Before packAndPrepareNew, "
            "exports_.dimension_0()=" << exports_.dimension_0 () << std::endl;
          std::cerr << os.str ();
        }
        // Ask the source to pack data.  Also ask it whether there are a
        // constant number of packets per element (constantNumPackets is
        // an output argument).  If there are, constantNumPackets will
        // come back nonzero.  Otherwise, the source will fill the
        // numExportPacketsPerLID_ array.

        packAndPrepareNew (src, exportLIDs, exports_, numExportPacketsPerLID_,
                           constantNumPackets, distor);
        if (debug) {
          std::ostringstream os;
          const int myRank = this->getMap ()->getComm ()->getRank ();
          os << ">>> (Proc " << myRank << ") 5.0. After packAndPrepareNew, "
            "exports_.dimension_0()=" << exports_.dimension_0 () << std::endl;
          std::cerr << os.str ();
        }
      }
    }

    if (debug) {
      std::cerr << ">>> 6. releaseViews" << std::endl;
    }

    // FIXME (mfh 17 Feb 2014) See comments above; there is no need to
    // create host views of the source object's data.

    // We don't need the source's data anymore, so it can let go of
    // its views.  On an accelerator device with a separate memory
    // space (like a GPU), this frees host memory, since device memory
    // has the "master" version of the data.
    if (srcDistObj != NULL) {
      srcDistObj->releaseViews ();
    }

    // We only need to send data if the combine mode is not ZERO.
    if (CM != ZERO) {
      if (constantNumPackets != 0) {
        if (debug) {
          std::cerr << ">>> 7. Realloc imports_" << std::endl;
        }
        // There are a constant number of packets per element.  We
        // already know (from the number of "remote" (incoming)
        // elements) how many incoming elements we expect, so we can
        // resize the buffer accordingly.
        const size_t rbufLen = remoteLIDs.dimension_0 () * constantNumPackets;
        if (debug) {
          std::ostringstream os;
          os << "*** doTransferNew: imports_.dimension_0() = "
             << imports_.dimension_0 () << ", rbufLen = " << rbufLen
             << " = " << remoteLIDs.dimension_0 ()
             << " * " << constantNumPackets << std::endl;
          std::cerr << os.str ();
        }
        reallocImportsIfNeeded (rbufLen, debug);
      }

      // FIXME (mfh 28 Mar 2016) Could this possibly help fix #227 ???
      execution_space::fence ();

      if (debug) {
        std::cerr << ">>> 8. Copy numExportPacketsPerLID to host" << std::endl;
      }

      // Do we need to do communication (via doPostsAndWaits)?
      bool needCommunication = true;
      if (revOp == DoReverse && ! isDistributed ()) {
        needCommunication = false;
      }
      // FIXME (mfh 30 Jun 2013): Checking whether the source object
      // is distributed requires a cast to DistObject.  If it's not a
      // DistObject, then I'm not quite sure what to do.  Perhaps it
      // would be more appropriate for SrcDistObject to have an
      // isDistributed() method.  For now, I'll just assume that we
      // need to do communication unless the cast succeeds and the
      // source is not distributed.
      else if (revOp == DoForward && srcDistObj != NULL &&
               ! srcDistObj->isDistributed ()) {
        needCommunication = false;
      }

      // FIXME (mfh 17 Feb 2014) Distributor doesn't actually inspect
      // the contents of the "exports" or "imports" arrays, other than
      // to do a deep copy in the (should be technically unnecessary,
      // but isn't for some odd reason) "self-message" case.
      // Distributor thus doesn't need host views; it could do just
      // fine with device views, assuming that MPI knows how to read
      // device memory (which doesn't even require UVM).

      if (needCommunication) {
        if (debug) {
          std::cerr << ">>> 9. Communicate" << std::endl;
        }

        if (revOp == DoReverse) {
          if (debug) {
            std::cerr << ">>> 9.0. Reverse mode" << std::endl;
          }

#ifdef HAVE_TPETRA_TRANSFER_TIMERS
          Teuchos::TimeMonitor doPostsAndWaitsMon (*doPostsAndWaitsTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
          if (constantNumPackets == 0) { //variable num-packets-per-LID:
            if (debug) {
              std::cerr << ">>> 9.1. Variable # packets / LID: first comm" << std::endl;
            }
            numExportPacketsPerLID_.template sync<Kokkos::HostSpace> ();
            numImportPacketsPerLID_.template sync<Kokkos::HostSpace> ();
            distor.doReversePostsAndWaits (create_const_view (numExportPacketsPerLID_.template view<Kokkos::HostSpace> ()),
                                           1,
                                           numImportPacketsPerLID_.template view<Kokkos::HostSpace> ());
            size_t totalImportPackets = 0;
            // FIXME (mfh 17 Feb 2014) This would be a good place for
            // a Kokkos reduction.  numImportPacketsPerLID_ has as
            // many entries as the number of LIDs on the calling
            // process.
            {
              typedef decltype (numImportPacketsPerLID_) dual_view_type;
              typedef typename dual_view_type::t_host host_view_type;
              typedef typename host_view_type::const_type const_host_view_type;

              const_host_view_type host_numImportPacketsPerLID =
                numImportPacketsPerLID_.template view<Kokkos::HostSpace> ();
              const view_size_type numLids = host_numImportPacketsPerLID.size ();
              for (view_size_type i = 0; i < numLids; ++i) {
                totalImportPackets += host_numImportPacketsPerLID[i];
              }
            }

            if (debug) {
              std::cerr << ">>> 9.2. Realloc" << std::endl;
            }

            reallocImportsIfNeeded (totalImportPackets, debug);

            if (debug) {
              std::cerr << ">>> 9.3. Second comm" << std::endl;
            }

            if (packOnHost) {
              numExportPacketsPerLID_.template sync<host_memory_space> ();
              numImportPacketsPerLID_.template sync<host_memory_space> ();
              // imports_ is for output only, so we don't need to sync it.
              imports_.template modify<host_memory_space> ();
              distor.doReversePostsAndWaits (create_const_view (exports_.template view<host_memory_space> ()),
                                             getArrayViewFromDualView (numExportPacketsPerLID_),
                                             imports_.template view<host_memory_space> (),
                                             getArrayViewFromDualView (numImportPacketsPerLID_));
            }
            else {
              // FIXME (mfh 25 Apr 2016) Once doReversePostsAndWaits
              // can take numExportPacketsPerLID and
              // numImportPacketsPerLID as View or DualView, rather
              // than as Teuchos::ArrayView, then we can use their
              // device versions.  For now, we'll use their host
              // versions.
              numExportPacketsPerLID_.template sync<Kokkos::HostSpace> ();
              numImportPacketsPerLID_.template sync<Kokkos::HostSpace> ();
              // imports_ is for output only, so we don't need to sync it.
              imports_.template modify<dev_memory_space> ();
              distor.doReversePostsAndWaits (create_const_view (exports_.template view<dev_memory_space> ()),
                                             getArrayViewFromDualView (numExportPacketsPerLID_),
                                             imports_.template view<dev_memory_space> (),
                                             getArrayViewFromDualView (numImportPacketsPerLID_));
            }
          }
          else {
            if (debug) {
              const int myRank = this->getMap ()->getComm ()->getRank ();
              std::ostringstream os;
              os << ">>> (Proc " << myRank << "): 9.1. Const # packets per "
                "LID: exports_.dimension_0() = " << exports_.dimension_0 ()
                 << ", imports_.dimension_0() = " << imports_.dimension_0 ()
                 << std::endl;
              std::cerr << os.str ();
            }
            if (packOnHost) {
              // imports_ is for output only, so we don't need to sync it.
              imports_.template modify<host_memory_space> ();
              distor.doReversePostsAndWaits (create_const_view (exports_.template view<host_memory_space> ()),
                                             constantNumPackets,
                                             imports_.template view<host_memory_space> ());
            }
            else { // pack on device
              // imports_ is for output only, so we don't need to sync it.
              imports_.template modify<dev_memory_space> ();
              distor.doReversePostsAndWaits (create_const_view (exports_.template view<dev_memory_space> ()),
                                             constantNumPackets,
                                             imports_.template view<dev_memory_space> ());
            }
          }
        }
        else { // revOp == DoForward
          if (debug) {
            std::cerr << ">>> 9.0. Forward mode" << std::endl;
          }

#ifdef HAVE_TPETRA_TRANSFER_TIMERS
          Teuchos::TimeMonitor doPostsAndWaitsMon (*doPostsAndWaitsTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
          if (constantNumPackets == 0) { //variable num-packets-per-LID:
            if (debug) {
              std::cerr << ">>> 9.1. Variable # packets / LID: first comm" << std::endl;
            }

            numExportPacketsPerLID_.template sync<Kokkos::HostSpace> ();
            numImportPacketsPerLID_.template sync<Kokkos::HostSpace> ();
            distor.doPostsAndWaits (create_const_view (numExportPacketsPerLID_.template view<Kokkos::HostSpace> ()), 1,
                                    numImportPacketsPerLID_.template view<Kokkos::HostSpace> ());
            size_t totalImportPackets = 0;
            {
              typedef decltype (numImportPacketsPerLID_) dual_view_type;
              typedef typename dual_view_type::t_host host_view_type;
              typedef typename host_view_type::const_type const_host_view_type;
              const_host_view_type host_numImportPacketsPerLID =
                numImportPacketsPerLID_.template view<Kokkos::HostSpace> ();

              // FIXME (mfh 17 Feb 2014) This would be a good place for
              // a Kokkos reduction.  numImportPacketsPerLID_ has as
              // many entries as the number of LIDs on the calling
              // process.
              const view_size_type numLids = host_numImportPacketsPerLID.size ();
              for (view_size_type i = 0; i < numLids; ++i) {
                totalImportPackets += host_numImportPacketsPerLID[i];
              }
            }

            if (debug) {
              std::cerr << ">>> 9.2. Realloc" << std::endl;
            }

            reallocImportsIfNeeded (totalImportPackets, debug);

            if (debug) {
              std::cerr << ">>> 9.3. Second comm" << std::endl;
            }

            if (packOnHost) {
              numExportPacketsPerLID_.template sync<host_memory_space> ();
              numImportPacketsPerLID_.template sync<host_memory_space> ();
              // imports_ is for output only, so we don't need to sync it.
              imports_.template modify<host_memory_space> ();
              distor.doPostsAndWaits (create_const_view (exports_.template view<host_memory_space> ()),
                                      getArrayViewFromDualView (numExportPacketsPerLID_),
                                      imports_.template view<host_memory_space> (),
                                      getArrayViewFromDualView (numImportPacketsPerLID_));
            }
            else { // pack on device
              // FIXME (mfh 25 Apr 2016) Once doReversePostsAndWaits
              // can take numExportPacketsPerLID and
              // numImportPacketsPerLID as a View or DualView, rather
              // than as a Teuchos::ArrayView, then we can use their
              // device version.  For now, we'll use the host version.
              numExportPacketsPerLID_.template sync<host_memory_space> ();
              numImportPacketsPerLID_.template sync<host_memory_space> ();
              // imports_ is for output only, so we don't need to sync it.
              imports_.template modify<dev_memory_space> ();
              distor.doPostsAndWaits (create_const_view (exports_.template view<dev_memory_space> ()),
                                      getArrayViewFromDualView (numExportPacketsPerLID_),
                                      imports_.template view<dev_memory_space> (),
                                      getArrayViewFromDualView (numImportPacketsPerLID_));
            }
          }
          else {
            if (debug) {
              const int myRank = this->getMap ()->getComm ()->getRank ();
              std::ostringstream os;
              os << ">>> (Proc " << myRank << "): 9.1. Const # packets per "
                "LID: exports_.dimension_0()=" << exports_.dimension_0 ()
                 << ", imports_.dimension_0() = " << imports_.dimension_0 ()
                 << std::endl;
              std::cerr << os.str ();
            }

            if (packOnHost) {
              // imports_ is for output only, so we don't need to sync it.
              imports_.template modify<host_memory_space> ();
              distor.doPostsAndWaits (create_const_view (exports_.template view<host_memory_space> ()),
                                      constantNumPackets,
                                      imports_.template view<host_memory_space> ());
            }
            else { // pack on device
              // imports_ is for output only, so we don't need to sync it.
              imports_.template modify<dev_memory_space> ();
              distor.doPostsAndWaits (create_const_view (exports_.template view<dev_memory_space> ()),
                                      constantNumPackets,
                                      imports_.template view<dev_memory_space> ());
            }
          }
        }

        if (debug) {
          std::cerr << ">>> 10. Copy numImportPacketsPerLID to device" << std::endl;
        }

        if (debug) {
          std::cerr << ">>> 11. unpackAndCombineNew" << std::endl;
        }

        {
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
          Teuchos::TimeMonitor unpackAndCombineMon (*unpackAndCombineTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS

          // NOTE (mfh 26 Apr 2016) We don't actually need to sync the
          // input DualViews, but they DO need to be most recently
          // updated in the same memory space.
          //
          // FIXME (mfh 26 Apr 2016) Check that all input DualViews
          // were most recently updated in the same memory space, and
          // sync them to the same place (based on packOnHost) if not.
          unpackAndCombineNew (remoteLIDs, imports_, numImportPacketsPerLID_,
                               constantNumPackets, distor, CM);
        }
      }
    } // if (CM != ZERO)

    if (debug) {
      std::cerr << ">>> 12. releaseViews" << std::endl;
    }

    // FIXME (mfh 17 Dec(??? probably Feb) 2014) We don't have to call
    // releaseViews() on the destination object any more, since MPI
    // knows how to read device memory.

    this->releaseViews ();

    if (debug) {
      std::cerr << ">>> 13. Done with doTransferNew" << std::endl;
    }
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  print (std::ostream &os) const
  {
    using Teuchos::FancyOStream;
    using Teuchos::getFancyOStream;
    using Teuchos::RCP;
    using Teuchos::rcpFromRef;
    using std::endl;

    RCP<FancyOStream> out = getFancyOStream (rcpFromRef (os));
    this->describe (*out, Teuchos::VERB_DEFAULT);
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  createViews () const
  {}

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  createViewsNonConst (KokkosClassic::ReadWriteOption /*rwo*/)
  {}

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node, const bool classic>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node, classic>::
  releaseViews () const
  {}

  template<class DistObjectType>
  void
  removeEmptyProcessesInPlace (Teuchos::RCP<DistObjectType>& input,
                               const Teuchos::RCP<const Map<typename DistObjectType::local_ordinal_type,
                                                            typename DistObjectType::global_ordinal_type,
                                                            typename DistObjectType::node_type> >& newMap)
  {
    input->removeEmptyProcessesInPlace (newMap);
    if (newMap.is_null ()) { // my process is excluded
      input = Teuchos::null;
    }
  }

  template<class DistObjectType>
  void
  removeEmptyProcessesInPlace (Teuchos::RCP<DistObjectType>& input)
  {
    using Teuchos::RCP;
    typedef typename DistObjectType::local_ordinal_type LO;
    typedef typename DistObjectType::global_ordinal_type GO;
    typedef typename DistObjectType::node_type NT;
    typedef Map<LO, GO, NT> map_type;

    RCP<const map_type> newMap = input->getMap ()->removeEmptyProcesses ();
    removeEmptyProcessesInPlace<DistObjectType> (input, newMap);
  }

// Explicit instantiation macro for general DistObject.
#define TPETRA_DISTOBJECT_INSTANT(SCALAR, LO, GO, NODE) \
  template class DistObject< SCALAR , LO , GO , NODE >;

// Explicit instantiation macro for DistObject<char, ...>.
// The "SLGN" stuff above doesn't work for Packet=char.
#define TPETRA_DISTOBJECT_INSTANT_CHAR(LO, GO, NODE) \
  template class DistObject< char , LO , GO , NODE >;

} // namespace Tpetra

#endif // TPETRA_DISTOBJECT_DEF_HPP
