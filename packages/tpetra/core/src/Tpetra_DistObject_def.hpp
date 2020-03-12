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
#include "Tpetra_Details_reallocDualViewIfNeeded.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_checkGlobalError.hpp"
#include "Tpetra_Details_Profiling.hpp"
#include "Tpetra_Util.hpp" // Details::createPrefix
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include <typeinfo>
#include <memory>
#include <sstream>

namespace Tpetra {

  namespace { // (anonymous)
    template<class DeviceType, class IndexType = size_t>
    struct SumFunctor {
      SumFunctor (const Kokkos::View<const size_t*, DeviceType>& viewToSum) :
        viewToSum_ (viewToSum) {}
      KOKKOS_INLINE_FUNCTION void operator() (const IndexType i, size_t& lclSum) const {
        lclSum += viewToSum_(i);
      }
      Kokkos::View<const size_t*, DeviceType> viewToSum_;
    };

    template<class DeviceType, class IndexType = size_t>
    size_t
    countTotalImportPackets (const Kokkos::View<const size_t*, DeviceType>& numImportPacketsPerLID)
    {
      using Kokkos::parallel_reduce;
      typedef DeviceType DT;
      typedef typename DT::execution_space DES;
      typedef Kokkos::RangePolicy<DES, IndexType> range_type;

      const IndexType numOut = numImportPacketsPerLID.extent (0);
      size_t totalImportPackets = 0;
      parallel_reduce ("Count import packets",
                       range_type (0, numOut),
                       SumFunctor<DeviceType, IndexType> (numImportPacketsPerLID),
                       totalImportPackets);
      return totalImportPackets;
    }
  } // namespace (anonymous)


  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
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

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
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

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
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
                << exports_.extent (0)
                << endl
                << "Import buffer size (in packets): "
                << imports_.extent (0)
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

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  removeEmptyProcessesInPlace (const Teuchos::RCP<const map_type>& /* newMap */)
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

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  doImport (const SrcDistObject& source,
            const Import<LocalOrdinal, GlobalOrdinal, Node>& importer,
            const CombineMode CM,
            const bool restrictedMode)
  {
    using Details::Behavior;
    using std::endl;
    const char modeString[] = "doImport (forward mode)";

    // mfh 18 Oct 2017: Set TPETRA_VERBOSE to true for copious debug
    // output to std::cerr on every MPI process.  This is unwise for
    // runs with large numbers of MPI processes.
    const bool verbose = Behavior::verbose("DistObject");
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("DistObject", modeString);
      std::ostringstream os;
      os << *prefix << "Start" << endl;
      std::cerr << os.str ();
    }
    this->doTransfer (source, importer, modeString, DoForward, CM,
                      restrictedMode);
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str ();
    }
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  doExport (const SrcDistObject& source,
            const Export<LocalOrdinal, GlobalOrdinal, Node>& exporter,
            const CombineMode CM,
            const bool restrictedMode)
  {
    using Details::Behavior;
    using std::endl;
    const char modeString[] = "doExport (forward mode)";

    // mfh 18 Oct 2017: Set TPETRA_VERBOSE to true for copious debug
    // output to std::cerr on every MPI process.  This is unwise for
    // runs with large numbers of MPI processes.
    const bool verbose = Behavior::verbose("DistObject");
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("DistObject", modeString);
      std::ostringstream os;
      os << *prefix << "Start" << endl;
      std::cerr << os.str ();
    }
    this->doTransfer (source, exporter, modeString, DoForward, CM,
                      restrictedMode);
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str ();
    }
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  doImport (const SrcDistObject& source,
            const Export<LocalOrdinal, GlobalOrdinal, Node>& exporter,
            const CombineMode CM,
            const bool restrictedMode)
  {
    using Details::Behavior;
    using std::endl;
    const char modeString[] = "doImport (reverse mode)";

    // mfh 18 Oct 2017: Set TPETRA_VERBOSE to true for copious debug
    // output to std::cerr on every MPI process.  This is unwise for
    // runs with large numbers of MPI processes.
    const bool verbose = Behavior::verbose("DistObject");
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("DistObject", modeString);
      std::ostringstream os;
      os << *prefix << "Start" << endl;
      std::cerr << os.str ();
    }
    this->doTransfer (source, exporter, modeString, DoReverse, CM,
                      restrictedMode);
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str ();
    }
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  doExport (const SrcDistObject& source,
            const Import<LocalOrdinal, GlobalOrdinal, Node> & importer,
            const CombineMode CM,
            const bool restrictedMode)
  {
    using Details::Behavior;
    using std::endl;
    const char modeString[] = "doExport (reverse mode)";

    // mfh 18 Oct 2017: Set TPETRA_VERBOSE to true for copious debug
    // output to std::cerr on every MPI process.  This is unwise for
    // runs with large numbers of MPI processes.
    const bool verbose = Behavior::verbose("DistObject");
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("DistObject", modeString);
      std::ostringstream os;
      os << *prefix << "Start" << endl;
      std::cerr << os.str ();
    }
    this->doTransfer (source, importer, modeString, DoReverse, CM,
                      restrictedMode);
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Done" << endl;
      std::cerr << os.str ();
    }
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  isDistributed () const {
    return map_->isDistributed ();
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  constantNumberOfPackets () const {
    return 0; // default implementation; subclasses may override
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  doTransfer (const SrcDistObject& src,
              const ::Tpetra::Details::Transfer<local_ordinal_type, global_ordinal_type, node_type>& transfer,
              const char modeString[],
              const ReverseOption revOp,
              const CombineMode CM,
              bool restrictedMode)
  {
    using Details::Behavior;
    using Details::getDualViewCopyFromArrayView;
    using Details::ProfilingRegion;
    using std::endl;
    const char funcName[] = "Tpetra::DistObject::doTransfer";

    ProfilingRegion region_doTransfer(funcName);
    const bool verbose = Behavior::verbose("DistObject");
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      std::ostringstream os;
      prefix = this->createPrefix("DistObject", "doTransfer");
      os << *prefix << "Source type: " << Teuchos::typeName(src)
         << ", Target type: " << Teuchos::typeName(*this) << endl;
      std::cerr << os.str();
    }

    // "Restricted Mode" does two things:
    // 1) Skips copyAndPermute
    // 2) Allows the "target" Map of the transfer to be a subset of
    //    the Map of *this, in a "locallyFitted" sense.
    //
    // This cannot be used if #2 is not true, OR there are permutes.
    // Source Maps still need to match

    // mfh 18 Oct 2017: Set TPETRA_DEBUG to true to enable extra debug
    // checks.  These may communicate more.
    const bool debug = Behavior::debug("DistObject");
    if (debug) {
      if (! restrictedMode && revOp == DoForward) {
        const bool myMapSameAsTransferTgtMap =
          this->getMap ()->isSameAs (* (transfer.getTargetMap ()));
        TEUCHOS_TEST_FOR_EXCEPTION
          (! myMapSameAsTransferTgtMap, std::invalid_argument,
           "Tpetra::DistObject::" << modeString << ": For forward-mode "
           "communication, the target DistObject's Map must be the same "
           "(in the sense of Tpetra::Map::isSameAs) as the input "
           "Export/Import object's target Map.");
      }
      else if (! restrictedMode && revOp == DoReverse) {
        const bool myMapSameAsTransferSrcMap =
          this->getMap ()->isSameAs (* (transfer.getSourceMap ()));
        TEUCHOS_TEST_FOR_EXCEPTION
          (! myMapSameAsTransferSrcMap, std::invalid_argument,
           "Tpetra::DistObject::" << modeString << ": For reverse-mode "
           "communication, the target DistObject's Map must be the same "
         "(in the sense of Tpetra::Map::isSameAs) as the input "
           "Export/Import object's source Map.");
      }
      else if (restrictedMode && revOp == DoForward) {
        const bool myMapLocallyFittedTransferTgtMap =
          this->getMap ()->isLocallyFitted (* (transfer.getTargetMap ()));
        TEUCHOS_TEST_FOR_EXCEPTION
          (! myMapLocallyFittedTransferTgtMap , std::invalid_argument,
           "Tpetra::DistObject::" << modeString << ": For forward-mode "
           "communication using restricted mode, Export/Import object's "
           "target Map must be locally fitted (in the sense of "
           "Tpetra::Map::isLocallyFitted) to target DistObject's Map.");
      }
      else { // if (restrictedMode && revOp == DoReverse) {
        const bool myMapLocallyFittedTransferSrcMap =
          this->getMap ()->isLocallyFitted (* (transfer.getSourceMap ()));
        TEUCHOS_TEST_FOR_EXCEPTION
          (! myMapLocallyFittedTransferSrcMap, std::invalid_argument,
           "Tpetra::DistObject::" << modeString << ": For reverse-mode "
           "communication using restricted mode, Export/Import object's "
           "source Map must be locally fitted (in the sense of "
           "Tpetra::Map::isLocallyFitted) to target DistObject's Map.");
      }

      // SrcDistObject need not even _have_ Maps.  However, if the
      // source object is a DistObject, it has a Map, and we may
      // compare that Map with the Transfer's Maps.
      const this_type* srcDistObj = dynamic_cast<const this_type*> (&src);
      if (srcDistObj != nullptr) {
        if (revOp == DoForward) {
          const bool srcMapSameAsImportSrcMap =
            srcDistObj->getMap ()->isSameAs (* (transfer.getSourceMap ()));
          TEUCHOS_TEST_FOR_EXCEPTION
            (! srcMapSameAsImportSrcMap, std::invalid_argument,
             "Tpetra::DistObject::" << modeString << ": For forward-mode "
             "communication, the source DistObject's Map must be the same "
             "as the input Export/Import object's source Map.");
        }
        else { // revOp == DoReverse
          const bool srcMapSameAsImportTgtMap =
            srcDistObj->getMap ()->isSameAs (* (transfer.getTargetMap ()));
          TEUCHOS_TEST_FOR_EXCEPTION
            (! srcMapSameAsImportTgtMap, std::invalid_argument,
             "Tpetra::DistObject::" << modeString << ": For reverse-mode "
             "communication, the source DistObject's Map must be the same "
             "as the input Export/Import object's target Map.");
        }
      }
    }

    const size_t numSameIDs = transfer.getNumSameIDs ();
    Distributor& distor = transfer.getDistributor ();

    TEUCHOS_TEST_FOR_EXCEPTION
      (debug && restrictedMode &&
       (transfer.getPermuteToLIDs_dv().extent(0) != 0 ||
        transfer.getPermuteFromLIDs_dv().extent(0) != 0),
       std::invalid_argument,
       "Tpetra::DistObject::" << modeString << ": Transfer object "
       "cannot have permutes in restricted mode.");

    // Do we need all communication buffers to live on host?
    const bool commOnHost = ! Behavior::assumeMpiIsCudaAware ();
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "doTransfer: Use new interface; "
        "commOnHost=" << (commOnHost ? "true" : "false") << endl;
      std::cerr << os.str ();
    }

    using const_lo_dv_type =
      Kokkos::DualView<const local_ordinal_type*, buffer_device_type>;
    const_lo_dv_type permToLIDs = (revOp == DoForward) ?
      transfer.getPermuteToLIDs_dv () :
      transfer.getPermuteFromLIDs_dv ();
    const_lo_dv_type permFromLIDs = (revOp == DoForward) ?
      transfer.getPermuteFromLIDs_dv () :
      transfer.getPermuteToLIDs_dv ();
    const_lo_dv_type remoteLIDs = (revOp == DoForward) ?
      transfer.getRemoteLIDs_dv () :
      transfer.getExportLIDs_dv ();
    const_lo_dv_type exportLIDs = (revOp == DoForward) ?
      transfer.getExportLIDs_dv () :
      transfer.getRemoteLIDs_dv ();
    doTransferNew (src, CM, numSameIDs, permToLIDs, permFromLIDs,
                   remoteLIDs, exportLIDs, distor, revOp, commOnHost,
                   restrictedMode);

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Tpetra::DistObject::doTransfer: Done!" << endl;
      std::cerr << os.str ();
    }
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  reallocImportsIfNeeded (const size_t newSize,
                          const bool verbose,
                          const std::string* prefix)
  {
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Realloc (if needed) imports_ from "
         << imports_.extent (0) << " to " << newSize << std::endl;
      std::cerr << os.str ();
    }
    using ::Tpetra::Details::reallocDualViewIfNeeded;
    const bool reallocated =
      reallocDualViewIfNeeded (this->imports_, newSize, "imports");
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Finished realloc'ing imports_" << std::endl;
      std::cerr << os.str ();
    }
    return reallocated;
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  reallocArraysForNumPacketsPerLid (const size_t numExportLIDs,
                                    const size_t numImportLIDs)
  {
    using Details::Behavior;
    using ::Tpetra::Details::dualViewStatusToString;
    using ::Tpetra::Details::reallocDualViewIfNeeded;
    using std::endl;
    // If an array is already allocated, and if is at least
    // tooBigFactor times bigger than it needs to be, free it and
    // reallocate to the size we need, in order to save space.
    // Otherwise, take subviews to reduce allocation size.
    constexpr size_t tooBigFactor = 10;

    const bool verbose = Behavior::verbose("DistObject");
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("DistObject",
        "reallocArraysForNumPacketsPerLid");
      std::ostringstream os;
      os << *prefix
         << "numExportLIDs: " << numExportLIDs
         << ", numImportLIDs: " << numImportLIDs
         << endl;
      os << *prefix << "DualView status before:" << endl
         << *prefix
         << dualViewStatusToString (this->numExportPacketsPerLID_,
                                    "numExportPacketsPerLID_")
         << endl
         << *prefix
         << dualViewStatusToString (this->numImportPacketsPerLID_,
                                    "numImportPacketsPerLID_")
         << endl;
      std::cerr << os.str ();
    }

    // Reallocate numExportPacketsPerLID_ if needed.
    const bool firstReallocated =
      reallocDualViewIfNeeded (this->numExportPacketsPerLID_,
                               numExportLIDs,
                               "numExportPacketsPerLID",
                               tooBigFactor,
                               true); // need fence before, if realloc'ing

    // If we reallocated above, then we fenced after that
    // reallocation.  This means that we don't need to fence again,
    // before the next reallocation.
    const bool needFenceBeforeNextAlloc = ! firstReallocated;
    const bool secondReallocated =
      reallocDualViewIfNeeded (this->numImportPacketsPerLID_,
                               numImportLIDs,
                               "numImportPacketsPerLID",
                               tooBigFactor,
                               needFenceBeforeNextAlloc);

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "DualView status after:" << endl
         << *prefix << dualViewStatusToString (this->numExportPacketsPerLID_,
                                               "numExportPacketsPerLID_")
         << endl
         << *prefix << dualViewStatusToString (this->numImportPacketsPerLID_,
                                               "numImportPacketsPerLID_")
         << endl;
      std::cerr << os.str ();
    }

    return firstReallocated || secondReallocated;
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  doTransferNew (const SrcDistObject& src,
                 const CombineMode CM,
                 const size_t numSameIDs,
                 const Kokkos::DualView<const local_ordinal_type*,
                   buffer_device_type>& permuteToLIDs,
                 const Kokkos::DualView<const local_ordinal_type*,
                   buffer_device_type>& permuteFromLIDs,
                 const Kokkos::DualView<const local_ordinal_type*,
                   buffer_device_type>& remoteLIDs,
                 const Kokkos::DualView<const local_ordinal_type*,
                   buffer_device_type>& exportLIDs,
                 Distributor& distor,
                 const ReverseOption revOp,
                 const bool commOnHost,
                 const bool restrictedMode)
  {
    using Details::Behavior;
    using ::Tpetra::Details::dualViewStatusToString;
    using ::Tpetra::Details::getArrayViewFromDualView;
    using Details::ProfilingRegion;
    using Kokkos::Compat::getArrayView;
    using Kokkos::Compat::getConstArrayView;
    using Kokkos::Compat::getKokkosViewDeepCopy;
    using Kokkos::Compat::create_const_view;
    using std::endl;
    using DT = device_type;
    using DES = typename DT::execution_space;
    const char funcName[] = "Tpetra::DistObject::doTransferNew";

    ProfilingRegion region_dTN(funcName);
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
    // FIXME (mfh 04 Feb 2019) Deprecate Teuchos::TimeMonitor in favor
    // of Kokkos profiling.
    Teuchos::TimeMonitor doXferMon (*doXferTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS

    const bool debug = Behavior::debug("DistObject");
    const bool verbose = Behavior::verbose("DistObject");
    // Prefix for verbose output.  Use a pointer, so we don't pay for
    // string construction unless needed.  We set this below.
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      prefix = this->createPrefix("DistObject", "doTransferNew");
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Input arguments:" << endl
         << *prefix << "  combineMode: " << combineModeToString (CM) << endl
         << *prefix << "  numSameIDs: " << numSameIDs << endl
         << *prefix << "  "
         << dualViewStatusToString (permuteToLIDs, "permuteToLIDs") << endl
         << *prefix << "  "
         << dualViewStatusToString (permuteFromLIDs, "permuteFromLIDs") << endl
         << *prefix << "  "
         << dualViewStatusToString (remoteLIDs, "remoteLIDs") << endl
         << *prefix << "  "
         << dualViewStatusToString (exportLIDs, "exportLIDs") << endl
         << *prefix << "  revOp: Do" << (revOp == DoReverse ? "Reverse" : "Forward") << endl
         << *prefix << "  commOnHost: " << (commOnHost ? "true" : "false") << endl;
      std::cerr << os.str ();
    }

    {
      ProfilingRegion region_cs ("Tpetra::DistObject::doTransferNew::checkSizes");
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "1. checkSizes" << endl;
        std::cerr << os.str ();
      }
      const bool checkSizesResult = this->checkSizes (src);
      TEUCHOS_TEST_FOR_EXCEPTION
        (! checkSizesResult, std::invalid_argument,
         "Tpetra::DistObject::doTransfer: checkSizes() indicates that the "
         "destination object is not a legal target for redistribution from the "
         "source object.  This probably means that they do not have the same "
         "dimensions.  For example, MultiVectors must have the same number of "
         "rows and columns.");
    }

    // NOTE (mfh 26 Apr 2016) Chris Baker's implementation understood
    // that if CM == INSERT || CM == REPLACE, the target object could
    // be write only.  We don't optimize for that here.

    if (!restrictedMode && numSameIDs + permuteToLIDs.extent (0) != 0) {
      // There is at least one GID to copy or permute.
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "2. copyAndPermute" << endl;
        std::cerr << os.str ();
      }
      ProfilingRegion region_cp
        ("Tpetra::DistObject::doTransferNew::copyAndPermute");
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
      // FIXME (mfh 04 Feb 2019) Deprecate Teuchos::TimeMonitor in favor
      // of Kokkos profiling.
      Teuchos::TimeMonitor copyAndPermuteMon (*copyAndPermuteTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS

      if (numSameIDs + permuteToLIDs.extent (0) != 0) {
        // There is at least one GID to copy or permute.
        if (verbose) {
          std::ostringstream os;
          os << *prefix << "2. copyAndPermute" << endl;
          std::cerr << os.str ();
        }
        this->copyAndPermute (src, numSameIDs, permuteToLIDs,
                              permuteFromLIDs);
        if (verbose) {
          std::ostringstream os;
          os << *prefix << "After copyAndPermute:" << endl
             << *prefix << "  "
             << dualViewStatusToString (permuteToLIDs, "permuteToLIDs")
             << endl
             << *prefix << "  "
             << dualViewStatusToString (permuteFromLIDs, "permuteFromLIDs")
             << endl;
          std::cerr << os.str ();
        }
      }
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
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "constantNumPackets=" << constantNumPackets << endl;
      std::cerr << os.str ();
    }

    // We only need to pack communication buffers if the combine mode
    // is not ZERO. A "ZERO combine mode" means that the results are
    // the same as if we had received all zeros, and added them to the
    // existing values. That means we don't need to communicate.
    if (CM != ZERO) {
      if (constantNumPackets == 0) {
        if (verbose) {
          std::ostringstream os;
          os << *prefix << "3. (Re)allocate num{Ex,Im}portPacketsPerLID"
             << endl;
          std::cerr << os.str ();
        }
        // This only reallocates if necessary, that is, if the sizes
        // don't match.
        this->reallocArraysForNumPacketsPerLid (exportLIDs.extent (0),
                                                remoteLIDs.extent (0));
      }

      if (verbose) {
        std::ostringstream os;
        os << *prefix << "4. packAndPrepare: before, "
           << dualViewStatusToString (this->exports_, "exports_")
           << endl;
        std::cerr << os.str ();
      }
      {
        ProfilingRegion region_pp
          ("Tpetra::DistObject::doTransferNew::packAndPrepare");
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
        // FIXME (mfh 04 Feb 2019) Deprecate Teuchos::TimeMonitor in
        // favor of Kokkos profiling.
        Teuchos::TimeMonitor packAndPrepareMon (*packAndPrepareTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS

        // Ask the source to pack data.  Also ask it whether there are
        // a constant number of packets per element
        // (constantNumPackets is an output argument).  If there are,
        // constantNumPackets will come back nonzero.  Otherwise, the
        // source will fill the numExportPacketsPerLID_ array.

        // FIXME (mfh 18 Oct 2017) if (! commOnHost), sync to device?
        // Alternately, make packAndPrepare take a "commOnHost"
        // argument to tell it where to leave the data?
        //
        // NOTE (mfh 04 Feb 2019) Subclasses of DistObject should have
        // the freedom to pack and unpack either on host or device.
        // We should prefer sync'ing only on demand.  Thus, we can
        // answer the above question: packAndPrepare should not
        // take a commOnHost argument, and doTransferNew should sync
        // where needed, if needed.
        this->packAndPrepare (src, exportLIDs, this->exports_,
                              this->numExportPacketsPerLID_,
                              constantNumPackets, distor);
        if (commOnHost) {
          if (this->exports_.need_sync_host ()) {
            this->exports_.sync_host ();
          }
        }
        else { // ! commOnHost
          if (this->exports_.need_sync_device ()) {
            this->exports_.sync_device ();
          }
        }
      }
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "5.1. After packAndPrepare, "
           << dualViewStatusToString (this->exports_, "exports_")
           << endl;
        std::cerr << os.str ();
      }
    } // if (CM != ZERO)

    // We only need to send data if the combine mode is not ZERO.
    if (CM != ZERO) {
      if (constantNumPackets != 0) {
        // There are a constant number of packets per element.  We
        // already know (from the number of "remote" (incoming)
        // elements) how many incoming elements we expect, so we can
        // resize the buffer accordingly.
        const size_t rbufLen = remoteLIDs.extent (0) * constantNumPackets;
        reallocImportsIfNeeded (rbufLen, verbose, prefix.get ());
      }

      // Do we need to do communication (via doPostsAndWaits)?
      bool needCommunication = true;

      // This may be NULL.  It will be used below.
      const this_type* srcDistObj = dynamic_cast<const this_type*> (&src);

      if (revOp == DoReverse && ! this->isDistributed ()) {
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

      if (! needCommunication) {
        if (verbose) {
          std::ostringstream os;
          os << *prefix << "Comm not needed; skipping" << endl;
          std::cerr << os.str ();
        }
      }
      else {
        ProfilingRegion region_dpw
          ("Tpetra::DistObject::doTransferNew::doPostsAndWaits");
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
        // FIXME (mfh 04 Feb 2019) Deprecate Teuchos::TimeMonitor in
        // favor of Kokkos profiling.
        Teuchos::TimeMonitor doPostsAndWaitsMon (*doPostsAndWaitsTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS

        if (verbose) {
          std::ostringstream os;
          os << *prefix << "7.0. "
             << (revOp == DoReverse ? "Reverse" : "Forward")
             << " mode" << endl;
          std::cerr << os.str ();
        }

        if (constantNumPackets == 0) { // variable num packets per LID
          if (verbose) {
            std::ostringstream os;
            os << *prefix << "7.1. Variable # packets / LID: first comm "
               << "(commOnHost = " << (commOnHost ? "true" : "false") << ")"
               << endl;
            std::cerr << os.str ();
          }
          size_t totalImportPackets = 0;
          if (commOnHost) {
            if (this->numExportPacketsPerLID_.need_sync_host ()) {
              this->numExportPacketsPerLID_.sync_host ();
            }
            if (this->numImportPacketsPerLID_.need_sync_host ()) {
              this->numImportPacketsPerLID_.sync_host ();
            }
            this->numImportPacketsPerLID_.modify_host (); // out arg
            auto numExp_h =
              create_const_view (this->numExportPacketsPerLID_.view_host ());
            auto numImp_h = this->numImportPacketsPerLID_.view_host ();

            // MPI communication happens here.
            if (verbose) {
              std::ostringstream os;
              os << *prefix << "Call do"
                 << (revOp == DoReverse ? "Reverse" : "") << "PostsAndWaits"
                 << endl;
              std::cerr << os.str ();
            }
            if (revOp == DoReverse) {
              distor.doReversePostsAndWaits (numExp_h, 1, numImp_h);
            }
            else {
              distor.doPostsAndWaits (numExp_h, 1, numImp_h);
            }
            DES().fence (); // just in case UVM doesn't behave right

            if (verbose) {
              std::ostringstream os;
              os << *prefix << "Count totalImportPackets" << std::endl;
              std::cerr << os.str ();
            }
            using the_dev_type = typename decltype (numImp_h)::device_type;
            totalImportPackets = countTotalImportPackets<the_dev_type> (numImp_h);
          }
          else { // ! commOnHost
            if (this->numExportPacketsPerLID_.need_sync_device ()) {
              this->numExportPacketsPerLID_.sync_device ();
            }
            if (this->numImportPacketsPerLID_.need_sync_device ()) {
              this->numImportPacketsPerLID_.sync_device ();
            }
            this->numImportPacketsPerLID_.modify_device (); // out arg
            auto numExp_d = create_const_view
              (this->numExportPacketsPerLID_.view_device ());
            auto numImp_d = this->numImportPacketsPerLID_.view_device ();

            // MPI communication happens here.
            if (verbose) {
              std::ostringstream os;
              os << *prefix << "Call do"
                 << (revOp == DoReverse ? "Reverse" : "") << "PostsAndWaits"
                 << endl;
              std::cerr << os.str ();
            }
            if (revOp == DoReverse) {
              distor.doReversePostsAndWaits (numExp_d, 1, numImp_d);
            }
            else {
              distor.doPostsAndWaits (numExp_d, 1, numImp_d);
            }
            DES().fence (); // just in case UVM doesn't behave right

            if (verbose) {
              std::ostringstream os;
              os << *prefix << "Count totalImportPackets" << std::endl;
              std::cerr << os.str ();
            }
            using the_dev_type = typename decltype (numImp_d)::device_type;
            totalImportPackets = countTotalImportPackets<the_dev_type> (numImp_d);
          }

          if (verbose) {
            std::ostringstream os;
            os << *prefix << "totalImportPackets=" << totalImportPackets << endl;
            std::cerr << os.str ();
          }
          this->reallocImportsIfNeeded (totalImportPackets, verbose,
                                        prefix.get ());
          if (verbose) {
            std::ostringstream os;
            os << *prefix << "7.3. Second comm" << std::endl;
            std::cerr << os.str ();
          }

          // mfh 04 Feb 2019: Distributor expects the "num packets per
          // LID" arrays on host, so that it can issue MPI sends and
          // receives correctly.
          if (this->numExportPacketsPerLID_.need_sync_host ()) {
            this->numExportPacketsPerLID_.sync_host ();
          }
          if (this->numImportPacketsPerLID_.need_sync_host ()) {
            this->numImportPacketsPerLID_.sync_host ();
          }

          // NOTE (mfh 25 Apr 2016, 01 Aug 2017) doPostsAndWaits and
          // doReversePostsAndWaits currently want
          // numExportPacketsPerLID and numImportPacketsPerLID as
          // Teuchos::ArrayView, rather than as Kokkos::View.
          //
          // NOTE (mfh 04 Feb 2019) This does NOT copy from host to
          // device.  The above syncs might.
          auto numExportPacketsPerLID_av =
            getArrayViewFromDualView (this->numExportPacketsPerLID_);
          auto numImportPacketsPerLID_av =
            getArrayViewFromDualView (this->numImportPacketsPerLID_);

          // imports_ is for output only, so we don't need to sync it
          // before marking it as modified.  However, in order to
          // prevent spurious debug-mode errors (e.g., "modified on
          // both device and host"), we first need to clear its
          // "modified" flags.
          this->imports_.clear_sync_state ();

          if (verbose) {
            std::ostringstream os;
            os << *prefix << "Comm on "
               << (commOnHost ? "host" : "device")
               << "; call do" << (revOp == DoReverse ? "Reverse" : "")
               << "PostsAndWaits" << endl;
            std::cerr << os.str ();
          }

          if (commOnHost) {
            this->imports_.modify_host ();
            if (revOp == DoReverse) {
              distor.doReversePostsAndWaits
                (create_const_view (this->exports_.view_host ()),
                 numExportPacketsPerLID_av,
                 this->imports_.view_host (),
                 numImportPacketsPerLID_av);
            }
            else {
              distor.doPostsAndWaits
                (create_const_view (this->exports_.view_host ()),
                 numExportPacketsPerLID_av,
                 this->imports_.view_host (),
                 numImportPacketsPerLID_av);
            }
          }
          else { // pack on device
            this->imports_.modify_device ();
            if (revOp == DoReverse) {
              distor.doReversePostsAndWaits
                (create_const_view (this->exports_.view_device ()),
                 numExportPacketsPerLID_av,
                 this->imports_.view_device (),
                 numImportPacketsPerLID_av);
            }
            else {
              distor.doPostsAndWaits
                (create_const_view (this->exports_.view_device ()),
                 numExportPacketsPerLID_av,
                 this->imports_.view_device (),
                 numImportPacketsPerLID_av);
            }
          }
        }
        else { // constant number of packets per LID
          if (verbose) {
            std::ostringstream os;
            os << *prefix << "7.1. Const # packets per LID: " << endl
               << *prefix << "  "
               << dualViewStatusToString (this->exports_, "exports_")
               << endl
               << *prefix << "  "
               << dualViewStatusToString (this->exports_, "imports_")
               << endl;
            std::cerr << os.str ();
          }
          // imports_ is for output only, so we don't need to sync it
          // before marking it as modified.  However, in order to
          // prevent spurious debug-mode errors (e.g., "modified on
          // both device and host"), we first need to clear its
          // "modified" flags.
          this->imports_.clear_sync_state ();

          if (verbose) {
            std::ostringstream os;
            os << *prefix << "7.2. Comm on "
               << (commOnHost ? "host" : "device")
               << "; call do" << (revOp == DoReverse ? "Reverse" : "")
               << "PostsAndWaits" << endl;
            std::cerr << os.str ();
          }
          if (commOnHost) {
            this->imports_.modify_host ();
            if (revOp == DoReverse) {
              distor.doReversePostsAndWaits
                (create_const_view (this->exports_.view_host ()),
                 constantNumPackets,
                 this->imports_.view_host ());
            }
            else {
              distor.doPostsAndWaits
                (create_const_view (this->exports_.view_host ()),
                 constantNumPackets,
                 this->imports_.view_host ());
            }
          }
          else { // pack on device
            this->imports_.modify_device ();
            if (revOp == DoReverse) {
              distor.doReversePostsAndWaits
                (create_const_view (this->exports_.view_device ()),
                 constantNumPackets,
                 this->imports_.view_device ());
            }
            else {
              distor.doPostsAndWaits
                (create_const_view (this->exports_.view_device ()),
                 constantNumPackets,
                 this->imports_.view_device ());
            }
          } // commOnHost
        } // constant or variable num packets per LID

        if (verbose) {
          std::ostringstream os;
          os << *prefix << "8. unpackAndCombine" << endl;
          std::cerr << os.str ();
        }
        ProfilingRegion region_uc
          ("Tpetra::DistObject::doTransferNew::unpackAndCombine");
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
        // FIXME (mfh 04 Feb 2019) Deprecate Teuchos::TimeMonitor in
        // favor of Kokkos profiling.
        Teuchos::TimeMonitor unpackAndCombineMon (*unpackAndCombineTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS

        if (debug) {
          std::ostringstream lclErrStrm;
          bool lclSuccess = false;
          try {
            this->unpackAndCombine (remoteLIDs, this->imports_,
                                    this->numImportPacketsPerLID_,
                                    constantNumPackets, distor, CM);
            lclSuccess = true;
          }
          catch (std::exception& e) {
            lclErrStrm << "unpackAndCombine threw an exception: "
                       << endl << e.what();
          }
          catch (...) {
            lclErrStrm << "unpackAndCombine threw an exception "
              "not a subclass of std::exception.";
          }
          const char gblErrMsgHeader[] = "Tpetra::DistObject::"
            "doTransferNew threw an exception in unpackAndCombine on "
            "one or more processes in the DistObject's communicator.";
          auto comm = getMap()->getComm();
          Details::checkGlobalError(std::cerr, lclSuccess,
                                    lclErrStrm.str().c_str(),
                                    gblErrMsgHeader, *comm);
        }
        else {
          this->unpackAndCombine (remoteLIDs, this->imports_,
                                  this->numImportPacketsPerLID_,
                                  constantNumPackets, distor, CM);
        }
      } // if (needCommunication)
    } // if (CM != ZERO)

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "9. Done!" << endl;
      std::cerr << os.str ();
    }
  }


  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  copyAndPermute
  (const SrcDistObject&,
   const size_t,
   const Kokkos::DualView<
     const local_ordinal_type*,
     buffer_device_type>&,
   const Kokkos::DualView<
     const local_ordinal_type*,
     buffer_device_type>&)
  {}

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  packAndPrepare
  (const SrcDistObject&,
   const Kokkos::DualView<
     const local_ordinal_type*,
     buffer_device_type>&,
   Kokkos::DualView<
     packet_type*,
     buffer_device_type>&,
   Kokkos::DualView<
     size_t*,
     buffer_device_type>,
   size_t&,
   Distributor&)
  {}

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  unpackAndCombine
  (const Kokkos::DualView<
     const local_ordinal_type*,
     buffer_device_type>& /* importLIDs */,
   Kokkos::DualView<
     packet_type*,
     buffer_device_type> /* imports */,
   Kokkos::DualView<
     size_t*,
     buffer_device_type> /* numPacketsPerLID */,
   const size_t /* constantNumPackets */,
   Distributor& /* distor */,
   const CombineMode /* combineMode */)
  {}


  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  print (std::ostream& os) const
  {
    using Teuchos::FancyOStream;
    using Teuchos::getFancyOStream;
    using Teuchos::RCP;
    using Teuchos::rcpFromRef;
    using std::endl;

    RCP<FancyOStream> out = getFancyOStream (rcpFromRef (os));
    this->describe (*out, Teuchos::VERB_DEFAULT);
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::unique_ptr<std::string>
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  createPrefix(const char className[],
               const char methodName[]) const
  {
    auto map = this->getMap();
    auto comm = map.is_null() ? Teuchos::null : map->getComm();
    return Details::createPrefix(
      comm.getRawPtr(), className, methodName);
  }

  template<class DistObjectType>
  void
  removeEmptyProcessesInPlace(
    Teuchos::RCP<DistObjectType>& input,
    const Teuchos::RCP<const Map<
      typename DistObjectType::local_ordinal_type,
      typename DistObjectType::global_ordinal_type,
      typename DistObjectType::node_type>>& newMap)
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
    auto newMap = input->getMap ()->removeEmptyProcesses ();
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
