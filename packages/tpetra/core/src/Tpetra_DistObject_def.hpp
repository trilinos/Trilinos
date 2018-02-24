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
#include "Tpetra_Details_reallocDualViewIfNeeded.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include <memory>

namespace Tpetra {

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
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  DistObject (const DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>& rhs) :
    map_ (rhs.map_)
  {}

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  ~DistObject ()
  {}

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

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
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

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  doImport (const SrcDistObject& source,
            const Import<LocalOrdinal, GlobalOrdinal, Node>& importer,
            CombineMode CM)
  {
    using std::endl;
    const char modeString[] = "doImport (forward mode)";

    // mfh 18 Oct 2017: Set TPETRA_VERBOSE to true for copious debug
    // output to std::cerr on every MPI process.  This is unwise for
    // runs with large numbers of MPI processes.
    const bool verbose = ::Tpetra::Details::Behavior::verbose ();
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      int myRank = 0;
      auto map = this->getMap ();
      if (! map.is_null ()) {
        auto comm = map->getComm ();
        if (! comm.is_null ()) {
          myRank = comm->getRank ();
        }
      }
      prefix = [myRank] () {
        std::ostringstream os;
        os << "(Proc " << myRank << ") ";
        return std::unique_ptr<std::string> (new std::string (os.str ()));
      } ();
      std::ostringstream os;
      os << *prefix << "Tpetra::DistObject::" << modeString << ":" << endl;
      std::cerr << os.str ();
    }
    this->doTransfer (source, importer, modeString, DoForward, CM);
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Tpetra::DistObject::" << modeString << ": Done!"
         << endl;
      std::cerr << os.str ();
    }
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  doExport (const SrcDistObject& source,
            const Export<LocalOrdinal, GlobalOrdinal, Node>& exporter,
            CombineMode CM)
  {
    using std::endl;
    const char modeString[] = "doExport (forward mode)";

    // mfh 18 Oct 2017: Set TPETRA_VERBOSE to true for copious debug
    // output to std::cerr on every MPI process.  This is unwise for
    // runs with large numbers of MPI processes.
    const bool verbose = ::Tpetra::Details::Behavior::verbose ();
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      int myRank = 0;
      auto map = this->getMap ();
      if (! map.is_null ()) {
        auto comm = map->getComm ();
        if (! comm.is_null ()) {
          myRank = comm->getRank ();
        }
      }
      prefix = [myRank] () {
        std::ostringstream os;
        os << "(Proc " << myRank << ") ";
        return std::unique_ptr<std::string> (new std::string (os.str ()));
      } ();
      std::ostringstream os;
      os << *prefix << "Tpetra::DistObject::" << modeString << ":" << endl;
      std::cerr << os.str ();
    }
    this->doTransfer (source, exporter, modeString, DoForward, CM);
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Tpetra::DistObject::" << modeString << ": Done!"
         << endl;
      std::cerr << os.str ();
    }
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  doImport (const SrcDistObject& source,
            const Export<LocalOrdinal, GlobalOrdinal, Node>& exporter,
            CombineMode CM)
  {
    using std::endl;
    const char modeString[] = "doImport (reverse mode)";

    // mfh 18 Oct 2017: Set TPETRA_VERBOSE to true for copious debug
    // output to std::cerr on every MPI process.  This is unwise for
    // runs with large numbers of MPI processes.
    const bool verbose = ::Tpetra::Details::Behavior::verbose ();
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      int myRank = 0;
      auto map = this->getMap ();
      if (! map.is_null ()) {
        auto comm = map->getComm ();
        if (! comm.is_null ()) {
          myRank = comm->getRank ();
        }
      }
      prefix = [myRank] () {
        std::ostringstream os;
        os << "(Proc " << myRank << ") ";
        return std::unique_ptr<std::string> (new std::string (os.str ()));
      } ();
      std::ostringstream os;
      os << *prefix << "Tpetra::DistObject::" << modeString << ":" << endl;
      std::cerr << os.str ();
    }
    this->doTransfer (source, exporter, modeString, DoReverse, CM);
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Tpetra::DistObject::" << modeString << ": Done!"
         << endl;
      std::cerr << os.str ();
    }
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  doExport (const SrcDistObject& source,
            const Import<LocalOrdinal, GlobalOrdinal, Node> & importer,
            CombineMode CM)
  {
    using std::endl;
    const char modeString[] = "doExport (reverse mode)";

    // mfh 18 Oct 2017: Set TPETRA_VERBOSE to true for copious debug
    // output to std::cerr on every MPI process.  This is unwise for
    // runs with large numbers of MPI processes.
    const bool verbose = ::Tpetra::Details::Behavior::verbose ();
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      int myRank = 0;
      auto map = this->getMap ();
      if (! map.is_null ()) {
        auto comm = map->getComm ();
        if (! comm.is_null ()) {
          myRank = comm->getRank ();
        }
      }
      prefix = [myRank] () {
        std::ostringstream os;
        os << "(Proc " << myRank << ") ";
        return std::unique_ptr<std::string> (new std::string (os.str ()));
      } ();
      std::ostringstream os;
      os << *prefix << "Tpetra::DistObject::" << modeString << ":" << endl;
      std::cerr << os.str ();
    }
    this->doTransfer (source, importer, modeString, DoReverse, CM);
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Tpetra::DistObject::" << modeString << ": Done!"
         << endl;
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
              const Details::Transfer<local_ordinal_type, global_ordinal_type, node_type>& transfer,
              const char modeString[],
              const ReverseOption revOp,
              const CombineMode CM)
  {
    using Tpetra::Details::getDualViewCopyFromArrayView;
    using std::endl;
    typedef LocalOrdinal LO;
    typedef device_type DT;

    // mfh 18 Oct 2017: Set TPETRA_DEBUG to true to enable extra debug
    // checks.  These may communicate more.
    const bool debug = ::Tpetra::Details::Behavior::debug ();
    if (debug) {
      if (revOp == DoForward) {
        const bool myMapSameAsTransferTgtMap =
          this->getMap ()->isSameAs (* (transfer.getTargetMap ()));
        TEUCHOS_TEST_FOR_EXCEPTION
          (! myMapSameAsTransferTgtMap, std::invalid_argument,
           "Tpetra::DistObject::" << modeString << ": For forward-mode "
           "communication, the target DistObject's Map must be the same "
           "(in the sense of Tpetra::Map::isSameAs) as the input "
           "Export/Import object's target Map.");
      }
      else { // revOp == DoReverse
        const bool myMapSameAsTransferSrcMap =
          this->getMap ()->isSameAs (* (transfer.getSourceMap ()));
        TEUCHOS_TEST_FOR_EXCEPTION
          (! myMapSameAsTransferSrcMap, std::invalid_argument,
           "Tpetra::DistObject::" << modeString << ": For reverse-mode "
           "communication, the target DistObject's Map must be the same "
           "(in the sense of Tpetra::Map::isSameAs) as the input "
           "Export/Import object's source Map.");
      }

      // SrcDistObject need not even _have_ Maps.  However, if the
      // source object is a DistObject, it has a Map, and we may
      // compare that Map with the Transfer's Maps.
      const this_type* srcDistObj = dynamic_cast<const this_type*> (&src);
      if (srcDistObj != NULL) {
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

    // mfh 03 Aug 2017, 17 Oct 2017: Set TPETRA_VERBOSE to true for
    // copious debug output to std::cerr on every MPI process.  This
    // is unwise for runs with large numbers of MPI processes.
    const bool verbose = ::Tpetra::Details::Behavior::verbose ();
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      int myRank = 0;
      auto map = this->getMap ();
      if (! map.is_null ()) {
        auto comm = map->getComm ();
        if (! comm.is_null ()) {
          myRank = comm->getRank ();
        }
      }
      prefix = [myRank] () {
        std::ostringstream os;
        os << "(Proc " << myRank << ") ";
        return std::unique_ptr<std::string> (new std::string (os.str ()));
      } ();
      std::ostringstream os;
      os << *prefix << "Tpetra::DistObject::doTransfer:" << endl;
      std::cerr << os.str ();
    }

    const size_t numSameIDs = transfer.getNumSameIDs ();
    typedef Teuchos::ArrayView<const LocalOrdinal> view_type;
    const view_type permuteToLIDs_ = (revOp == DoForward) ?
      transfer.getPermuteToLIDs () : transfer.getPermuteFromLIDs ();
    const view_type permuteFromLIDs_ = (revOp == DoForward) ?
      transfer.getPermuteFromLIDs () : transfer.getPermuteToLIDs ();
    const view_type exportLIDs_ = (revOp == DoForward) ?
      transfer.getExportLIDs () : transfer.getRemoteLIDs ();
    const view_type remoteLIDs_ = (revOp == DoForward) ?
      transfer.getRemoteLIDs () : transfer.getExportLIDs ();
    Distributor& distor = transfer.getDistributor ();

    if (this->useNewInterface ()) {
      using ::Tpetra::Details::Behavior;
      // Do we need all communication buffers to live on host?
      const bool commOnHost = ! Behavior::assumeMpiIsCudaAware ();
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "doTransfer: Use new interface; "
          "commOnHost=" << (commOnHost ? "true" : "false") << endl;
        std::cerr << os.str ();
      }

      // Convert arguments to Kokkos::DualView.  This currently
      // involves deep copy, either to host or to device (depending on
      // commOnHost).  At some point, we need to change the interface
      // of doTransfer so it takes DualView (or just View) rather than
      // Teuchos::ArrayView, so that we won't need this deep copy.
      //
      // We don't need to sync the arguments.  commOnHost determines
      // where the most recent version lives.
      Kokkos::DualView<LO*, DT> permuteToLIDs =
        getDualViewCopyFromArrayView<LO, DT> (permuteToLIDs_,
                                              "permuteToLIDs",
                                              commOnHost);
      Kokkos::DualView<LO*, DT> permuteFromLIDs =
        getDualViewCopyFromArrayView<LO, DT> (permuteFromLIDs_,
                                              "permuteFromLIDs",
                                              commOnHost);
      // No need to sync this.  packAndPrepareNew will use it to
      // determine where to pack (in host or device memory).
      Kokkos::DualView<LO*, DT> remoteLIDs =
        getDualViewCopyFromArrayView<LO, DT> (remoteLIDs_,
                                              "remoteLIDs",
                                              commOnHost);
      Kokkos::DualView<LO*, DT> exportLIDs =
        getDualViewCopyFromArrayView<LO, DT> (exportLIDs_,
                                              "exportLIDs",
                                              commOnHost);
      doTransferNew (src, CM, numSameIDs, permuteToLIDs, permuteFromLIDs,
                     remoteLIDs, exportLIDs, distor, revOp, commOnHost);
    }
    else {
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "doTransfer: Use old interface" << endl;
        std::cerr << os.str ();
      }
      doTransferOld (src, CM, numSameIDs, permuteToLIDs_, permuteFromLIDs_,
                     remoteLIDs_, exportLIDs_, distor, revOp);
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Tpetra::DistObject::doTransfer: Done!" << endl;
      std::cerr << os.str ();
    }
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  reallocImportsIfNeeded (const size_t newSize, const bool verbose)
  {
    if (verbose) {
      const int myRank = this->getMap ()->getComm ()->getRank ();
      std::ostringstream os;
      os << "(Proc " << myRank << ") Reallocate (if needed) imports_ from "
         << imports_.dimension_0 () << " to " << newSize << std::endl;
      std::cerr << os.str ();
    }
    using Details::reallocDualViewIfNeeded;
    const bool reallocated =
      reallocDualViewIfNeeded (this->imports_, newSize, "imports");
    if (verbose) {
      const int myRank = this->getMap ()->getComm ()->getRank ();
      std::ostringstream os;
      os << "(Proc " << myRank << ") Finished reallocating imports_"
         << std::endl;
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
    using Details::reallocDualViewIfNeeded;
    using Details::dualViewStatusToString;
    using std::endl;
    // If an array is already allocated, and if is at least
    // tooBigFactor times bigger than it needs to be, free it and
    // reallocate to the size we need, in order to save space.
    // Otherwise, take subviews to reduce allocation size.
    constexpr size_t tooBigFactor = 10;

    const bool verbose = ::Tpetra::Details::Behavior::verbose ();
    if (verbose) {
      const int myRank = this->getMap ()->getComm ()->getRank ();
      std::ostringstream os;
      os << "(Proc " << myRank << ") reallocArraysForNumPacketsPerLid before:"
         << endl
         << "(Proc " << myRank << ")   "
         << dualViewStatusToString (this->numExportPacketsPerLID_, "numExportPacketsPerLID_")
         << endl
         << "(Proc " << myRank << ")   "
         << dualViewStatusToString (this->numImportPacketsPerLID_, "numImportPacketsPerLID_")
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
      const int myRank = this->getMap ()->getComm ()->getRank ();
      std::ostringstream os;
      os << "(Proc " << myRank << ") reallocArraysForNumPacketsPerLid before:"
         << endl
         << "(Proc " << myRank << ")   "
         << dualViewStatusToString (this->numExportPacketsPerLID_, "numExportPacketsPerLID_")
         << endl
         << "(Proc " << myRank << ")   "
         << dualViewStatusToString (this->numImportPacketsPerLID_, "numImportPacketsPerLID_")
         << endl;
      std::cerr << os.str ();
    }

    return firstReallocated || secondReallocated;
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
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
    using Tpetra::Details::getArrayViewFromDualView;
    using Tpetra::Details::reallocDualViewIfNeeded;

    // mfh 03 Aug 2017: Set this to true for copious debug output to
    // std::cerr on every MPI process.  This is unwise for runs with
    // large numbers of MPI processes.
    constexpr bool debug = false;

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
        this->reallocArraysForNumPacketsPerLid (exportLIDs.size (),
                                                remoteLIDs.size ());
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
        reallocDualViewIfNeeded (this->exports_, exportsLen, "exports");
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
            {
              typedef typename Kokkos::DualView<size_t*,
                device_type>::t_host::execution_space host_exec_space;
              typedef Kokkos::RangePolicy<host_exec_space, Array_size_type> range_type;
              const size_t* const arrayToSum = numImportPacketsPerLID.getRawPtr ();
              Kokkos::parallel_reduce ("Count import packets",
                                       range_type (0, numImportPacketsPerLID.size ()),
                                       [=] (const Array_size_type& i, size_t& lclSum) {
                                         lclSum += arrayToSum[i];
                                       }, totalImportPackets);
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
            {
              typedef typename Kokkos::DualView<size_t*,
                device_type>::t_host::execution_space host_exec_space;
              typedef Kokkos::RangePolicy<host_exec_space, Array_size_type> range_type;
              const size_t* const arrayToSum = numImportPacketsPerLID.getRawPtr ();
              Kokkos::parallel_reduce ("Count import packets",
                                       range_type (0, numImportPacketsPerLID.size ()),
                                       [=] (const Array_size_type& i, size_t& lclSum) {
                                         lclSum += arrayToSum[i];
                                       }, totalImportPackets);
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

  namespace { // (anonymous)
    template<class DeviceType, class IndexType = size_t>
    struct SumFunctor {
      SumFunctor (const Kokkos::View<const size_t*, DeviceType>& viewToSum) :
        viewToSum_ (viewToSum) {}
      KOKKOS_FUNCTION void operator() (const IndexType& i, size_t& lclSum) const {
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

      const IndexType numOut = numImportPacketsPerLID.dimension_0 ();
      size_t totalImportPackets = 0;
      parallel_reduce ("Count import packets",
                       range_type (0, numOut),
                       SumFunctor<DeviceType, IndexType> (numImportPacketsPerLID),
                       totalImportPackets);
      return totalImportPackets;
    }
  } // namespace (anonymous)

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  doTransferNew (const SrcDistObject& src,
                 const CombineMode CM,
                 const size_t numSameIDs,
                 const Kokkos::DualView<const local_ordinal_type*,
                   device_type>& permuteToLIDs,
                 const Kokkos::DualView<const local_ordinal_type*,
                   device_type>& permuteFromLIDs,
                 const Kokkos::DualView<const local_ordinal_type*,
                   device_type>& remoteLIDs,
                 const Kokkos::DualView<const local_ordinal_type*,
                   device_type>& exportLIDs,
                 Distributor& distor,
                 const ReverseOption revOp,
                 const bool commOnHost)
  {
    using Details::dualViewStatusToString;
    using Tpetra::Details::getArrayViewFromDualView;
    using Kokkos::Compat::getArrayView;
    using Kokkos::Compat::getConstArrayView;
    using Kokkos::Compat::getKokkosViewDeepCopy;
    using Kokkos::Compat::create_const_view;
    using std::endl;
    typedef LocalOrdinal LO;
    typedef device_type DT;

    typedef typename Kokkos::DualView<LO*, DT>::t_dev::execution_space DES;
    //typedef typename Kokkos::DualView<LO*, DT>::t_dev::memory_space DMS; // unused
    //typedef typename Kokkos::DualView<LO*, DT>::t_dev::memory_space HMS; // unused

    // DistObject's communication buffers (exports_,
    // numExportPacketsPerLID_, imports_, and numImportPacketsPerLID_)
    // may have different memory spaces than device_type would
    // indicate.  See GitHub issue #1088.  Abbreviations: "communication
    // host memory space" and "communication device memory space."
    typedef typename Kokkos::DualView<size_t*,
      buffer_device_type>::t_dev::memory_space CDMS;
    typedef typename Kokkos::DualView<size_t*,
      buffer_device_type>::t_host::memory_space CHMS;

    // mfh 03 Aug 2017, 17 Oct 2017: Set TPETRA_VERBOSE to true for
    // copious debug output to std::cerr on every MPI process.  This
    // is unwise for runs with large numbers of MPI processes.
    const bool verbose = ::Tpetra::Details::Behavior::verbose ();
    // Prefix for verbose output.  Use a pointer, so we don't pay for
    // string construction unless needed.  We set this below.
    std::unique_ptr<std::string> prefix;
    if (verbose) {
      auto map = this->getMap ();
      auto comm = map.is_null () ? Teuchos::null : map->getComm ();
      const int myRank = comm.is_null () ? 0 : comm->getRank ();
      std::ostringstream os;
      os << "(Proc " << myRank << ") ";
      prefix = std::unique_ptr<std::string> (new std::string (os.str ()));
    }

    if (verbose) {
      std::ostringstream os;
      os << *prefix << "Tpetra::CrsMatrix::doTransferNew: Input arguments:" << endl
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

#ifdef HAVE_TPETRA_TRANSFER_TIMERS
    Teuchos::TimeMonitor doXferMon (*doXferTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS

    {
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

    if (numSameIDs + permuteToLIDs.dimension_0 () != 0) {
      // There is at least one GID to copy or permute.
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "2. copyAndPermuteNew" << endl;
        std::cerr << os.str ();
      }
      {
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
        Teuchos::TimeMonitor copyAndPermuteMon (*copyAndPermuteTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
        this->copyAndPermuteNew (src, numSameIDs, permuteToLIDs,
                                 permuteFromLIDs);
      }
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "After copyAndPermuteNew:" << endl
           << *prefix << "  "
           << dualViewStatusToString (permuteToLIDs, "permuteToLIDs")
           << endl
           << *prefix << "  "
           << dualViewStatusToString (permuteFromLIDs, "permuteFromLIDs")
           << endl;
        std::cerr << os.str ();
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
        this->reallocArraysForNumPacketsPerLid (exportLIDs.dimension_0 (),
                                                remoteLIDs.dimension_0 ());
      }

      if (verbose) {
        std::ostringstream os;
        os << *prefix << "4. packAndPrepareNew: before, "
           << dualViewStatusToString (this->exports_, "exports_")
           << endl;
        std::cerr << os.str ();
      }
      {
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
        Teuchos::TimeMonitor packAndPrepareMon (*packAndPrepareTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
        // Ask the source to pack data.  Also ask it whether there are
        // a constant number of packets per element
        // (constantNumPackets is an output argument).  If there are,
        // constantNumPackets will come back nonzero.  Otherwise, the
        // source will fill the numExportPacketsPerLID_ array.
        this->packAndPrepareNew (src, exportLIDs, this->exports_,
                                 this->numExportPacketsPerLID_,
                                 constantNumPackets, distor);
        // FIXME (mfh 18 Oct 2017) if (! commOnHost), sync to device?
        // Alternately, make packAndPrepareNew take a "commOnHost"
        // argument to tell it where to leave the data?
        if (commOnHost) {
          typedef typename Kokkos::View<char*, buffer_device_type>::HostMirror::device_type
            buffer_host_device_type;
          typedef typename buffer_host_device_type::memory_space
            buffer_host_memory_space;
          this->exports_.template sync<buffer_host_memory_space> ();
        }
        else { // ! commOnHost
          typedef typename buffer_device_type::memory_space buffer_dev_memory_space;
          this->exports_.template sync<buffer_dev_memory_space> ();
        }
      }
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "5.1. After packAndPrepareNew, "
           << dualViewStatusToString (this->exports_, "exports_")
           << endl;
        std::cerr << os.str ();
      }
    } // if (CM != ZERO)

    // We only need to send data if the combine mode is not ZERO.
    if (CM != ZERO) {
      if (constantNumPackets != 0) {
        if (verbose) {
          std::ostringstream os;
          os << *prefix << "6. Realloc imports_" << std::endl;
          std::cerr << os.str ();
        }
        // There are a constant number of packets per element.  We
        // already know (from the number of "remote" (incoming)
        // elements) how many incoming elements we expect, so we can
        // resize the buffer accordingly.
        const size_t rbufLen = remoteLIDs.dimension_0 () * constantNumPackets;
        reallocImportsIfNeeded (rbufLen, verbose);
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

      if (verbose) {
        std::ostringstream os;
        os << *prefix << "needCommunication="
           << (needCommunication ? "true" : "false") << endl;
        std::cerr << os.str ();
      }

      // FIXME (mfh 17 Feb 2014) Distributor doesn't actually inspect
      // the contents of the "exports" or "imports" arrays, other than
      // to do a deep copy in the (should be technically unnecessary,
      // but isn't for some odd reason) "self-message" case.
      // Distributor thus doesn't need host views; it could do just
      // fine with device views, assuming that MPI knows how to read
      // device memory (which doesn't even require UVM).

      if (needCommunication) {
        if (revOp == DoReverse) {
          if (verbose) {
            std::ostringstream os;
            os << *prefix << "7.0. Reverse mode" << endl;
            std::cerr << os.str ();
          }
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
          Teuchos::TimeMonitor doPostsAndWaitsMon (*doPostsAndWaitsTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
          if (constantNumPackets == 0) { //variable num-packets-per-LID:
            if (verbose) {
              std::ostringstream os;
              os << *prefix << "7.1. Variable # packets / LID: first comm "
                 << "(commOnHost = " << (commOnHost ? "true" : "false") << ")"
                 << endl;
              std::cerr << os.str ();
            }
            size_t totalImportPackets = 0;
            if (commOnHost) {
              this->numExportPacketsPerLID_.template sync<CHMS> ();
              this->numImportPacketsPerLID_.template sync<CHMS> ();
              this->numImportPacketsPerLID_.template modify<CHMS> (); // output argument
              auto numExp_h = create_const_view (this->numExportPacketsPerLID_.template view<CHMS> ());
              auto numImp_h = this->numImportPacketsPerLID_.template view<CHMS> ();

              // MPI communication happens here.
              distor.doReversePostsAndWaits (numExp_h, 1, numImp_h);

              DES::fence (); // just in case UVM doesn't behave right
              typedef typename decltype (numImp_h)::device_type the_dev_type;
              totalImportPackets = countTotalImportPackets<the_dev_type> (numImp_h);
            }
            else {
              this->numExportPacketsPerLID_.template sync<CDMS> ();
              this->numImportPacketsPerLID_.template sync<CDMS> ();
              this->numImportPacketsPerLID_.template modify<CDMS> (); // output argument
              auto numExp_d = create_const_view (this->numExportPacketsPerLID_.template view<CDMS> ());
              auto numImp_d = this->numImportPacketsPerLID_.template view<CDMS> ();

              // MPI communication happens here.
              distor.doReversePostsAndWaits (numExp_d, 1, numImp_d);

              DES::fence (); // just in case UVM doesn't behave right
              typedef typename decltype (numImp_d)::device_type the_dev_type;
              totalImportPackets = countTotalImportPackets<the_dev_type> (numImp_d);
            }

            if (verbose) {
              std::ostringstream os;
              os << *prefix << "totalImportPackets=" << totalImportPackets
                 << endl;
              std::cerr << os.str ();
            }
            this->reallocImportsIfNeeded (totalImportPackets, verbose);
            if (verbose) {
              std::ostringstream os;
              os << *prefix << "7.3. Second comm" << std::endl;
              std::cerr << os.str ();
            }

            // NOTE (mfh 25 Apr 2016, 01 Aug 2017) Since we need to
            // launch MPI communication on host, we will need
            // numExportPacketsPerLID and numImportPacketsPerLID on
            // host.
            this->numExportPacketsPerLID_.template sync<CHMS> ();
            this->numImportPacketsPerLID_.template sync<CHMS> ();

            // NOTE (mfh 25 Apr 2016, 01 Aug 2017) doPostsAndWaits and
            // doReversePostsAndWaits currently want
            // numExportPacketsPerLID and numImportPacketsPerLID as
            // Teuchos::ArrayView, rather than as Kokkos::View.
            auto numExportPacketsPerLID_av =
              getArrayViewFromDualView (this->numExportPacketsPerLID_);
            auto numImportPacketsPerLID_av =
              getArrayViewFromDualView (this->numImportPacketsPerLID_);

            // imports_ is for output only, so we don't need to sync
            // it before marking it as modified.  However, in order to
            // prevent spurious debug-mode errors (e.g., "modified on
            // both device and host"), we first need to clear its
            // "modified" flags.
            this->imports_.modified_device() = 0;
            this->imports_.modified_host() = 0;

            if (commOnHost) {
              this->imports_.template modify<CHMS> ();
              distor.doReversePostsAndWaits (create_const_view (this->exports_.template view<CHMS> ()),
                                             numExportPacketsPerLID_av,
                                             this->imports_.template view<CHMS> (),
                                             numImportPacketsPerLID_av);
            }
            else {
              this->imports_.template modify<CDMS> ();
              distor.doReversePostsAndWaits (create_const_view (this->exports_.template view<CDMS> ()),
                                             numExportPacketsPerLID_av,
                                             this->imports_.template view<CDMS> (),
                                             numImportPacketsPerLID_av);
            }
          }
          else {
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

            // imports_ is for output only, so we don't need to sync
            // it before marking it as modified.  However, in order to
            // prevent spurious debug-mode errors (e.g., "modified on
            // both device and host"), we first need to clear its
            // "modified" flags.
            this->imports_.modified_device() = 0;
            this->imports_.modified_host() = 0;

            if (commOnHost) {
              this->imports_.template modify<CHMS> ();
              distor.doReversePostsAndWaits (create_const_view (this->exports_.template view<CHMS> ()),
                                             constantNumPackets,
                                             this->imports_.template view<CHMS> ());
            }
            else { // pack on device
              this->imports_.template modify<CDMS> ();
              distor.doReversePostsAndWaits (create_const_view (this->exports_.template view<CDMS> ()),
                                             constantNumPackets,
                                             this->imports_.template view<CDMS> ());
            }
          }
        }
        else { // revOp == DoForward
          if (verbose) {
            std::cerr << ">>> 7.0. Forward mode" << std::endl;
          }

#ifdef HAVE_TPETRA_TRANSFER_TIMERS
          Teuchos::TimeMonitor doPostsAndWaitsMon (*doPostsAndWaitsTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
          if (constantNumPackets == 0) { //variable num-packets-per-LID:
            if (verbose) {
              std::cerr << ">>> 7.1. Variable # packets / LID: first comm" << std::endl;
            }

            size_t totalImportPackets = 0;
            if (commOnHost) {
              this->numExportPacketsPerLID_.template sync<CHMS> ();
              this->numImportPacketsPerLID_.template sync<CHMS> ();
              this->numImportPacketsPerLID_.template modify<CHMS> (); // output argument
              auto numExp_h = create_const_view (this->numExportPacketsPerLID_.template view<CHMS> ());
              auto numImp_h = this->numImportPacketsPerLID_.template view<CHMS> ();

              // MPI communication happens here.
              distor.doPostsAndWaits (numExp_h, 1, numImp_h);

              DES::fence (); // just in case UVM doesn't behave right
              typedef typename decltype (numImp_h)::device_type the_dev_type;
              totalImportPackets = countTotalImportPackets<the_dev_type> (numImp_h);
            }
            else {
              this->numExportPacketsPerLID_.template sync<CDMS> ();
              this->numImportPacketsPerLID_.template sync<CDMS> ();
              this->numImportPacketsPerLID_.template modify<CDMS> (); // output argument
              auto numExp_d = create_const_view (this->numExportPacketsPerLID_.template view<CDMS> ());
              auto numImp_d = this->numImportPacketsPerLID_.template view<CDMS> ();

              // MPI communication happens here.
              distor.doPostsAndWaits (numExp_d, 1, numImp_d);

              DES::fence (); // just in case UVM doesn't behave right
              typedef typename decltype (numImp_d)::device_type the_dev_type;
              totalImportPackets = countTotalImportPackets<the_dev_type> (numImp_d);
            }

            this->reallocImportsIfNeeded (totalImportPackets, verbose);

            if (verbose) {
              std::cerr << ">>> 7.3. Second comm" << std::endl;
            }

            // NOTE (mfh 25 Apr 2016, 01 Aug 2017) Since we need to
            // launch MPI communication on host, we will need
            // numExportPacketsPerLID and numImportPacketsPerLID on
            // host.
            this->numExportPacketsPerLID_.template sync<CHMS> ();
            this->numImportPacketsPerLID_.template sync<CHMS> ();

            // NOTE (mfh 25 Apr 2016, 01 Aug 2017) doPostsAndWaits and
            // doReversePostsAndWaits currently want
            // numExportPacketsPerLID and numImportPacketsPerLID as
            // Teuchos::ArrayView, rather than as Kokkos::View.
            auto numExportPacketsPerLID_av =
              getArrayViewFromDualView (this->numExportPacketsPerLID_);
            auto numImportPacketsPerLID_av =
              getArrayViewFromDualView (this->numImportPacketsPerLID_);

            // imports_ is for output only, so we don't need to sync
            // it before marking it as modified.  However, in order to
            // prevent spurious debug-mode errors (e.g., "modified on
            // both device and host"), we first need to clear its
            // "modified" flags.
            this->imports_.modified_device() = 0;
            this->imports_.modified_host() = 0;

            if (commOnHost) {
              this->imports_.template modify<CHMS> ();
              distor.doPostsAndWaits (create_const_view (this->exports_.template view<CHMS> ()),
                                      numExportPacketsPerLID_av,
                                      this->imports_.template view<CHMS> (),
                                      numImportPacketsPerLID_av);
            }
            else { // pack on device
              this->imports_.template modify<CDMS> ();
              distor.doPostsAndWaits (create_const_view (this->exports_.template view<CDMS> ()),
                                      numExportPacketsPerLID_av,
                                      this->imports_.template view<CDMS> (),
                                      numImportPacketsPerLID_av);
            }
          }
          else { // constant number of packets per LID
            if (verbose) {
              std::ostringstream os;
              os << *prefix << "7.1. Const # packets per LID: "
                 << "exports_.dimension_0()=" << exports_.dimension_0 ()
                 << ", imports_.dimension_0() = " << imports_.dimension_0 ()
                 << endl;
              std::cerr << os.str ();
            }
            // imports_ is for output only, so we don't need to sync
            // it before marking it as modified.  However, in order to
            // prevent spurious debug-mode errors (e.g., "modified on
            // both device and host"), we first need to clear its
            // "modified" flags.
            this->imports_.modified_device() = 0;
            this->imports_.modified_host() = 0;

            if (commOnHost) {
              if (verbose) {
                std::ostringstream os;
                os << *prefix << "7.2. Comm buffers on host" << endl;
                std::cerr << os.str ();
              }
              this->imports_.template modify<CHMS> ();
              distor.doPostsAndWaits (create_const_view (this->exports_.template view<CHMS> ()),
                                      constantNumPackets,
                                      this->imports_.template view<CHMS> ());
            }
            else { // pack on device
              if (verbose) {
                std::ostringstream os;
                os << *prefix << "7.2. Comm buffers on device" << endl;
                std::cerr << os.str ();
              }
              this->imports_.template modify<CDMS> ();
              distor.doPostsAndWaits (create_const_view (this->exports_.template view<CDMS> ()),
                                      constantNumPackets,
                                      this->imports_.template view<CDMS> ());
            }
          }
        }

        {
          if (verbose) {
            std::ostringstream os;
            os << *prefix << "8. unpackAndCombineNew" << endl;
            std::cerr << os.str ();
          }
#ifdef HAVE_TPETRA_TRANSFER_TIMERS
          Teuchos::TimeMonitor unpackAndCombineMon (*unpackAndCombineTimer_);
#endif // HAVE_TPETRA_TRANSFER_TIMERS
          // NOTE (mfh 26 Apr 2016) We don't actually need to sync the
          // input DualViews, but they DO need to be most recently
          // updated in the same memory space.
          //
          // FIXME (mfh 26 Apr 2016) Check that all input DualViews
          // were most recently updated in the same memory space, and
          // sync them to the same place (based on commOnHost) if not.
          this->unpackAndCombineNew (remoteLIDs, this->imports_,
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

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  createViews () const
  {}

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  createViewsNonConst (KokkosClassic::ReadWriteOption /*rwo*/)
  {}

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>::
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
