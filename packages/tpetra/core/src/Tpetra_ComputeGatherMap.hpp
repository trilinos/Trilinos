// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Tpetra_ComputeGatherMap_hpp
#define __Tpetra_ComputeGatherMap_hpp

/// \file Tpetra_ComputeGatherMap.hpp
/// \brief From a distributed map build a map with all GIDs on the root node
/// \author Mark Hoemmen
///
#include "Tpetra_Map.hpp"
#include <numeric>

// Macro that marks a function as "possibly unused," in order to
// suppress build warnings.
#if ! defined(TRILINOS_UNUSED_FUNCTION)
#  if defined(__GNUC__) || (defined(__INTEL_COMPILER) && !defined(_MSC_VER))
#    define TRILINOS_UNUSED_FUNCTION __attribute__((__unused__))
#  elif defined(__clang__)
#    if __has_attribute(unused)
#      define TRILINOS_UNUSED_FUNCTION __attribute__((__unused__))
#    else
#      define TRILINOS_UNUSED_FUNCTION
#    endif // Clang has 'unused' attribute
#  elif defined(__IBMCPP__)
// IBM's C++ compiler for Blue Gene/Q (V12.1) implements 'used' but not 'unused'.
//
// http://pic.dhe.ibm.com/infocenter/compbg/v121v141/index.jsp
#    define TRILINOS_UNUSED_FUNCTION
#  else // some other compiler
#    define TRILINOS_UNUSED_FUNCTION
#  endif
#endif // ! defined(TRILINOS_UNUSED_FUNCTION)


namespace Tpetra {
  namespace Details {

    namespace {
#ifdef HAVE_MPI
      // We're communicating integers of type IntType.  Figure
      // out the right MPI_Datatype for IntType.  Usually it
      // is int or long, so these are good enough for now.
      template<class IntType> TRILINOS_UNUSED_FUNCTION MPI_Datatype
      getMpiDatatype () {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Not implemented for IntType != int, long, or long long");
      }

      template<> TRILINOS_UNUSED_FUNCTION MPI_Datatype
      getMpiDatatype<int> () { return MPI_INT; }

      template<> TRILINOS_UNUSED_FUNCTION MPI_Datatype
      getMpiDatatype<long> () { return MPI_LONG; }

      template<> TRILINOS_UNUSED_FUNCTION MPI_Datatype
      getMpiDatatype<long long> () { return MPI_LONG_LONG; }
#endif // HAVE_MPI

      template<class IntType>
      void
      gather (const IntType sendbuf[], const int sendcnt,
              IntType recvbuf[], const int recvcnt,
              const int root, const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
      {
        using Teuchos::RCP;
        using Teuchos::rcp_dynamic_cast;
#ifdef HAVE_MPI
        using Teuchos::MpiComm;

        // Get the raw MPI communicator, so we don't have to wrap MPI_Gather in Teuchos.
        RCP<const MpiComm<int> > mpiComm = rcp_dynamic_cast<const MpiComm<int> > (comm);
        if (! mpiComm.is_null ()) {
          MPI_Datatype sendtype = getMpiDatatype<IntType> ();
          MPI_Datatype recvtype = sendtype;
          MPI_Comm rawMpiComm = * (mpiComm->getRawMpiComm ());
          const int err = MPI_Gather (reinterpret_cast<void*> (const_cast<IntType*> (sendbuf)),
                                      sendcnt,
                                      sendtype,
                                      reinterpret_cast<void*> (recvbuf),
                                      recvcnt,
                                      recvtype,
                                      root,
                                      rawMpiComm);
          TEUCHOS_TEST_FOR_EXCEPTION(
            err != MPI_SUCCESS, std::runtime_error, "MPI_Gather failed");
          return;
        }
#endif // HAVE_MPI

        TEUCHOS_TEST_FOR_EXCEPTION(
          recvcnt > sendcnt, std::invalid_argument,
          "gather: If the input communicator contains only one process, "
          "then you cannot receive more entries than you send.  "
          "You aim to receive " << recvcnt << " entries, "
          "but to send " << sendcnt << ".");
        // Serial communicator case: just copy.  recvcnt is the
        // amount to receive, so it's the amount to copy.
        std::copy (sendbuf, sendbuf + recvcnt, recvbuf);
      }

      template<class IntType>
      void
      gatherv (const IntType sendbuf[], const int sendcnt,
               IntType recvbuf[], const int recvcnts[], const int displs[],
               const int root, const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
      {
        using Teuchos::RCP;
        using Teuchos::rcp_dynamic_cast;
#ifdef HAVE_MPI
        using Teuchos::MpiComm;

        // Get the raw MPI communicator, so we don't have to wrap
        // MPI_Gather in Teuchos.
        RCP<const MpiComm<int> > mpiComm =
          rcp_dynamic_cast<const MpiComm<int> > (comm);
        if (! mpiComm.is_null ()) {
          MPI_Datatype sendtype = getMpiDatatype<IntType> ();
          MPI_Datatype recvtype = sendtype;
          MPI_Comm rawMpiComm = * (mpiComm->getRawMpiComm ());
          const int err = MPI_Gatherv (reinterpret_cast<void*> (const_cast<IntType*> (sendbuf)),
                                       sendcnt,
                                       sendtype,
                                       reinterpret_cast<void*> (recvbuf),
                                       const_cast<int*> (recvcnts),
                                       const_cast<int*> (displs),
                                       recvtype,
                                       root,
                                       rawMpiComm);
          TEUCHOS_TEST_FOR_EXCEPTION(
            err != MPI_SUCCESS, std::runtime_error, "MPI_Gatherv failed");
          return;
        }
#endif // HAVE_MPI
        TEUCHOS_TEST_FOR_EXCEPTION(
          recvcnts[0] > sendcnt,
          std::invalid_argument,
          "gatherv: If the input communicator contains only one process, "
          "then you cannot receive more entries than you send.  "
          "You aim to receive " << recvcnts[0] << " entries, "
          "but to send " << sendcnt << ".");
        // Serial communicator case: just copy.  recvcnts[0] is the
        // amount to receive, so it's the amount to copy.  Start
        // writing to recvbuf at the offset displs[0].
        std::copy (sendbuf, sendbuf + recvcnts[0], recvbuf + displs[0]);
      }
    } // namespace (anonymous)


    // Given an arbitrary Map, compute a Map containing all the GIDs
    // in the same order as in (the one-to-one version of) map, but
    // all owned exclusively by Proc 0.
    template<class MapType>
    Teuchos::RCP<const MapType>
    computeGatherMap (Teuchos::RCP<const MapType> map,
                      const Teuchos::RCP<Teuchos::FancyOStream>& err,
                      const bool dbg=false)
    {
      using Tpetra::createOneToOne;
      using Tpetra::global_size_t;
      using Teuchos::Array;
      using Teuchos::ArrayRCP;
      using Teuchos::ArrayView;
      using Teuchos::as;
      using Teuchos::Comm;
      using Teuchos::RCP;
      using std::endl;
      typedef typename MapType::local_ordinal_type LO;
      typedef typename MapType::global_ordinal_type GO;
      typedef typename MapType::node_type NT;

      RCP<const Comm<int> > comm = map->getComm ();
      const int numProcs = comm->getSize ();
      const int myRank = comm->getRank ();

      if (! err.is_null ()) {
        err->pushTab ();
      }
      if (dbg) {
        *err << myRank << ": computeGatherMap:" << endl;
      }
      if (! err.is_null ()) {
        err->pushTab ();
      }

      RCP<const MapType> oneToOneMap;
      if (map->isContiguous ()) {
        oneToOneMap = map; // contiguous Maps are always 1-to-1
      } else {
        if (dbg) {
          *err << myRank << ": computeGatherMap: Calling createOneToOne" << endl;
        }
        // It could be that Map is one-to-one, but the class doesn't
        // give us a way to test this, other than to create the
        // one-to-one Map.
        oneToOneMap = createOneToOne<LO, GO, NT> (map);
      }

      RCP<const MapType> gatherMap;
      if (numProcs == 1) {
        gatherMap = oneToOneMap;
      } else {
        if (dbg) {
          *err << myRank << ": computeGatherMap: Gathering Map counts" << endl;
        }
        // Gather each process' count of Map elements to Proc 0,
        // into the recvCounts array.  This will tell Proc 0 how
        // many GIDs to expect from each process when calling
        // MPI_Gatherv.  Counts and offsets are all int, because
        // that's what MPI uses.  Teuchos::as will at least prevent
        // bad casts to int in a debug build.
        //
        // Yeah, it's not memory scalable.  Guess what: We're going
        // to bring ALL the data to Proc 0, not just the receive
        // counts.  The first MPI_Gather is only the beginning...
        // Seriously, if you want to make this memory scalable, the
        // right thing to do (after the Export to the one-to-one
        // Map) is to go round robin through the processes, having
        // each send a chunk of data (with its GIDs, to get the
        // order of rows right) at a time, and overlapping writing
        // to the file (resp. reading from it) with receiving (resp.
        // sending) the next chunk.
        const int myEltCount = as<int> (oneToOneMap->getLocalNumElements ());
        Array<int> recvCounts (numProcs);
        const int rootProc = 0;
        gather<int> (&myEltCount, 1, recvCounts.getRawPtr (), 1, rootProc, comm);

        ArrayView<const GO> myGlobalElts = oneToOneMap->getLocalElementList ();
        const int numMyGlobalElts = as<int> (myGlobalElts.size ());
        // Only Proc 0 needs to receive and store all the GIDs (from
        // all processes).
        ArrayRCP<GO> allGlobalElts;
        if (myRank == 0) {
          allGlobalElts = Teuchos::arcp<GO> (oneToOneMap->getGlobalNumElements ());
          std::fill (allGlobalElts.begin (), allGlobalElts.end (), 0);
        }

        if (dbg) {
          *err << myRank << ": computeGatherMap: Computing MPI_Gatherv "
            "displacements" << endl;
        }
        // Displacements for gatherv() in this case (where we are
        // gathering into a contiguous array) are an exclusive
        // partial sum (first entry is zero, second starts the
        // partial sum) of recvCounts.
        Array<int> displs (numProcs, 0);
        std::partial_sum (recvCounts.begin (), recvCounts.end () - 1,
                          displs.begin () + 1);
        if (dbg) {
          *err << myRank << ": computeGatherMap: Calling MPI_Gatherv" << endl;
        }
        gatherv<GO> (myGlobalElts.getRawPtr (), numMyGlobalElts,
                     allGlobalElts.getRawPtr (), recvCounts.getRawPtr (),
                     displs.getRawPtr (), rootProc, comm);

        if (dbg) {
          *err << myRank << ": computeGatherMap: Creating gather Map" << endl;
        }
        // Create a Map with all the GIDs, in the same order as in
        // the one-to-one Map, but owned by Proc 0.
        ArrayView<const GO> allElts (NULL, 0);
        if (myRank == 0) {
          allElts = allGlobalElts ();
        }
        const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid ();
        gatherMap = rcp (new MapType (INVALID, allElts,
                                      oneToOneMap->getIndexBase (),
                                      comm));
      }
      if (! err.is_null ()) {
        err->popTab ();
      }
      if (dbg) {
        *err << myRank << ": computeGatherMap: done" << endl;
      }
      if (! err.is_null ()) {
        err->popTab ();
      }
      return gatherMap;
    }

  } // namespace Details
} // namespace Tpetra

#endif // __Tpetra_ComputeGatherMap_hpp
