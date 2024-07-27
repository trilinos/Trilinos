// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DIRECTORYIMPL_DEF_HPP
#define TPETRA_DIRECTORYIMPL_DEF_HPP

/// \file Tpetra_DirectoryImpl_def.hpp
/// \brief Definition of implementation details of Tpetra::Directory.

#include "Tpetra_Distributor.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_TieBreak.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_Details_FixedHashTable.hpp"
#include "Teuchos_Comm.hpp"
#include <memory>
#include <sstream>

// FIXME (mfh 16 Apr 2013) GIANT HACK BELOW
#ifdef HAVE_TPETRACORE_MPI
#  include <mpi.h>
#endif // HAVE_TPETRACORE_MPI
// FIXME (mfh 16 Apr 2013) GIANT HACK ABOVE

namespace Tpetra {
  namespace Details {
    template<class LO, class GO, class NT>
    LookupStatus
    Directory<LO, GO, NT>::
    getEntries (const map_type& map,
                const Teuchos::ArrayView<const GO> &globalIDs,
                const Teuchos::ArrayView<int> &nodeIDs,
                const Teuchos::ArrayView<LO> &localIDs,
                const bool computeLIDs) const
    {
      // Ensure that globalIDs, nodeIDs, and localIDs (if applicable)
      // all have the same size, before modifying any output arguments.
      TEUCHOS_TEST_FOR_EXCEPTION(nodeIDs.size() != globalIDs.size(),
        std::invalid_argument, Teuchos::typeName(*this) << "::getEntries(): "
        "Output arrays do not have the right sizes.  nodeIDs.size() = "
        << nodeIDs.size() << " != globalIDs.size() = " << globalIDs.size()
        << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(
        computeLIDs && localIDs.size() != globalIDs.size(),
        std::invalid_argument, Teuchos::typeName(*this) << "::getEntries(): "
        "Output array do not have the right sizes.  localIDs.size() = "
        << localIDs.size() << " != globalIDs.size() = " << globalIDs.size()
        << ".");

      // Initially, fill nodeIDs and localIDs (if applicable) with
      // invalid values.  The "invalid" process ID is -1 (this means
      // the same thing as MPI_ANY_SOURCE to Teuchos, so it's an
      // "invalid" process ID); the invalid local ID comes from
      // OrdinalTraits.
      std::fill (nodeIDs.begin(), nodeIDs.end(), -1);
      if (computeLIDs) {
        std::fill (localIDs.begin(), localIDs.end(),
                   Teuchos::OrdinalTraits<LO>::invalid ());
      }
      // Actually do the work.
      return this->getEntriesImpl (map, globalIDs, nodeIDs, localIDs, computeLIDs);
    }


    template<class LO, class GO, class NT>
    ReplicatedDirectory<LO, GO, NT>::
    ReplicatedDirectory (const map_type& map) :
      numProcs_ (map.getComm ()->getSize ())
    {}


    template<class LO, class GO, class NT>
    bool
    ReplicatedDirectory<LO, GO, NT>::
    isOneToOne (const Teuchos::Comm<int>& /* comm */) const
    {
      // A locally replicated Map is one-to-one only if there is no
      // replication, that is, only if the Map's communicator only has
      // one process.
      return (numProcs_ == 1);
    }


    template<class LO, class GO, class NT>
    std::string
    ReplicatedDirectory<LO, GO, NT>::description () const
    {
      std::ostringstream os;
      os << "ReplicatedDirectory"
         << "<" << Teuchos::TypeNameTraits<LO>::name ()
         << ", " << Teuchos::TypeNameTraits<GO>::name ()
         << ", " << Teuchos::TypeNameTraits<NT>::name () << ">";
      return os.str ();
    }


    template<class LO, class GO, class NT>
    ContiguousUniformDirectory<LO, GO, NT>::
    ContiguousUniformDirectory (const map_type& map)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(! map.isContiguous (), std::invalid_argument,
        Teuchos::typeName (*this) << " constructor: Map is not contiguous.");
      TEUCHOS_TEST_FOR_EXCEPTION(! map.isUniform (), std::invalid_argument,
        Teuchos::typeName (*this) << " constructor: Map is not uniform.");
    }


    template<class LO, class GO, class NT>
    std::string
    ContiguousUniformDirectory<LO, GO, NT>::description () const
    {
      std::ostringstream os;
      os << "ContiguousUniformDirectory"
         << "<" << Teuchos::TypeNameTraits<LO>::name ()
         << ", " << Teuchos::TypeNameTraits<GO>::name ()
         << ", " << Teuchos::TypeNameTraits<NT>::name () << ">";
      return os.str ();
    }


    template<class LO, class GO, class NT>
    LookupStatus
    ContiguousUniformDirectory<LO, GO, NT>::
    getEntriesImpl (const map_type& map,
                    const Teuchos::ArrayView<const GO> &globalIDs,
                    const Teuchos::ArrayView<int> &nodeIDs,
                    const Teuchos::ArrayView<LO> &localIDs,
                    const bool computeLIDs) const
    {
      using Teuchos::Comm;
      using Teuchos::RCP;
      typedef typename Teuchos::ArrayView<const GO>::size_type size_type;
      const LO invalidLid = Teuchos::OrdinalTraits<LO>::invalid ();
      LookupStatus res = AllIDsPresent;

      RCP<const Comm<int> > comm = map.getComm ();
      const GO g_min = map.getMinAllGlobalIndex ();

      // Let N_G be the global number of elements in the Map,
      // and P be the number of processes in its communicator.
      // Then, N_G = P * N_L + R = R*(N_L + 1) + (P - R)*N_L.
      //
      // The first R processes own N_L+1 elements.
      // The remaining P-R processes own N_L elements.
      //
      // Let g be the current GID, g_min be the global minimum GID,
      // and g_0 = g - g_min.  If g is a valid GID in this Map, then
      // g_0 is in [0, N_G - 1].
      //
      // If g is a valid GID in this Map and g_0 < R*(N_L + 1), then
      // the rank of the process that owns g is floor(g_0 / (N_L +
      // 1)), and its corresponding local index on that process is g_0
      // mod (N_L + 1).
      //
      // Let g_R = g_0 - R*(N_L + 1).  If g is a valid GID in this Map
      // and g_0 >= R*(N_L + 1), then the rank of the process that
      // owns g is then R + floor(g_R / N_L), and its corresponding
      // local index on that process is g_R mod N_L.

      const size_type N_G =
        static_cast<size_type> (map.getGlobalNumElements ());
      const size_type P = static_cast<size_type> (comm->getSize ());
      const size_type N_L  = N_G / P;
      const size_type R = N_G - N_L * P; // N_G mod P
      const size_type N_R = R * (N_L + static_cast<size_type> (1));

#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(
        N_G != P*N_L + R, std::logic_error,
        "Tpetra::ContiguousUniformDirectory::getEntriesImpl: "
        "N_G = " << N_G << " != P*N_L + R = " << P << "*" << N_L << " + " << R
        << " = " << P*N_L + R << ".  "
        "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

      const size_type numGids = globalIDs.size (); // for const loop bound
      // Avoid signed/unsigned comparisons below, in case GO is
      // unsigned.  (Integer literals are generally signed.)
      const GO ONE = static_cast<GO> (1);

      if (computeLIDs) {
        for (size_type k = 0; k < numGids; ++k) {
          const GO g_0 = globalIDs[k] - g_min;

          // The first test is a little strange just in case GO is
          // unsigned.  Compilers raise a warning on tests like "x <
          // 0" if x is unsigned, but don't usually raise a warning if
          // the expression is a bit more complicated than that.
          if (g_0 + ONE < ONE || g_0 >= static_cast<GO> (N_G)) {
            nodeIDs[k] = -1;
            localIDs[k] = invalidLid;
            res = IDNotPresent;
          }
          else if (g_0 < static_cast<GO> (N_R)) {
            // The GID comes from the initial sequence of R processes.
            nodeIDs[k] = static_cast<int> (g_0 / static_cast<GO> (N_L + 1));
            localIDs[k] = static_cast<LO> (g_0 % static_cast<GO> (N_L + 1));
          }
          else if (g_0 >= static_cast<GO> (N_R)) {
            // The GID comes from the remaining P-R processes.
            const GO g_R = g_0 - static_cast<GO> (N_R);
            nodeIDs[k] = static_cast<int> (R + g_R / N_L);
            localIDs[k] = static_cast<int> (g_R % N_L);
          }
#ifdef HAVE_TPETRA_DEBUG
          else {
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
              "Tpetra::ContiguousUniformDirectory::getEntriesImpl: "
              "should never get here.  "
              "Please report this bug to the Tpetra developers.");
          }
#endif // HAVE_TPETRA_DEBUG
        }
      }
      else { // don't compute local indices
        for (size_type k = 0; k < numGids; ++k) {
          const GO g_0 = globalIDs[k] - g_min;
          // The first test is a little strange just in case GO is
          // unsigned.  Compilers raise a warning on tests like "x <
          // 0" if x is unsigned, but don't usually raise a warning if
          // the expression is a bit more complicated than that.
          if (g_0 + ONE < ONE || g_0 >= static_cast<GO> (N_G)) {
            nodeIDs[k] = -1;
            res = IDNotPresent;
          }
          else if (g_0 < static_cast<GO> (N_R)) {
            // The GID comes from the initial sequence of R processes.
            nodeIDs[k] = static_cast<int> (g_0 / static_cast<GO> (N_L + 1));
          }
          else if (g_0 >= static_cast<GO> (N_R)) {
            // The GID comes from the remaining P-R processes.
            const GO g_R = g_0 - static_cast<GO> (N_R);
            nodeIDs[k] = static_cast<int> (R + g_R / N_L);
          }
#ifdef HAVE_TPETRA_DEBUG
          else {
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
              "Tpetra::ContiguousUniformDirectory::getEntriesImpl: "
              "should never get here.  "
              "Please report this bug to the Tpetra developers.");
          }
#endif // HAVE_TPETRA_DEBUG
        }
      }
      return res;
    }

    template<class LO, class GO, class NT>
    DistributedContiguousDirectory<LO, GO, NT>::
    DistributedContiguousDirectory (const map_type& map)
    {
      using Teuchos::arcp;
      using Teuchos::gatherAll;
      using Teuchos::RCP;

      RCP<const Teuchos::Comm<int> > comm = map.getComm ();

      TEUCHOS_TEST_FOR_EXCEPTION(! map.isDistributed (), std::invalid_argument,
        Teuchos::typeName (*this) << " constructor: Map is not distributed.");
      TEUCHOS_TEST_FOR_EXCEPTION(! map.isContiguous (), std::invalid_argument,
        Teuchos::typeName (*this) << " constructor: Map is not contiguous.");

      const int numProcs = comm->getSize ();

      // Make room for the min global ID on each process, plus one
      // entry at the end for the "max cap."
      allMinGIDs_ = arcp<GO> (numProcs + 1);
      // Get my process' min global ID.
      GO minMyGID = map.getMinGlobalIndex ();
      // Gather all of the min global IDs into the first numProcs
      // entries of allMinGIDs_.

      // FIXME (mfh 16 Apr 2013) GIANT HACK BELOW
      //
      // The purpose of this giant hack is that gatherAll appears to
      // interpret the "receive count" argument differently than
      // MPI_Allgather does.  Matt Bettencourt reports Valgrind issues
      // (memcpy with overlapping data) with MpiComm<int>::gatherAll,
      // which could relate either to this, or to OpenMPI.
#ifdef HAVE_TPETRACORE_MPI
      MPI_Datatype rawMpiType = MPI_INT;
      bool useRawMpi = true;
      if (typeid (GO) == typeid (int)) {
        rawMpiType = MPI_INT;
      } else if (typeid (GO) == typeid (long)) {
        rawMpiType = MPI_LONG;
      } else {
        useRawMpi = false;
      }
      if (useRawMpi) {
        using Teuchos::rcp_dynamic_cast;
        using Teuchos::MpiComm;
        RCP<const MpiComm<int> > mpiComm =
          rcp_dynamic_cast<const MpiComm<int> > (comm);
        // It could be a SerialComm instead, even in an MPI build, so
        // be sure to check.
        if (! comm.is_null ()) {
          MPI_Comm rawMpiComm = * (mpiComm->getRawMpiComm ());
          const int err =
            MPI_Allgather (&minMyGID, 1, rawMpiType,
                           allMinGIDs_.getRawPtr (), 1, rawMpiType,
                           rawMpiComm);
          TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
            "Tpetra::DistributedContiguousDirectory: MPI_Allgather failed");
        } else {
          gatherAll<int, GO> (*comm, 1, &minMyGID, numProcs, allMinGIDs_.getRawPtr ());
        }
      } else {
        gatherAll<int, GO> (*comm, 1, &minMyGID, numProcs, allMinGIDs_.getRawPtr ());
      }
#else // NOT HAVE_TPETRACORE_MPI
      gatherAll<int, GO> (*comm, 1, &minMyGID, numProcs, allMinGIDs_.getRawPtr ());
#endif // HAVE_TPETRACORE_MPI
      // FIXME (mfh 16 Apr 2013) GIANT HACK ABOVE

      //gatherAll<int, GO> (*comm, 1, &minMyGID, numProcs, allMinGIDs_.getRawPtr ());

      // Put the max cap at the end.  Adding one lets us write loops
      // over the global IDs with the usual strict less-than bound.
      allMinGIDs_[numProcs] = map.getMaxAllGlobalIndex ()
                            + Teuchos::OrdinalTraits<GO>::one ();
    }

    template<class LO, class GO, class NT>
    std::string
    DistributedContiguousDirectory<LO, GO, NT>::description () const
    {
      std::ostringstream os;
      os << "DistributedContiguousDirectory"
         << "<" << Teuchos::TypeNameTraits<LO>::name ()
         << ", " << Teuchos::TypeNameTraits<GO>::name ()
         << ", " << Teuchos::TypeNameTraits<NT>::name () << ">";
      return os.str ();
    }

    template<class LO, class GO, class NT>
    LookupStatus
    ReplicatedDirectory<LO, GO, NT>::
    getEntriesImpl (const map_type& map,
                    const Teuchos::ArrayView<const GO> &globalIDs,
                    const Teuchos::ArrayView<int> &nodeIDs,
                    const Teuchos::ArrayView<LO> &localIDs,
                    const bool computeLIDs) const
    {
      using Teuchos::Array;
      using Teuchos::ArrayRCP;
      using Teuchos::ArrayView;
      using Teuchos::as;
      using Teuchos::Comm;
      using Teuchos::RCP;

      LookupStatus res = AllIDsPresent;
      RCP<const Teuchos::Comm<int> > comm = map.getComm ();
      const int myRank = comm->getRank ();

      // Map is on one process or is locally replicated.
      typename ArrayView<int>::iterator procIter = nodeIDs.begin();
      typename ArrayView<LO>::iterator lidIter = localIDs.begin();
      typename ArrayView<const GO>::iterator gidIter;
      for (gidIter = globalIDs.begin(); gidIter != globalIDs.end(); ++gidIter) {
        if (map.isNodeGlobalElement (*gidIter)) {
          *procIter++ = myRank;
          if (computeLIDs) {
            *lidIter++ = map.getLocalElement (*gidIter);
          }
        }
        else {
          // Advance the pointers, leaving these values set to invalid
          procIter++;
          if (computeLIDs) {
            lidIter++;
          }
          res = IDNotPresent;
        }
      }
      return res;
    }

    template<class LO, class GO, class NT>
    LookupStatus
    DistributedContiguousDirectory<LO, GO, NT>::
    getEntriesImpl (const map_type& map,
                    const Teuchos::ArrayView<const GO> &globalIDs,
                    const Teuchos::ArrayView<int> &nodeIDs,
                    const Teuchos::ArrayView<LO> &localIDs,
                    const bool computeLIDs) const
    {
      using Teuchos::Array;
      using Teuchos::ArrayRCP;
      using Teuchos::ArrayView;
      using Teuchos::as;
      using Teuchos::Comm;
      using Teuchos::RCP;

      RCP<const Teuchos::Comm<int> > comm = map.getComm ();
      const int numProcs = comm->getSize ();
      const LO LINVALID = Teuchos::OrdinalTraits<LO>::invalid();
      const GO noGIDsOnProc = std::numeric_limits<GO>::max();
      LookupStatus res = AllIDsPresent;

      // Find the first initialized GID for search below
      int firstProcWithGIDs;
      for (firstProcWithGIDs = 0; firstProcWithGIDs < numProcs;
           firstProcWithGIDs++) {
        if (allMinGIDs_[firstProcWithGIDs] != noGIDsOnProc) break;
      }

      // If Map is empty, return invalid values for all requested lookups
      // This case should not happen, as an empty Map is not considered
      // Distributed.
      if (firstProcWithGIDs == numProcs) {
        // No entries in Map
        res = (globalIDs.size() > 0) ? IDNotPresent : AllIDsPresent;
        std::fill(nodeIDs.begin(), nodeIDs.end(), -1);
        if (computeLIDs)
          std::fill(localIDs.begin(), localIDs.end(), LINVALID);
        return res;
      }

      const GO one = as<GO> (1);
      const GO nOverP = as<GO> (map.getGlobalNumElements ()
                      / as<global_size_t>(numProcs - firstProcWithGIDs));

      // Map is distributed but contiguous.
      typename ArrayView<int>::iterator procIter = nodeIDs.begin();
      typename ArrayView<LO>::iterator lidIter = localIDs.begin();
      typename ArrayView<const GO>::iterator gidIter;
      for (gidIter = globalIDs.begin(); gidIter != globalIDs.end(); ++gidIter) {
        LO LID = LINVALID; // Assume not found until proven otherwise
        int image = -1;
        GO GID = *gidIter;
        // Guess uniform distribution (TODO: maybe replace by a binary search)
        // We go through all this trouble to avoid overflow and
        // signed / unsigned casting mistakes (that were made in
        // previous versions of this code).
        int curRank;
        const GO firstGuess = firstProcWithGIDs + GID / std::max(nOverP, one);
        curRank = as<int>(std::min(firstGuess, as<GO>(numProcs - 1)));

        // This while loop will stop because 
        // allMinGIDs_[np] == global num elements
        while (allMinGIDs_[curRank] == noGIDsOnProc) curRank++;

        bool badGID = false;
        while (curRank >= firstProcWithGIDs && GID < allMinGIDs_[curRank]) {
          curRank--;
          while (curRank >= firstProcWithGIDs && 
                 allMinGIDs_[curRank] == noGIDsOnProc) curRank--;
        }
        if (curRank < firstProcWithGIDs) {
          // GID is lower than all GIDs in map
          badGID = true;
        }
        else if (curRank == numProcs) {
          // GID is higher than all GIDs in map
          badGID = true;
        }
        else {
          // we have the lower bound; now limit from above
          int above = curRank + 1;
          while (allMinGIDs_[above] == noGIDsOnProc) above++;

          while (GID >= allMinGIDs_[above]) {
            curRank = above;
            if (curRank == numProcs) {
              // GID is higher than all GIDs in map
              badGID = true;
              break;
            }
            above++;
            while (allMinGIDs_[above] == noGIDsOnProc) above++;
          }
        }

        if (!badGID) {
          image = curRank;
          LID = as<LO> (GID - allMinGIDs_[image]);
        }
        else {
          res = IDNotPresent;
        }
        *procIter++ = image;
        if (computeLIDs) {
          *lidIter++ = LID;
        }
      }
      return res;
    }

    template<class LO, class GO, class NT>
    DistributedNoncontiguousDirectory<LO, GO, NT>::
    DistributedNoncontiguousDirectory (const map_type& map) :
      oneToOneResult_ (ONE_TO_ONE_NOT_CALLED_YET), // to be revised below
      locallyOneToOne_ (true), // to be revised below
      useHashTables_ (false) // to be revised below
    {
      initialize (map, Teuchos::null);
    }

    template<class LO, class GO, class NT>
    DistributedNoncontiguousDirectory<LO, GO, NT>::
    DistributedNoncontiguousDirectory (const map_type& map,
                                       const tie_break_type& tie_break) :
      oneToOneResult_ (ONE_TO_ONE_NOT_CALLED_YET), // to be revised below
      locallyOneToOne_ (true), // to be revised below
      useHashTables_ (false) // to be revised below
    {
      initialize (map, Teuchos::ptrFromRef (tie_break));
    }

    template<class LO, class GO, class NT>
    void
    DistributedNoncontiguousDirectory<LO, GO, NT>::
    initialize (const map_type& map,
                Teuchos::Ptr<const tie_break_type> tie_break)
    {
      using Teuchos::arcp;
      using Teuchos::Array;
      using Teuchos::ArrayRCP;
      using Teuchos::ArrayView;
      using Teuchos::as;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::typeName;
      using Teuchos::TypeNameTraits;
      using std::cerr;
      using std::endl;
      typedef Array<int>::size_type size_type;

      // This class' implementation of getEntriesImpl() currently
      // encodes the following assumptions:
      //
      // 1. global_size_t >= GO
      // 2. global_size_t >= int
      // 3. global_size_t >= LO
      //
      // We check these assumptions here.
      TEUCHOS_TEST_FOR_EXCEPTION(sizeof(global_size_t) < sizeof(GO),
        std::logic_error, typeName (*this) << ": sizeof(Tpetra::"
        "global_size_t) = " << sizeof(global_size_t) << " < sizeof(Global"
        "Ordinal = " << TypeNameTraits<LO>::name () << ") = " << sizeof(GO)
        << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(sizeof(global_size_t) < sizeof(int),
        std::logic_error, typeName (*this) << ": sizeof(Tpetra::"
        "global_size_t) = " << sizeof(global_size_t) << " < sizeof(int) = "
        << sizeof(int) << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(sizeof(global_size_t) < sizeof(LO),
        std::logic_error, typeName (*this) << ": sizeof(Tpetra::"
        "global_size_t) = " << sizeof(global_size_t) << " < sizeof(Local"
        "Ordinal = " << TypeNameTraits<LO>::name () << ") = " << sizeof(LO)
        << ".");

      RCP<const Teuchos::Comm<int> > comm = map.getComm ();
      const LO LINVALID = Teuchos::OrdinalTraits<LO>::invalid ();
      const GO minAllGID = map.getMinAllGlobalIndex ();
      const GO maxAllGID = map.getMaxAllGlobalIndex ();

      // The "Directory Map" (see below) will have a range of elements
      // from the minimum to the maximum GID of the user Map, and a
      // minimum GID of minAllGID from the user Map.  It doesn't
      // actually have to store all those entries, though do beware of
      // calling getLocalElementList on it (see Bug 5822).
      const global_size_t numGlobalEntries = maxAllGID - minAllGID + 1;

      // We can't afford to replicate the whole directory on each
      // process, so create the "Directory Map", a uniform contiguous
      // Map that describes how we will distribute the directory over
      // processes.
      //
      // FIXME (mfh 08 May 2012) Here we're setting minAllGID to be
      // the index base.  The index base should be separate from the
      // minimum GID.
      directoryMap_ = rcp (new map_type (numGlobalEntries, minAllGID, comm,
                                         GloballyDistributed));
      // The number of Directory elements that my process owns.
      const size_t dir_numMyEntries = directoryMap_->getLocalNumElements ();

      // Fix for Bug 5822: If the input Map is "sparse," that is if
      // the difference between the global min and global max GID is
      // much larger than the global number of elements in the input
      // Map, then it's possible that the Directory Map might have
      // many more entries than the input Map on this process.  This
      // can cause memory scalability issues.  In that case, we switch
      // from the array-based implementation of Directory storage to
      // the hash table - based implementation.  We don't use hash
      // tables all the time, because they are slower in the common
      // case of a nonsparse Map.
      //
      // NOTE: This is a per-process decision.  Some processes may use
      // array-based storage, whereas others may use hash table -
      // based storage.

      // A hash table takes a constant factor more space, more or
      // less, than an array.  Thus, it's not worthwhile, even in
      // terms of memory usage, always to use a hash table.
      // Furthermore, array lookups are faster than hash table
      // lookups, so it may be worthwhile to use an array even if it
      // takes more space.  The "sparsity threshold" governs when to
      // switch to a hash table - based implementation.
      const size_t inverseSparsityThreshold = 10;
      useHashTables_ =
        (dir_numMyEntries >= inverseSparsityThreshold * map.getLocalNumElements());

      // Get list of process IDs that own the directory entries for the
      // Map GIDs.  These will be the targets of the sends that the
      // Distributor will do.
      const int myRank = comm->getRank ();
      const size_t numMyEntries = map.getLocalNumElements ();
      Array<int> sendImageIDs (numMyEntries);
      ArrayView<const GO> myGlobalEntries = map.getLocalElementList ();
      // An ID not present in this lookup indicates that it lies outside
      // of the range [minAllGID,maxAllGID] (from map_).  this means
      // something is wrong with map_, our fault.
      const LookupStatus lookupStatus =
        directoryMap_->getRemoteIndexList (myGlobalEntries, sendImageIDs);
      TEUCHOS_TEST_FOR_EXCEPTION(
        lookupStatus == IDNotPresent, std::logic_error, Teuchos::typeName(*this)
        << " constructor: the Directory Map could not find out where one or "
        "more of my Map's indices should go.  The input to getRemoteIndexList "
        "is " << Teuchos::toString (myGlobalEntries) << ", and the output is "
        << Teuchos::toString (sendImageIDs ()) << ".  The input Map itself has "
        "the following entries on the calling process " <<
        map.getComm ()->getRank () << ": " <<
        Teuchos::toString (map.getLocalElementList ()) << ", and has "
        << map.getGlobalNumElements () << " total global indices in ["
        << map.getMinAllGlobalIndex () << "," << map.getMaxAllGlobalIndex ()
        << "].  The Directory Map has "
        << directoryMap_->getGlobalNumElements () << " total global indices in "
        "[" << directoryMap_->getMinAllGlobalIndex () << "," <<
        directoryMap_->getMaxAllGlobalIndex () << "], and the calling process "
        "has GIDs [" << directoryMap_->getMinGlobalIndex () << "," <<
        directoryMap_->getMaxGlobalIndex () << "].  "
        "This probably means there is a bug in Map or Directory.  "
        "Please report this bug to the Tpetra developers.");

      // Initialize the distributor using the list of process IDs to
      // which to send.  We'll use the distributor to send out triples
      // of (GID, process ID, LID).  We're sending the entries to the
      // processes that the Directory Map says should own them, which is
      // why we called directoryMap_->getRemoteIndexList() above.
      Distributor distor (comm);
      const size_t numReceives = distor.createFromSends (sendImageIDs);

      // NOTE (mfh 21 Mar 2012) The following code assumes that
      // sizeof(GO) >= sizeof(int) and sizeof(GO) >= sizeof(LO).
      //
      // Create and fill buffer of (GID, PID, LID) triples to send
      // out.  We pack the (GID, PID, LID) triples into a single Array
      // of GO, casting the PID from int to GO and the LID from LO to
      // GO as we do so.
      //
      // FIXME (mfh 23 Mar 2014) This assumes that sizeof(LO) <=
      // sizeof(GO) and sizeof(int) <= sizeof(GO).  The former is
      // required, and the latter is generally the case, but we should
      // still check for this.
      const int packetSize = 3; // We're sending triples, so packet size is 3.
      Kokkos::View<GO*, Kokkos::HostSpace> exportEntries("exportEntries", packetSize * numMyEntries);
      {
        size_t exportIndex = 0;
        for (size_t i = 0; i < numMyEntries; ++i) {
          exportEntries[exportIndex++] = myGlobalEntries[i];
          exportEntries[exportIndex++] = as<GO> (myRank);
          exportEntries[exportIndex++] = as<GO> (i);
        }
      }
      // Buffer of data to receive.  The Distributor figured out for
      // us how many packets we're receiving, when we called its
      // createFromSends() method to set up the distribution plan.
      Kokkos::View<GO*, Kokkos::HostSpace> importElements("importElements", packetSize * distor.getTotalReceiveLength());

      // Distribute the triples of (GID, process ID, LID).
      distor.doPostsAndWaits(exportEntries, packetSize, importElements);

      // Unpack the redistributed data.  Both implementations of
      // Directory storage map from an LID in the Directory Map (which
      // is the LID of the GID to store) to either a PID or an LID in
      // the input Map.  Each "packet" (contiguous chunk of
      // importElements) contains a triple: (GID, PID, LID).
      if (useHashTables_) {
        // Create the hash tables.  We know exactly how many elements
        // to expect in each hash table.  FixedHashTable's constructor
        // currently requires all the keys and values at once, so we
        // have to extract them in temporary arrays.  It may be
        // possible to rewrite FixedHashTable to use a "start fill" /
        // "end fill" approach that avoids the temporary arrays, but
        // we won't try that for now.

        // The constructors of Array and ArrayRCP that take a number
        // of elements all initialize the arrays.  Instead, allocate
        // raw arrays, then hand them off to ArrayRCP, to avoid the
        // initial unnecessary initialization without losing the
        // benefit of exception safety (and bounds checking, in a
        // debug build).
        LO* tableKeysRaw = NULL;
        LO* tableLidsRaw = NULL;
        int* tablePidsRaw = NULL;
        try {
          tableKeysRaw = new LO [numReceives];
          tableLidsRaw = new LO [numReceives];
          tablePidsRaw = new int [numReceives];
        } catch (...) {
          if (tableKeysRaw != NULL) {
            delete [] tableKeysRaw;
          }
          if (tableLidsRaw != NULL) {
            delete [] tableLidsRaw;
          }
          if (tablePidsRaw != NULL) {
            delete [] tablePidsRaw;
          }
          throw;
        }
        ArrayRCP<LO> tableKeys (tableKeysRaw, 0, numReceives, true);
        ArrayRCP<LO> tableLids (tableLidsRaw, 0, numReceives, true);
        ArrayRCP<int> tablePids (tablePidsRaw, 0, numReceives, true);

        if (tie_break.is_null ()) {
          // Fill the temporary arrays of keys and values.
          size_type importIndex = 0;
          for (size_type i = 0; i < static_cast<size_type> (numReceives); ++i) {
            const GO curGID = importElements[importIndex++];
            const LO curLID = directoryMap_->getLocalElement (curGID);
            TEUCHOS_TEST_FOR_EXCEPTION(
              curLID == LINVALID, std::logic_error,
              Teuchos::typeName(*this) << " constructor: Incoming global index "
              << curGID << " does not have a corresponding local index in the "
              "Directory Map.  Please report this bug to the Tpetra developers.");
            tableKeys[i] = curLID;
            tablePids[i] = importElements[importIndex++];
            tableLids[i] = importElements[importIndex++];
          }
          // Set up the hash tables.  The hash tables' constructor
          // detects whether there are duplicates, so that we can set
          // locallyOneToOne_.
          lidToPidTable_ =
            rcp (new lidToPidTable_type (tableKeys (), tablePids ()));
          locallyOneToOne_ = ! (lidToPidTable_->hasDuplicateKeys ());
          lidToLidTable_ =
            rcp (new lidToLidTable_type (tableKeys (), tableLids ()));
        }
        else { // tie_break is NOT null

          // For each directory Map LID received, collect all the
          // corresponding (PID,LID) pairs.  If the input Map is not
          // one-to-one, corresponding directory Map LIDs will have
          // more than one pair.  In that case, we will use the
          // TieBreak object to pick exactly one pair.
          typedef std::map<LO, std::vector<std::pair<int, LO> > > pair_table_type;
          pair_table_type ownedPidLidPairs;

          // For each directory Map LID received, collect the zero or
          // more input Map (PID,LID) pairs into ownedPidLidPairs.
          size_type importIndex = 0;
          for (size_type i = 0; i < static_cast<size_type> (numReceives); ++i) {
            const GO curGID = importElements[importIndex++];
            const LO dirMapLid = directoryMap_->getLocalElement (curGID);
            TEUCHOS_TEST_FOR_EXCEPTION(
              dirMapLid == LINVALID, std::logic_error,
              Teuchos::typeName(*this) << " constructor: Incoming global index "
              << curGID << " does not have a corresponding local index in the "
              "Directory Map.  Please report this bug to the Tpetra developers.");
            tableKeys[i] = dirMapLid;
            const int PID = importElements[importIndex++];
            const int LID = importElements[importIndex++];

            // These may change below.  We fill them in just to ensure
            // that they won't have invalid values.
            tablePids[i] = PID;
            tableLids[i] = LID;

            // For every directory Map LID, we have to remember all
            // (PID, LID) pairs.  The TieBreak object will arbitrate
            // between them in the loop below.
            ownedPidLidPairs[dirMapLid].push_back (std::make_pair (PID, LID));
          }

          // Use TieBreak to arbitrate between (PID,LID) pairs
          // corresponding to each directory Map LID.
          //
          // FIXME (mfh 23 Mar 2014) How do I know that i is the same
          // as the directory Map LID?
          // KDD 21 Mar 2018:  It isn't, especially if the user's IDs are not
          // contiguous, but the directory map is.  Need to set tablePids[i]
          // and tableLids[i], so need to loop over numReceives (as that is
          // how those arrays are allocated).  FIXED

          for (size_type i = 0; i < static_cast<size_type> (numReceives); ++i) {
            const LO dirMapLid = tableKeys[i];
            const std::vector<std::pair<int, LO> >& pidLidList =
                                                    ownedPidLidPairs[dirMapLid];
            const size_t listLen = pidLidList.size();
            if (listLen == 0) continue;  // KDD This will never happen
            const GO dirMapGid = directoryMap_->getGlobalElement (dirMapLid);
            if (listLen > 1) {
              locallyOneToOne_ = false;
            }
            // If there is some (PID,LID) pair for the current input
            // Map LID, then it makes sense to invoke the TieBreak
            // object to arbitrate between the options.  Even if
            // there is only one (PID,LID) pair, we still want to
            // give the TieBreak object a chance to do whatever it
            // likes to do, in terms of side effects (e.g., track
            // (PID,LID) pairs).
            const size_type index =
              static_cast<size_type> (tie_break->selectedIndex (dirMapGid,
                                                                pidLidList));
            tablePids[i] = pidLidList[index].first;
            tableLids[i] = pidLidList[index].second;
          }

          // Set up the hash tables.
          lidToPidTable_ =
            rcp (new lidToPidTable_type (tableKeys (), tablePids ()));
          lidToLidTable_ =
            rcp (new lidToLidTable_type (tableKeys (), tableLids ()));
        }
      }
      else {
        if (tie_break.is_null ()) {
          // Use array-based implementation of Directory storage.
          // Allocate these arrays and fill them with invalid values,
          // in case the input Map's GID list is sparse (i.e., does
          // not populate all GIDs from minAllGID to maxAllGID).
          PIDs_ = arcp<int> (dir_numMyEntries);
          std::fill (PIDs_.begin (), PIDs_.end (), -1);
          LIDs_ = arcp<LO> (dir_numMyEntries);
          std::fill (LIDs_.begin (), LIDs_.end (), LINVALID);
          // Fill in the arrays with PIDs resp. LIDs.
          size_type importIndex = 0;
          for (size_type i = 0; i < static_cast<size_type> (numReceives); ++i) {
            const GO curGID = importElements[importIndex++];
            const LO curLID = directoryMap_->getLocalElement (curGID);
            TEUCHOS_TEST_FOR_EXCEPTION(curLID == LINVALID, std::logic_error,
              Teuchos::typeName(*this) << " constructor: Incoming global index "
              << curGID << " does not have a corresponding local index in the "
              "Directory Map.  Please report this bug to the Tpetra developers.");

            // If PIDs_[curLID] is not -1, then curGID is a duplicate
            // on the calling process, so the Directory is not locally
            // one-to-one.
            if (PIDs_[curLID] != -1) {
              locallyOneToOne_ = false;
            }
            PIDs_[curLID] = importElements[importIndex++];
            LIDs_[curLID] = importElements[importIndex++];
          }
        }
        else {
          PIDs_ = arcp<int> (dir_numMyEntries);
          LIDs_ = arcp<LO> (dir_numMyEntries);
          std::fill (PIDs_.begin (), PIDs_.end (), -1);

          // All received (PID, LID) pairs go into ownedPidLidPairs.
          // This is a map from the directory Map's LID to the (PID,
          // LID) pair (where the latter LID comes from the input Map,
          // not the directory Map).  If the input Map is not
          // one-to-one, corresponding LIDs will have
          // ownedPidLidPairs[curLID].size() > 1.  In that case, we
          // will use the TieBreak object to pick exactly one pair.
          Array<std::vector<std::pair<int, LO> > > ownedPidLidPairs (dir_numMyEntries);
          size_t importIndex = 0;
          for (size_t i = 0; i < numReceives; ++i) {
            const GO  GID = importElements[importIndex++];
            const int PID = importElements[importIndex++];
            const LO  LID = importElements[importIndex++];

            const LO dirMapLid = directoryMap_->getLocalElement (GID);
            TEUCHOS_TEST_FOR_EXCEPTION(
              dirMapLid == LINVALID, std::logic_error,
              Teuchos::typeName(*this) << " constructor: Incoming global index "
              << GID << " does not have a corresponding local index in the "
              "Directory Map.  Please report this bug to the Tpetra developers.");
            ownedPidLidPairs[dirMapLid].push_back (std::make_pair (PID, LID));
          }

          // Use TieBreak to arbitrate between (PID,LID) pairs
          // corresponding to each directory Map LID.
          //
          // FIXME (mfh 23 Mar 2014) How do I know that i is the same
          // as the directory Map LID?
          // KDD 21 Mar 2018:  It isn't, especially if the user's IDs are not
          // contiguous.  Loop over all ownedPidLidPairs; skip those that have
          // empty lists.  FIXED

          for (size_t i = 0; i < dir_numMyEntries; ++i) {
            const std::vector<std::pair<int, LO> >& pidLidList =
              ownedPidLidPairs[i];
            const size_t listLen = pidLidList.size();
            if (listLen == 0) continue;  // KDD will happen for GIDs not in
                                         // KDD the user's source map
            const LO dirMapLid = static_cast<LO> (i);
            const GO dirMapGid = directoryMap_->getGlobalElement (dirMapLid);
            if (listLen > 1) {
              locallyOneToOne_ = false;
            }
            // If there is some (PID,LID) pair for the current input
            // Map LID, then it makes sense to invoke the TieBreak
            // object to arbitrate between the options.  Even if
            // there is only one (PID,LID) pair, we still want to
            // give the TieBreak object a chance to do whatever it
            // likes to do, in terms of side effects (e.g., track
            // (PID,LID) pairs).
            const size_type index =
              static_cast<size_type> (tie_break->selectedIndex (dirMapGid,
                                                                pidLidList));
            PIDs_[i] = pidLidList[index].first;
            LIDs_[i] = pidLidList[index].second;
          }
        }
      }
    }

    template<class LO, class GO, class NT>
    std::string
    DistributedNoncontiguousDirectory<LO, GO, NT>::description () const
    {
      std::ostringstream os;
      os << "DistributedNoncontiguousDirectory"
         << "<" << Teuchos::TypeNameTraits<LO>::name ()
         << ", " << Teuchos::TypeNameTraits<GO>::name ()
         << ", " << Teuchos::TypeNameTraits<NT>::name () << ">";
      return os.str ();
    }

    template<class LO, class GO, class NT>
    LookupStatus
    DistributedNoncontiguousDirectory<LO, GO, NT>::
    getEntriesImpl (const map_type& map,
                    const Teuchos::ArrayView<const GO> &globalIDs,
                    const Teuchos::ArrayView<int> &nodeIDs,
                    const Teuchos::ArrayView<LO> &localIDs,
                    const bool computeLIDs) const
    {
      using Teuchos::Array;
      using Teuchos::ArrayRCP;
      using Teuchos::ArrayView;
      using Teuchos::as;
      using Teuchos::RCP;
      using Teuchos::toString;
      using Details::Behavior;
      using Details::verbosePrintArray;
      using std::cerr;
      using std::endl;
      using size_type = typename Array<GO>::size_type;
      const char funcPrefix[] = "Tpetra::"
        "DistributedNoncontiguousDirectory::getEntriesImpl: ";
      const char errSuffix[] =
        "  Please report this bug to the Tpetra developers.";

      RCP<const Teuchos::Comm<int> > comm = map.getComm ();
      const bool verbose = Behavior::verbose ("Directory") ||
        Behavior::verbose ("Tpetra::Directory");
      const size_t maxNumToPrint = verbose ?
        Behavior::verbosePrintCountThreshold() : size_t(0);

      std::unique_ptr<std::string> procPrefix;
      if (verbose) {
        std::ostringstream os;
        os << "Proc ";
        if (map.getComm ().is_null ()) {
          os << "?";
        }
        else {
          os << map.getComm ()->getRank ();
        }
        os << ": ";
        procPrefix = std::unique_ptr<std::string>(
          new std::string(os.str()));
        os << funcPrefix << "{";
        verbosePrintArray(os, globalIDs, "GIDs", maxNumToPrint);
        os << ", ";
        verbosePrintArray(os, nodeIDs, "PIDs", maxNumToPrint);
        os << ", ";
        verbosePrintArray(os, localIDs, "LIDs", maxNumToPrint);
        os << ", computeLIDs: "
           << (computeLIDs ? "true" : "false") << "}" << endl;
        cerr << os.str ();
      }

      const size_t numEntries = globalIDs.size ();
      const LO LINVALID = Teuchos::OrdinalTraits<LO>::invalid();
      LookupStatus res = AllIDsPresent;

      //
      // Set up directory structure.
      //

      // If we're computing LIDs, we also have to include them in each
      // packet, along with the GID and process ID.
      const int packetSize = computeLIDs ? 3 : 2;

      // For data distribution, we use: Surprise!  A Distributor!
      Distributor distor (comm);

      // Get directory locations for the requested list of entries.
      Array<int> dirImages (numEntries);
      if (verbose) {
        std::ostringstream os;
        os << *procPrefix << "Call directoryMap_->getRemoteIndexList"
           << endl;
        cerr << os.str ();
      }
      res = directoryMap_->getRemoteIndexList (globalIDs, dirImages ());
      if (verbose) {
        std::ostringstream os;
        os << *procPrefix << "Director Map getRemoteIndexList out ";
        verbosePrintArray(os, dirImages, "PIDs", maxNumToPrint);
        os << endl;
        cerr << os.str ();
      }

      // Check for unfound globalIDs and set corresponding nodeIDs to -1
      size_t numMissing = 0;
      if (res == IDNotPresent) {
        for (size_t i=0; i < numEntries; ++i) {
          if (dirImages[i] == -1) {
            nodeIDs[i] = -1;
            if (computeLIDs) {
              localIDs[i] = LINVALID;
            }
            numMissing++;
          }
        }
      }

      Array<GO> sendGIDs;
      Array<int> sendImages;
      if (verbose) {
        std::ostringstream os;
        os << *procPrefix << "Call Distributor::createFromRecvs"
           << endl;
        cerr << os.str ();
      }
      distor.createFromRecvs (globalIDs, dirImages (), sendGIDs, sendImages);
      if (verbose) {
        std::ostringstream os;
        os << *procPrefix << "Distributor::createFromRecvs result: "
           << "{";
        verbosePrintArray(os, sendGIDs, "sendGIDs", maxNumToPrint);
        os << ", ";
        verbosePrintArray(os, sendImages, "sendPIDs", maxNumToPrint);
        os << "}" << endl;
        cerr << os.str ();
      }
      const size_type numSends = sendGIDs.size ();

      //
      // mfh 13 Nov 2012:
      //
      // The code below temporarily stores LO, GO, and int values in
      // an array of global_size_t.  If one of the signed types (LO
      // and GO should both be signed) happened to be -1 (or some
      // negative number, but -1 is the one that came up today), then
      // conversion to global_size_t will result in a huge
      // global_size_t value, and thus conversion back may overflow.
      // (Teuchos::as doesn't know that we meant it to be an LO or GO
      // all along.)
      //
      // The overflow normally would not be a problem, since it would
      // just go back to -1 again.  However, Teuchos::as does range
      // checking on conversions in a debug build, so it throws an
      // exception (std::range_error) in this case.  Range checking is
      // generally useful in debug mode, so we don't want to disable
      // this behavior globally.
      //
      // We solve this problem by forgoing use of Teuchos::as for the
      // conversions below from LO, GO, or int to global_size_t, and
      // the later conversions back from global_size_t to LO, GO, or
      // int.
      //
      // I've recorded this discussion as Bug 5760.
      //

      //    global_size_t >= GO
      //    global_size_t >= size_t >= int
      //    global_size_t >= size_t >= LO
      // Therefore, we can safely store all of these in a global_size_t
      Kokkos::View<global_size_t*, Kokkos::HostSpace> exports("exports", packetSize * numSends);
      {
        // Packet format:
        // - If computing LIDs: (GID, PID, LID)
        // - Otherwise:         (GID, PID)
        //
        // "PID" means "process ID" (a.k.a. "node ID," a.k.a. "rank").

        // Current position to which to write in exports array.  If
        // sending pairs, we pack the (GID, PID) pair for gid =
        // sendGIDs[k] in exports[2*k], exports[2*k+1].  If sending
        // triples, we pack the (GID, PID, LID) pair for gid =
        // sendGIDs[k] in exports[3*k, 3*k+1, 3*k+2].
        size_t exportsIndex = 0;

        if (useHashTables_) {
          if (verbose) {
            std::ostringstream os;
            os << *procPrefix << "Pack exports (useHashTables_ true)"
               << endl;
            cerr << os.str ();
          }
          for (size_type gidIndex = 0; gidIndex < numSends; ++gidIndex) {
            const GO curGID = sendGIDs[gidIndex];
            // Don't use as() here (see above note).
            exports[exportsIndex++] = static_cast<global_size_t> (curGID);
            const LO curLID = directoryMap_->getLocalElement (curGID);
            TEUCHOS_TEST_FOR_EXCEPTION
              (curLID == LINVALID, std::logic_error, funcPrefix <<
               "Directory Map's global index " << curGID << " lacks "
               "a corresponding local index." << errSuffix);
            // Don't use as() here (see above note).
            exports[exportsIndex++] =
              static_cast<global_size_t> (lidToPidTable_->get (curLID));
            if (computeLIDs) {
              // Don't use as() here (see above note).
              exports[exportsIndex++] =
                static_cast<global_size_t> (lidToLidTable_->get (curLID));
            }
          }
        }
        else {
          if (verbose) {
            std::ostringstream os;
            os << *procPrefix << "Pack exports (useHashTables_ false)"
               << endl;
            cerr << os.str ();
          }
          for (size_type gidIndex = 0; gidIndex < numSends; ++gidIndex) {
            const GO curGID = sendGIDs[gidIndex];
            // Don't use as() here (see above note).
            exports[exportsIndex++] = static_cast<global_size_t> (curGID);
            const LO curLID = directoryMap_->getLocalElement (curGID);
            TEUCHOS_TEST_FOR_EXCEPTION
              (curLID == LINVALID, std::logic_error, funcPrefix <<
               "Directory Map's global index " << curGID << " lacks "
               "a corresponding local index." << errSuffix);
            // Don't use as() here (see above note).
            exports[exportsIndex++] =
              static_cast<global_size_t> (PIDs_[curLID]);
            if (computeLIDs) {
              // Don't use as() here (see above note).
              exports[exportsIndex++] =
                static_cast<global_size_t> (LIDs_[curLID]);
            }
          }
        }

        TEUCHOS_TEST_FOR_EXCEPTION
          (exportsIndex > exports.size(), std::logic_error,
           funcPrefix << "On Process " << comm->getRank () << ", "
           "exportsIndex = " << exportsIndex << " > exports.size() = "
           << exports.size () << "." << errSuffix);
      }

      TEUCHOS_TEST_FOR_EXCEPTION
        (numEntries < numMissing, std::logic_error, funcPrefix <<
         "On Process " << comm->getRank () << ", numEntries = " <<
         numEntries << " < numMissing = " << numMissing << "." <<
         errSuffix);

      //
      // mfh 13 Nov 2012: See note above on conversions between
      // global_size_t and LO, GO, or int.
      //
      const size_t numRecv = numEntries - numMissing;

      {
        const size_t importLen = packetSize * distor.getTotalReceiveLength ();
        const size_t requiredImportLen = numRecv * packetSize;
        const int myRank = comm->getRank ();
        TEUCHOS_TEST_FOR_EXCEPTION(
          importLen < requiredImportLen, std::logic_error,
          "Tpetra::Details::DistributedNoncontiguousDirectory::getEntriesImpl: "
          "On Process " << myRank << ": The 'imports' array must have length "
          "at least " << requiredImportLen << ", but its actual length is " <<
          importLen << ".  numRecv: " << numRecv << ", packetSize: " <<
          packetSize << ", numEntries (# GIDs): " << numEntries <<
          ", numMissing: " << numMissing << ": distor.getTotalReceiveLength(): "
          << distor.getTotalReceiveLength () << ".  " << std::endl <<
          "Distributor description: " << distor.description () << ".  "
          << std::endl <<
          "Please report this bug to the Tpetra developers.");
      }

      Kokkos::View<global_size_t*, Kokkos::HostSpace> imports("imports", packetSize * distor.getTotalReceiveLength());
      // FIXME (mfh 20 Mar 2014) One could overlap the sort2() below
      // with communication, by splitting this call into doPosts and
      // doWaits.  The code is still correct in this form, however.
      if (verbose) {
        std::ostringstream os;
        os << *procPrefix << "Call doPostsAndWaits: {"
           << "packetSize: " << packetSize << ", ";
        verbosePrintArray(os, exports, "exports", maxNumToPrint);
        os << "}" << endl;
        cerr << os.str ();
      }
      distor.doPostsAndWaits(exports, packetSize, imports);
      if (verbose) {
        std::ostringstream os;
        os << *procPrefix << "doPostsAndWaits result: ";
        verbosePrintArray(os, imports, "imports", maxNumToPrint);
        os << endl;
        cerr << os.str ();
      }

      Array<GO> sortedIDs (globalIDs); // deep copy (for later sorting)
      Array<GO> offset (numEntries); // permutation array (sort2 output)
      for (GO ii = 0; ii < static_cast<GO> (numEntries); ++ii) {
        offset[ii] = ii;
      }
      sort2 (sortedIDs.begin(), sortedIDs.begin() + numEntries, offset.begin());
      if (verbose) {
        std::ostringstream os;
        os << *procPrefix;
        verbosePrintArray(os, sortedIDs, "sortedIDs", maxNumToPrint);
        os << ", ";
        verbosePrintArray(os, offset, "offset", maxNumToPrint);
        os << endl;
        cerr << os.str ();
      }

      size_t importsIndex = 0;

      // we know these conversions are in range, because we loaded this data
      //
      // Don't use as() for conversions here; we know they are in range.
      for (size_t i = 0; i < numRecv; ++i) {
        const GO curGID = static_cast<GO> (imports[importsIndex++]);
        auto p1 = std::equal_range (sortedIDs.begin(),
                                    sortedIDs.end(), curGID);
        if (p1.first != p1.second) {
          const size_t j = p1.first - sortedIDs.begin();
          nodeIDs[offset[j]] =
            static_cast<int> (imports[importsIndex++]);
          if (computeLIDs) {
            localIDs[offset[j]] =
              static_cast<LO> (imports[importsIndex++]);
          }
          if (nodeIDs[offset[j]] == -1) {
            res = IDNotPresent;
          }
        }
      }

      TEUCHOS_TEST_FOR_EXCEPTION
        (size_t (importsIndex) > size_t (imports.size ()),
         std::logic_error, funcPrefix << "On Process " <<
         comm->getRank () << ": importsIndex = " << importsIndex
         << " > imports.size() = " << imports.size () << ".  "
         "numRecv: " << numRecv
         << ", packetSize: " << packetSize << ", "
         "numEntries (# GIDs): " << numEntries
         << ", numMissing: " << numMissing
         << ": distor.getTotalReceiveLength(): "
         << distor.getTotalReceiveLength () << "." << errSuffix);
      if (verbose) {
        std::ostringstream os;
        os << *procPrefix << funcPrefix << "Done!" << endl;
        cerr << os.str ();
      }
      return res;
    }


    template<class LO, class GO, class NT>
    bool
    DistributedNoncontiguousDirectory<LO, GO, NT>::
    isOneToOne (const Teuchos::Comm<int>& comm) const
    {
      if (oneToOneResult_ == ONE_TO_ONE_NOT_CALLED_YET) {
        const int lcl121 = isLocallyOneToOne () ? 1 : 0;
        int gbl121 = 0;
        Teuchos::reduceAll<int, int> (comm, Teuchos::REDUCE_MIN, lcl121,
                                      Teuchos::outArg (gbl121));
        oneToOneResult_ = (gbl121 == 1) ? ONE_TO_ONE_TRUE : ONE_TO_ONE_FALSE;
      }
      return (oneToOneResult_ == ONE_TO_ONE_TRUE);
    }
  } // namespace Details
} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//
#define TPETRA_DIRECTORYIMPL_INSTANT(LO,GO,NODE) \
  namespace Details { \
    template class Directory< LO , GO , NODE >; \
    template class ReplicatedDirectory< LO , GO , NODE >; \
    template class ContiguousUniformDirectory< LO, GO, NODE >; \
    template class DistributedContiguousDirectory< LO , GO , NODE >; \
    template class DistributedNoncontiguousDirectory< LO , GO , NODE >; \
  }

#endif // TPETRA_DIRECTORYIMPL_DEF_HPP
