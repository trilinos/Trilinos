// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DIRECTORY_HPP
#define TPETRA_DIRECTORY_HPP

#include "Tpetra_Distributor.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_DirectoryImpl.hpp"
#include "Tpetra_Directory_decl.hpp"

namespace Tpetra {

  template<class LO, class GO, class NT>
  Directory<LO, GO, NT>::Directory () :
    impl_ (NULL)
  {}

  template<class LO, class GO, class NT>
  Directory<LO, GO, NT>::~Directory () {
    if (impl_ != NULL) {
      delete impl_;
      impl_ = NULL;
    }
  }

  template<class LO, class GO, class NT>
  bool
  Directory<LO, GO, NT>::initialized () const {
    return impl_ != NULL;
  }


  template<class LO, class GO, class NT>
  void
  Directory<LO, GO, NT>::
  initialize (const Map<LO, GO, NT>& map,
              const Tpetra::Details::TieBreak<LO,GO>& tieBreak)
  {
    if (initialized ()) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        impl_ == NULL, std::logic_error, "Tpetra::Directory::initialize: "
        "The Directory claims that it has been initialized, "
        "but its implementation object has not yet been created.  "
        "Please report this bug to the Tpetra developers.");
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        impl_ != NULL, std::logic_error, "Tpetra::Directory::initialize: "
        "Directory implementation has already been initialized, "
        "but initialized() returns false.  "
        "Please report this bug to the Tpetra developers.");

      // Create an implementation object of the appropriate type,
      // depending on whether the Map is distributed or replicated,
      // and contiguous or noncontiguous.
      //
      // mfh 06 Apr 2014: When a distributed noncontiguous Directory
      // takes a TieBreak, all the entries (local indices and process
      // ranks) owned by the Directory on the calling process pass
      // through the TieBreak object.  This may have side effects,
      // such as the TieBreak object remembering whether there were
      // any duplicates on the calling process.  We want to extend use
      // of a TieBreak object to other kinds of Directories.  For a
      // distributed contiguous Directory, the calling process owns
      // all of the (PID,LID) pairs in the input Map.  For a locally
      // replicated contiguous Directory, Process 0 owns all of the
      // (PID,LID) pairs in the input Map.
      //
      // It may seem silly to pass in a TieBreak when there are no
      // ties to break.  However, the TieBreak object gets to see all
      // (PID,LID) pairs that the Directory owns on the calling
      // process, and interface of TieBreak allows side effects.
      // Users may wish to exploit them regardless of the kind of Map
      // they pass in.
      const ::Tpetra::Details::Directory<LO, GO, NT>* dir = NULL;
      bool usedTieBreak = false;
      if (map.isDistributed ()) {
        if (map.isUniform ()) {
          dir = new ::Tpetra::Details::ContiguousUniformDirectory<LO, GO, NT> (map);
        }
        else if (map.isContiguous ()) {
          dir = new ::Tpetra::Details::DistributedContiguousDirectory<LO, GO, NT> (map);
        }
        else {
          dir = new ::Tpetra::Details::DistributedNoncontiguousDirectory<LO, GO, NT> (map, tieBreak);
          usedTieBreak = true;
        }
      }
      else {
        dir = new ::Tpetra::Details::ReplicatedDirectory<LO, GO, NT> (map);

        if (tieBreak.mayHaveSideEffects () && map.getLocalNumElements () != 0) {
          // We need the second clause in the above test because Map's
          // interface provides an inclusive range of local indices.
          const int myRank = map.getComm ()->getRank ();
          // In a replicated Directory, Process 0 owns all the
          // Directory's entries.  This is an arbitrary assignment; any
          // one process would do.
          if (myRank == 0) {
            std::vector<std::pair<int, LO> > pidLidList (1);
            const LO minLocInd = map.getMinLocalIndex ();
            const LO maxLocInd = map.getMaxLocalIndex ();
            for (LO locInd = minLocInd; locInd <= maxLocInd; ++locInd) {
              pidLidList[0] = std::make_pair (myRank, locInd);
              const GO globInd = map.getGlobalElement (locInd);
              // We don't care about the return value; we just want to
              // invoke the side effects.
              (void) tieBreak.selectedIndex (globInd, pidLidList);
            }
          }
        }
        usedTieBreak = true;
      } // done with all different Map cases

      // If we haven't already used the TieBreak object, use it now.
      // This code appears twice because ReplicatedDirectory is a
      // special case: we already know what gets replicated.
      if (! usedTieBreak && tieBreak.mayHaveSideEffects () &&
          map.getLocalNumElements () != 0) {
        // We need the third clause in the above test because Map's
        // interface provides an inclusive range of local indices.
        std::vector<std::pair<int, LO> > pidLidList (1);
        const LO minLocInd = map.getMinLocalIndex ();
        const LO maxLocInd = map.getMaxLocalIndex ();
        const int myRank = map.getComm ()->getRank ();
        for (LO locInd = minLocInd; locInd <= maxLocInd; ++locInd) {
          pidLidList[0] = std::make_pair (myRank, locInd);
          const GO globInd = map.getGlobalElement (locInd);
          // We don't care about the return value; we just want to
          // invoke the side effects.
          (void) tieBreak.selectedIndex (globInd, pidLidList);
        }
      }

      impl_ = dir;
    }
  }

  template<class LO, class GO, class NT>
  void
  Directory<LO, GO, NT>::initialize (const Map<LO, GO, NT>& map)
  {
    if (initialized ()) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        impl_ == NULL, std::logic_error, "Tpetra::Directory::initialize: "
        "The Directory claims that it has been initialized, "
        "but its implementation object has not yet been created.  "
        "Please report this bug to the Tpetra developers.");
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        impl_ != NULL, std::logic_error, "Tpetra::Directory::initialize: "
        "Directory implementation has already been initialized, "
        "but initialized() returns false.  "
        "Please report this bug to the Tpetra developers.");

      // Create an implementation object of the appropriate type,
      // depending on whether the Map is distributed or replicated,
      // and contiguous or noncontiguous.
      const ::Tpetra::Details::Directory<LO, GO, NT>* dir = NULL;
      if (map.isDistributed ()) {
        if (map.isUniform ()) {
          dir = new ::Tpetra::Details::ContiguousUniformDirectory<LO, GO, NT> (map);
        }
        else if (map.isContiguous ()) {
          dir = new ::Tpetra::Details::DistributedContiguousDirectory<LO, GO, NT> (map);
        }
        else {
          dir = new ::Tpetra::Details::DistributedNoncontiguousDirectory<LO, GO, NT> (map);
        }
      }
      else {
        dir = new ::Tpetra::Details::ReplicatedDirectory<LO, GO, NT> (map);
      }
      TEUCHOS_TEST_FOR_EXCEPTION(
        dir == NULL, std::logic_error, "Tpetra::Directory::initialize: "
        "Failed to create Directory implementation.  "
        "Please report this bug to the Tpetra developers.");
      impl_ = dir;
    }
  }

  template<class LO, class GO, class NT>
  LookupStatus
  Directory<LO, GO, NT>::
  getDirectoryEntries (const Map<LO, GO, NT>& map,
                       const Teuchos::ArrayView<const GO>& globalIDs,
                       const Teuchos::ArrayView<int>& nodeIDs) const
  {
    if (! initialized ()) {
      // This const_cast is super wrong, but "mutable" is also a lie,
      // and Map's interface needs this method to be marked const for
      // some reason.
      const_cast<Directory<LO, GO, NT>* > (this)->initialize (map);
    }
    const bool computeLIDs = false;
    return impl_->getEntries (map, globalIDs, nodeIDs, Teuchos::null, computeLIDs);
  }

  template<class LO, class GO, class NT>
  LookupStatus
  Directory<LO, GO, NT>::
  getDirectoryEntries (const Map<LO, GO, NT>& map,
                       const Teuchos::ArrayView<const GO>& globalIDs,
                       const Teuchos::ArrayView<int>& nodeIDs,
                       const Teuchos::ArrayView<LO>& localIDs) const
  {
    if (! initialized ()) {
      // This const_cast is super wrong, but "mutable" is also a lie,
      // and Map's interface needs this method to be marked const for
      // some reason.
      const_cast<Directory<LO, GO, NT>* > (this)->initialize (map);
    }
    const bool computeLIDs = true;
    return impl_->getEntries (map, globalIDs, nodeIDs, localIDs, computeLIDs);
  }

  template<class LO, class GO, class NT>
  bool Directory<LO, GO, NT>::isOneToOne (const Map<LO, GO, NT>& map) const {
    if (! initialized ()) {
      // This const_cast is super wrong, but "mutable" is also a lie,
      // and Map's interface needs this method to be marked const for
      // some reason.
      const_cast<Directory<LO, GO, NT>* > (this)->initialize (map);
    }
    return impl_->isOneToOne (* (map.getComm ()));
  }

  template<class LO, class GO, class NT>
  std::string
  Directory<LO, GO, NT>::description () const
  {
    using Teuchos::TypeNameTraits;

    std::ostringstream os;
    os << "Directory"
       << "<" << TypeNameTraits<LO>::name ()
       << ", " << TypeNameTraits<GO>::name ()
       << ", " << TypeNameTraits<NT>::name () << ">";
    return os.str ();
  }

} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_DIRECTORY_INSTANT(LO,GO,NODE) \
  template class Directory< LO , GO , NODE >;

#endif // TPETRA_DIRECTORY_HPP
