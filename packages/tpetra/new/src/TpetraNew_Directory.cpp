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

#include "TpetraNew_Directory.hpp"
#include "Tpetra_Distributor.hpp"
#include "TpetraNew_Map.hpp"
#include "TpetraNew_DirectoryImpl.hpp"

namespace TpetraNew {

  Directory::Directory () :
    impl_ (nullptr)
  {}

  Directory::~Directory () {
    if (impl_ != nullptr) {
      delete impl_;
      impl_ = nullptr;
    }
  }

  bool
  Directory::initialized () const {
    return impl_ != nullptr;
  }

  void
  Directory::
  initialize (const map_type& map,
              const ::Tpetra::Details::TieBreak<local_ordinal_type, global_ordinal_type>& tieBreak)
  {
    using LO = local_ordinal_type;
    using GO = global_ordinal_type;
    const char prefix[] = "TpetraNew::Directory::initialize: ";
    const char suffix[] = ".  Please report this bug to the Tpetra developers.";
    
    if (initialized ()) {
      TEUCHOS_TEST_FOR_EXCEPTION
	(impl_ == nullptr, std::logic_error,
	 prefix << "The Directory claims that it has been initialized, "
	 "but its implementation object has not yet been created" << suffix);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION
	(impl_ != nullptr, std::logic_error, prefix <<
	 "Directory implementation has already been initialized, "
	 "but initialized() returns false" << suffix);

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
      const Details::Directory* dir = nullptr;
      bool usedTieBreak = false;
      if (map.isDistributed ()) {
        if (map.isUniform ()) {
          dir = new Details::ContiguousUniformDirectory (map);
        }
        else if (map.isContiguous ()) {
          dir = new Details::DistributedContiguousDirectory (map);
        }
        else {
          dir = new Details::DistributedNoncontiguousDirectory (map, tieBreak);
          usedTieBreak = true;
        }
      }
      else {
        dir = new Details::ReplicatedDirectory (map);

        if (tieBreak.mayHaveSideEffects () && map.getNodeNumElements () != 0) {
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
          map.getNodeNumElements () != 0) {
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

  void
  Directory::initialize (const map_type& map)
  {
    if (initialized ()) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        impl_ == nullptr, std::logic_error, "Tpetra::Directory::initialize: "
        "The Directory claims that it has been initialized, "
        "but its implementation object has not yet been created.  "
        "Please report this bug to the Tpetra developers.");
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        impl_ != nullptr, std::logic_error, "Tpetra::Directory::initialize: "
        "Directory implementation has already been initialized, "
        "but initialized() returns false.  "
        "Please report this bug to the Tpetra developers.");

      // Create an implementation object of the appropriate type,
      // depending on whether the Map is distributed or replicated,
      // and contiguous or noncontiguous.
      const Details::Directory* dir = nullptr;
      if (map.isDistributed ()) {
        if (map.isUniform ()) {
          dir = new Details::ContiguousUniformDirectory (map);
        }
        else if (map.isContiguous ()) {
          dir = new Details::DistributedContiguousDirectory (map);
        }
        else {
          dir = new Details::DistributedNoncontiguousDirectory (map);
        }
      }
      else {
        dir = new Details::ReplicatedDirectory (map);
      }
      TEUCHOS_TEST_FOR_EXCEPTION(
        dir == nullptr, std::logic_error, "Tpetra::Directory::initialize: "
        "Failed to create Directory implementation.  "
        "Please report this bug to the Tpetra developers.");
      impl_ = dir;
    }
  }

  ::Tpetra::LookupStatus
  Directory::
  getDirectoryEntries (const map_type& map,
                       const Teuchos::ArrayView<const global_ordinal_type>& globalIDs,
                       const Teuchos::ArrayView<int>& nodeIDs) const
  {
    if (! initialized ()) {
      // This const_cast is super wrong, but "mutable" is also a lie,
      // and Map's interface needs this method to be marked const for
      // some reason.
      const_cast<Directory*> (this)->initialize (map);
    }
    const bool computeLIDs = false;
    return impl_->getEntries (map, globalIDs, nodeIDs, Teuchos::null, computeLIDs);
  }

  ::Tpetra::LookupStatus
  Directory::
  getDirectoryEntries (const map_type& map,
                       const Teuchos::ArrayView<const global_ordinal_type>& globalIDs,
                       const Teuchos::ArrayView<int>& nodeIDs,
                       const Teuchos::ArrayView<local_ordinal_type>& localIDs) const
  {
    if (! initialized ()) {
      // This const_cast is super wrong, but "mutable" is also a lie,
      // and Map's interface needs this method to be marked const for
      // some reason.
      const_cast<Directory*> (this)->initialize (map);
    }
    const bool computeLIDs = true;
    return impl_->getEntries (map, globalIDs, nodeIDs, localIDs, computeLIDs);
  }

  bool Directory::isOneToOne (const map_type& map) const {
    if (! initialized ()) {
      // This const_cast is super wrong, but "mutable" is also a lie,
      // and Map's interface needs this method to be marked const for
      // some reason.
      const_cast<Directory*> (this)->initialize (map);
    }
    return impl_->isOneToOne (* (map.getComm ()));
  }

  std::string
  Directory::description () const
  {
    // std::ostringstream os;
    // os << "Directory";
    // return os.str ();
    return "Directory";
  }

} // namespace TpetraNew

