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

#ifndef TPETRA_DIRECTORY_HPP
#define TPETRA_DIRECTORY_HPP

#include <Teuchos_as.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_Distributor.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_DirectoryImpl.hpp>
#include <Tpetra_DirectoryImpl_def.hpp>

#ifdef DOXYGEN_USE_ONLY
#  include "Tpetra_Directory_decl.hpp"
#endif

namespace Tpetra {

  template<class LO, class GO, class NT>
  Directory<LO, GO, NT>::
  Directory (const Teuchos::RCP<const Map<LO, GO, NT> >& map)
  {
    // Create an implementation object of the appropriate type,
    // depending on whether the Map is distributed or replicated, and
    // contiguous or noncontiguous.
    RCP<const Details::Directory<LO, GO, NT> > dir;
    if (map->isDistributed ()) {
      if (map->isContiguous ()) {
        dir = rcp (new Details::DistributedContiguousDirectory<LO, GO, NT> (map));
      }
      else {
        dir = rcp (new Details::DistributedNoncontiguousDirectory<LO, GO, NT> (map));
      }
    }
    else {
      dir = rcp (new Details::ReplicatedDirectory<LO, GO, NT> (map));
    }
    TEUCHOS_TEST_FOR_EXCEPTION(dir.is_null (), std::logic_error, "Tpetra::"
      "Directory constructor failed to create Directory implementation.  "
      "Please report this bug to the Tpetra developers.");
    impl_ = dir;
  }

  template<class LO, class GO, class NT>
  Directory<LO, GO, NT>::Directory() {}

  template<class LO, class GO, class NT>
  Directory<LO, GO, NT>::~Directory() {}

  template<class LO, class GO, class NT>
  LookupStatus
  Directory<LO, GO, NT>::
  getDirectoryEntries (const Teuchos::ArrayView<const GO>& globalIDs,
                       const Teuchos::ArrayView<int>& nodeIDs) const
  {
    const bool computeLIDs = false;
    return impl_->getEntries (globalIDs, nodeIDs, Teuchos::null, computeLIDs);
  }

  template<class LO, class GO, class NT>
  LookupStatus
  Directory<LO, GO, NT>::
  getDirectoryEntries (const Teuchos::ArrayView<const GO>& globalIDs,
                       const Teuchos::ArrayView<int>& nodeIDs,
                       const Teuchos::ArrayView<LO>& localIDs) const
  {
    const bool computeLIDs = true;
    return impl_->getEntries (globalIDs, nodeIDs, localIDs, computeLIDs);
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
  \
  template class Directory< LO , GO , NODE >; \

#endif // TPETRA_DIRECTORY_HPP
