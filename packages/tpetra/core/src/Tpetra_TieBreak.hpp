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

#ifndef TPETRA_DETAILS_TIEBREAK_HPP
#define TPETRA_DETAILS_TIEBREAK_HPP

/// \file Tpetra_TieBreak.hpp
/// \brief Interface for breaking ties in ownership.

#include <Teuchos_RCP.hpp>

namespace Tpetra {
namespace Details {

  /// \class TieBreak
  /// \brief Interface for breaking ties in ownership.
  /// \tparam LocalOrdinal The type of local indices.
  /// \tparam GlobalOrdinal The type of global indices.
  ///
  /// This class provides an abstract interface to a way for Directory
  /// to "break ties" in ownership of global indices.  To "break a
  /// tie" for a global index GID means: Given a set of one or more
  /// (PID, LID) pairs, where each PID is a process that owns GID and
  /// LID is the corresponding local index of GID, choose exactly one
  /// (PID, LID) pair from that set.  Furthermore, this choice must be
  /// a <i>global</i> choice, that is, the same on all participating
  /// processes.
  template <typename LocalOrdinal, typename GlobalOrdinal>
  class TieBreak {
  public:
    /// \brief Representation of a global index on a process.
    ///
    /// This struct holds a global index (GID), a process that owns it
    /// (PID), and its local index (LID) on that process.
    ///
    /// FIXME (mfh 15 Sep 2013) This should be an implementation
    /// detail of subclasses; there should be no need to expose it in
    /// the public interface.
    ///
    /// FIXME (mfh 15 Sep 2013) PID should go last, so this struct
    /// would pack into 128 bits if <tt>LocalOrdinal</tt> is 32 bits
    /// and <tt>GlobalOrdinal</tt> is 64 bits.  I would fix the order
    /// myself, but I'm not sure if any downstream code depends on it.
    struct Triplet {
      LocalOrdinal LID;
      GlobalOrdinal GID;
      int PID;
    };

    /// \brief Whether selectedIndex() may have side effects.
    ///
    /// If you are defining your own implementation (i.e., subclass)
    /// of this class, and if you know that your implementation of
    /// selectedIndex() does <i>not</i> have side effects, then you
    /// may redefine this method to return false.  This may skip using
    /// the TieBreak object in certain cases, to save local work.
    /// Otherwise, we must use the TieBreak object in case the
    /// implementation has desired side effects.
    virtual bool mayHaveSideEffects () const {
      return true;
    }

    //! Virtual destructor (for memory safety of derived classes).
    virtual ~TieBreak () {}

    /// \brief Break any ties in ownership of the given global index GID.
    ///
    /// Given a global index GID, and a set of (PID, LID) pairs (of
    /// processes that own GID, and the local index of GID on that
    /// process), return the index of one pair in the list.
    ///
    /// This method must always be called collectively over all the
    /// processes in the Directory's communicator.  Subclasses reserve
    /// the right to use communication (either point-to-point or
    /// collective) over that communicator, but are not required to do
    /// so.  However, their decisions are required to be
    /// <i>consistent</i> over those processes.  This means that if
    /// multiple processes call this method with the same GID and the
    /// same list of pairs, all these processes must return the same
    /// index.  (It would also be a good idea for subclasses not to be
    /// sensitive to the order of pairs.)
    ///
    /// FIXME (mfh 17 Sep 2013) I'm not a fan of the name of this
    /// method.  We should call it something more indicative of its
    /// function, like "arbitrateOwnership" or "breakTie".
    virtual std::size_t
    selectedIndex (GlobalOrdinal GID,
                   const std::vector<std::pair<int, LocalOrdinal> >& pid_and_lid) const = 0;
  };

} // namespace Tpetra
} // namespace Details

#endif // TPETRA_DETAILS_TIEBREAK_HPP
