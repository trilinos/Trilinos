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

#ifndef TPETRA_DETAILS_TRANSFER_DECL_HPP
#define TPETRA_DETAILS_TRANSFER_DECL_HPP

#include "Tpetra_Details_Transfer_fwd.hpp"
#include "Tpetra_Map_decl.hpp"
#include "Teuchos_Describable.hpp"

namespace Tpetra {
  //
  // Forward declaration.  The "doxygen" bit simply tells Doxygen (our
  // automatic documentation generation system) to skip forward
  // declarations.
  //
#ifndef DOXYGEN_SHOULD_SKIP_THIS
class Distributor;
#endif // DOXYGEN_SHOULD_SKIP_THIS

namespace Details {

/// \class Transfer
/// \brief Common base class of Import and Export
/// \warning This is an implementation detail of Tpetra.  We make no
///   promises of backwards compatibility with this class.
template<class LO,
         class GO,
         class NT>
class Transfer : public Teuchos::Describable {
public:
  //! Destructor (declared virtual for memory safety of derived classes).
  virtual ~Transfer () {}

  /// \brief The specialization of Map used by this class and subclasses.
  ///
  /// The initial two colons avoid confusion between Tpetra::Map and
  /// Tpetra::Detaills::Map.
  typedef ::Tpetra::Map<LO, GO, NT> map_type;

  /// \brief Number of initial identical IDs.
  ///
  /// The number of IDs that are identical between the source and
  /// target Maps, up to the first different ID.
  virtual size_t getNumSameIDs() const = 0;

  /// \brief Number of IDs to permute but not to communicate.
  ///
  /// The number of IDs that are local to the calling process, but not
  /// part of the first getNumSameIDs() entries.  The Import or Export
  /// will permute these entries locally (without distributed-memory
  /// communication).
  virtual size_t getNumPermuteIDs () const = 0;

  //! List of local IDs in the source Map that are permuted.
  virtual Teuchos::ArrayView<const LO> getPermuteFromLIDs () const = 0;

  //! List of local IDs in the target Map that are permuted.
  virtual Teuchos::ArrayView<const LO> getPermuteToLIDs () const = 0;

  //! Number of entries not on the calling process.
  virtual size_t getNumRemoteIDs () const = 0;

  //! List of entries in the target Map to receive from other processes.
  virtual Teuchos::ArrayView<const LO> getRemoteLIDs () const = 0;

  //! Number of entries that must be sent by the calling process to other processes.
  virtual size_t getNumExportIDs () const = 0;

  //! List of entries in the source Map that will be sent to other processes.
  virtual Teuchos::ArrayView<const LO> getExportLIDs () const = 0;

  /// \brief List of processes to which entries will be sent.
  ///
  /// The entry with local ID getExportLIDs()[i] will be sent to
  /// process getExportPiDs()[i].
  virtual Teuchos::ArrayView<const int> getExportPIDs () const = 0;

  //! The source Map used to construct this Export or Import.
  virtual Teuchos::RCP<const map_type> getSourceMap () const = 0;

  //! The target Map used to construct this Export or Import.
  virtual Teuchos::RCP<const map_type> getTargetMap () const = 0;

  //! The Distributor that this Export or Import object uses to move data.
  virtual ::Tpetra::Distributor& getDistributor () const = 0;

  /// \brief Is this Export or Import locally complete?
  ///
  /// If this is an Export, then do all source Map indices on the
  /// calling process exist on at least one process (not necessarily
  /// this one) in the target Map?
  ///
  /// If this is an Import, then do all target Map indices on the
  /// calling process exist on at least one process (not necessarily
  /// this one) in the source Map?
  ///
  /// It's not necessarily an error for an Export or Import not to be
  /// locally complete on one or more processes.  For example, this
  /// may happen in the common use case of "restriction" -- that is,
  /// taking a subset of a large object.  Nevertheless, you may find
  /// this predicate useful for figuring out whether you set up your
  /// Maps in the way that you expect.
  virtual bool isLocallyComplete () const = 0;

  /// \brief Describe this object in a human-readable way to the given
  ///   output stream.
  ///
  /// You must call this method as a collective over all processes in
  /// the communicator of the source and target Map of this object.
  ///
  /// \param out [out] Output stream to which to write.  Only Process
  ///   0 in this object's communicator may write to the output
  ///   stream.
  ///
  /// \param verbLevel [in] Verbosity level.  This also controls
  ///   whether this method does any communication.  At verbosity
  ///   levels higher (greater) than Teuchos::VERB_LOW, this method
  ///   behaves as a collective over the object's communicator.
  ///
  /// Teuchos::FancyOStream wraps std::ostream.  It adds features like
  /// tab levels.  If you just want to wrap std::cout, try this:
  /// \code
  /// auto out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::out));
  /// \endcode
  virtual void
  describe (Teuchos::FancyOStream& out,
            const Teuchos::EVerbosityLevel verbLevel =
              Teuchos::Describable::verbLevel_default) const;

protected:
  /// \brief Implementation of describe() for subclasses
  ///   (Tpetra::Import and Tpetra::Export).
  ///
  /// \param out [out] Output stream to which to write.  Only Process
  ///   0 in this object's communicator may write to the output stream.
  ///
  /// \param className [in] Name of the subclass of Transfer calling
  ///   this method.
  ///
  /// \param verbLevel [in] Verbosity level.  This also controls
  ///   whether this method does any communication.  At verbosity
  ///   levels higher (greater) than Teuchos::VERB_LOW, this method
  ///   behaves as a collective over the object's communicator.
  void
  describeImpl (Teuchos::FancyOStream& out,
                const std::string& className,
                const Teuchos::EVerbosityLevel verbLevel =
                  Teuchos::Describable::verbLevel_default) const;

private:
  /// \brief Print "global" (not necessarily just on Process 0)
  ///   information for describe().
  ///
  /// This is an implementation detail of describe().
  ///
  /// \param out [out] The output stream argument given to describe().
  void
  globalDescribe (Teuchos::FancyOStream& out,
                  const Teuchos::EVerbosityLevel vl) const;

  /// \brief Print the calling process' verbose describe() information
  ///   to the given output string.
  ///
  /// This is an implementation detail of describe().
  std::string
  localDescribeToString (const Teuchos::EVerbosityLevel vl) const;
};

} // namespace Details
} // namespace Tpetra

// Explicit instantiation macro.
// Only invoke this when in the Tpetra::Details namespace.
// Most users do not need to use this.
//
// LO: The local ordinal type.
// GO: The global ordinal type.
// NODE: The Kokkos Node type.
#define TPETRA_DETAILS_TRANSFER_INSTANT(LO, GO, NODE) \
  template class Transfer< LO , GO , NODE >;

#endif // TPETRA_DETAILS_TRANSFER_DECL_HPP
