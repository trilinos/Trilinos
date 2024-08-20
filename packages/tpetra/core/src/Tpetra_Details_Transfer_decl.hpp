// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_TRANSFER_DECL_HPP
#define TPETRA_DETAILS_TRANSFER_DECL_HPP

#include "Tpetra_Details_Transfer_fwd.hpp"
#include "Tpetra_ImportExportData_fwd.hpp"
#include "Tpetra_Map_decl.hpp"
#include "Teuchos_Describable.hpp"

namespace Tpetra {

// Forward declaration.  The macro tells Doxygen (our automatic
// documentation generation system) to skip forward declarations.
#ifndef DOXYGEN_SHOULD_SKIP_THIS
class Distributor;
#endif // DOXYGEN_SHOULD_SKIP_THIS

namespace Details {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/// \brief Expert function that messes with internals of Transfer.
/// Used for testing. Do not use unless you know what you are doing.
template <class LO, class GO, class NT>
void
expertSetRemoteLIDsContiguous(Transfer<LO, GO, NT> transfer, bool contig);


/// \brief Expert function that messes with internals of Transfer.
/// Used for testing. Do not use unless you know what you are doing.
template <class LO, class GO, class NT>
void
expertSetExportLIDsContiguous(Transfer<LO, GO, NT> transfer, bool contig);
#endif // DOXYGEN_SHOULD_SKIP_THIS

/// \class Transfer
/// \brief Common base class of Import and Export
/// \warning This is an implementation detail of Tpetra.  We make no
///   promises of backwards compatibility with this class.
template<class LO,
         class GO,
         class NT>
class Transfer : public Teuchos::Describable {
public:
  /// \brief Map specialization used by this class and subclasses.
  ///
  /// The initial two colons avoid confusion between Tpetra::Map and
  /// Tpetra::Details::Map.
  using map_type = ::Tpetra::Map<LO, GO, NT>;

private:
  using execution_space = typename NT::device_type::execution_space;
  using memory_space =
    ::Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename NT::device_type>;
  using device_type = Kokkos::Device<execution_space, memory_space>;

public:
  // Don't be tempted to comment this out if code doesn't build.
  // The point is to ensure correct initialization of TransferData_.
  Transfer () = delete;

  /// \brief Four-argument constructor (most often used).
  ///
  /// \pre </tt> ! source.is_null() </tt>
  ///
  /// \param source [in] Source Map of the Export or Import.
  /// \param target [in] Target Map of the Export or Import.
  ///   May be null only if Export or Import is using one of
  ///   the special source-Map-only constructors.
  /// \param out [in/out] Stream for verbose debugging output.
  ///   If null, Transfer will wrap and use std::cerr.
  /// \param plist [in/out] Parameters; may be null.
  ///
  /// \param className [in] Either "Export" or "Import".  Used to
  ///   control verbose debugging output (sometimes you might want it
  ///   only for one or the other class).
  Transfer (const Teuchos::RCP<const map_type>& source,
            const Teuchos::RCP<const map_type>& target,
            const Teuchos::RCP<Teuchos::FancyOStream>& out,
            const Teuchos::RCP<Teuchos::ParameterList>& plist,
            const std::string& className);

  Transfer (const Transfer<LO, GO, NT>& rhs) = default;

  struct reverse_tag {};
  /// \brief Reverse-mode "copy" constructor.
  ///
  /// Use this for constructing an Export from an Import, or an Import
  /// from an Export.
  Transfer (const Transfer<LO, GO, NT>& rhs, reverse_tag tag);

  Transfer<LO, GO, NT>& operator= (const Transfer<LO, GO, NT>&) = default;

  //! Destructor (declared virtual for memory safety of derived classes).
  virtual ~Transfer () = default;

  /// \brief Number of initial identical IDs.
  ///
  /// The number of IDs that are identical between the source and
  /// target Maps, up to the first different ID.
  size_t getNumSameIDs() const;

  /// \brief Number of IDs to permute but not to communicate.
  ///
  /// The number of IDs that are local to the calling process, but not
  /// part of the first getNumSameIDs() entries.  The Import or Export
  /// will permute these entries locally (without distributed-memory
  /// communication).
  size_t getNumPermuteIDs () const;

  /// \brief List of local IDs in the source Map that are permuted, as
  ///   a const DualView (that is sync'd to both host and device).
  Kokkos::DualView<const LO*, device_type> getPermuteFromLIDs_dv () const;

  //! List of local IDs in the source Map that are permuted.
  Teuchos::ArrayView<const LO> getPermuteFromLIDs () const;

  /// \brief List of local IDs in the target Map that are permuted, as
  ///   a const DualView (that is sync'd to both host and device).
  Kokkos::DualView<const LO*, device_type> getPermuteToLIDs_dv () const;

  //! List of local IDs in the target Map that are permuted.
  Teuchos::ArrayView<const LO> getPermuteToLIDs () const;

  //! Number of entries not on the calling process.
  size_t getNumRemoteIDs () const;

  /// \brief List of entries in the target Map to receive from other
  ///   processes, as a const DualView (that is sync'd to both host
  ///   and device).
  Kokkos::DualView<const LO*, device_type> getRemoteLIDs_dv () const;

  //! List of entries in the target Map to receive from other processes.
  Teuchos::ArrayView<const LO> getRemoteLIDs () const;

  //! Number of entries that must be sent by the calling process to other processes.
  size_t getNumExportIDs () const;

  /// \brief List of entries in the source Map that will be sent to
  ///   other processes, as a const DualView (that is sync'd to both
  ///   host and device).
  Kokkos::DualView<const LO*, device_type> getExportLIDs_dv () const;

  //! List of entries in the source Map that will be sent to other processes.
  Teuchos::ArrayView<const LO> getExportLIDs () const;

  /// \brief List of processes to which entries will be sent.
  ///
  /// The entry with local ID getExportLIDs()[i] will be sent to
  /// process getExportPiDs()[i].
  Teuchos::ArrayView<const int> getExportPIDs () const;

  //! The source Map used to construct this Export or Import.
  Teuchos::RCP<const map_type> getSourceMap () const;

  //! The target Map used to construct this Export or Import.
  Teuchos::RCP<const map_type> getTargetMap () const;

  //! The Distributor that this Export or Import object uses to move data.
  ::Tpetra::Distributor& getDistributor () const;

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
  bool isLocallyComplete () const;


  void detectRemoteExportLIDsContiguous() const;

  bool areRemoteLIDsContiguous() const;

  bool areExportLIDsContiguous() const;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  //! Friend declaration for non-member function that allows setting a flag.
  friend
  void expertSetRemoteLIDsContiguous<LO,GO,NT>(Transfer<LO, GO, NT> transfer, bool contig);

  //! Friend declaration for non-member function that allows setting a flag.
  friend
  void expertSetExportLIDsContiguous<LO,GO,NT>(Transfer<LO, GO, NT> transfer, bool contig);
#endif // DOXYGEN_SHOULD_SKIP_THIS

  /// \brief Are source and target map locally fitted?
  ///
  /// Returns whether source and target map are locally fitted on the
  /// calling rank. This is can be more efficient that calling
  /// isLocallyFitted() on the maps directly, since no indices need to
  /// be compared.
  bool isLocallyFitted () const;

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
  //! All the data needed for executing the Export communication plan.
  Teuchos::RCP<ImportExportData<LO, GO, NT> > TransferData_;

  //! Valid (nonnull) output stream for verbose output.
  Teuchos::FancyOStream& verboseOutputStream () const;

  //! Whether to print verbose debugging output.
  bool verbose () const;

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
  /// \brief Set parameters.
  ///
  /// Please see the Export or Import class documentation for a list
  /// of all accepted parameters and their default values.
  ///
  /// \param plist [in/out] Parameters; may be null.
  /// \param className [in] "Export" or "Import".
  void
  setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& plist,
                    const std::string& className);

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
  template class Transfer< LO , GO , NODE >;          \
  template void expertSetRemoteLIDsContiguous< LO , GO , NODE >(Transfer<LO, GO, NODE> transfer, bool contig);    \
  template void expertSetExportLIDsContiguous< LO , GO , NODE >(Transfer<LO, GO, NODE> transfer, bool contig);

#endif // TPETRA_DETAILS_TRANSFER_DECL_HPP
