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

#ifndef TPETRA_EXPORT_DECL_HPP
#define TPETRA_EXPORT_DECL_HPP

#include "Tpetra_Details_Transfer.hpp"
#include "Tpetra_Export_fwd.hpp"
#include "Tpetra_Import_fwd.hpp"
#include "Tpetra_ImportExportData_fwd.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_RCP.hpp"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Teuchos {
template<class T> class Array; // forward declaration
class ParameterList; // forward declaration
} // namespace Teuchos
#endif // DOXYGEN_SHOULD_SKIP_THIS

namespace Tpetra {
  /// \brief Communication plan for data redistribution from a
  ///   (possibly) multiply-owned to a uniquely-owned distribution.
  ///
  /// \tparam LocalOrdinal The type of local indices.  See the
  ///   documentation of Map for requirements.
  /// \tparam GlobalOrdinal The type of global indices.  See the
  ///   documentation of Map for requirements.
  /// \tparam Node The Kokkos Node type.  See the documentation of Map
  ///   for requirements.
  ///
  /// Tpetra users should use this class to construct a communication
  /// plan between two data distributions (i.e., two Map objects).
  /// The plan can be called repeatedly by computational classes to
  /// perform communication according to the same pattern.
  /// Constructing the plan may be expensive, both in terms of
  /// communication and computation.  However, it can be reused
  /// inexpensively.
  ///
  /// Tpetra has two classes for data redistribution: Import and
  /// Export.  Import is for redistributing data from a uniquely-owned
  /// distribution to a possibly multiply-owned distribution.  Export
  /// is for redistributing data from a possibly multiply-owned
  /// distribution to a uniquely-owned distribution.
  ///
  /// The names "Import" and "Export" have nothing to do with the
  /// direction in which data moves relative to the calling process;
  /// any process may do both receives and sends in an Import or
  /// Export.  Rather, the names suggest what happens in their most
  /// common use case, the communication pattern for sparse
  /// matrix-vector multiply.  Import "brings in" remote source vector
  /// data (from the domain Map to the column Map) for local
  /// computation, and Export "pushes" the result back (from the row
  /// Map to the range Map).  Import and Export have other uses as
  /// well.
  ///
  /// One use case of Export is finite element assembly.  For example,
  /// one way to compute a distributed forcing term vector is to use
  /// an overlapping distribution for the basis functions' domains.
  /// An Export with the SUM combine mode combines each process'
  /// contribution to the integration into a single nonoverlapping
  /// distribution.
  ///
  /// Epetra separated Import and Export for performance reasons.  The
  /// implementation is different, depending on which direction is the
  /// uniquely-owned Map.  Tpetra retains this convention.
  ///
  /// This class is templated on the same template arguments as Map:
  /// the local ordinal type <tt>LocalOrdinal</tt>, the global ordinal
  /// type <tt>GlobalOrdinal</tt>, and the Kokkos <tt>Node</tt> type.
  ///
  /// This method accepts an optional list of parameters, either
  /// through the constructor or through the setParameterList()
  /// method.  Most users do not need to worry about these parameters;
  /// the default values are fine.  However, for expert users, we
  /// expose the following parameter:
  /// - "Barrier between receives and sends" (\c bool): Whether to
  ///   execute a barrier between receives and sends, when executing
  ///   the Import (i.e., when calling DistObject's doImport()
  ///   (forward mode) or doExport() (reverse mode)).
  ///
  template<class LocalOrdinal,
           class GlobalOrdinal,
           class Node>
  class Export:
    public ::Tpetra::Details::Transfer<LocalOrdinal, GlobalOrdinal, Node>
  {
  private:
    friend class Import<LocalOrdinal,GlobalOrdinal,Node>;
  public:
    //! The specialization of Map used by this class.
    typedef ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> map_type;

    //! @name Constructor/Destructor Methods
    //@{

    /// \brief Construct a Export object from the source and target Map.
    ///
    /// \param source [in] The source distribution.  This may be a
    ///   multiply owned (overlapping) distribution.
    ///
    /// \param target [in] The target distribution.  This <i>must</i>
    ///   be a uniquely owned (nonoverlapping) distribution.
    Export (const Teuchos::RCP<const map_type>& source,
            const Teuchos::RCP<const map_type>& target);

    /// \brief Construct an Export from the source and target Maps,
    ///   with an output stream for debugging output.
    ///
    /// \param source [in] The source distribution.  This may be a
    ///   multiply owned (overlapping) distribution.
    ///
    /// \param target [in] The target distribution.  This <i>must</i>
    ///   be a uniquely owned (nonoverlapping) distribution.
    ///
    /// \param out [in/out] Output stream for debugging output.
    Export (const Teuchos::RCP<const map_type>& source,
            const Teuchos::RCP<const map_type>& target,
            const Teuchos::RCP<Teuchos::FancyOStream>& out);

    /// \brief Constructor (with list of parameters)
    ///
    /// \param source [in] The source distribution.  This may be a
    ///   multiply owned (overlapping) distribution.
    ///
    /// \param target [in] The target distribution.  This <i>must</i>
    ///   be a uniquely owned (nonoverlapping) distribution.
    ///
    /// \param plist [in/out] List of parameters.  Currently passed
    ///   directly to the Distributor that implements communication.
    ///   If you don't know what this should be, you should use the
    ///   two-argument constructor, listed above.
    Export (const Teuchos::RCP<const map_type>& source,
            const Teuchos::RCP<const map_type>& target,
            const Teuchos::RCP<Teuchos::ParameterList>& plist);

    /// \brief Constructor (with list of parameters and debugging
    ///   output stream)
    ///
    /// \param source [in] The source distribution.  This may be a
    ///   multiply owned (overlapping) distribution.
    ///
    /// \param target [in] The target distribution.  This <i>must</i>
    ///   be a uniquely owned (nonoverlapping) distribution.
    ///
    /// \param out [in/out] Output stream for debugging output.
    ///
    /// \param plist [in/out] List of parameters.  Currently passed
    ///   directly to the Distributor that implements communication.
    ///   If you don't know what this should be, you should use the
    ///   two-argument constructor, listed above.
    Export (const Teuchos::RCP<const map_type>& source,
            const Teuchos::RCP<const map_type>& target,
            const Teuchos::RCP<Teuchos::FancyOStream>& out,
            const Teuchos::RCP<Teuchos::ParameterList>& plist);

    /// \brief Copy constructor.
    ///
    /// \note Currently this only makes a shallow copy of the Export's
    ///   underlying data.
    Export (const Export<LocalOrdinal,GlobalOrdinal,Node>& rhs);

    /// \brief "Copy" constructor from an Export object.
    ///
    /// This constructor creates an Export object from the "reverse"
    /// of the given Import object.  This method is mainly useful for
    /// Tpetra developers, for example when building the explicit
    /// transpose of a sparse matrix.
    Export (const Import<LocalOrdinal,GlobalOrdinal,Node> & importer);

    //! Destructor.
    virtual ~Export ();

    /// \brief Set parameters.
    ///
    /// Please see the class documentation for a list of all accepted
    /// parameters and their default values.
    void setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& plist);

    //@}
    //! @name Export Attribute Methods
    //@{

    /// \brief Number of initial identical IDs.
    ///
    /// The number of IDs that are identical between the source and
    /// target Maps, up to the first different ID.
    size_t getNumSameIDs() const;

    /// \brief Number of IDs to permute but not to communicate.
    ///
    /// The number of IDs that are local to the calling process, but
    /// not part of the first getNumSameIDs() entries.  The Import
    /// will permute these entries locally (without distributed-memory
    /// communication).
    size_t getNumPermuteIDs() const;

    //! List of local IDs in the source Map that are permuted.
    Teuchos::ArrayView<const LocalOrdinal> getPermuteFromLIDs() const;

    //! List of local IDs in the target Map that are permuted.
    Teuchos::ArrayView<const LocalOrdinal> getPermuteToLIDs() const;

    //! Number of entries not on the calling process.
    size_t getNumRemoteIDs() const;

    //! List of entries in the target Map to receive from other processes.
    Teuchos::ArrayView<const LocalOrdinal> getRemoteLIDs() const;

    //! Number of entries that must be sent by the calling process to other processes.
    size_t getNumExportIDs() const;

    //! List of entries in the source Map that will be sent to other processes.
    Teuchos::ArrayView<const LocalOrdinal> getExportLIDs() const;

    /// \brief List of processes to which entries will be sent.
    ///
    /// The entry with local ID getExportLIDs()[i] will be sent to
    /// process getExportPiDs()[i].
    Teuchos::ArrayView<const int> getExportPIDs() const;

    //! The source Map used to construct this Export.
    Teuchos::RCP<const map_type> getSourceMap () const;

    //! The target Map used to construct this Export.
    Teuchos::RCP<const map_type> getTargetMap () const;

    //! The Distributor that this Export object uses to move data.
    Distributor & getDistributor() const;

    /// \brief Do all source Map indices on the calling process exist
    ///   on at least one process (not necessarily this one) in the
    ///   target Map?
    ///
    /// It's not necessarily an error for an Export not to be locally
    /// complete on one or more processes.  Nevertheless, you may find
    /// this predicate useful for figuring out whether you set up your
    /// Maps in the way that you expect.
    bool isLocallyComplete () const;

    //! Assignment operator
    Export<LocalOrdinal,GlobalOrdinal,Node>&
    operator= (const Export<LocalOrdinal,GlobalOrdinal,Node>& rhs);

    //@}
    //! @name I/O Methods
    //@{

    /// \brief Describe this object in a human-readable way to the
    ///   given output stream.
    ///
    /// You must call this method as a collective over all processes
    /// in the communicator of the source and target Map of this
    /// object.
    ///
    /// \param out [out] Output stream to which to write.  Only
    ///   Process 0 in this object's communicator may write to the
    ///   output stream.
    ///
    /// \param verbLevel [in] Verbosity level.  This also controls
    ///   whether this method does any communication.  At verbosity
    ///   levels higher (greater) than Teuchos::VERB_LOW, this method
    ///   behaves as a collective over the object's communicator.
    ///
    /// Teuchos::FancyOStream wraps std::ostream.  It adds features
    /// like tab levels.  If you just want to wrap std::cout, try
    /// this:
    /// \code
    /// auto out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::out));
    /// \endcode
    virtual void
    describe (Teuchos::FancyOStream& out,
              const Teuchos::EVerbosityLevel verbLevel =
                Teuchos::Describable::verbLevel_default) const;

    /// \brief Print the Export's data to the given output stream.
    ///
    /// This method assumes that the given output stream can be
    /// written on all process(es) in the Export's communicator.  The
    /// resulting output is useful mainly for debugging.
    ///
    /// \note This method tries its best (by using barriers at the end
    ///   of each iteration of a for loop over all communicator ranks)
    ///   to ensure ordered deterministic output.  However, the
    ///   assumption that all processes can write to the stream means
    ///   that there are no ordering guarantees other than what the
    ///   operating and run-time system provide.  (MPI synchronization
    ///   may be separate from output stream synchronization, so the
    ///   barriers only improve the chances that output can complete
    ///   before the next process starts writing.)
    virtual void print (std::ostream& os) const;

    //@}

  private:
    //! All the data needed for executing the Export communication plan.
    Teuchos::RCP<ImportExportData<LocalOrdinal,GlobalOrdinal,Node> > ExportData_;
    //! Output stream for debugging output.
    Teuchos::RCP<Teuchos::FancyOStream> out_;
    //! Whether to print copious debugging output on all processes.
    bool debug_;

    //! @name Initialization helper functions (called by the constructor)
    //@{

    //==============================================================================
    // sets up numSameIDs_, numPermuteIDs_, and the export IDs
    // these variables are already initialized to 0 by the ImportExportData ctr.
    // also sets up permuteToLIDs_, permuteFromLIDs_, and exportLIDs_
    void setupSamePermuteExport(Teuchos::Array<GlobalOrdinal> & exportGIDs);
    void setupRemote(Teuchos::Array<GlobalOrdinal> & exportGIDs);
    //@}
  }; // class Export

  /** \brief Non-member constructor for Export objects.

      Creates a Export object from the given source and target maps.
      \pre <tt>src != null</tt>
      \pre <tt>tgt != null</tt>
      \return The Export object. If <tt>src == tgt</tt>, returns \c null.
        (Debug mode: throws std::runtime_error if one of \c src or \c tgt is \c null.)

      \relatesalso Export
    */
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Export<LocalOrdinal, GlobalOrdinal, Node> >
  createExport (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& src,
                const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& tgt)
  {
    if (src == tgt) {
      return Teuchos::null;
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      src == Teuchos::null || tgt == Teuchos::null, std::runtime_error,
      "Tpetra::createExport(): neither source nor target map may be null:"
      << std::endl << "source: " << src << std::endl << "target: " << tgt
      << std::endl);
#endif // HAVE_TPETRA_DEBUG
    return Teuchos::rcp (new Export<LocalOrdinal, GlobalOrdinal, Node> (src, tgt));
  }

} // namespace Tpetra

#endif // TPETRA_EXPORT_DECL_HPP
