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

#ifndef TPETRA_IMPORT_DECL_HPP
#define TPETRA_IMPORT_DECL_HPP

#include <Tpetra_ConfigDefs.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>

namespace Tpetra {
  //
  // Forward declarations.  The "doxygen" bit simply tells Doxygen
  // (our automatic documentation generation system) to skip forward
  // declarations.
  //
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  class Distributor;

  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  class ImportExportData;

  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  class Map;
#endif // DOXYGEN_SHOULD_SKIP_THIS

  /// \brief Communication plan for data redistribution from a
  ///   uniquely-owned to a (possibly) multiply-owned distribution.
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
  /// As mentioned above, one use case of Import is bringing in remote
  /// source vector data for a distributed sparse matrix-vector
  /// multiply.  The source vector itself is uniquely owned, but must
  /// be brought in into an overlapping distribution so that each
  /// process can compute its part of the target vector without
  /// further communication.
  ///
  /// Epetra separated Import and Export for performance reasons.  The
  /// implementation is different, depending on which direction is the
  /// uniquely-owned Map.  Tpetra retains this convention.
  ///
  /// This class is templated on the same template arguments as Map:
  /// the local ordinal type <tt>LocalOrdinal</tt>, the global ordinal
  /// type <tt>GlobalOrdinal</tt>, and the Kokkos <tt>Node</tt> type.
  template <class LocalOrdinal,
            class GlobalOrdinal = LocalOrdinal,
            class Node = Kokkos::DefaultNode::DefaultNodeType>
  class Import: public Teuchos::Describable {
  public:
    //! The specialization of Map used by this class.
    typedef Map<LocalOrdinal,GlobalOrdinal,Node> map_type;

    //! @name Constructor/Destructor Methods
    //@{

    /// \brief Construct an Import from the source and target Maps.
    ///
    /// \param source [in] The source distribution.  This <i>must</i>
    ///   be a uniquely owned (nonoverlapping) distribution.
    ///
    /// \param target [in] The target distribution.  This may be a
    ///   multiply owned (overlapping) distribution.
    Import (const Teuchos::RCP<const map_type>& source,
            const Teuchos::RCP<const map_type>& target);

    /// \brief Construct an Import from the source and target Maps,
    ///   with an output stream for debugging output.
    ///
    /// \param source [in] The source distribution.  This <i>must</i>
    ///   be a uniquely owned (nonoverlapping) distribution.
    ///
    /// \param target [in] The target distribution.  This may be a
    ///   multiply owned (overlapping) distribution.
    ///
    /// \param out [in/out] Output stream for debugging output.
    Import (const Teuchos::RCP<const map_type>& source,
            const Teuchos::RCP<const map_type>& target,
            const RCP<Teuchos::FancyOStream>& out);

    /// \brief Constructor (with list of parameters)
    ///
    /// \param source [in] The source distribution.  This <i>must</i>
    ///   be a uniquely owned (nonoverlapping) distribution.
    ///
    /// \param target [in] The target distribution.  This may be a
    ///   multiply owned (overlapping) distribution.
    ///
    /// \param plist [in/out] List of parameters.  Currently passed
    ///   directly to the Distributor that implements communication.
    ///   If you don't know what this should be, you should use the
    ///   two-argument constructor, listed above.
    Import (const Teuchos::RCP<const map_type>& source,
            const Teuchos::RCP<const map_type>& target,
            const Teuchos::RCP<Teuchos::ParameterList>& plist);

    /// \brief Constructor (with list of parameters and debug output stream)
    ///
    /// \param source [in] The source distribution.  This <i>must</i>
    ///   be a uniquely owned (nonoverlapping) distribution.
    ///
    /// \param target [in] The target distribution.  This may be a
    ///   multiply owned (overlapping) distribution.
    ///
    /// \param out [in/out] Output stream (for printing copious debug
    ///   output on all processes, if that option is enabled).
    ///
    /// \param plist [in/out] List of parameters.  Currently passed
    ///   directly to the Distributor that implements communication.
    ///   If you don't know what this should be, you should use the
    ///   two-argument constructor, listed above.
    Import (const Teuchos::RCP<const map_type>& source,
            const Teuchos::RCP<const map_type>& target,
            const RCP<Teuchos::FancyOStream>& out,
            const Teuchos::RCP<Teuchos::ParameterList>& plist);

    /// \brief Copy constructor.
    ///
    /// \note Currently this only makes a shallow copy of the Import's
    ///   underlying data.
    Import(const Import<LocalOrdinal,GlobalOrdinal,Node> & import);

    //! Destructor.
    ~Import();

    //@}
    //! @name Import Attribute Methods
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
    ArrayView<const LocalOrdinal> getPermuteFromLIDs() const;

    //! List of local IDs in the target Map that are permuted.
    ArrayView<const LocalOrdinal> getPermuteToLIDs() const;

    //! Number of entries not on the calling process.
    size_t getNumRemoteIDs() const;

    //! List of entries in the target Map to receive from other processes.
    ArrayView<const LocalOrdinal> getRemoteLIDs() const;

    //! Number of entries that must be sent by the calling process to other processes.
    size_t getNumExportIDs() const;

    //! List of entries in the source Map that will be sent to other processes.
    ArrayView<const LocalOrdinal> getExportLIDs() const;

    /// \brief List of processes to which entries will be sent.
    ///
    /// The entry with Local ID <tt>getExportLIDs()[i]</tt> will be
    /// sent to process <tt>getExportImageIDs()[i]</tt>.
    ArrayView<const int> getExportImageIDs() const;

    //! The Source Map used to construct this Import object.
    const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& getSourceMap() const;

    //! The Target Map used to construct this Import object.
    const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& getTargetMap() const;

    //! The Distributor that this Import object uses to move data.
    Distributor & getDistributor() const;

    //! Assignment operator.
    Import<LocalOrdinal,GlobalOrdinal,Node>&
    operator= (const Import<LocalOrdinal,GlobalOrdinal,Node>& Source);

    //@}
    //! @name I/O Methods
    //@{

    /// \brief Print the Import's data to the given output stream.
    ///
    /// This method assumes that the given output stream can be
    /// written on all process(es) in the Import's communicator.  The
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
    //! All the data needed for executing the Import communication plan.
    RCP<ImportExportData<LocalOrdinal,GlobalOrdinal,Node> > ImportData_;
    //! Temporary array used for initialization.
    RCP<Array<GlobalOrdinal> > remoteGIDs_;
    //! Output stream for debug output.
    RCP<Teuchos::FancyOStream> out_;
    //! Whether to print copious debug output on each process.
    bool debug_;

    //! @name Initialization helper functions (called by the constructor)
    //@{

    //==============================================================================
    // sets up numSameIDs_, numPermuteIDs_, and numRemoteIDs_
    // these variables are already initialized to 0 by the ImportExportData ctr.
    // also sets up permuteToLIDs_, permuteFromLIDs_, and remoteLIDs_

    /// \brief Compute the necessary receives for the Import.
    ///
    /// This routine fills in the following fields of ImportData_:
    ///
    ///   - numSameIDs_ (the number of consecutive initial GIDs owned
    ///     by both the source and target Maps)
    ///   - permuteToLIDs_ (for each of the remaining GIDs g in the
    ///     target Map, if the source Map also owns g, then
    ///     permuteToLIDs_ gets the corresponding LID in the target,
    ///     and permuteFromLIDs_ gets the corresponding LID in the
    ///     source)
    ///   - permuteFromLIDs_ (see permuteToLIDs_)
    ///   - remoteLIDs_ (the LIDs of the GIDs that are owned by the
    ///     target Map, but not by the source Map)
    ///
    /// It also fills in the temporary remoteGIDs_ array with the GIDs
    /// that are owned by the target Map but not by the source Map.
    ///
    /// The name for this routine comes from what it does.  It first
    /// finds the GIDs that are the same (representing elements which
    /// require neither communication nor permutation).  Then it finds
    /// permutation IDs (which require permutation, but no
    /// communication, because they are in a possibly different order
    /// in the source and target Maps, but owned by the same process)
    /// and remote IDs (which require communication, because they are
    /// owned by the target Map but not by the source Map).
    ///
    /// This routine does not communicate, except perhaps for the
    /// TPETRA_ABUSE_WARNING (that is only triggered if there are
    /// remote IDs but the source is not distributed).
    void setupSamePermuteRemote();

    /// \brief Compute the send communication plan from the receives.
    ///
    /// This routine is called after setupSamePermuteRemote(), if the
    /// source Map is distributed.  It uses the <tt>remoteGIDs_</tt>
    /// temporary array that was allocated by that routine.  After
    /// this routine completes, the <tt>remoteGIDs_</tt> array is no
    /// longer needed.
    ///
    /// Algorithm:
    ///
    /// 1. Identify which GIDs are in the target Map but not in the
    ///    source Map.  These correspond to required receives.  Store
    ///    them for now in <tt>remoteGIDs_</tt>.  Find the process IDs
    ///    of the remote GIDs to receive.
    ///
    /// 2. Invoke Distributor's createFromRecvs() using the above
    ///    remote GIDs and remote process IDs as input.  This sets up
    ///    the Distributor and computes the send GIDs and process IDs.
    ///
    /// 3. Use the source Map to compute the send LIDs from the send
    ///    GIDs.
    ///
    /// This routine fills in the <tt>remoteLIDs_</tt> field of
    /// <tt>ImportData_</tt>.
    void setupExport();
    //@}
  }; // class Import

  /** \brief Nonmember constructor for Import.

      Create a Import object from the given source and target Maps.
      \pre <tt>src != null</tt>
      \pre <tt>tgt != null</tt>
      \return The Import object. If <tt>src == tgt</tt>, returns \c null.
        (Debug mode: throws std::runtime_error if one of \c src or \c tgt is \c null.)

      \relatesalso Import
    */
  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> >
  createImport (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& src,
                const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& tgt)
  {
    if (src == tgt) {
      return Teuchos::null;
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      src == null || tgt == null, std::runtime_error,
      "Tpetra::createImport(): neither source nor target map may be null:"
      << std::endl << "source: " << src << std::endl << "target: " << tgt
      << std::endl);
#endif // HAVE_TPETRA_DEBUG
    return Teuchos::rcp (new Import<LocalOrdinal, GlobalOrdinal, Node> (src, tgt));
  }

  /** \brief Nonmember constructor for Import that takes a ParameterList.

      Create a Import object from the given source and target Maps,
      using the given list of parameters.
      \pre <tt>src != null</tt>
      \pre <tt>tgt != null</tt>
      \return The Import object. If <tt>src == tgt</tt>, returns \c null.
        (Debug mode: throws std::runtime_error if one of \c src or \c tgt is \c null.)

      \relatesalso Import
    */
  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> >
  createImport (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& src,
                const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& tgt,
                const Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    if (src == tgt) {
      return Teuchos::null;
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      src == null || tgt == null, std::runtime_error,
      "Tpetra::createImport(): neither source nor target map may be null:"
      << std::endl << "source: " << src << std::endl << "target: " << tgt
      << std::endl);
#endif // HAVE_TPETRA_DEBUG
    typedef Import<LocalOrdinal, GlobalOrdinal, Node> import_type;
    return Teuchos::rcp (new import_type (src, tgt, plist));
  }
} // namespace Tpetra

#endif // TPETRA_IMPORT_DECL_HPP
