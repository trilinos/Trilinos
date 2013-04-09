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
  ///   (possibly) multiply-owned to a uniquely-owned distribution.
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
  template <class LocalOrdinal,
            class GlobalOrdinal = LocalOrdinal,
            class Node = Kokkos::DefaultNode::DefaultNodeType>
  class Export: public Teuchos::Describable {
  public:
    //! The specialization of Map used by this class.
    typedef Map<LocalOrdinal,GlobalOrdinal,Node> map_type;

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

    /// \brief Constructor (with list of parameters).
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

    /// \brief Copy constructor.
    ///
    /// \note Currently this only makes a shallow copy of the Export's
    ///   underlying data.
    Export (const Export<LocalOrdinal,GlobalOrdinal,Node>& rhs);

    //! Destructor.
    ~Export();

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
    /// The entry with local ID getExportLIDs()[i] will be sent to
    /// process getExportImageIDs()[i].
    ArrayView<const int> getExportImageIDs() const;

    //! The source Map used to construct this Export.
    const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getSourceMap() const;

    //! The target Map used to construct this Export.
    const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getTargetMap() const;

    //! The Distributor that this Export object uses to move data.
    Distributor & getDistributor() const;

    //! Assignment operator
    Export<LocalOrdinal,GlobalOrdinal,Node>&
    operator= (const Export<LocalOrdinal,GlobalOrdinal,Node>& rhs);

    //@}
    //! @name I/O Methods
    //@{

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
    RCP<ImportExportData<LocalOrdinal,GlobalOrdinal,Node> > ExportData_;

    //! @name Initialization helper functions (called by the constructor)
    //@{

    //==============================================================================
    // sets up numSameIDs_, numPermuteIDs_, and the export IDs
    // these variables are already initialized to 0 by the ImportExportData ctr.
    // also sets up permuteToLIDs_, permuteFromLIDs_, exportGIDs_, and exportLIDs_
    void setupSamePermuteExport();
    void setupRemote();
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
      src == null || tgt == null, std::runtime_error,
      "Tpetra::createExport(): neither source nor target map may be null:"
      << std::endl << "source: " << src << std::endl << "target: " << tgt
      << std::endl);
#endif // HAVE_TPETRA_DEBUG
    return Teuchos::rcp (new Export<LocalOrdinal, GlobalOrdinal, Node> (src, tgt));
  }
} // namespace Tpetra

#endif // TPETRA_EXPORT_DECL_HPP
