// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_EXPORT_DECL_HPP
#define TPETRA_EXPORT_DECL_HPP

#include "Tpetra_Details_Transfer.hpp"
#include "Tpetra_Export_fwd.hpp"
#include "Tpetra_Import_fwd.hpp"
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
  /// the default values are fine.  
  ///
  template<class LocalOrdinal,
           class GlobalOrdinal,
           class Node>
  class Export:
    public ::Tpetra::Details::Transfer<LocalOrdinal, GlobalOrdinal, Node>
  {
  private:
    friend class Import<LocalOrdinal,GlobalOrdinal,Node>;
    using base_type =
      ::Tpetra::Details::Transfer<LocalOrdinal, GlobalOrdinal, Node>;
  public:
    //! The specialization of Map used by this class.
    using map_type = ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;

    //! @name Constructors, assignment, and destructor
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

    //! Assignment operator
    Export<LocalOrdinal,GlobalOrdinal,Node>&
    operator= (const Export<LocalOrdinal,GlobalOrdinal,Node>& rhs) = default;

    //! Destructor.
    virtual ~Export () = default;

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
    //! @name Initialization helper functions (called by the constructor)
    //@{
    //! Set up same, permute, and export IDs.
    void setupSamePermuteExport(Teuchos::Array<GlobalOrdinal> & exportGIDs);

    //! Set up remote IDs.
    void setupRemote(Teuchos::Array<GlobalOrdinal> & exportGIDs);
    //@}
  }; // class Export

  /// \brief Nonmember "constructor" for Export objects.
  ///
  /// Create an Export object from the given source and target Maps.
  ///
  /// \pre <tt>src != null</tt>
  /// \pre <tt>tgt != null</tt>
  ///
  /// \return The Export object. If <tt>src == tgt</tt>, returns \c null.
  ///
  /// (Debug mode: throws std::runtime_error if one of \c src or \c tgt is \c null.)
  ///
  /// \relatesalso Export
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Export<LocalOrdinal, GlobalOrdinal, Node> >
  createExport (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& src,
                const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& tgt)
  {
    if (src == tgt) {
      return Teuchos::null;
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION
      (src == Teuchos::null || tgt == Teuchos::null, std::runtime_error,
       "Tpetra::createExport: Neither source nor target map may be null.");
#endif // HAVE_TPETRA_DEBUG
    using export_type = Export<LocalOrdinal, GlobalOrdinal, Node>;
    return Teuchos::rcp (new export_type (src, tgt));
  }

} // namespace Tpetra

#endif // TPETRA_EXPORT_DECL_HPP
