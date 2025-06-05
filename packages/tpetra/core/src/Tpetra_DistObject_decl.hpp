// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// clang-format off
#ifndef TPETRA_DISTOBJECT_DECL_HPP
#define TPETRA_DISTOBJECT_DECL_HPP

/// \file Tpetra_DistObject_decl.hpp
/// \brief Declaration of the Tpetra::DistObject class

#include "Tpetra_Details_DistributorActor.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_SrcDistObject.hpp"
#include "Tpetra_DistObject_fwd.hpp"
#include "Kokkos_ArithTraits.hpp"
#include <memory>
#include <type_traits>

// #ifndef HAVE_TPETRA_TRANSFER_TIMERS
// #  define HAVE_TPETRA_TRANSFER_TIMERS 1
// #endif // HAVE_TPETRA_TRANSFER_TIMERS

#ifdef HAVE_TPETRA_TRANSFER_TIMERS
#  undef HAVE_TPETRA_TRANSFER_TIMERS
#endif // HAVE_TPETRA_TRANSFER_TIMERS

namespace KokkosClassic {
  /// \brief Read/write options for non-const views.
  ///
  /// \warning This is NOT for users!  This only exists for backwards
  ///   compatibility.
  enum ReadWriteOption {
    ReadWrite = 0, /*!< Indicates that the view may be safely read and written. */
    OverwriteAll = 1 /*!< Indicates that the contents of the view are undefined until set on the host. */
  };
} // namespace KokkosClassic

namespace Tpetra {

  /// \brief Remove processes which contain no elements in this object's Map.
  ///
  /// \tparam DistObjectType A specialization of DistObject.
  ///
  /// \warning This method is ONLY for use by experts.  The fact that
  ///   the documentation of this method starts with a "Vocabulary"
  ///   section should give you proper respect for the complicated
  ///   semantics of this method in a parallel MPI run.
  /// \warning We make NO promises of backwards compatibility.
  ///   This method may change or disappear at any time.
  ///
  /// Vocabulary:
  /// - The Map returned by <tt>input->getMap() on input to this
  ///   method is the "original Map."
  /// - The communicator returned by <tt>input->getComm() on
  ///   input to this method is the "original communicator."
  /// - All processes in the original communicator which contain zero
  ///   elements in the original Map are "excluded processes."
  /// - All other processes in the original communicator are "included
  ///   processes."
  ///
  /// Preconditions:
  /// - The nonnull object \c input is distributed over the
  ///   original Map.
  /// - The input Map <tt>newMap</tt> <i>must</i> be the same as the
  ///   result of calling removeEmptyProcesses() on the original Map.
  /// - On excluded processes, <tt>newMap</tt> must be
  ///   <tt>Teuchos::null</tt>.  (This is what
  ///   <tt>getMap()->removeEmptyProcesses()</tt> returns anyway on
  ///   excluded processes.)
  ///
  /// This method has collective semantics over the original
  /// communicator.  On included processes, reassign this object's Map
  /// (that would be returned by getMap()) to the input \c newMap, and
  /// do any work that needs to be done to restore correct semantics.
  /// The input DistObject \c input will be nonnull on return.  On
  /// excluded processes, free any data in \c input that need freeing,
  /// do any other work that needs to be done to restore correct
  /// semantics, and set \c input to null before returning.
  ///
  /// The two-argument version of this function is useful if you have
  /// already precomputed the new Map that excludes processes with
  /// zero elements.  For example, you might want to apply this Map to
  /// several different MultiVector instances.  The one-argument
  /// version of this function is useful if you want the DistObject to
  /// compute the new Map itself, because you only plan to use it for
  /// that one DistObject instance.
  ///
  /// Here is a sample use case.  Suppose that \c input is some
  /// subclass of DistObject, like MultiVector, CrsGraph, or
  /// CrsMatrix.  Suppose also that \c map_type is the corresponding
  /// specialization of Map.
  /// \code
  /// RCP<const map_type> origRowMap = input->getMap ();
  /// RCP<const map_type> newRowMap = origRowMap->removeEmptyProcesses ();
  /// removeEmptyProcessesInPlace (input, newRowMap);
  /// // Either (both the new Map and input are null), or
  /// // (both the new Map and input are not null).
  /// assert ((newRowMap.is_null () && input.is_null ()) ||
  ///         (! newRowMap.is_null () && ! input.is_null ()));
  /// \endcode
  ///
  /// \warning On excluded processes, calling this function
  ///   invalidates any other references to the input DistObject
  ///   <tt>input</tt>.  Calling any methods (other than the
  ///   destructor) on the input on excluded processes has undefined
  ///   behavior in that case, and may result in deadlock.
  ///
  /// \note The name differs from Map's method
  ///   removeEmptyProcesses(), in order to emphasize that the
  ///   operation on DistObject happens in place, modifying the
  ///   input, whereas the operation removeEmptyProcess() on Map
  ///   does not modify the input.
  ///
  /// \note To implementers of DistObject subclasses: The default
  ///   implementation of this class throws std::logic_error.
  ///
  /// \note To implementers of DistObject subclasses: On exit, the
  ///   only method of this object which is safe to call on excluded
  ///   processes is the destructor, or this method with the original
  ///   Map.  This implies that subclasses' destructors must not
  ///   contain communication operations.
  template<class DistObjectType>
  void
  removeEmptyProcessesInPlace (Teuchos::RCP<DistObjectType>& input,
                               const Teuchos::RCP<const Map<typename DistObjectType::local_ordinal_type,
                                                            typename DistObjectType::global_ordinal_type,
                                                            typename DistObjectType::node_type> >& newMap);

  /// \brief Remove processes which contain no elements in this object's Map.
  ///
  /// \tparam DistObjectType A specialization of DistObject.
  ///
  /// \warning This method is ONLY for use by experts.
  /// \warning We make NO promises of backwards compatibility.
  ///   This method may change or disappear at any time.
  ///
  /// This method behaves just like the two-argument version of
  /// removeEmptyProcessesInPlace(), except that it first calls
  /// removeEmptyProcesses() on the input DistObject's Map to compute
  /// the new Map.
  ///
  /// The two-argument version of this function is useful if you have
  /// already precomputed the new Map that excludes processes with
  /// zero elements.  For example, you might want to apply this Map to
  /// several different MultiVector instances.  The one-argument
  /// version of this function is useful if you want the DistObject to
  /// compute the new Map itself, because you only plan to use it for
  /// that one DistObject instance.
  ///
  /// Here is a sample use case.  Suppose that \c input is some
  /// subclass of DistObject, like MultiVector, CrsGraph, or
  /// CrsMatrix.  Suppose also that \c map_type is the corresponding
  /// specialization of Map.
  /// \code
  /// removeEmptyProcessesInPlace (input);
  /// RCP<const map_type> newRowMap;
  /// if (! input.is_null ()) {
  ///   newRowMap = input->getMap ();
  /// }
  /// // Either (both the new Map and input are null), or
  /// // (both the new Map and input are not null).
  /// assert ((newRowMap.is_null () && input.is_null ()) ||
  ///         (! newRowMap.is_null () && ! input.is_null ()));
  /// \endcode
  template<class DistObjectType>
  void
  removeEmptyProcessesInPlace (Teuchos::RCP<DistObjectType>& input);

  /// \class DistObject
  /// \brief Base class for distributed Tpetra objects that support
  ///   data redistribution.
  ///
  /// DistObject is a base class for all Tpetra distributed global
  /// objects, including CrsGraph, CrsMatrix, MultiVector, and Vector.
  /// Classes that inherit from DistObject work as either the source
  /// or the target of an Export or Import (parallel redistribution)
  /// operation.  If you want to write a class that can work as the
  /// source or target of an Export or Import operation, that class
  /// must inherit from DistObject.
  ///
  /// \tparam LocalOrdinal The type of local indices.  Same as Map's
  ///   \c LocalOrdinal template parameter.  See Map's documentation
  ///   for a discussion of the types that are valid to use here.
  ///
  /// \tparam GlobalOrdinal The type of global indices.  Same as Map's
  ///   \c GlobalOrdinal template parameter.  See Map's documentation
  ///   for a discussion of the types that are valid to use here.
  ///
  /// \tparam Node Same as Map's \c Node template parameter.  See
  ///   Map's documentation for a discussion of the types that are
  ///   valid to use here.
  ///
  /// \section Tpetra_DistObject_Summary Summary
  ///
  /// Most Tpetra users will only use this class' methods to perform
  /// data redistribution for subclasses such as CrsGraph, CrsMatrix,
  /// MultiVector, and Vector.  DistObject provides four methods for
  /// redistributing data: two versions of <tt>doImport()</tt>, and
  /// two versions of <tt>doExport()</tt>.  Import operations
  /// redistribute data from a nonoverlapping (one-to-one)
  /// distribution to a possibly overlapping distribution.  Export
  /// operations redistribute data from a possibly overlapping
  /// distribution to a nonoverlapping (one-to-one) distribution.
  /// Once you have precomputed a data redistribution plan (an Import
  /// or Export object), you may use the plan to redistribute an input
  /// object's data into this object, by calling one of these methods.
  /// The input object of <tt>doImport()</tt> or <tt>doExport()</tt>
  /// is always the "source" of the redistribution operation, which
  /// sends the data.  The <tt>*this</tt> object is the target, which
  /// receives and combines the data.  It has the distribution given
  /// by <tt>this->getMap()</tt>.
  ///
  /// \section Tpetra_DistObject_FwdRev Forward or reverse redistribution modes
  ///
  /// Both Import and Export operations occur in two modes: forward
  /// and reverse.  Forward mode is the usual case, where you are
  /// calling a method with its matching plan type
  /// (<tt>doImport()</tt> for an Import plan, or <tt>doExport()</tt>
  /// for an Export plan).  In that case, the input DistObject must
  /// have the same Map as the source Map of the plan, and the target
  /// DistObject must have the same Map as the target Map of the plan.
  /// Reverse mode is also possible, where you call a method with the
  /// opposite plan type (<tt>doImport()</tt> for an Export plan, or
  /// <tt>doExport()</tt> for an Import plan).  In that case, the
  /// source DistObject's Map must be the same as the target Map of
  /// the plan, and the target DistObject's Map must be the same as
  /// the source Map of the plan.  If you call <tt>doImport()</tt>, we
  /// still call this an Import operation, even if you are using an
  /// Export plan in reverse.  Similarly, if you call
  /// <tt>doExport()</tt>, we call this an Export operation.
  ///
  /// Most users will want to use forward mode.  However, reverse mode
  /// is useful for some applications.  For example, suppose you are
  /// solving a nonlinear partial differential equation using the
  /// finite element method, with Newton's method for the nonlinear
  /// equation.  When assembling into a vector, it is convenient and
  /// efficient to do local assembly first into a vector with an
  /// overlapping distribution, then do global assembly via forward
  /// mode Export into a vector with a nonoverlapping distribution.
  /// After the linear solve, you may want to bring the resulting
  /// nonoverlapping distribution vector back to the overlapping
  /// distribution for another update phase.  This would be a reverse
  /// mode Import, using the precomputed Export object.
  ///
  /// Another use case for reverse mode is in CrsMatrix, for the
  /// transpose version of distributed sparse matrix-vector multiply
  /// ("mat-vec").  Non-transpose mat-vec (a function from the domain
  /// Map to the range Map) does an Import to bring in the source
  /// vector's data from the domain Map to the column Map of the
  /// sparse matrix, and an Export (if necessary) to bring the results
  /// from the row Map of the sparse matrix to the range Map.
  /// Transpose mat-vec (a function from the range Map to the domain
  /// Map) uses these precomputed Import and Export objects in reverse
  /// mode: first the Export in reverse mode to Import the source
  /// vector's data to the row Map, and then the Import in reverse
  /// mode to Export the results to the domain Map.  Reverse mode lets
  /// us reuse the precomputed data redistribution plans for the
  /// transpose case.
  ///
  /// \section Tpetra_DistObject_ImplSubclass How to implement a subclass
  ///
  /// If you want to implement your own DistObject subclass, you
  /// <i>must</i> implement at least the following methods:
  /// <ul>
  /// <li> checkSizes() </li>
  /// <li> copyAndPermute() </li>
  /// <li> packAndPrepare() </li>
  /// <li> unpackAndCombine() </li>
  /// </ul>
  /// You <i>may</i> also implement constantNumberOfPackets(), if
  /// appropriate.
  ///
  /// There is also an "old" interface, which is deprecated (see
  /// GitHub Issue #4853) and will be removed soon.  Do not implement
  /// the "old" interface.
  ///
  /// DistObject implements SrcDistObject, because we presume that if
  /// an object can be the target of an Import or Export, it can also
  /// be the source of an Import or Export.
  template <class Packet,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  class DistObject :
    virtual public SrcDistObject,
    virtual public Teuchos::Describable
  {
  public:
    //! @name Typedefs
    //@{

    /// \brief The type of each datum being sent or received in an Import or Export.
    ///
    /// Note that this type does not always correspond to the
    /// <tt>Scalar</tt> template parameter of subclasses.
    using packet_type = typename ::Kokkos::ArithTraits<Packet>::val_type;
    //! The type of local indices.
    using local_ordinal_type = LocalOrdinal;
    //! The type of global indices.
    using global_ordinal_type = GlobalOrdinal;
    //! The Node type.  If you don't know what this is, don't use it.
    using node_type = Node;

    //! The Kokkos Device type.
    using device_type = typename Node::device_type;
    //! The Kokkos execution space.
    using execution_space = typename device_type::execution_space;

    //! The type of the Map specialization to use with this class.
    using map_type = Map<local_ordinal_type, global_ordinal_type, node_type>;

    //@}
    //! @name Constructors and destructor
    //@{

    /// \brief Constructor
    ///
    /// \param map [in] Map over which the object is distributed.
    explicit DistObject (const Teuchos::RCP<const map_type>& map);

    //! Copy constructor (default).
    DistObject (const DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>&) = default;

    //! Assignment operator (default).
    DistObject& operator= (const DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>&) = default;

    //! Move constructor (default).
    DistObject (DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>&&) = default;

    //! Move assignment (default).
    DistObject& operator= (DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>&&) = default;

    /// \brief Destructor (virtual for memory safety of derived classes).
    ///
    /// \note To Tpetra developers: See the C++ Core Guidelines C.21
    ///   ("If you define or <tt>=delete</tt> any default operation,
    ///   define or <tt>=delete</tt> them all"), in particular the
    ///   AbstractBase example, for why this destructor declaration
    ///   implies that we need the above four <tt>=default</tt>
    ///   declarations for copy construction, move construction, copy
    ///   assignment, and move assignment.
    virtual ~DistObject () = default;

    //@}
    //! @name Public methods for redistributing data
    //@{

    /// \brief Import data into this object using an Import object
    ///   ("forward mode").
    ///
    /// The input DistObject is always the source of the data
    /// redistribution operation, and the <tt>*this</tt> object is
    /// always the target.
    ///
    /// If you don't know the difference between forward and reverse
    /// mode, then you probably want forward mode.  Use this method
    /// with your precomputed Import object if you want to do an
    /// Import, else use doExport() with a precomputed Export object.
    ///
    /// "Restricted Mode" does two things:
    /// <ol>
    /// <li> Skips copyAndPermute </li>
    /// <li> Allows the "target" Map of the transfer to be a subset of
    ///      the Map of <tt>*this</tt>, in a "locallyFitted" sense. </li>
    /// </ol>
    /// This cannot be used if (2) is not true, OR there are permutes.
    /// The "source" maps still need to match.
    ///
    /// \param source [in] The "source" object for redistribution.
    /// \param importer [in] Precomputed data redistribution plan.
    ///   Its source Map must be the same as the input DistObject's Map,
    ///   and its target Map must be the same as <tt>this->getMap()</tt>.
    /// \param CM [in] How to combine incoming data with the same
    ///   global index.
    void
    doImport (const SrcDistObject& source,
              const Import<LocalOrdinal, GlobalOrdinal, Node>& importer,
              const CombineMode CM,
              const bool restrictedMode = false);

    /// \brief Export data into this object using an Export object
    ///   ("forward mode").
    ///
    /// The input DistObject is always the source of the data
    /// redistribution operation, and the <tt>*this</tt> object is
    /// always the target.
    ///
    /// If you don't know the difference between forward and reverse
    /// mode, then you probably want forward mode.  Use this method
    /// with your precomputed Export object if you want to do an
    /// Export, else use doImport() with a precomputed Import object.
    ///
    /// "Restricted Mode" does two things:
    /// <ol>
    /// <li> Skips copyAndPermute </li>
    /// <li> Allows the "target" Map of the transfer to be a subset of
    ///      the Map of <tt>*this</tt>, in a "locallyFitted" sense. </li>
    /// </ol>
    /// This cannot be used if (2) is not true, OR there are permutes.
    /// The "source" maps still need to match.
    ///
    /// \param source [in] The "source" object for redistribution.
    /// \param exporter [in] Precomputed data redistribution plan.
    ///   Its source Map must be the same as the input DistObject's Map,
    ///   and its target Map must be the same as <tt>this->getMap()</tt>.
    /// \param CM [in] How to combine incoming data with the same
    ///   global index.
    void
    doExport (const SrcDistObject& source,
              const Export<LocalOrdinal, GlobalOrdinal, Node>& exporter,
              const CombineMode CM,
              const bool restrictedMode = false);

    /// \brief Import data into this object using an Export object
    ///   ("reverse mode").
    ///
    /// The input DistObject is always the source of the data
    /// redistribution operation, and the <tt>*this</tt> object is
    /// always the target.
    ///
    /// If you don't know the difference between forward and reverse
    /// mode, then you probably want forward mode.  Use the version of
    /// doImport() that takes a precomputed Import object in that
    /// case.
    ///
    /// "Restricted Mode" does two things:
    /// <ol>
    /// <li> Skips copyAndPermute </li>
    /// <li> Allows the "target" Map of the transfer to be a subset of
    ///      the Map of <tt>*this</tt>, in a "locallyFitted" sense. </li>
    /// </ol>
    /// This cannot be used if (2) is not true, OR there are permutes.
    /// The "source" maps still need to match.
    ///
    /// \param source [in] The "source" object for redistribution.
    /// \param exporter [in] Precomputed data redistribution plan.
    ///   Its <i>target</i> Map must be the same as the input DistObject's Map,
    ///   and its <i>source</i> Map must be the same as <tt>this->getMap()</tt>.
    ///   (Note the difference from forward mode.)
    /// \param CM [in] How to combine incoming data with the same
    ///   global index.
    void
    doImport (const SrcDistObject& source,
              const Export<LocalOrdinal, GlobalOrdinal, Node>& exporter,
              const CombineMode CM,
              const bool restrictedMode = false);

    /// \brief Export data into this object using an Import object
    ///   ("reverse mode").
    ///
    /// The input DistObject is always the source of the data
    /// redistribution operation, and the <tt>*this</tt> object is
    /// always the target.
    ///
    /// If you don't know the difference between forward and reverse
    /// mode, then you probably want forward mode.  Use the version of
    /// doExport() that takes a precomputed Export object in that
    /// case.
    ///
    /// "Restricted Mode" does two things:
    /// <ol>
    /// <li> Skips copyAndPermute </li>
    /// <li> Allows the "target" Map of the transfer to be a subset of
    ///      the Map of <tt>*this</tt>, in a "locallyFitted" sense. </li>
    /// </ol>
    /// This cannot be used if (2) is not true, OR there are permutes.
    /// The "source" maps still need to match.
    ///
    /// \param source [in] The "source" object for redistribution.
    /// \param importer [in] Precomputed data redistribution plan.
    ///   Its <i>target</i> Map must be the same as the input DistObject's Map,
    ///   and its <i>source</i> Map must be the same as <tt>this->getMap()</tt>.
    ///   (Note the difference from forward mode.)
    /// \param CM [in] How to combine incoming data with the same
    ///   global index.
    void
    doExport (const SrcDistObject& source,
              const Import<LocalOrdinal, GlobalOrdinal, Node>& importer,
              const CombineMode CM,
              const bool restrictedMode = false);

    void
    beginImport(const SrcDistObject& source,
                const Import<LocalOrdinal, GlobalOrdinal, Node>& importer,
                const CombineMode CM,
                const bool restrictedMode = false);

    void
    beginExport(const SrcDistObject& source,
                const Export<LocalOrdinal, GlobalOrdinal, Node>& exporter,
                const CombineMode CM,
                const bool restrictedMode = false);

    void
    beginImport(const SrcDistObject& source,
                const Export<LocalOrdinal, GlobalOrdinal, Node>& exporter,
                const CombineMode CM,
                const bool restrictedMode = false);

    void
    beginExport(const SrcDistObject& source,
                const Import<LocalOrdinal, GlobalOrdinal, Node>& importer,
                const CombineMode CM,
                const bool restrictedMode = false);

    void
    endImport(const SrcDistObject& source,
              const Import<LocalOrdinal, GlobalOrdinal, Node>& importer,
              const CombineMode CM,
              const bool restrictedMode = false);

    void
    endExport(const SrcDistObject& source,
              const Export<LocalOrdinal, GlobalOrdinal, Node>& exporter,
              const CombineMode CM,
              const bool restrictedMode = false);

    void
    endImport(const SrcDistObject& source,
              const Export<LocalOrdinal, GlobalOrdinal, Node>& exporter,
              const CombineMode CM,
              const bool restrictedMode = false);

    void
    endExport(const SrcDistObject& source,
              const Import<LocalOrdinal, GlobalOrdinal, Node>& importer,
              const CombineMode CM,
              const bool restrictedMode = false);

    /// \brief Whether the data from an import/export operation has
    ///        arrived, and is ready for the unpack and combine step.
    bool transferArrived() const;

    //@}
    //! @name Attribute accessor methods
    //@{

    /// \brief Whether this is a globally distributed object.
    ///
    /// For a definition of "globally distributed" (and its opposite,
    /// "locally replicated"), see the documentation of Map's
    /// isDistributed() method.
    bool isDistributed () const;

    /// \brief The Map describing the parallel distribution of this object.
    ///
    /// Note that some Tpetra objects might be distributed using
    /// multiple Map objects.  For example, CrsMatrix has both a row
    /// Map and a column Map.  It is up to the subclass to decide
    /// which Map to use when invoking the DistObject constructor.
    virtual Teuchos::RCP<const map_type> getMap () const { return map_; }

    //@}
    //! @name I/O methods
    //@{

    /// \brief Print this object to the given output stream.
    ///
    /// We generally assume that all MPI processes can print to the
    /// given stream.
    void print (std::ostream& os) const;

    //@}
    //! @name Implementation of Teuchos::Describable
    //@{

    /// \brief One-line descriptiion of this object.
    ///
    /// We declare this method virtual so that subclasses of
    /// DistObject may override it.
    virtual std::string description () const;

    /// \brief Print a descriptiion of this object to the given output stream.
    ///
    /// We declare this method virtual so that subclasses of
    /// Distobject may override it.
    virtual void
    describe (Teuchos::FancyOStream &out,
              const Teuchos::EVerbosityLevel verbLevel =
              Teuchos::Describable::verbLevel_default) const;

    //@}
    //! @name Methods for use only by experts
    //@{

    /// \brief Remove processes which contain no entries in this
    ///   object's Map.
    ///
    /// \warning This method is ONLY for use by experts.  We highly
    ///   recommend using the nonmember function of the same name
    ///   defined in this file.
    ///
    /// \warning We make NO promises of backwards compatibility.
    ///   This method may change or disappear at any time.
    ///
    /// On input, this object is distributed over the Map returned by
    /// getMap() (the "original Map," with its communicator, the
    /// "original communicator").  The input \c newMap of this method
    /// <i>must</i> be the same as the result of calling
    /// <tt>getMap()->removeEmptyProcesses()</tt>.  On processes in
    /// the original communicator which contain zero entries
    /// ("excluded processes," as opposed to "included processes"),
    /// the input \c newMap must be \c Teuchos::null (which is what
    /// <tt>getMap()->removeEmptyProcesses()</tt> returns anyway).
    ///
    /// On included processes, reassign this object's Map (that would
    /// be returned by getMap()) to the input \c newMap, and do any
    /// work that needs to be done to restore correct semantics.  On
    /// excluded processes, free any data that needs freeing, and do
    /// any other work that needs to be done to restore correct
    /// semantics.
    ///
    /// This method has collective semantics over the original
    /// communicator.  On exit, the only method of this object which
    /// is safe to call on excluded processes is the destructor.  This
    /// implies that subclasses' destructors must not contain
    /// communication operations.
    ///
    /// \return The object's new Map.  Its communicator is a new
    ///   communicator, distinct from the old Map's communicator,
    ///   which contains a subset of the processes in the old
    ///   communicator.
    ///
    /// \note The name differs from Map's method
    ///   removeEmptyProcesses(), in order to emphasize that the
    ///   operation on DistObject happens in place, modifying the
    ///   input, whereas the operation removeEmptyProcess() on Map
    ///   does not modify the input.
    ///
    /// \note To implementers of DistObject subclasses: The default
    ///   implementation of this class throws std::logic_error.
    virtual void
    removeEmptyProcessesInPlace (const Teuchos::RCP<const map_type>& newMap);

    //@}

  protected:
    /// \enum ReverseOption
    /// \brief Whether the data transfer should be performed in
    ///   forward or reverse mode.
    ///
    /// "Reverse mode" means calling doExport() with an Import object,
    /// or calling doImport() with an Export object.  "Forward mode"
    /// means calling doExport() with an Export object, or calling
    /// doImport() with an Import object.
    enum ReverseOption {
      DoForward, //*!< Perform the transfer in forward mode.
      DoReverse  //*!< Perform the transfer in reverse mode.
    };

    /// \brief Whether the implementation's instance promises always
    ///   to have a constant number of packets per LID (local index),
    ///   and if so, how many packets per LID there are.
    ///
    /// If this method returns zero, the instance says that it might
    /// possibly have a different number of packets for each LID
    /// (local index) to send or receive.  If it returns nonzero, the
    /// instance promises that the number of packets is the same for
    /// all LIDs, and that the return value is this number of packets
    /// per LID.
    ///
    /// The default implementation of this method returns zero.  This
    /// does not affect the behavior of doTransfer() in any way.  If a
    /// nondefault implementation returns nonzero, doTransfer() will
    /// use this information to avoid unnecessary allocation and / or
    /// resizing of arrays.
    virtual size_t constantNumberOfPackets () const;

    /// \brief Redistribute data across (MPI) processes.
    ///
    /// \param src [in] The source object, to redistribute into
    ///   the target object, which is <tt>*this</tt> object.
    ///
    /// \param transfer [in] The Export or Import object representing
    ///   the communication pattern.  (Details::Transfer is the common
    ///   base class of these two objects.)
    ///
    /// \param modeString [in] Human-readable string, for verbose
    ///   debugging output and error output, explaining what function
    ///   called this method.  Example: "doImport (forward)",
    ///   "doExport (reverse)".
    ///
    /// \param revOp [in] Whether to do a forward or reverse mode
    ///   redistribution.
    ///
    /// \param CM [in] The combine mode that describes how to combine
    ///   values that map to the same global ID on the same process.
    virtual void
    doTransfer (const SrcDistObject& src,
                const ::Tpetra::Details::Transfer<local_ordinal_type, global_ordinal_type, node_type>& transfer,
                const char modeString[],
                const ReverseOption revOp,
                const CombineMode CM,
                const bool restrictedMode);

    /// \brief Reallocate numExportPacketsPerLID_ and/or
    ///   numImportPacketsPerLID_, if necessary.
    ///
    /// \param numExportLIDs [in] Number of entries in the exportLIDs
    ///   input array argument of doTransfer().
    ///
    /// \param numImportLIDs [in] Number of entries in the remoteLIDs
    ///   input array argument of doTransfer().
    ///
    /// \return Whether we actually reallocated either of the arrays.
    ///
    /// \warning This is an implementation detail of doTransferNew().
    ///   This needs to be protected, but that doesn't mean users
    ///   should call this method.
    virtual bool
    reallocArraysForNumPacketsPerLid (const size_t numExportLIDs,
                                      const size_t numImportLIDs);


    /// \typedef buffer_memory_space
    /// \brief Kokkos memory space for communication buffers.
    using buffer_memory_space =
      ::Tpetra::Details::DefaultTypes::comm_buffer_memory_space<device_type>;

  public:
    /// \typedef buffer_device_type
    /// \brief Kokkos::Device specialization for communication buffers.
    ///
    /// See #1088 for why this is not just \c device_type.
    ///
    /// This needs to be public so that I can declare functions like
    /// packAndPrepareWithOwningPIDs.
    ///
    /// \warning This is an implementation detail.  DO NOT DEPEND ON
    ///   IT.  It may disappear or change at any time.
    using buffer_device_type =
      Kokkos::Device<typename device_type::execution_space,
                     buffer_memory_space>;
  protected:
    /// \brief Implementation detail of doTransfer.
    ///
    /// LID DualViews come from the Transfer object given to
    /// doTransfer.  They are <i>always</i> sync'd on both host and
    /// device.  Users must never attempt to modify or sync them.
    void beginTransfer(const SrcDistObject& src,
                       const ::Tpetra::Details::Transfer<local_ordinal_type, global_ordinal_type, node_type>& transfer,
                       const char modeString[],
                       const ReverseOption revOp,
                       const CombineMode CM,
                       const bool restrictedMode);

    void endTransfer(const SrcDistObject& src,
                     const ::Tpetra::Details::Transfer<local_ordinal_type, global_ordinal_type, node_type>& transfer,
                     const char modeString[],
                     const ReverseOption revOp,
                     const CombineMode CM,
                     const bool restrictedMode);

    void doPosts(const Details::DistributorPlan& distributorPlan,
                 size_t constantNumPackets,
                 bool commOnHost,
                 std::shared_ptr<std::string> prefix,
                 const bool canTryAliasing,
                 const CombineMode CM);

    void doPackAndPrepare(const SrcDistObject& src,
                          const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& exportLIDs,
                          size_t& constantNumPackets,
                          const execution_space &space);

    void doUnpackAndCombine(const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& remoteLIDs,
                            size_t constantNumPackets,
                            CombineMode CM,
                            const execution_space &space);

    /// \name Methods implemented by subclasses and used by doTransfer().
    ///
    /// The doTransfer() method uses the subclass' implementations of
    /// these methods to implement data transfer.  Subclasses of
    /// DistObject must implement these methods.  This is an instance
    /// of the <a
    /// href="http://en.wikipedia.org/wiki/Template_method_pattern">Template
    /// Method Pattern</a>.  ("Template" here doesn't mean "C++
    /// template"; it means "pattern with holes that are filled in by
    /// the subclass' method implementations.")
    //@{

    /// \brief Compare the source and target (\e this) objects for compatibility.
    ///
    /// \return True if they are compatible, else false.
    virtual bool
    checkSizes (const SrcDistObject& source) = 0;

    /// \brief Perform copies and permutations that are local to the
    ///   calling (MPI) process.
    ///
    /// Subclasses <i>must</i> reimplement this function.  Its default
    /// implementation does nothing.  Note that the <t>target</i>
    /// object of the Export or Import, namely <tt>*this</tt>, packs
    /// the <i>source</i> object's data.
    ///
    /// \pre permuteToLIDs and permuteFromLIDs are sync'd to both host
    ///   and device.  That is,
    ///   <tt>permuteToLIDs.need_sync_host()</tt>,
    ///   <tt>permuteToLIDs.need_sync_device()</tt>,
    ///   <tt>permuteFromLIDs.need_sync_host()</tt>, and
    ///   <tt>permuteFromLIDs.need_sync_device()</tt> are all false.
    ///
    /// \param source [in] On entry, the source object of the Export
    ///   or Import operation.
    /// \param numSameIDs [in] The number of elements that are the
    ///   same on the source and target objects.  These elements live
    ///   on the same process in both the source and target objects.
    /// \param permuteToLIDs [in] List of the elements that are
    ///   permuted.  They are listed by their local index (LID) in the
    ///   destination object.
    /// \param permuteFromLIDs [in] List of the elements that are
    ///   permuted.  They are listed by their local index (LID) in the
    ///   source object.
    /// \param CM [in] CombineMode to be used during copyAndPermute; 
    ///   may or may not be used by the particular object being called;
    ///   behavior with respect to CombineMode may differ by object.
    virtual void
    copyAndPermute (const SrcDistObject& source,
                    const size_t numSameIDs,
                    const Kokkos::DualView<const local_ordinal_type*,
                      buffer_device_type>& permuteToLIDs,
                    const Kokkos::DualView<const local_ordinal_type*,
                      buffer_device_type>& permuteFromLIDs,
                    const CombineMode CM);

  // clang-format on
  /*! \brief Same as copyAndPermute, but do operations in \c space
   */
  virtual void copyAndPermute(
      const SrcDistObject &source, const size_t numSameIDs,
      const Kokkos::DualView<const local_ordinal_type *, buffer_device_type>
          &permuteToLIDs,
      const Kokkos::DualView<const local_ordinal_type *, buffer_device_type>
          &permuteFromLIDs,
      const CombineMode CM, const execution_space &space);
  // clang-format off

    /// \brief Pack data and metadata for communication (sends).
    ///
    /// Subclasses <i>must</i> reimplement this function.  Its default
    /// implementation does nothing.  Note that the <t>target</i>
    /// object of the Export or Import, namely <tt>*this</tt>, packs
    /// the <i>source</i> object's data.
    ///
    /// \pre exportLIDs is sync'd to both host and device.  That is,
    ///   <tt>exportLIDs.need_sync_host ()</tt> and
    ///   <tt>exportLIDs.need_sync_device()</tt> are both false.
    ///
    /// \param source [in] Source object for the redistribution.
    ///
    /// \param exportLIDs [in] List of the entries (as local IDs in
    ///   the source object) that Tpetra will send to other processes.
    ///
    /// \param exports [out] On exit, the packed data to send.
    ///   Implementations must reallocate this as needed (prefer
    ///   reusing the existing allocation if possible), and may modify
    ///   and/or sync this wherever they like.
    ///
    /// \param numPacketsPerLID [out] On exit, the implementation of
    ///   this method must do one of two things: either set
    ///   <tt>numPacketsPerLID[i]</tt> to the number of packets to be
    ///   packed for <tt>exportLIDs[i]</tt> and set
    ///   <tt>constantNumPackets</tt> to zero, or set
    ///   <tt>constantNumPackets</tt> to a nonzero value.  If the
    ///   latter, the implementation must not modify the entries of
    ///   <tt>numPacketsPerLID</tt>.  If the former, the
    ///   implementation may sync <tt>numPacketsPerLID</tt> this
    ///   wherever it likes, either to host or to device.  The
    ///   allocation belongs to DistObject, not to subclasses; don't
    ///   be tempted to change this to pass by reference.
    ///
    /// \param constantNumPackets [out] On exit, 0 if the number of
    ///   packets per LID could differ, else (if nonzero) the number
    ///   of packets per LID (which must be constant).
    virtual void
    packAndPrepare (const SrcDistObject& source,
                    const Kokkos::DualView<const local_ordinal_type*,
                      buffer_device_type>& exportLIDs,
                    Kokkos::DualView<packet_type*,
                      buffer_device_type>& exports,
                    Kokkos::DualView<size_t*,
                      buffer_device_type> numPacketsPerLID,
                    size_t& constantNumPackets);

    /*! \brief Same as packAndPrepare, but in an execution space instance
    */
    virtual void
    packAndPrepare (const SrcDistObject& source,
                    const Kokkos::DualView<const local_ordinal_type*,
                      buffer_device_type>& exportLIDs,
                    Kokkos::DualView<packet_type*,
                      buffer_device_type>& exports,
                    Kokkos::DualView<size_t*,
                      buffer_device_type> numPacketsPerLID,
                    size_t& constantNumPackets,
                    const execution_space &space);

    /// \brief Perform any unpacking and combining after
    ///   communication.
    ///
    /// Subclasses <i>must</i> reimplement this function.  Its default
    /// implementation does nothing.  Note that the <t>target</i>
    /// object of the Export or Import, namely <tt>*this</tt>, unpacks
    /// the received data into itself, possibly modifying its entries.
    ///
    /// \pre importLIDs is sync'd to both host and device.  That is,
    ///   <tt>importLIDs.need_sync_host ()</tt> and
    ///   <tt>importLIDs.need_sync_device()</tt> are both false.
    ///
    /// \param importLIDs [in] List of the entries (as LIDs in the
    ///   destination object) we received from other processes.
    ///
    /// \param imports [in/out] On input: Buffer of received data to
    ///   unpack.  DistObject promises nothing about where this is
    ///   sync'd.  Implementations may sync this wherever they like,
    ///   either to host or to device.  The allocation belongs to
    ///   DistObject, not to subclasses; don't be tempted to change
    ///   this to pass by reference.
    ///
    /// \param numPacketsPerLID [in/out] On input: If
    ///   <tt>constantNumPackets</tt> is zero, then
    ///   <tt>numPacketsPerLID[i]</tt> contains the number of packets
    ///   imported for </tt>importLIDs[i]</tt>.  DistObject promises
    ///   nothing about where this is sync'd.  Implementations may
    ///   sync this wherever they like, either to host or to device.
    ///   The allocation belongs to DistObject, not to subclasses;
    ///   don't be tempted to change this to pass by reference.
    ///
    /// \param constantNumPackets [in] If nonzero, then the number of
    ///   packets per LID is the same for all entries ("constant") and
    ///   <tt>constantNumPackets</tt> is that number.  If zero, then
    ///   <tt>numPacketsPerLID[i]</tt> is the number of packets to
    ///   unpack for LID <tt>importLIDs[i]</tt>.
    ///
    /// \param combineMode [in] The CombineMode to use when combining
    ///   the imported entries with existing entries.
    virtual void
    unpackAndCombine (const Kokkos::DualView<const local_ordinal_type*,
                        buffer_device_type>& importLIDs,
                      Kokkos::DualView<packet_type*,
                        buffer_device_type> imports,
                      Kokkos::DualView<size_t*,
                        buffer_device_type> numPacketsPerLID,
                      const size_t constantNumPackets,
                      const CombineMode combineMode);

    virtual void
    unpackAndCombine (const Kokkos::DualView<const local_ordinal_type*,
                        buffer_device_type>& importLIDs,
                      Kokkos::DualView<packet_type*,
                        buffer_device_type> imports,
                      Kokkos::DualView<size_t*,
                        buffer_device_type> numPacketsPerLID,
                      const size_t constantNumPackets,
                      const CombineMode combineMode,
                      const execution_space &space);

    //! The Map over which this object is distributed.
    Teuchos::RCP<const map_type> map_;

  protected:
    std::unique_ptr<std::string>
    createPrefix(const char className[],
                 const char methodName[]) const;

    /// \brief Buffer into which packed data are imported (received
    ///   from other processes).
    ///
    /// Unfortunately, I had to declare these protected, because
    /// CrsMatrix uses them at one point.  Please, nobody else use
    /// them.
    Kokkos::DualView<packet_type*, buffer_device_type> imports_;

    /// \brief Reallocate imports_ if needed.
    ///
    /// This unfortunately must be declared protected, for the same
    /// reason that imports_ is declared protected.
    ///
    /// \param newSize [in] New size of imports_.
    /// \param verbose [in] Whether to print verbose debugging output
    ///   to stderr on every (MPI) process in the communicator.
    /// \param prefix [in] If <tt>verbose</tt> is <tt>true</tt>, then
    ///   this is a nonnull prefix to print at the beginning of each
    ///   line of verbose debugging output.  Otherwise, not used.
    ///
    /// \return Whether we actually reallocated.
    ///
    /// We don't need a "reallocExportsIfNeeded" method, because
    /// <tt>exports_</tt> always gets passed into packAndPrepare()
    /// by nonconst reference.  Thus, that method can resize the
    /// DualView without needing to call other DistObject methods.
    virtual bool
    reallocImportsIfNeeded (const size_t newSize,
                            const bool verbose,
                            const std::string* prefix,
                            const bool remoteLIDsContiguous=false,
                            const CombineMode CM=INSERT);

    /// \brief Number of packets to receive for each receive operation.
    ///
    /// This array is used in Distributor::doPosts() (and
    /// doReversePosts()) when starting the ireceive operation.
    ///
    /// This may be ignored in doTransfer() if constantNumPackets
    /// is nonzero, indicating a constant number of packets per LID.
    /// (For example, MultiVector sets the constantNumPackets output
    /// argument of packAndPrepare() to the number of columns in
    /// the multivector.)
    ///
    /// Unfortunately, I had to declare this protected, because
    /// CrsMatrix uses it at one point.  Please, nobody else use it.
    Kokkos::DualView<size_t*, buffer_device_type> numImportPacketsPerLID_;

    /// \brief Buffer from which packed data are exported (sent to
    ///   other processes).
    ///
    /// Unfortunately, I had to declare this protected, because
    /// CrsMatrix uses it at one point.  Please, nobody else use it.
    Kokkos::DualView<packet_type*, buffer_device_type> exports_;

    /// \brief Number of packets to send for each send operation.
    ///
    /// This array is used in Distributor::doPosts() (and
    /// doReversePosts()) for preparing for the send operation.
    ///
    /// This may be ignored in doTransfer() if constantNumPackets is
    /// nonzero, indicating a constant number of packets per LID.
    /// (For example, MultiVector sets the constantNumPackets output
    /// argument of packAndPrepare() to the number of columns in the
    /// multivector.)
    ///
    /// Unfortunately, I had to declare this protected, because
    /// CrsMatrix uses them at one point.  Please, nobody else use it.
    Kokkos::DualView<size_t*, buffer_device_type> numExportPacketsPerLID_;

  private:
    using this_type = DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>;

    Details::DistributorActor distributorActor_;

#ifdef HAVE_TPETRA_TRANSFER_TIMERS
    Teuchos::RCP<Teuchos::Time> doXferTimer_;
    Teuchos::RCP<Teuchos::Time> copyAndPermuteTimer_;
    Teuchos::RCP<Teuchos::Time> packAndPrepareTimer_;
    Teuchos::RCP<Teuchos::Time> doPostsAndWaitsTimer_;
    Teuchos::RCP<Teuchos::Time> unpackAndCombineTimer_;
#endif // HAVE_TPETRA_TRANSFER_TIMERS
  }; // class DistObject
} // namespace Tpetra

#endif // TPETRA_DISTOBJECT_DECL_HPP
