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

#ifndef TPETRA_DISTOBJECT_KA_DECL_HPP
#define TPETRA_DISTOBJECT_KA_DECL_HPP

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_SrcDistObject.hpp"
#include <Kokkos_NodeAPIConfigDefs.hpp> // enum KokkosClassic::ReadWriteOption

#if TPETRA_USE_KOKKOS_DISTOBJECT

#include "KokkosCompat_View.hpp"
#include "Kokkos_Core.hpp"

// #ifndef HAVE_TPETRA_TRANSFER_TIMERS
// #  define HAVE_TPETRA_TRANSFER_TIMERS 1
// #endif // HAVE_TPETRA_TRANSFER_TIMERS

#ifdef HAVE_TPETRA_TRANSFER_TIMERS
#  undef HAVE_TPETRA_TRANSFER_TIMERS
#endif // HAVE_TPETRA_TRANSFER_TIMERS


namespace Tpetra {
  /// \class DistObjectKA
  /// \brief Base class for distributed Tpetra objects that support
  ///   data redistribution.
  ///
  /// This class is a variant of DistObject that uses Kokkos::View to
  /// pack and unpack data.  This improves performance on GPU
  /// architectures.  TPETRA_USE_KOKKOS_DISTOBJECT must be defined in
  /// order to use this class.
  template <class Packet,
            class LocalOrdinal = Map<>::local_ordinal_type,
            class GlobalOrdinal = typename Map<LocalOrdinal>::global_ordinal_type,
            class Node = typename Map<LocalOrdinal, GlobalOrdinal>::node_type>
  class DistObjectKA :
    virtual public SrcDistObject,
    virtual public Teuchos::Describable {
  public:
    //! @name Typedefs
    //@{

    /// \brief The type of each datum being sent or received in an Import or Export.
    ///
    /// Note that this type does not always correspond to the
    /// <tt>Scalar</tt> template parameter of subclasses.
    typedef Packet packet_type;
    //! The type of local indices.
    typedef LocalOrdinal local_ordinal_type;
    //! The type of global indices.
    typedef GlobalOrdinal global_ordinal_type;
    //! The Kokkos Node type.
    typedef Node node_type;

    //! The Kokkos Device type
    typedef typename Kokkos::Compat::NodeDevice<node_type>::type device_type;
    typedef typename device_type::size_type view_size_type;

    //@}
    //! @name Constructors and destructor
    //@{

    //! Constructor.
    explicit DistObjectKA (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map);

    //! Copy constructor.
    DistObjectKA (const DistObjectKA<Packet,LocalOrdinal,GlobalOrdinal,Node>& rhs);

    //! Destructor (virtual for memory safety of derived classes).
    virtual ~DistObjectKA ();

    //@}
    //! @name Public methods for redistributing data
    //@{

    /// \brief Import data into this object using an Import object ("forward mode").
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
    /// \param source [in] The "source" object for redistribution.
    /// \param importer [in] Precomputed data redistribution plan.
    ///   Its source Map must be the same as the input DistObject's Map,
    ///   and its target Map must be the same as <tt>this->getMap()</tt>.
    /// \param CM [in] How to combine incoming data with the same
    ///   global index.
    void
    doImport (const SrcDistObject& source,
              const Import<LocalOrdinal,GlobalOrdinal,Node>& importer,
              CombineMode CM);

    /// \brief Export data into this object using an Export object ("forward mode").
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
    /// \param source [in] The "source" object for redistribution.
    /// \param exporter [in] Precomputed data redistribution plan.
    ///   Its source Map must be the same as the input DistObject's Map,
    ///   and its target Map must be the same as <tt>this->getMap()</tt>.
    /// \param CM [in] How to combine incoming data with the same
    ///   global index.
    void
    doExport (const SrcDistObject& source,
              const Export<LocalOrdinal,GlobalOrdinal,Node>& exporter,
              CombineMode CM);

    /// \brief Import data into this object using an Export object ("reverse mode").
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
    /// \param source [in] The "source" object for redistribution.
    /// \param exporter [in] Precomputed data redistribution plan.
    ///   Its <i>target</i> Map must be the same as the input DistObject's Map,
    ///   and its <i>source</i> Map must be the same as <tt>this->getMap()</tt>.
    ///   (Note the difference from forward mode.)
    /// \param CM [in] How to combine incoming data with the same
    ///   global index.
    void
    doImport (const SrcDistObject& source,
              const Export<LocalOrdinal,GlobalOrdinal,Node>& exporter,
              CombineMode CM);

    /// \brief Export data into this object using an Import object ("reverse mode").
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
    /// \param source [in] The "source" object for redistribution.
    /// \param importer [in] Precomputed data redistribution plan.
    ///   Its <i>target</i> Map must be the same as the input DistObject's Map,
    ///   and its <i>source</i> Map must be the same as <tt>this->getMap()</tt>.
    ///   (Note the difference from forward mode.)
    /// \param CM [in] How to combine incoming data with the same
    ///   global index.
    void
    doExport (const SrcDistObject& source,
              const Import<LocalOrdinal,GlobalOrdinal,Node>& importer,
              CombineMode CM);

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
    virtual Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >
    getMap() const { return map_; }

    //@}
    //! @name I/O methods
    //@{

    /// \brief Print this object to the given output stream.
    ///
    /// We generally assume that all MPI processes can print to the
    /// given stream.
    void print (std::ostream &os) const;

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
              const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;
    //@}
    //! @name Methods for use only by experts
    //@{

    /// \brief Remove processes which contain no elements in this object's Map.
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
    /// the original communicator which contain zero elements
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
    removeEmptyProcessesInPlace (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& newMap);

    // Forward declaration of nonmember function.
    template<class PT, class LO, class GO, class NT>
    friend void
    removeEmptyProcessesInPlace (Teuchos::RCP<Tpetra::DistObjectKA<PT, LO, GO, NT> >& input,
                                 const Teuchos::RCP<const Map<LO, GO, NT> >& newMap);
    // Forward declaration of nonmember function.
    template<class PT, class LO, class GO, class NT>
    friend void
    removeEmptyProcessesInPlace (Teuchos::RCP<Tpetra::DistObjectKA<PT, LO, GO, NT> >& input);
    //@}

  protected:
    /// \enum ReverseOption
    /// \brief Whether the data transfer should be performed in forward or reverse mode.
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
    ///   to have a constant number of packets per LID, and if so, how
    ///   many packets per LID there are.
    ///
    /// If this method returns zero, the instance says that it might
    /// possibly have a different number of packets for each LID to
    /// send or receive.  If it returns nonzero, the instance promises
    /// that the number of packets is the same for all LIDs, and that
    /// the return value is this number of packets per LID.
    ///
    /// The default implementation of this method returns zero.  This
    /// does not affect the behavior of doTransfer() in any way.  If a
    /// nondefault implementation returns nonzero, doTransfer() will
    /// use this information to avoid unnecessary allocation and / or
    /// resizing of arrays.
    virtual size_t constantNumberOfPackets () const;

    /// \brief Redistribute data across memory images.
    ///
    /// \param src [in] The source object, to redistribute into
    ///   the target object, which is <tt>*this</tt> object.
    ///
    /// \param CM [in] The combine mode that describes how to combine
    ///   values that map to the same global ID on the same process.
    ///
    /// \param permuteToLIDs [in] See copyAndPermute().
    ///
    /// \param permuteFromLIDs [in] See copyAndPermute().
    ///
    /// \param remoteLIDs [in] List of entries (as local IDs) in the
    ///   destination object to receive from other processes.
    ///
    /// \param exportLIDs [in] See packAndPrepare().
    ///
    /// \param distor [in/out] The Distributor object that knows how
    ///   to redistribute data.
    ///
    /// \param revOp [in] Whether to do a forward or reverse mode
    ///   redistribution.
    virtual void
    doTransfer (const SrcDistObject& src,
                CombineMode CM,
                size_t numSameIDs,
                const Teuchos::ArrayView<const LocalOrdinal> &permuteToLIDs,
                const Teuchos::ArrayView<const LocalOrdinal> &permuteFromLIDs,
                const Teuchos::ArrayView<const LocalOrdinal> &remoteLIDs,
                const Teuchos::ArrayView<const LocalOrdinal> &exportLIDs,
                Distributor &distor,
                ReverseOption revOp);

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

    /// \brief Perform copies and permutations that are local to this process.
    ///
    /// \param source [in] On entry, the source object, from which we
    ///   are distributing.  We distribute to the destination object,
    ///   which is <tt>*this</tt> object.
    /// \param numSameIDs [in] The umber of elements that
    ///   are the same on the source and destination (this) objects.
    ///   These elements are owned by the same process in both the
    ///   source and destination objects.  No permutation occurs.
    /// \param numPermuteIDs [in] The number of elements that are
    ///   locally permuted between the source and destination objects.
    /// \param permuteToLIDs [in] List of the elements that are
    ///   permuted.  They are listed by their LID in the destination
    ///   object.
    /// \param permuteFromLIDs [in] List of the elements that are
    ///   permuted.  They are listed by their LID in the source
    ///   object.
    virtual void
    copyAndPermute (
      const SrcDistObject& source,
      size_t numSameIDs,
      const Kokkos::View<const LocalOrdinal*, device_type> &permuteToLIDs,
      const Kokkos::View<const LocalOrdinal*, device_type> &permuteFromLIDs) = 0;

    /// \brief Perform any packing or preparation required for communication.
    ///
    /// \param source [in] Source object for the redistribution.
    ///
    /// \param exportLIDs [in] List of the entries (as local IDs in
    ///   the source object) we will be sending to other images.
    ///
    /// \param exports [out] On exit, the buffer for data to send.
    ///
    /// \param numPacketsPerLID [out] On exit, the implementation of
    ///   this method must do one of two things: set
    ///   numPacketsPerLID[i] to contain the number of packets to be
    ///   exported for exportLIDs[i] and set constantNumPackets to
    ///   zero, or set constantNumPackets to a nonzero value.  If the
    ///   latter, the implementation need not fill numPacketsPerLID.
    ///
    /// \param constantNumPackets [out] On exit, 0 if numPacketsPerLID
    ///   has variable contents (different size for each LID).  If
    ///   nonzero, then it is expected that the number of packets per
    ///   LID is constant, and that constantNumPackets is that value.
    ///
    /// \param distor [in] The Distributor object we are using.
    virtual void
    packAndPrepare (
      const SrcDistObject& source,
      const Kokkos::View<const LocalOrdinal*, device_type> &exportLIDs,
      Kokkos::View<Packet*, device_type> &exports,
      const Kokkos::View<size_t*, device_type> &numPacketsPerLID,
      size_t& constantNumPackets,
      Distributor &distor) = 0;

    /// \brief Perform any unpacking and combining after communication.
    ///
    /// \param importLIDs [in] List of the entries (as LIDs in the
    ///   destination object) we received from other images.
    ///
    /// \param imports [in] Buffer containing data we received.
    ///
    /// \param numPacketsPerLID [in] If constantNumPackets is zero,
    ///   then numPacketsPerLID[i] contains the number of packets
    ///   imported for importLIDs[i].
    ///
    /// \param constantNumPackets [in] If nonzero, then
    ///   numPacketsPerLID is constant (same value in all entries) and
    ///   constantNumPackets is that value.  If zero, then
    ///   numPacketsPerLID[i] is the number of packets imported for
    ///   importLIDs[i].
    ///
    /// \param distor [in] The Distributor object we are using.
    ///
    /// \param CM [in] The combine mode to use when combining the
    ///   imported entries with existing entries.
    virtual void
    unpackAndCombine (
      const Kokkos::View<const LocalOrdinal*, device_type> &importLIDs,
      const Kokkos::View<const Packet*, device_type> &imports,
      const Kokkos::View<size_t*, device_type> &numPacketsPerLID,
      size_t constantNumPackets,
      Distributor &distor,
      CombineMode CM) = 0;
    //@}

    /// \brief Hook for creating a const view.
    ///
    /// doTransfer() calls this on the source object.  By default,
    /// it does nothing, but the source object can use this as a hint
    /// to fetch data from a compute buffer on an off-CPU device (such
    /// as a GPU) into host memory.
    virtual void createViews () const;

    /// \brief Hook for creating a nonconst view.
    ///
    /// doTransfer() calls this on the destination (<tt>*this</tt>)
    /// object.  By default, it does nothing, but the destination
    /// object can use this as a hint to fetch data from a compute
    /// buffer on an off-CPU device (such as a GPU) into host memory.
    ///
    /// \param rwo [in] Whether to create a write-only or a
    ///   read-and-write view.  For Kokkos Node types where compute
    ///   buffers live in a separate memory space (e.g., in the device
    ///   memory of a discrete accelerator like a GPU), a write-only
    ///   view only requires copying from host memory to the compute
    ///   buffer, whereas a read-and-write view requires copying both
    ///   ways (once to read, from the compute buffer to host memory,
    ///   and once to write, back to the compute buffer).
    virtual void createViewsNonConst (KokkosClassic::ReadWriteOption rwo);

    /// \brief Hook for releasing views.
    ///
    /// doTransfer() calls this on both the source and destination
    /// objects, once it no longer needs to access that object's data.
    /// By default, this method does nothing.  Implementations may use
    /// this as a hint to free host memory which is a view of a
    /// compute buffer, once the host memory view is no longer needed.
    /// Some implementations may prefer to mirror compute buffers in
    /// host memory; for these implementations, releaseViews() may do
    /// nothing.
    virtual void releaseViews () const;

    //! The Map over which this object is distributed.
    Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > map_;

  private:
    //! Buffer into which packed data are imported (received from other processes).
    Kokkos::View<packet_type*,device_type> imports_;

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
    Kokkos::View<size_t*,device_type> numImportPacketsPerLID_;
    typename Kokkos::View<size_t*,device_type>::HostMirror host_numImportPacketsPerLID_;

    //! Buffer from which packed data are exported (sent to other processes).
    Kokkos::View<packet_type*,device_type> exports_;

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
    Kokkos::View<size_t*,device_type> numExportPacketsPerLID_;

#ifdef HAVE_TPETRA_TRANSFER_TIMERS
    Teuchos::RCP<Teuchos::Time> doXferTimer_;
    Teuchos::RCP<Teuchos::Time> copyAndPermuteTimer_;
    Teuchos::RCP<Teuchos::Time> packAndPrepareTimer_;
    Teuchos::RCP<Teuchos::Time> doPostsAndWaitsTimer_;
    Teuchos::RCP<Teuchos::Time> unpackAndCombineTimer_;
#endif // HAVE_TPETRA_TRANSFER_TIMERS

  }; // class DistObject


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

} // namespace Tpetra

#endif /* TPETRA_USE_KOKKOS_DISTOBJECT */

#endif /* TPETRA_DISTOBJECT_DECL_HPP */
