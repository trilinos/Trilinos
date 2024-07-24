// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// clang-format off
#ifndef TPETRA_MULTIVECTOR_DECL_HPP
#define TPETRA_MULTIVECTOR_DECL_HPP

/// \file Tpetra_MultiVector_decl.hpp
/// \brief Declaration of the Tpetra::MultiVector class

#include "Tpetra_MultiVector_fwd.hpp"
#include "Tpetra_Vector_fwd.hpp"
#include "Tpetra_FEMultiVector_fwd.hpp"
#include "Tpetra_DistObject.hpp"
#include "Tpetra_Map_fwd.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Kokkos_DualView.hpp"
#include "Teuchos_BLAS_types.hpp"
#include "Teuchos_DataAccess.hpp"
#include "Teuchos_Range1D.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "Kokkos_InnerProductSpaceTraits.hpp"
#include "Tpetra_KokkosRefactor_Details_MultiVectorLocalDeepCopy.hpp"
#include "Tpetra_Access.hpp"
#include "Tpetra_Details_WrappedDualView.hpp"
#include <type_traits>

#ifdef HAVE_TPETRACORE_TEUCHOSNUMERICS
#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Teuchos {
  template<class OrdinalType, class ScalarType>
  class SerialDenseMatrix; // forward declaration
}
#endif // DOXYGEN_SHOULD_SKIP_THIS
#endif // HAVE_TPETRACORE_TEUCHOSNUMERICS


namespace Tpetra {


  /// \brief Copy the contents of the MultiVector \c src into \c dst.
  /// \relatesalso MultiVector
  ///
  /// \pre The two inputs must have the same communicator.
  /// \pre The Map of \c src must be compatible with the Map of \c dst.
  /// \pre The two inputs must have the same number of columns.
  ///
  /// Copy the contents of the MultiVector \c src into the MultiVector
  /// \c dst.  ("Copy the contents" means the same thing as "deep
  /// copy.")  The two MultiVectors need not necessarily have the same
  /// template parameters, but the assignment of their entries must
  /// make sense.  Furthermore, their Maps must be compatible, that
  /// is, the MultiVectors' local dimensions must be the same on all
  /// processes.
  ///
  /// This method must always be called as a collective operation on
  /// all processes over which the multivector is distributed.  This
  /// is because the method reserves the right to check for
  /// compatibility of the two Maps, at least in debug mode, and throw
  /// if they are not compatible.
  template <class DS, class DL, class DG, class DN,
            class SS, class SL, class SG, class SN>
  void
  deep_copy (MultiVector<DS, DL, DG, DN>& dst,
             const MultiVector<SS, SL, SG, SN>& src);

#ifdef HAVE_TPETRACORE_TEUCHOSNUMERICS
  /// \brief Copy the contents of a Teuchos::SerialDenseMatrix into
  ///   the local part of the given Tpetra::MultiVector.
  /// \relatesalso MultiVector
  ///
  /// \pre <tt>src.numRows() == dst.getLocalLength()</tt>
  /// \pre <tt>src.numCols() == dst.getNumVectors()</tt>
  template <class ST, class LO, class GO, class NT>
  void
  deep_copy (MultiVector<ST, LO, GO, NT>& dst,
             const Teuchos::SerialDenseMatrix<int, ST>& src);

  /// \brief Copy the local part of the Tpetra::MultiVector into the
  ///   Teuchos::SerialDenseMatrix.
  /// \relatesalso MultiVector
  ///
  /// \pre <tt>src.numRows() == dst.getLocalLength()</tt>
  /// \pre <tt>src.numCols() == dst.getNumVectors()</tt>
  template <class ST, class LO, class GO, class NT>
  void
  deep_copy (Teuchos::SerialDenseMatrix<int, ST>& dst,
             const MultiVector<ST, LO, GO, NT>& src);
#endif // HAVE_TPETRACORE_TEUCHOSNUMERICS

  /// \brief Return a deep copy of the given MultiVector.
  /// \relatesalso MultiVector
  ///
  /// \note MultiVector's constructor returns a <i>shallow</i> copy of
  ///   its input, by default.  If you want a deep copy, use the
  ///   two-argument copy constructor with Teuchos::Copy as the second
  ///   argument, or call this function (createCopy).
  template <class ST, class LO, class GO, class NT>
  MultiVector<ST, LO, GO, NT>
  createCopy (const MultiVector<ST, LO, GO, NT>& src);

  /// \brief Nonmember MultiVector "constructor": Create a MultiVector
  ///   from a given Map.
  /// \relatesalso MultiVector
  /// \relatesalso Vector
  ///
  /// \param map [in] Map describing the distribution of rows of the
  ///   resulting MultiVector.
  /// \param numVectors [in] Number of columns of the resulting
  ///   MultiVector.
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  createMultiVector (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& map,
                     const size_t numVectors);

  // WARNING NOT FOR USERS
  // This means we don't need to make MultiVector a friend of
  // Vector or of itself (with different template parameters).
  template<class SC, class LO, class GO, class NT>
  Teuchos::ArrayView<const size_t>
  getMultiVectorWhichVectors (const MultiVector<SC, LO, GO, NT>& X);


  /// \brief One or more distributed dense vectors.
  ///
  /// A "multivector" contains one or more dense vectors.  All the
  /// vectors in a multivector have the same distribution of rows in
  /// parallel over the communicator used to create the multivector.
  /// Multivectors containing more than one vector are useful for
  /// algorithms that solve multiple linear systems at once, or that
  /// solve for a cluster of eigenvalues and their corresponding
  /// eigenvectors at once.  These "block" algorithms often have
  /// accuracy or performance advantages over corresponding algorithms
  /// that solve for only one vector at a time.  For example, working
  /// with multiple vectors at a time allows Tpetra to use faster BLAS
  /// 3 routines for local computations.  It may also reduce the
  /// number of parallel reductions.
  ///
  /// The Vector class implements the MultiVector interface, so if you
  /// only wish to work with a single vector at a time, you may simply
  /// use Vector instead of MultiVector.  However, if you are writing
  /// solvers or preconditioners, you would do better to write to the
  /// MultiVector interface and always assume that each MultiVector
  /// contains more than one vector.  This will make your solver or
  /// preconditioner more compatible with other Trilinos packages, and
  /// it will also let you exploit the performance optimizations
  /// mentioned above.
  ///
  /// \tparam Scalar The type of each entry of the multivector.  (You
  ///   can use real-valued or complex-valued types here, unlike in
  ///   Epetra, where the scalar type is always \c double.)
  /// \tparam LocalOrdinal The type of local indices.  See the
  ///   documentation of Map for requirements.
  /// \tparam GlobalOrdinal The type of global indices.  See the
  ///   documentation of Map for requirements.
  /// \tparam Node The Kokkos Node type.
  ///
  /// \section Kokkos_KR_MV_prereq Prerequisites
  ///
  /// Before reading the rest of this documentation, it helps to know
  /// a little bit about Kokkos.  In particular, you should know about
  /// execution spaces, memory spaces, and shallow copy semantics.
  /// You should also know something about the Teuchos memory
  /// management classes, in particular Teuchos::RCP, though it helps
  /// to know a bit about Teuchos::ArrayRCP and Teuchos::ArrayView as
  /// well.  You may also want to know about the differences between
  /// BLAS 1, 2, and 3 operations, and learn a little bit about MPI
  /// (the Message Passing Interface for distributed-memory
  /// programming).  You won't have to use MPI directly to use
  /// MultiVector, but it helps to be familiar with the general idea
  /// of distributed storage of data over a communicator.
  ///
  /// \section Kokkos_KR_MV_view A MultiVector is a view of data
  ///
  /// A MultiVector is a view of data.  A <i>view</i> behaves like a
  /// pointer; it provides access to the original multivector's data
  /// without copying the data.  This means that the copy constructor
  /// and assignment operator (<tt>operator=</tt>) do <i>shallow</i>
  /// copies.  They do <i>not</i> copy the data; they just copy
  /// pointers and other "metadata."  If you would like to copy a
  /// MultiVector into an existing MultiVector, call the nonmember
  /// function deep_copy().  If you would like to create a new
  /// MultiVector which is a deep copy of an existing MultiVector,
  /// call the nonmember function createCopy(), or use the
  /// two-argument copy constructor with Teuchos::Copy as the second
  /// argument.
  ///
  /// Views have the additional property that they automatically
  /// handle deallocation.  They use <i>reference counting</i> for
  /// this, much like how std::shared_ptr works.  That means you do
  /// not have to worry about "freeing" a MultiVector after it has
  /// been created.  Furthermore, you may pass shallow copies around
  /// without needing to worry about which is the "master" view of the
  /// data.  There is no "master" view of the data; when the last view
  /// falls out of scope, the data will be deallocated.
  ///
  /// This is what the documentation means when it speaks of <i>view
  /// semantics</i>.  The opposite of that is <i>copy</i> or
  /// <i>container</i> semantics, where the copy constructor and
  /// <tt>operator=</tt> do deep copies (of the data).  We say that
  /// "std::vector has container semantics," and "MultiVector has view
  /// semantics."
  ///
  /// MultiVector also has "subview" methods that give results
  /// analogous to the Kokkos::subview() function.  That is, they
  /// return a MultiVector which views some subset of another
  /// MultiVector's rows and columns.  The subset of columns in a view
  /// need not be contiguous.  For example, given a multivector X with
  /// 43 columns, it is possible to have a multivector Y which is a
  /// view of columns 1, 3, and 42 (zero-based indices) of X.  We call
  /// such multivectors <i>noncontiguous</i>.  They have the the
  /// property that isConstantStride() returns false.
  ///
  /// Noncontiguous multivectors lose some performance advantages.
  /// For example, local computations may be slower, since Tpetra
  /// cannot use BLAS 3 routines (e.g., matrix-matrix multiply) on a
  /// noncontiguous multivectors without copying into temporary
  /// contiguous storage.  For performance reasons, if you get a
  /// Kokkos::View of a noncontiguous MultiVector's local data, it
  /// does <i>not</i> correspond to the columns that the MultiVector
  /// views.
  ///
  /// \section Kokkos_KR_MV_dual DualView semantics
  ///
  /// Tpetra was designed to perform well on many different kinds of
  /// computers.  Some computers have different <i>memory spaces</i>.
  /// For example, GPUs (Graphics Processing Units) by NVIDIA have
  /// "device memory" and "host memory."  The GPU has faster access to
  /// device memory than host memory, but usually there is less device
  /// memory than host memory.  Intel's "Knights Landing" architecture
  /// has two different memory spaces, also with different capacity
  /// and performance characteristics.  Some architectures let the
  /// processor address memory in any space, possibly with a
  /// performance penalty.  Others can only access data in certain
  /// spaces, and require a special system call to copy memory between
  /// spaces.
  ///
  /// The Kokkos package provides abstractions for handling multiple
  /// memory spaces.  In particular, Kokkos::DualView lets users
  /// "mirror" data that live in one space, with data in another
  /// space.  It also lets users manually mark data in one space as
  /// modified (modify()), and synchronize (sync()) data from one
  /// space to another.  The latter only actually copies if the data
  /// have been marked as modified.  Users can access data in a
  /// particular space by calling view().  All three of these methods
  /// -- modify(), sync(), and view() -- are templated on the memory
  /// space.  This is how users select the memory space on which they
  /// want the method to act.
  ///
  /// MultiVector implements "DualView semantics."  This means that it
  /// implements the above three operations:
  /// <ul>
  /// <li> modify(): Mark data in a memory space as modified (or about
  ///      to be modified) </li>
  /// <li> sync(): If data in the target memory space are least
  ///      recently modified compared with the other space, copy data
  ///      to the target memory space </li>
  /// <li> getLocalView(): Return a Kokkos::View of the data in a
  ///      given memory space </li>
  /// </ul>
  ///
  /// If your computer only has one memory space, as with conventional
  /// single-core or multicore processors, you don't have to worry
  /// about this.  You can ignore the modify() and sync() methods in
  /// that case.
  ///
  /// \section Kokkos_KR_MV_access How to access the local data
  ///
  /// The getLocalView() method for getting a Kokkos::View is the main
  /// way to access a MultiVector's local data.  If you want to read
  /// or write the actual values in a multivector, this is what you
  /// want.  The resulting Kokkos::View behaves like a 2-D array.  You
  /// can address it using an index pair (i,j), where i is the local
  /// row index, and j is the column index.
  ///
  /// MultiVector also has methods that return an
  /// Teuchos::ArrayRCP<Scalar> ("1-D view"), or a
  /// Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > ("2-D view").
  /// These exist only for backwards compatibility, and also give
  /// access to the local data.
  ///
  /// All of these views only view <i>local</i> data.  This means that
  /// the corresponding rows of the multivector are owned by the
  /// calling (MPI) process.  You may <i>not</i> use these methods to
  /// access <i>remote</i> data, that is, rows that do not belong to
  /// the calling process.
  ///
  /// MultiVector's public interface also has methods for modifying
  /// local data, like sumIntoLocalValue() and replaceGlobalValue().
  /// These methods act on host data <i>only</i>.  To access or modify
  /// device data, you must get the Kokkos::View and work with it
  /// directly.
  ///
  /// \section Kokkos_KR_MV_why Why won't you give me a raw pointer?
  ///
  /// Tpetra was designed to allow different data representations
  /// underneath the same interface.  This lets Tpetra run correctly
  /// and efficiently on many different kinds of hardware.  These
  /// different kinds of hardware all have in common the following:
  /// <ul>
  /// <li> Data layout matters a lot for performance </li>
  /// <li> The right layout for your data depends on the hardware </li>
  /// <li> Data may be distributed over different memory spaces in
  ///      hardware, and efficient code must respect this, whether or
  ///      not the programming model presents the different memories
  ///      as a single address space </li>
  /// <li> Copying between different data layouts or memory spaces is
  ///      expensive and should be avoided whenever possible </li>
  /// <li> Optimal data layout may require control over initialization
  ///      of storage </li>
  /// </ul>
  /// These conclusions have practical consequences for the
  /// MultiVector interface.  In particular, we have deliberately made
  /// it difficult for you to access data directly by raw pointer.
  /// This is because the underlying layout may not be what you
  /// expect.  The memory might not even be accessible from the host
  /// CPU.  Instead, we give access through a Kokkos::View, which
  /// behaves like a 2-D array.  You can ask the Kokkos::View for a
  /// raw pointer by calling its <tt>data()</tt> method, but
  /// then you are responsible for understanding its layout in memory.
  ///
  /// \section Kokkos_KR_MV_dist Parallel distribution of data
  ///
  /// A MultiVector's rows are distributed over processes in its (row)
  /// Map's communicator.  A MultiVector is a DistObject; the Map of
  /// the DistObject tells which process in the communicator owns
  /// which rows.  This means that you may use Import and Export
  /// operations to migrate between different distributions.  Please
  /// refer to the documentation of Map, Import, and Export for more
  /// information.
  ///
  /// MultiVector includes methods that perform parallel all-reduces.
  /// These include inner products and various kinds of norms.  All of
  /// these methods have the same blocking semantics as
  /// <tt>MPI_Allreduce</tt>.
  ///
  /// \warning Some computational methods, such as inner products and
  ///   norms, may return incorrect results if the MultiVector's Map
  ///   is overlapping (not one-to-one) but not locally replicated.
  ///   That is, if some but not all rows are shared by more than one
  ///   process in the communicator, then inner products and norms may
  ///   be wrong.  This behavior may change in future releases.
  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  class MultiVector :
    public DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node>
  {
  public:
    //! @name Typedefs to facilitate template metaprogramming.
    //@{

    /// \brief The type of each entry in the MultiVector.
    using scalar_type = Scalar;
    /// \brief The type used internally in place of \c Scalar.
    ///
    /// Some \c Scalar types might not work with Kokkos on all
    /// execution spaces, due to missing CUDA device macros or missing
    /// volatile overloads of some methods.  The C++ standard type
    /// <tt>std::complex<T></tt> has this problem.  To fix this, we
    /// replace <tt>std::complex<T></tt> values internally with the
    /// bitwise identical type <tt>Kokkos::complex<T></tt>.  The
    /// latter is the <tt>impl_scalar_type</tt> corresponding to
    /// <tt>Scalar = std::complex<T></tt>.
    ///
    /// Most users don't need to know about this.  Just be aware that
    /// if you ask for a Kokkos::View or Kokkos::DualView of the
    /// MultiVector's data, its entries have type \c impl_scalar_type,
    /// not \c scalar_type.
    using impl_scalar_type =
      typename Kokkos::ArithTraits<Scalar>::val_type;

    //! The type of the Map specialization used by this class.
    using map_type = Map<LocalOrdinal, GlobalOrdinal, Node>;
    //! The type of local indices that this class uses.
    using local_ordinal_type = typename map_type::local_ordinal_type;
    //! The type of global indices that this class uses.
    using global_ordinal_type = typename map_type::global_ordinal_type;
    //! This class' preferred Kokkos device type.
    using device_type = typename map_type::device_type;
    //! Legacy thing that you should not use any more.
    using node_type = typename map_type::node_type;

    /// \brief Type of an inner ("dot") product result.
    ///
    /// This is usually the same as <tt>impl_scalar_type</tt>, but may
    /// differ if <tt>impl_scalar_type</tt> is e.g., an uncertainty
    /// quantification type from the Stokhos package.
    using dot_type =
      typename Kokkos::Details::InnerProductSpaceTraits<impl_scalar_type>::dot_type;

    /// \brief Type of a norm result.
    ///
    /// This is usually the same as the type of the magnitude
    /// (absolute value) of <tt>impl_scalar_type</tt>, but may differ if
    /// <tt>impl_scalar_type</tt> is e.g., an uncertainty quantification
    /// type from the Stokhos package.
    using mag_type = typename Kokkos::ArithTraits<impl_scalar_type>::mag_type;

    /// \brief Type of the (new) Kokkos execution space.
    ///
    /// The execution space implements parallel operations, like
    /// parallel_for, parallel_reduce, and parallel_scan.
    using execution_space = typename device_type::execution_space;

    /// \brief Kokkos::DualView specialization used by this class.
    ///
    /// This is of interest to users who already have a
    /// Kokkos::DualView, and want the MultiVector to view it.  By
    /// "view" it, we mean that the MultiVector doesn't copy the data
    /// in the DualView; it just hangs on to the pointer.
    ///
    /// We take particular care to template the DualView on an
    /// execution space, rather than a memory space.  This ensures
    /// that Tpetra will use exactly the specified execution space(s)
    /// and no others.  This matters because View (and DualView)
    /// initialization is a parallel Kokkos kernel.  If the View is
    /// templated on an execution space, Kokkos uses that execution
    /// space (and only that execution space) to initialize the View.
    /// This is what we want.  If the View is templated on a
    /// <i>memory</i> space, Kokkos uses the memory space's default
    /// <i>execution</i> space to initialize.  This is not necessarily
    /// what we want.  For example, if building with OpenMP enabled,
    /// the default execution space for host memory is Kokkos::OpenMP,
    /// even if the user-specified DeviceType is Kokkos::Serial.  That
    /// is why we go through the trouble of asking for the
    /// execution_space's execution space.
    using dual_view_type = Kokkos::DualView<impl_scalar_type**,
                                            Kokkos::LayoutLeft,
                                            device_type>;
    using wrapped_dual_view_type = Details::WrappedDualView<dual_view_type>;

    using host_view_type = typename dual_view_type::t_host;
    using device_view_type = typename dual_view_type::t_dev;

    //@}
    //! @name Constructors and destructor
    //@{

    //! Default constructor: makes a MultiVector with no rows or columns.
    MultiVector ();

    /// \brief Basic constuctor.
    ///
    /// \param map [in] Map describing the distribution of rows.
    /// \param numVecs [in] Number of vectors (columns).
    /// \param zeroOut [in] Whether to initialize all the entries of
    ///   the MultiVector to zero.
    MultiVector (const Teuchos::RCP<const map_type>& map,
                 const size_t numVecs,
                 const bool zeroOut = true);

    /// \brief Copy constructor, with option to do deep or shallow copy.
    ///
    /// The current (so-called "Kokkos refactor," circa >= 2014/5)
    /// version of Tpetra, unlike the previous "classic" version,
    /// always has view semantics.  Thus, copyOrView = Teuchos::View
    /// has no effect, and copyOrView = Teuchos::Copy does not mark
    /// this MultiVector as having copy semantics.  However,
    /// copyOrView = Teuchos::Copy will make the resulting MultiVector
    /// a deep copy of the input MultiVector.
    ///
    MultiVector (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& source,
                 const Teuchos::DataAccess copyOrView);

    /// \brief Create multivector by copying two-dimensional array of local data.
    ///
    /// \param map [in] The Map describing the distribution of rows of
    ///   the multivector.
    /// \param view [in] A view of column-major dense matrix data.
    ///   The calling process will make a deep copy of this data.
    /// \param LDA [in] The leading dimension (a.k.a. "stride") of the
    ///   column-major input data.
    /// \param NumVectors [in] The number of columns in the input data.
    ///   This will be the number of vectors in the returned
    ///   multivector.
    ///
    /// \pre LDA >= A.size()
    /// \pre NumVectors > 0
    /// \post isConstantStride() == true
    MultiVector (const Teuchos::RCP<const map_type>& map,
                 const Teuchos::ArrayView<const Scalar>& A,
                 const size_t LDA,
                 const size_t NumVectors);

    /// \brief Create multivector by copying array of views of local data.
    ///
    /// \param map [in] The Map describing the distribution of rows of
    ///   the multivector.
    /// \param ArrayOfPtrs [in/out] Array of views of each column's data.
    ///   The calling process will make a deep copy of this data.
    /// \param NumVectors [in] The number of columns in the input
    ///   data, and the number of elements in ArrayOfPtrs.  This will
    ///   be the number of vectors in the returned multivector.
    ///
    /// \pre NumVectors > 0
    /// \pre NumVectors == ArrayOfPtrs.size()
    /// \post constantStride() == true
    MultiVector (const Teuchos::RCP<const map_type>& map,
                 const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> >&ArrayOfPtrs,
                 const size_t NumVectors);

    /// \brief Constructor, that takes a Kokkos::DualView of the
    ///   MultiVector's data, and returns a MultiVector that views
    ///   those data.
    ///
    /// To "view those data" means that this MultiVector and the input
    /// Kokkos::DualView point to the same data, just like two "raw"
    /// pointers (e.g., <tt>double*</tt>) can point to the same data.
    /// If you modify one, the other sees it (subject to the
    /// limitations of cache coherence).
    ///
    /// \param map [in] Map describing the distribution of rows.
    /// \param view [in] Kokkos::DualView of the data to view.
    MultiVector (const Teuchos::RCP<const map_type>& map,
                 const dual_view_type& view);

    /// \brief Constructor, that takes a Kokkos::View of the
    ///   MultiVector's data (living in the Device's memory space),
    ///   and returns a MultiVector that views those data.
    ///
    /// \param map [in] Map describing the distribution of rows.
    /// \param view [in] Kokkos::View of the data to view.
    ///
    /// Q: What's the difference between this constructor (that takes
    /// a Kokkos::View), and the constructor above that takes a
    /// Kokkos::DualView?
    ///
    /// A: Suppose that for the MultiVector's device type, there are
    /// actually two memory spaces (e.g., for Kokkos::Cuda with UVM
    /// off, assuming that this is allowed).  In order for MultiVector
    /// to implement DualView semantics correctly, this constructor
    /// must allocate a Kokkos::View of host memory (or lazily
    /// allocate it on modify() or sync()).
    ///
    /// Now suppose that you pass in the same Kokkos::View of device
    /// memory to two different MultiVector instances, X and Y.  Each
    /// would allocate its own Kokkos::View of host memory.  That
    /// means that X and Y would have different DualView instances,
    /// but their DualView instances would have the same device View.
    ///
    /// Suppose that you do the following:
    /// <ol>
    /// <li> Modify X on host (calling modify() correctly) </li>
    /// <li> Modify Y on host (calling modify() correctly) </li>
    /// <li> Sync Y to device (calling sync() correctly) </li>
    /// <li> Sync X to device (calling sync() correctly) </li>
    /// </ol>
    /// This would clobber Y's later changes on host, in favor of X's
    /// earlier changes on host.  That could be very confusing.  We
    /// allow this behavior because Kokkos::DualView allows it.  That
    /// is, Kokkos::DualView also lets you get the device View, and
    /// hand it off to another Kokkos::DualView.  It's confusing, but
    /// users need to know what they are doing if they start messing
    /// around with multiple memory spaces.
    MultiVector (const Teuchos::RCP<const map_type>& map,
                 const typename dual_view_type::t_dev& d_view);
   
    /// \brief Expert mode constructor, that takes a Kokkos::DualView
    ///   of the MultiVector's data and the "original"
    ///   Kokkos::DualView of the data, and returns a MultiVector that
    ///   views those data.
    ///
    /// \warning This constructor is only for expert users.  We make
    ///   no promises about backwards compatibility for this
    ///   interface.  It may change or go away at any time.  It is
    ///   mainly useful for Tpetra developers and we do not expect it
    ///   to be useful for anyone else.
    ///
    /// \param map [in] Map describing the distribution of rows.
    /// \param view [in] View of the data (shallow copy).
    /// \param origView [in] The </i>original</i> view of the data.
    ///
    /// The original view keeps the "original" dimensions.  Doing so
    /// lets us safely construct a column Map view of a (domain Map
    /// view of a (column Map MultiVector)).  The result of a
    /// Kokkos::subview does not remember the original dimensions of
    /// the view, and does not allow constructing a view with a
    /// superset of rows or columns, so we have to keep the original
    /// view.
    MultiVector (const Teuchos::RCP<const map_type>& map,
                 const dual_view_type& view,
                 const dual_view_type& origView);

    /// \brief Expert mode constructor, that takes a WrappedDualView
    ///   of the MultiVector's data.
    ///
    /// \warning This constructor is only for expert users.  We make
    ///   no promises about backwards compatibility for this
    ///   interface.  It may change or go away at any time.  It is
    ///   mainly useful for Tpetra developers and we do not expect it
    ///   to be useful for anyone else.
    MultiVector (const Teuchos::RCP<const map_type>& map,
                 const wrapped_dual_view_type& d_view);


  protected:

    /// \brief Single-column subview constructor, for derived classes ONLY.
    ///
    /// \param X [in] Input MultiVector to view (in possibly nonconst fashion).
    /// \param j [in] The column of X to view.
    MultiVector (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
                 const size_t j);

  public:

    /// \brief Expert mode constructor for noncontiguous views.
    ///
    /// \warning This constructor is only for expert users.  We make
    ///   no promises about backwards compatibility for this
    ///   interface.  It may change or go away at any time.  It is
    ///   mainly useful for Tpetra developers and we do not expect it
    ///   to be useful for anyone else.
    ///
    /// This constructor takes a Kokkos::DualView for the MultiVector
    /// to view, and a list of the columns to view, and returns a
    /// MultiVector that views those data.  The resulting MultiVector
    /// does <i>not</i> have constant stride, that is,
    /// isConstantStride() returns false.
    ///
    /// \param map [in] Map describing the distribution of rows.
    /// \param view [in] Device view to the data (shallow copy).
    /// \param whichVectors [in] Which columns (vectors) to view.
    MultiVector (const Teuchos::RCP<const map_type>& map,
                 const dual_view_type& view,
                 const Teuchos::ArrayView<const size_t>& whichVectors);

    /// \brief Expert mode constructor for noncontiguous views.
    ///
    /// \warning This constructor is only for expert users.  We make
    ///   no promises about backwards compatibility for this
    ///   interface.  It may change or go away at any time.  It is
    ///   mainly useful for Tpetra developers and we do not expect it
    ///   to be useful for anyone else.
    ///
    /// This constructor takes a Kokkos::DualView for the MultiVector
    /// to view, and a list of the columns to view, and returns a
    /// MultiVector that views those data.  The resulting MultiVector
    /// does <i>not</i> have constant stride, that is,
    /// isConstantStride() returns false.
    ///
    /// \param map [in] Map describing the distribution of rows.
    /// \param view [in] WrappedDualView to the data (shallow copy).
    /// \param whichVectors [in] Which columns (vectors) to view.
    MultiVector (const Teuchos::RCP<const map_type>& map,
                 const wrapped_dual_view_type& view,
                 const Teuchos::ArrayView<const size_t>& whichVectors);


    /// \brief Expert mode constructor for noncontiguous views, with
    ///   original view.
    ///
    /// \warning This constructor is only for expert users.  We make
    ///   no promises about backwards compatibility for this
    ///   interface.  It may change or go away at any time.  It is
    ///   mainly useful for Tpetra developers and we do not expect it
    ///   to be useful for anyone else.
    ///
    /// This constructor takes a Kokkos::DualView for the MultiVector
    /// to view, a view of the "original" data, and a list of the
    /// columns to view, and returns a MultiVector that views those
    /// data.  The resulting MultiVector does <i>not</i> have constant
    /// stride, that is, isConstantStride() returns false.
    ///
    /// \param map [in] Map describing the distribution of rows.
    /// \param view [in] View of the data (shallow copy).
    /// \param origView [in] The </i>original</i> view of the data.
    /// \param whichVectors [in] Which columns (vectors) to view.
    ///
    /// The original view keeps the "original" dimensions.  Doing so
    /// lets us safely construct a column Map view of a (domain Map
    /// view of a (column Map MultiVector)).  The result of a
    /// Kokkos::subview does not remember the original dimensions of
    /// the view, and does not allow constructing a view with a
    /// superset of rows or columns, so we have to keep the original
    /// view.
    MultiVector (const Teuchos::RCP<const map_type>& map,
                 const dual_view_type& view,
                 const dual_view_type& origView,
                 const Teuchos::ArrayView<const size_t>& whichVectors);

    /// \brief "Offset view" constructor; make a view of a contiguous
    ///   subset of rows on each process.
    ///
    /// Return a view of the MultiVector \c X, which views a subset of
    /// the rows of \c X.  Specify the subset by a subset Map of this
    /// MultiVector's current row Map, and an optional (local) offset.
    /// "View" means "alias": if the original (this) MultiVector's
    /// data change, the view will see the changed data.
    ///
    /// \param X [in] The MultiVector to view.
    /// \param subMap [in] The row Map for the new MultiVector.  This
    ///   must be a subset Map of the input MultiVector's row Map.
    /// \param offset [in] The local row offset at which to start the view.
    ///
    /// Suppose that you have a MultiVector X, and you want to view X,
    /// on all processes in X's (MPI) communicator, as split into two
    /// row blocks X1 and X2.  One could express this in Matlab
    /// notation as X = [X1; X2], except that here, X1 and X2 are
    /// views into X, rather than copies of X's data.  This method
    /// assumes that the <i>local</i> indices of X1 and X2 are each
    /// contiguous, and that the local indices of X2 follow those of
    /// X1.  If that is not the case, you cannot use views to divide X
    /// into blocks like this; you must instead use the Import or
    /// Export functionality, which copies the relevant rows of X.
    ///
    /// Here is how you would construct the views X1 and X2.
    /// \code
    /// MultiVector<> X (...); // the input MultiVector
    /// // ... fill X with data ...
    ///
    /// using Teuchos::RCP;
    ///
    /// // Map that on each process in X's communicator,
    /// // contains the global indices of the rows of X1.
    /// RCP<const Map<>> map1 (new Map<> (...));
    /// // Map that on each process in X's communicator,
    /// // contains the global indices of the rows of X2.
    /// RCP<const Map<>> map2 (new Map<> (...));
    ///
    /// // Create the first view X1.  The second argument, the offset,
    /// // is the index of the local row at which to start the view.
    /// // X1 is the topmost block of X, so the offset is zero.
    /// MultiVector<> X1 (X, map1, 0);
    ///
    /// // Create the second view X2.  X2 is directly below X1 in X,
    /// // so the offset is the local number of rows in X1.  This is
    /// // the same as the local number of entries in map1.
    /// MultiVector<> X1 (X, map2, X1.getLocalLength ());
    /// \endcode
    ///
    /// It is legal, in the above example, for X1 or X2 to have zero
    /// local rows on any or all process(es).  In that case, the
    /// corresponding Map must have zero local entries on that / those
    /// process(es).  In particular, if X2 has zero local rows on a
    /// process, then the corresponding offset on that process would
    /// be the number of local rows in X (and therefore in X1) on that
    /// process.  This is the only case in which the sum of the local
    /// number of entries in \c subMap (in this case, zero) and the
    /// offset may equal the number of local entries in
    /// <tt>*this</tt>.
    MultiVector (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
                 const Teuchos::RCP<const map_type>& subMap,
                 const local_ordinal_type rowOffset = 0);

    /// \brief "Offset view" constructor, that takes the new Map as a
    ///   <tt>const Map&</tt> rather than by RCP.
    ///
    /// This constructor exists for backwards compatibility.  It
    /// invokes the input Map's copy constructor, which is a shallow
    /// copy.  Maps are immutable anyway, so the copy is harmless.
    MultiVector (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
                 const map_type& subMap,
                 const size_t offset = 0);

    /// \brief Copy constructor (shallow copy).
    ///
    /// MultiVector's copy constructor always does a shallow copy.
    /// Use the nonmember function <tt>Tpetra::deep_copy</tt> (see
    /// below) to deep-copy one existing MultiVector to another, and
    /// use the two-argument "copy constructor" (in this file, with
    /// <tt>copyOrView=Teuchos::Copy</tt>) to create a MultiVector
    /// that is a deep copy of an existing MultiVector.
    MultiVector (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&) = default;

    //! Move constructor (shallow move).
    MultiVector (MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&&) = default;

    /// \brief Copy assigment (shallow copy).
    ///
    /// MultiVector's copy constructor always does a shallow copy.
    /// Use the nonmember function <tt>Tpetra::deep_copy</tt> (see
    /// below) to deep-copy one existing MultiVector to another.
    MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&
    operator= (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&) = default;

    //! Move assigment (shallow move).
    MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&
    operator= (MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&&) = default;

    /// \brief Destructor (virtual for memory safety of derived classes).
    ///
    /// \note To Tpetra developers: See the C++ Core Guidelines C.21
    ///   ("If you define or <tt>=delete</tt> any default operation,
    ///   define or <tt>=delete</tt> them all"), in particular the
    ///   AbstractBase example, for why this destructor declaration
    ///   implies that we need the above four <tt>=default</tt>
    ///   declarations for copy construction, move construction, copy
    ///   assignment, and move assignment.
    virtual ~MultiVector () = default;

    //! Swap contents of \c mv with contents of \c *this.
    void swap (MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& mv);

    //@}
    //! @name Post-construction modification routines
    //@{

  protected:
    /// \brief Whether sumIntoLocalValue and sumIntoGlobalValue should
    ///   use atomic updates by default.
    ///
    /// \warning This is an implementation detail.
    static const bool useAtomicUpdatesByDefault =
#ifdef KOKKOS_ENABLE_SERIAL
      ! std::is_same<execution_space, Kokkos::Serial>::value;
#else
      true;
#endif // KOKKOS_ENABLE_SERIAL

  public:
    /// \brief Replace value in host memory, using global row index.
    ///
    /// Replace the current value at row \c gblRow (a global index)
    /// and column \c col with the given value.  The column index is
    /// zero based.
    ///
    /// This method affects the host memory version of the data.  If
    /// \c device_type is a Kokkos device that has two memory spaces,
    /// and you want to modify the non-host version of the data, you
    /// must access the device View directly by calling
    /// getLocalView().  Please see modify(), sync(), and the
    /// discussion of DualView semantics elsewhere in the
    /// documentation.
    /// This method calls sync_host() before modifying
    /// host data, and modify_host() afterwards.
    ///
    /// This method does not have an "atomic" option like
    /// sumIntoGlobalValue.  This is deliberate.  Replacement is not
    /// commutative, unlike += (modulo rounding error).  Concurrent
    /// calls to replaceGlobalValue on different threads that modify
    /// the same entry/ies have undefined results.  (It's not just
    /// that one thread might win; it's that the value might get
    /// messed up.)
    ///
    /// \param gblRow [in] Global row index of the entry to modify.
    ///   This <i>must</i> be a valid global row index on the calling
    ///   process with respect to the MultiVector's Map.
    /// \param col [in] Column index of the entry to modify.
    /// \param value [in] Incoming value to add to the entry.
    void
    replaceGlobalValue (const GlobalOrdinal gblRow,
                        const size_t col,
                        const impl_scalar_type& value);

    /// \brief Like the above replaceGlobalValue, but only enabled if
    ///   T differs from impl_scalar_type.
    ///
    /// This method only exists if its template parameter \c T and
    /// impl_scalar_type differ, and if it is syntactically possible
    /// to convert \c T to impl_scalar_type.  This method is mainly
    /// useful for backwards compatibility, when the Scalar template
    /// parameter differs from impl_scalar_type.  That is commonly
    /// only the case when Scalar is std::complex<U> for some type U.
    ///
    /// This method affects the host memory version of the data.  If
    /// \c device_type is a Kokkos device that has two memory spaces,
    /// and you want to modify the non-host version of the data, you
    /// must access the device View directly by calling getLocalView().
    /// Please see modify(), sync(), and the discussion of DualView
    /// semantics elsewhere in the documentation.
    /// This method calls sync_host() before modifying
    /// host data, and modify_host() afterwards.
    ///
    /// This method does not have an "atomic" option like
    /// sumIntoGlobalValue.  This is deliberate.  Replacement is not
    /// commutative, unlike += (modulo rounding error).  Concurrent
    /// calls to replaceGlobalValue on different threads that modify
    /// the same entry/ies have undefined results.  (It's not just
    /// that one thread might win; it's that the value might get
    /// messed up.)
    ///
    /// \param gblRow [in] Global row index of the entry to modify.
    ///   This <i>must</i> be a valid global row index on the calling
    ///   process with respect to the MultiVector's Map.
    /// \param col [in] Column index of the entry to modify.
    /// \param value [in] Incoming value to add to the entry.
    template<typename T>
    typename std::enable_if<! std::is_same<T, impl_scalar_type>::value && std::is_convertible<T, impl_scalar_type>::value, void>::type
    replaceGlobalValue (GlobalOrdinal globalRow,
                        size_t col,
                        const T& value)
    {
      replaceGlobalValue (globalRow, col, static_cast<impl_scalar_type> (value));
    }

    /// \brief Update (+=) a value in host memory, using global row index.
    ///
    /// Add the given value to the existing value at row \c gblRow (a
    /// global index) and column \c col.  The column index is zero
    /// based.
    ///
    /// This method affects the host memory version of the data.  If
    /// \c device_type is a Kokkos device that has two memory spaces,
    /// and you want to modify the non-host version of the data, you
    /// must access the device View directly by calling
    /// getLocalView().  Please see modify(), sync(), and the
    /// discussion of DualView semantics elsewhere in the
    /// documentation.
    /// This method calls sync_host() before modifying
    /// host data, and modify_host() afterwards.
    ///
    /// \param gblRow [in] Global row index of the entry to modify.
    ///   This <i>must</i> be a valid global row index on the calling
    ///   process with respect to the MultiVector's Map.
    /// \param col [in] Column index of the entry to modify.
    /// \param value [in] Incoming value to add to the entry.
    /// \param atomic [in] Whether to use an atomic update.  If this
    ///   class' execution space is not Kokkos::Serial, then this is
    ///   true by default, else it is false by default.
    void
    sumIntoGlobalValue (const GlobalOrdinal gblRow,
                        const size_t col,
                        const impl_scalar_type& value,
                        const bool atomic = useAtomicUpdatesByDefault);

    /// \brief Like the above sumIntoGlobalValue, but only enabled if
    ///   T differs from impl_scalar_type.
    ///
    /// This method only exists if its template parameter \c T and
    /// impl_scalar_type differ, and if it is syntactically possible
    /// to convert \c T to impl_scalar_type.  This method is mainly
    /// useful for backwards compatibility, when the Scalar template
    /// parameter differs from impl_scalar_type.  That is commonly
    /// only the case when Scalar is std::complex<U> for some type U.
    ///
    /// This method affects the host memory version of the data.  If
    /// \c device_type is a Kokkos device that has two memory spaces,
    /// and you want to modify the non-host version of the data, you
    /// must access the device View directly by calling
    /// getLocalView().  Please see modify(), sync(), and the
    /// discussion of DualView semantics elsewhere in the
    /// documentation.
    /// This method calls sync_host() before modifying
    /// host data, and modify_host() afterwards.
    ///
    /// \param gblRow [in] Global row index of the entry to modify.
    ///   This <i>must</i> be a valid global row index on the calling
    ///   process with respect to the MultiVector's Map.
    /// \param col [in] Column index of the entry to modify.
    /// \param val [in] Incoming value to add to the entry.
    /// \param atomic [in] Whether to use an atomic update.  If this
    ///   class' execution space is not Kokkos::Serial, then this is
    ///   true by default, else it is false by default.
    template<typename T>
    typename std::enable_if<! std::is_same<T, impl_scalar_type>::value && std::is_convertible<T, impl_scalar_type>::value, void>::type
    sumIntoGlobalValue (const GlobalOrdinal gblRow,
                        const size_t col,
                        const T& val,
                        const bool atomic = useAtomicUpdatesByDefault)
    {
      sumIntoGlobalValue (gblRow, col, static_cast<impl_scalar_type> (val), atomic);
    }

    /// \brief Replace value in host memory, using local (row) index.
    ///
    /// Replace the current value at row \c lclRow (a local index) and
    /// column \c col with the given value.  The column index is zero
    /// based.
    ///
    /// This method affects the host memory version of the data.  If
    /// \c device_type is a Kokkos device that has two memory spaces,
    /// and you want to modify the non-host version of the data, you
    /// must access the device View directly by calling
    /// getLocalView().  Please see modify(), sync(), and the
    /// discussion of DualView semantics elsewhere in the
    /// documentation.
    /// This method calls sync_host() before modifying
    /// host data, and modify_host() afterwards.
    ///
    /// This method does not have an "atomic" option like
    /// sumIntoLocalValue.  This is deliberate.  Replacement is not
    /// commutative, unlike += (modulo rounding error).  Concurrent
    /// calls to replaceLocalValue on different threads that modify
    /// the same entry/ies have undefined results.  (It's not just
    /// that one thread might win; it's that the value might get
    /// messed up.)
    ///
    /// \param lclRow [in] Local row index of the entry to modify.
    ///   Must be a valid local index in this MultiVector's Map on the
    ///   calling process.
    /// \param col [in] Column index of the entry to modify.
    /// \param value [in] Incoming value to add to the entry.
    void
    replaceLocalValue (const LocalOrdinal lclRow,
                       const size_t col,
                       const impl_scalar_type& value);

    /// \brief Like the above replaceLocalValue, but only enabled if
    ///   T differs from impl_scalar_type.
    ///
    /// This method only exists if its template parameter \c T and
    /// impl_scalar_type differ, and if it is syntactically possible
    /// to convert \c T to impl_scalar_type.  This method is mainly
    /// useful for backwards compatibility, when the Scalar template
    /// parameter differs from impl_scalar_type.  That is commonly
    /// only the case when Scalar is std::complex<U> for some type U.
    ///
    /// This method affects the host memory version of the data.  If
    /// \c device_type is a Kokkos device that has two memory spaces,
    /// and you want to modify the non-host version of the data, you
    /// must access the device View directly by calling
    /// getLocalView().  Please see modify(), sync(), and the
    /// discussion of DualView semantics elsewhere in the
    /// documentation.
    /// This method calls sync_host() before modifying
    /// host data, and modify_host() afterwards.
    ///
    /// This method does not have an "atomic" option like
    /// sumIntoLocalValue.  This is deliberate.  Replacement is not
    /// commutative, unlike += (modulo rounding error).  Concurrent
    /// calls to replaceLocalValue on different threads that modify
    /// the same entry/ies have undefined results.  (It's not just
    /// that one thread might win; it's that the value might get
    /// messed up.)
    ///
    /// \param lclRow [in] Local row index of the entry to modify.
    ///   Must be a valid local index in this MultiVector's Map on the
    ///   calling process.
    /// \param col [in] Column index of the entry to modify.
    /// \param val [in] Incoming value to add to the entry.
    template<typename T>
    typename std::enable_if<! std::is_same<T, impl_scalar_type>::value && std::is_convertible<T, impl_scalar_type>::value, void>::type
    replaceLocalValue (const LocalOrdinal lclRow,
                       const size_t col,
                       const T& val)
    {
      replaceLocalValue (lclRow, col, static_cast<impl_scalar_type> (val));
    }

    /// \brief Update (+=) a value in host memory, using local row index.
    ///
    /// Add the given value to the existing value at row \c localRow
    /// (a local index) and column \c col.  The column index is zero
    /// based.
    ///
    /// This method affects the host memory version of the data.  If
    /// \c device_type is a Kokkos device that has two memory spaces,
    /// and you want to modify the non-host version of the data, you
    /// must access the device View directly by calling
    /// getLocalView().  Please see modify(), sync(), and the
    /// discussion of DualView semantics elsewhere in the
    /// documentation.
    /// This method calls sync_host() before modifying
    /// host data, and modify_host() afterwards.
    ///
    /// \param lclRow [in] Local row index of the entry to modify.
    ///   Must be a valid local index in this MultiVector's Map on the
    ///   calling process.
    /// \param col [in] Column index of the entry to modify.
    /// \param val [in] Incoming value to add to the entry.
    /// \param atomic [in] Whether to use an atomic update.  If this
    ///   class' execution space is not Kokkos::Serial, then this is
    ///   true by default, else it is false by default.
    void
    sumIntoLocalValue (const LocalOrdinal lclRow,
                       const size_t col,
                       const impl_scalar_type& val,
                       const bool atomic = useAtomicUpdatesByDefault);

    /// \brief Like the above sumIntoLocalValue, but only enabled if
    ///   T differs from impl_scalar_type.
    ///
    /// This method only exists if its template parameter \c T and
    /// impl_scalar_type differ, and if it is syntactically possible
    /// to convert \c T to impl_scalar_type.  This method is mainly
    /// useful for backwards compatibility, when the Scalar template
    /// parameter differs from impl_scalar_type.  That is commonly
    /// only the case when Scalar is std::complex<U> for some type U.
    ///
    /// This method affects the host memory version of the data.  If
    /// \c device_type is a Kokkos device that has two memory spaces,
    /// and you want to modify the non-host version of the data, you
    /// must access the device View directly by calling
    /// getLocalView().  Please see modify(), sync(), and the
    /// discussion of DualView semantics elsewhere in the
    /// documentation.
    /// This method calls sync_host() before modifying
    /// host data, and modify_host() afterwards.
    ///
    /// \param lclRow [in] Local row index of the entry to modify.
    /// \param col [in] Column index of the entry to modify.
    /// \param val [in] Incoming value to add to the entry.
    /// \param atomic [in] Whether to use an atomic update.  If this
    ///   class' execution space is not Kokkos::Serial, then this is
    ///   true by default, else it is false by default.
    template<typename T>
    typename std::enable_if<! std::is_same<T, impl_scalar_type>::value && std::is_convertible<T, impl_scalar_type>::value, void>::type
    sumIntoLocalValue (const LocalOrdinal lclRow,
                       const size_t col,
                       const T& val,
                       const bool atomic = useAtomicUpdatesByDefault)
    {
      sumIntoLocalValue (lclRow, col, static_cast<impl_scalar_type> (val), atomic);
    }

    //! Set all values in the multivector with the given value.
    void putScalar (const Scalar& value);

    /// \brief Set all values in the multivector with the given value.
    ///
    /// This method only exists if its template parameter \c T and
    /// impl_scalar_type differ, and if it is syntactically possible
    /// to convert \c T to impl_scalar_type.  This method is mainly
    /// useful for backwards compatibility, when the Scalar template
    /// parameter differs from impl_scalar_type.  That is commonly
    /// only the case when Scalar is std::complex<U> for some type U.
    template<typename T>
    typename std::enable_if<! std::is_same<T, impl_scalar_type>::value && std::is_convertible<T, impl_scalar_type>::value, void>::type
    putScalar (const T& value)
    {
      putScalar (static_cast<impl_scalar_type> (value));
    }

    /// \brief Set all values in the multivector to pseudorandom numbers.
    ///
    /// \note Do not expect repeatable results.
    /// \note Behavior of this method may or may not depend on
    ///   external use of the C library routines srand() and rand().
    ///   In particular, setting the seed there may not affect it
    ///   here.
    /// \warning This method does <i>not</i> promise to use a
    ///   distributed-memory parallel pseudorandom number generator.
    ///   Corresponding values on different processes might be
    ///   correlated.  It also does not promise to use a high-quality
    ///   pseudorandom number generator within each process.
    void randomize();

    /// \brief Set all values in the multivector to pseudorandom
    ///   numbers in the given range.
    ///
    /// \note Do not expect repeatable results.
    /// \note Behavior of this method may or may not depend on
    ///   external use of the C library routines srand() and rand().
    ///   In particular, setting the seed there may not affect it
    ///   here.
    /// \warning This method does <i>not</i> promise to use a
    ///   distributed-memory parallel pseudorandom number generator.
    ///   Corresponding values on different processes might be
    ///   correlated.  It also does not promise to use a high-quality
    ///   pseudorandom number generator within each process.
    void randomize (const Scalar& minVal, const Scalar& maxVal);

    /// \brief Replace the underlying Map in place.
    ///
    /// \warning The normal use case of this method, with an input Map
    ///   that is compatible with the object's current Map and has the
    ///   same communicator, is safe.  However, if the input Map has a
    ///   different communicator (with a different number of
    ///   processes, in particular) than this object's current Map,
    ///   the semantics of this method are tricky.  We recommend that
    ///   only experts try the latter use case.
    ///
    /// \pre If the new Map's communicator is similar to the original
    ///   Map's communicator, then the original Map and new Map must
    ///   be compatible: <tt>map->isCompatible (this->getMap ())</tt>.
    ///   "Similar" means that the communicators have the same number
    ///   of processes, though these need not be in the same order
    ///   (have the same assignments of ranks) or represent the same
    ///   communication contexts.  It means the same thing as the
    ///   MPI_SIMILAR return value of MPI_COMM_COMPARE.  See MPI 3.0
    ///   Standard, Section 6.4.1.
    ///
    /// \pre If the new Map's communicator contains more processes
    ///   than the original Map's communicator, then the projection of
    ///   the original Map onto the new communicator must be
    ///   compatible with the new Map.
    ///
    /// \pre If the new Map's communicator contains fewer processes
    ///   than the original Map's communicator, then the projection of
    ///   the new Map onto the original communicator must be
    ///   compatible with the original Map.
    ///
    /// This method replaces this object's Map with the given Map.
    /// This relabels the rows of the multivector using the global IDs
    /// in the input Map.  Thus, it implicitly applies a permutation,
    /// without actually moving data.  If the new Map's communicator
    /// has more processes than the original Map's communicator, it
    /// "projects" the MultiVector onto the new Map by filling in
    /// missing rows with zeros.  If the new Map's communicator has
    /// fewer processes than the original Map's communicator, the
    /// method "forgets about" any rows that do not exist in the new
    /// Map.  (It mathematical terms, if one considers a MultiVector
    /// as a function from one vector space to another, this operation
    /// <i>restricts</i> the range.)
    ///
    /// This method must always be called collectively on the
    /// communicator with the largest number of processes: either this
    /// object's current communicator
    /// (<tt>this->getMap()->getComm()</tt>), or the new Map's
    /// communicator (<tt>map->getComm()</tt>).  If the new Map's
    /// communicator has fewer processes, then the new Map must be
    /// null on processes excluded from the original communicator, and
    /// the current Map must be nonnull on all processes.  If the new
    /// Map has more processes, then it must be nonnull on all those
    /// processes, and the original Map must be null on those
    /// processes which are not in the new Map's communicator.  (The
    /// latter case can only happen to a MultiVector to which a
    /// replaceMap() operation has happened before.)
    ///
    /// \warning This method must always be called as a collective
    ///   operation on all processes in the original communicator
    ///   (<tt>this->getMap ()->getComm ()</tt>).  We reserve the
    ///   right to do checking in debug mode that requires this method
    ///   to be called collectively in order not to deadlock.
    ///
    /// \note This method does <i>not</i> do data redistribution.  If
    ///   you need to move data around, use Import or Export.
    void replaceMap (const Teuchos::RCP<const map_type>& map);

    /// \brief Sum values of a locally replicated multivector across all processes.
    ///
    /// \warning This method may only be called for locally replicated
    ///   MultiVectors.
    ///
    /// \pre isDistributed() == false
    void reduce();

    //@}

    //! @name Get a copy or view of a subset of rows and/or columns
    ///
    /// The following methods get either a (deep) copy or a view
    /// (shallow copy) of a subset of rows and/or columns of the
    /// MultiVector.  They return one of the following:
    ///
    /// <ul>
    /// <li> Another MultiVector </li>
    /// <li> A Kokkos::View or Kokkos::DualView </li>
    /// <li> A Teuchos::ArrayRCP (see the Teuchos Memory Management Classes) </li>
    /// </ul>
    ///
    /// We prefer use of Kokkos classes to Teuchos Memory Management
    /// Classes.  In particular, Teuchos::ArrayRCP reference counts
    /// are not thread safe, while Kokkos::View (and Kokkos::DualView)
    /// reference counts are thread safe.
    ///
    /// Not all of these methods are valid for a particular
    /// MultiVector. For instance, calling a method that accesses a
    /// view of the data in a 1-D format (i.e., get1dView) requires
    /// that the target MultiVector have constant stride.
    ///
    /// This category of methods also includes sync(), modify(), and
    /// getLocalView(), which help MultiVector implement DualView
    /// semantics.
    //@{

    //! Return a MultiVector with copies of selected columns.
    Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    subCopy (const Teuchos::Range1D& colRng) const;

    //! Return a MultiVector with copies of selected columns.
    Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    subCopy (const Teuchos::ArrayView<const size_t>& cols) const;

    //! Return a const MultiVector with const views of selected columns.
    Teuchos::RCP<const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    subView (const Teuchos::Range1D& colRng) const;

    //! Return a const MultiVector with const views of selected columns.
    Teuchos::RCP<const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    subView (const Teuchos::ArrayView<const size_t>& cols) const;

    //! Return a MultiVector with views of selected columns.
    Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    subViewNonConst (const Teuchos::Range1D& colRng);

    //! Return a MultiVector with views of selected columns.
    Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    subViewNonConst (const Teuchos::ArrayView<const size_t>& cols);

    /// \brief Return a const view of a subset of rows.
    ///
    /// Return a const (nonmodifiable) view of this MultiVector
    /// consisting of a subset of the rows, as specified by an offset
    /// and a subset Map of this MultiVector's current row Map.  If
    /// you want X1 or X2 to be nonconst (modifiable) views, use
    /// offsetViewNonConst() with the same arguments.  "View" means
    /// "alias": if the original (this) MultiVector's data change, the
    /// view will see the changed data.
    ///
    /// \param subMap [in] The row Map for the new MultiVector.  This
    ///   must be a subset Map of this MultiVector's row Map.
    /// \param offset [in] The local row offset at which to start the view.
    ///
    /// Suppose that you have a MultiVector X, and you want to view X,
    /// on all processes in X's (MPI) communicator, as split into two
    /// row blocks X1 and X2.  One could express this in Matlab
    /// notation as X = [X1; X2], except that here, X1 and X2 are
    /// views into X, rather than copies of X's data.  This method
    /// assumes that the <i>local</i> indices of X1 and X2 are each
    /// contiguous, and that the local indices of X2 follow those of
    /// X1.  If that is not the case, you cannot use views to divide X
    /// into blocks like this; you must instead use the Import or
    /// Export functionality, which copies the relevant rows of X.
    ///
    /// Here is how you would construct the views X1 and X2.
    /// \code
    /// using Teuchos::RCP;
    /// typedef Tpetra::Map<LO, GO, Node> map_type;
    /// typedef Tpetra::MultiVector<Scalar, LO, GO, Node> MV;
    ///
    /// MV X (...); // the input MultiVector
    /// // ... fill X with data ...
    ///
    /// // Map that on each process in X's communicator,
    /// // contains the global indices of the rows of X1.
    /// RCP<const map_type> map1 (new map_type (...));
    /// // Map that on each process in X's communicator,
    /// // contains the global indices of the rows of X2.
    /// RCP<const map_type> map2 (new map_type (...));
    ///
    /// // Create the first view X1.  The second argument, the offset,
    /// // is the index of the local row at which to start the view.
    /// // X1 is the topmost block of X, so the offset is zero.
    /// RCP<const MV> X1 = X.offsetView (map1, 0);
    ///
    /// // Create the second view X2.  X2 is directly below X1 in X,
    /// // so the offset is the local number of rows in X1.  This is
    /// // the same as the local number of entries in map1.
    /// RCP<const MV> X2 = X.offsetView (map2, X1->getLocalLength ());
    /// \endcode
    ///
    /// It is legal, in the above example, for X1 or X2 to have zero
    /// local rows on any or all process(es).  In that case, the
    /// corresponding Map must have zero local entries on that / those
    /// process(es).  In particular, if X2 has zero local rows on a
    /// process, then the corresponding offset on that process would
    /// be the number of local rows in X (and therefore in X1) on that
    /// process.  This is the only case in which the sum of the local
    /// number of entries in \c subMap (in this case, zero) and the \c
    /// offset may equal the number of local entries in
    /// <tt>*this</tt>.
    Teuchos::RCP<const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    offsetView (const Teuchos::RCP<const map_type>& subMap,
                const size_t offset) const;

    /// \brief Return a nonconst view of a subset of rows.
    ///
    /// Return a nonconst (modifiable) view of this MultiVector
    /// consisting of a subset of the rows, as specified by an offset
    /// and a subset Map of this MultiVector's current row Map.  If
    /// you want X1 or X2 to be const (nonmodifiable) views, use
    /// offsetView() with the same arguments.  "View" means "alias":
    /// if the original (this) MultiVector's data change, the view
    /// will see the changed data, and if the view's data change, the
    /// original MultiVector will see the changed data.
    ///
    /// \param subMap [in] The row Map for the new MultiVector.  This
    ///   must be a subset Map of this MultiVector's row Map.
    /// \param offset [in] The local row offset at which to start the view.
    ///
    /// See the documentation of offsetView() for a code example and
    /// an explanation of edge cases.
    Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    offsetViewNonConst (const Teuchos::RCP<const map_type>& subMap,
                        const size_t offset);

    //! Return a Vector which is a const view of column j.
    Teuchos::RCP<const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    getVector (const size_t j) const;

    //! Return a Vector which is a nonconst view of column j.
    Teuchos::RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    getVectorNonConst (const size_t j);

    //! Const view of the local values in a particular vector of this multivector.
    Teuchos::ArrayRCP<const Scalar> getData (size_t j) const;

    //! View of the local values in a particular vector of this multivector.
    Teuchos::ArrayRCP<Scalar> getDataNonConst (size_t j);

    /// \brief Fill the given array with a copy of this multivector's
    ///   local values.
    ///
    /// \param A [out] View of the array to fill.  We consider A as a
    ///   matrix with column-major storage.
    ///
    /// \param LDA [in] Leading dimension of the matrix A.
    void
    get1dCopy (const Teuchos::ArrayView<Scalar>& A,
               const size_t LDA) const;

    /// \brief Fill the given array with a copy of this multivector's
    ///   local values.
    ///
    /// \param ArrayOfPtrs [out] Array of arrays, one for each column
    ///   of the multivector.  On output, we fill ArrayOfPtrs[j] with
    ///   the data for column j of this multivector.
    void
    get2dCopy (const Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> >& ArrayOfPtrs) const;

    /// \brief Const persisting (1-D) view of this multivector's local values.
    ///
    /// This method assumes that the columns of the multivector are
    /// stored contiguously.  If not, this method throws
    /// std::runtime_error.
    Teuchos::ArrayRCP<const Scalar> get1dView () const;

    //! Return const persisting pointers to values.
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > get2dView () const;

    /// \brief Nonconst persisting (1-D) view of this multivector's local values.
    ///
    /// This method assumes that the columns of the multivector are
    /// stored contiguously.  If not, this method throws
    /// std::runtime_error.
    Teuchos::ArrayRCP<Scalar> get1dViewNonConst ();

    //! Return non-const persisting pointers to values.
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > get2dViewNonConst ();

    /// \brief Return a read-only, up-to-date view of this MultiVector's local data on host.
    /// This requires that there are no live device-space views.
    typename dual_view_type::t_host::const_type getLocalViewHost(Access::ReadOnlyStruct) const;

    /// \brief Return a mutable, up-to-date view of this MultiVector's local data on host.
    /// This requires that there are no live device-space views.
    typename dual_view_type::t_host getLocalViewHost(Access::ReadWriteStruct);

    /// \brief Return a mutable view of this MultiVector's local data on host, assuming all existing data will be overwritten.
    /// This requires that there are no live device-space views.
    typename dual_view_type::t_host getLocalViewHost(Access::OverwriteAllStruct);

    /// \brief Return a read-only, up-to-date view of this MultiVector's local data on device.
    /// This requires that there are no live host-space views.
    typename dual_view_type::t_dev::const_type getLocalViewDevice(Access::ReadOnlyStruct) const;

    /// \brief Return a mutable, up-to-date view of this MultiVector's local data on device.
    /// This requires that there are no live host-space views.
    typename dual_view_type::t_dev getLocalViewDevice(Access::ReadWriteStruct);

    /// \brief Return a mutable view of this MultiVector's local data on device, assuming all existing data will be overwritten.
    /// This requires that there are no live host-space views.
    typename dual_view_type::t_dev getLocalViewDevice(Access::OverwriteAllStruct);

    /// \brief Return the wrapped dual view holding this MultiVector's local data.
    ///
    /// \warning This method is ONLY for use by experts. We highly recommend accessing the local data
    /// by using the member functions getLocalViewHost and getLocalViewDevice.
    wrapped_dual_view_type getWrappedDualView() const;

    //! Whether this MultiVector needs synchronization to the given space.
    template<class TargetDeviceType>
    bool need_sync () const {
      return view_.getDualView().template need_sync<TargetDeviceType> ();
    }

    //! Whether this MultiVector needs synchronization to the host.
    bool need_sync_host () const;

    //! Whether this MultiVector needs synchronization to the device.
    bool need_sync_device () const;

    /// \brief Return a view of the local data on a specific device, with the given access mode.
    ///   The return type is either dual_view_type::t_dev, dual_view_type::t_host, or the const_type of
    ///   one of those.
    ///
    /// \tparam TargetDeviceType The Kokkos Device type whose data to return.
    ///
    /// For example, suppose you create a Tpetra::MultiVector for the
    /// Kokkos::Cuda device, like this:
    /// \code
    /// typedef Tpetra::KokkosCompat::KokkosDeviceWrapperNode<Kokkos::Cuda> > node_type;
    /// typedef Tpetra::Map<int, int, node_type> map_type;
    /// typedef Tpetra::MultiVector<float, int, int, node_type> mv_type;
    ///
    /// RCP<const map_type> map = ...;
    /// mv_type DV (map, 3);
    /// \endcode
    /// If you want to get the CUDA device Kokkos::View as read-write, do this:
    /// \code
    /// typedef typename mv_type::dual_view_type dual_view_type;
    /// typedef typename dual_view_type::t_dev device_view_type;
    /// device_view_type cudaView = DV.getLocalView<Kokkos::Cuda> (Access::ReadWrite);
    /// \endcode
    /// and if you want to get the host mirror of that View, do this:
    /// \code
    /// typedef typename dual_view_type::host_mirror_space host_execution_space;
    /// typedef typename dual_view_type::t_host host_view_type;
    /// host_view_type hostView = DV.getLocalView<host_execution_space> (Access::ReadWrite);
    /// \endcode
    template<class TargetDeviceType>
    typename std::remove_reference<decltype(std::declval<dual_view_type>().template view<TargetDeviceType>())>::type::const_type
    getLocalView (Access::ReadOnlyStruct s) const
    {
      return view_.template getView<TargetDeviceType>(s);
    }


    template<class TargetDeviceType>
    typename std::remove_reference<decltype(std::declval<dual_view_type>().template view<TargetDeviceType>())>::type
    getLocalView (Access::ReadWriteStruct s)
    {
      return view_.template getView<TargetDeviceType>(s);
    }

    template<class TargetDeviceType>
    typename std::remove_reference<decltype(std::declval<dual_view_type>().template view<TargetDeviceType>())>::type
    getLocalView (Access::OverwriteAllStruct s)
    {
      return view_.template getView<TargetDeviceType>(s);
    }

    //@}
    //! @name Mathematical methods
    //@{

    /// \brief Compute the dot product of each corresponding pair of
    ///   vectors (columns) in A and B.
    ///
    /// The "dot product" is the standard Euclidean inner product.  If
    /// the type of entries of the vectors (impl_scalar_type) is complex,
    /// then A is transposed, not <tt>*this</tt>.  For example, if x
    /// and y each have one column, then <tt>x.dot (y, dots)</tt>
    /// computes \f$y^* x = \bar{y}^T x = \sum_i \bar{y}_i \cdot x_i\f$.
    ///
    /// \pre <tt>*this</tt> and A have the same number of columns (vectors).
    /// \pre \c dots has at least as many entries as the number of columns in A.
    ///
    /// \post <tt>dots[j] == (this->getVector[j])->dot (* (A.getVector[j]))</tt>
    void
    dot (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
         const Teuchos::ArrayView<dot_type>& dots) const;

    /// \brief Compute the dot product of each corresponding pair of
    ///   vectors (columns) in A and B.
    /// \tparam T The output type of the dot products.
    ///
    /// This method only exists if dot_type and T are different types.
    /// For example, if impl_scalar_type and dot_type differ, then this
    /// method ensures backwards compatibility with the previous
    /// interface (that returned dot products as impl_scalar_type rather
    /// than as dot_type).  The complicated \c enable_if expression
    /// just ensures that the method only exists if dot_type and T are
    /// different types; the method still returns \c void, as above.
    template <typename T>
    typename std::enable_if< ! (std::is_same<dot_type, T>::value), void >::type
    dot (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
         const Teuchos::ArrayView<T> &dots) const
    {
      const size_t sz = static_cast<size_t> (dots.size ());
      Teuchos::Array<dot_type> dts (sz);
      this->dot (A, dts);
      for (size_t i = 0; i < sz; ++i) {
        // If T and dot_type differ, this does an implicit conversion.
        dots[i] = dts[i];
      }
    }

    //! Like the above dot() overload, but for std::vector output.
    template <typename T>
    typename std::enable_if< ! (std::is_same<dot_type, T>::value), void >::type
    dot (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
         std::vector<T>& dots) const
    {
      const size_t sz = dots.size ();
      Teuchos::Array<dot_type> dts (sz);
      this->dot (A, dts);
      for (size_t i = 0; i < sz; ++i) {
        // If T and dot_type differ, this does an implicit conversion.
        dots[i] = dts[i];
      }
    }

    /// \brief Compute the dot product of each corresponding pair of
    ///   vectors (columns) in A and B, storing the result in a device
    ///   View.
    ///
    /// The "dot product" is the standard Euclidean inner product.  If
    /// the type of entries of the vectors (impl_scalar_type) is complex,
    /// then A is transposed, not <tt>*this</tt>.  For example, if x
    /// and y each have one column, then <tt>x.dot (y, dots)</tt>
    /// computes \f$y^* x = \bar{y}^T x = \sum_i \bar{y}_i \cdot x_i\f$.
    ///
    /// \param A [in] MultiVector with which to dot \c *this.
    /// \param dots [out] Device View with getNumVectors() entries.
    ///
    /// \pre <tt>this->getNumVectors () == A.getNumVectors ()</tt>
    /// \pre <tt>dots.extent (0) == A.getNumVectors ()</tt>
    ///
    /// \post <tt>dots(j) == (this->getVector[j])->dot (* (A.getVector[j]))</tt>
    void
    dot (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
         const Kokkos::View<dot_type*, Kokkos::HostSpace>& norms) const;

    template<class ViewType>
    void
    dot (typename std::enable_if<std::is_same<typename ViewType::value_type,dot_type>::value &&
                                 std::is_same<typename ViewType::memory_space,typename device_type::memory_space>::value,
                                 const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>::type& A,
         const ViewType& dots) const {
      const Kokkos::View<dot_type*, Kokkos::HostSpace> h_dots("Tpetra::Dots",dots.extent(0));
      this->dot (A, h_dots);
      // DEEP_COPY REVIEW - NOT TESTED
      Kokkos::deep_copy(dots,h_dots);
    }

    /// \brief Compute the dot product of each corresponding pair of
    ///   vectors (columns) in A and B, storing the result in a device
    ///   view.
    /// \tparam T The output type of the dot products.
    ///
    /// This method only exists if dot_type and T are different types.
    /// For example, if Scalar and dot_type differ, then this method
    /// ensures backwards compatibility with the previous interface
    /// (that returned dot products as Scalar rather than as
    /// dot_type).  The complicated \c enable_if expression just
    /// ensures that the method only exists if dot_type and T are
    /// different types; the method still returns \c void, as above.
    template <typename T>
    typename std::enable_if< ! (std::is_same<dot_type, T>::value), void >::type
    dot (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
         const Kokkos::View<T*, device_type>& dots) const
    {
      const size_t numDots = dots.extent (0);
      Kokkos::View<dot_type*, device_type> dts ("MV::dot tmp", numDots);
      // Call overload that takes a Kokkos::View<dot_type*, device_type>.
      this->dot (A, dts);
      // FIXME (mfh 14 Jul 2014) Does this actually work if dot_type
      // and T differ?  We would need a test for this, but only the
      // Sacado and Stokhos packages are likely to care about this use
      // case.  It could also come up for Kokkos::complex ->
      // std::complex conversions, but those two implementations
      // should generally be bitwise compatible.
      // CT: no this can't possible work .....
      // DEEP_COPY REVIEW - NOT TESTED
      Kokkos::deep_copy (dots, dts);
    }

    //! Put element-wise absolute values of input Multi-vector in target: A = abs(this)
    void abs (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A);

    //! Put element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
    void reciprocal (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A);

    /// \brief Scale in place: <tt>this = alpha*this</tt>.
    ///
    /// Replace this MultiVector with alpha times this MultiVector.
    /// This method will always multiply, even if alpha is zero.  That
    /// means, for example, that if \c *this contains NaN entries
    /// before calling this method, the NaN entries will remain after
    /// this method finishes.
    void scale (const Scalar& alpha);

    /// \brief Scale each column in place: <tt>this[j] = alpha[j]*this[j]</tt>.
    ///
    /// Replace each column j of this MultiVector with
    /// <tt>alpha[j]</tt> times the current column j of this
    /// MultiVector.  This method will always multiply, even if all
    /// the entries of alpha are zero.  That means, for example, that
    /// if \c *this contains NaN entries before calling this method,
    /// the NaN entries will remain after this method finishes.
    void scale (const Teuchos::ArrayView<const Scalar>& alpha);

    /// \brief Scale each column in place: <tt>this[j] = alpha[j]*this[j]</tt>.
    ///
    /// Replace each column j of this MultiVector with
    /// <tt>alpha[j]</tt> times the current column j of this
    /// MultiVector.  This method will always multiply, even if all
    /// the entries of alpha are zero.  That means, for example, that
    /// if \c *this contains NaN entries before calling this method,
    /// the NaN entries will remain after this method finishes.
    void scale (const Kokkos::View<const impl_scalar_type*, device_type>& alpha);

    /// \brief Scale in place: <tt>this = alpha * A</tt>.
    ///
    /// Replace this MultiVector with scaled values of A.  This method
    /// will always multiply, even if alpha is zero.  That means, for
    /// example, that if \c *this contains NaN entries before calling
    /// this method, the NaN entries will remain after this method
    /// finishes.  It is legal for the input A to alias this
    /// MultiVector.
    void
    scale (const Scalar& alpha,
           const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A);

    /// \brief Update: <tt>this = beta*this + alpha*A</tt>.
    ///
    /// Update this MultiVector with scaled values of A.  If beta is
    /// zero, overwrite \c *this unconditionally, even if it contains
    /// NaN entries.  It is legal for the input A to alias this
    /// MultiVector.
    void
    update (const Scalar& alpha,
            const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
            const Scalar& beta);

    /// \brief Update: <tt>this = gamma*this + alpha*A + beta*B</tt>.
    ///
    /// Update this MultiVector with scaled values of A and B.  If
    /// gamma is zero, overwrite \c *this unconditionally, even if it
    /// contains NaN entries.  It is legal for the inputs A or B to
    /// alias this MultiVector.
    void
    update (const Scalar& alpha,
            const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
            const Scalar& beta,
            const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
            const Scalar& gamma);

    /// \brief Compute the one-norm of each vector (column), storing
    ///   the result in a host view.
    ///
    /// \param norms [out] Host View with getNumVectors() entries.
    ///
    /// \pre <tt>norms.extent (0) == this->getNumVectors ()</tt>
    /// \post <tt>norms(j) == (this->getVector[j])->norm1 (* (A.getVector[j]))</tt>
    ///
    /// The one-norm of a vector is the sum of the magnitudes of the
    /// vector's entries.  On exit, norms(j) is the one-norm of column
    /// j of this MultiVector.
    void
    norm1 (const Kokkos::View<mag_type*, Kokkos::HostSpace>& norms) const;

    template<class ViewType>
      typename std::enable_if<std::is_same<typename ViewType::value_type,mag_type>::value &&
                              std::is_same<typename ViewType::memory_space,typename device_type::memory_space>::value>::type
    norm1 (const ViewType& norms) const {
      // FIXME (mfh 11 Apr 2019) The enable_ifs make it useless for
      // this method to be templated.  (It only exists in case
      // HostSpace = device_type::memory_space.)
      using host_norms_view_type = Kokkos::View<mag_type*, Kokkos::HostSpace>;
      host_norms_view_type h_norms ("Tpetra::MV::h_norms", norms.extent (0));
      this->norm1 (h_norms);
      // DEEP_COPY REVIEW - HOST-TO-DEVICE
      Kokkos::deep_copy (execution_space(), norms, h_norms);
    }

    /// \brief Compute the one-norm of each vector (column), storing
    ///   the result in a device view.
    /// \tparam T The output type of the dot products.
    ///
    /// See the above norm1() method for documentation.
    ///
    /// This method only exists if mag_type and T are different types.
    /// For example, if Teuchos::ScalarTraits<Scalar>::magnitudeType
    /// and mag_type differ, then this method ensures backwards
    /// compatibility with the previous interface (that returned norms
    /// products as Teuchos::ScalarTraits<Scalar>::magnitudeType
    /// rather than as mag_type).  The complicated \c enable_if
    /// expression just ensures that the method only exists if
    /// mag_type and T are different types; the method still returns
    /// \c void, as above.
    template <typename T>
    typename std::enable_if< ! (std::is_same<mag_type, T>::value), void >::type
    norm1 (const Kokkos::View<T*, device_type>& norms) const
    {
      const size_t numNorms = norms.extent (0);
      Kokkos::View<mag_type*, device_type> tmpNorms ("MV::norm1 tmp", numNorms);
      // Call overload that takes a Kokkos::View<mag_type*, device_type>.
      this->norm1 (tmpNorms);
      // FIXME (mfh 15 Jul 2014) Does this actually work if mag_type
      // and T differ?  We would need a test for this, but only the
      // Sacado and Stokhos packages are likely to care about this use
      // case.  It could also come up with Kokkos::complex ->
      // std::complex conversion.
      // DEEP_COPY REVIEW - NOT TESTED
      Kokkos::deep_copy (norms, tmpNorms);
    }

    /// \brief Compute the one-norm of each vector (column).
    ///
    /// See the uppermost norm1() method above for documentation.
    void norm1 (const Teuchos::ArrayView<mag_type>& norms) const;

    /// \brief Compute the one-norm of each vector (column).
    /// \tparam T The output type of the norms.
    ///
    /// See the uppermost norm1() method above for documentation.
    ///
    /// This method only exists if mag_type and T are different types.
    /// For example, if Teuchos::ScalarTraits<Scalar>::magnitudeType
    /// and mag_type differ, then this method ensures backwards
    /// compatibility with the previous interface (that returned norms
    /// as Teuchos::ScalarTraits<Scalar>::magnitudeType
    /// rather than as mag_type).  The complicated \c enable_if
    /// expression just ensures that the method only exists if
    /// mag_type and T are different types; the method still returns
    /// \c void, as above.
    template <typename T>
    typename std::enable_if< ! (std::is_same<mag_type,T>::value), void >::type
    norm1 (const Teuchos::ArrayView<T>& norms) const
    {
      typedef typename Teuchos::ArrayView<T>::size_type size_type;
      const size_type sz = norms.size ();
      Teuchos::Array<mag_type> theNorms (sz);
      this->norm1 (theNorms);
      for (size_type i = 0; i < sz; ++i) {
        // If T and mag_type differ, this does an implicit conversion.
        norms[i] = theNorms[i];
      }
    }

    /// \brief Compute the two-norm of each vector (column), storing
    ///   the result in a host View.
    ///
    /// \param norms [out] Host View with getNumVectors() entries.
    ///
    /// \pre <tt>norms.extent (0) == this->getNumVectors ()</tt>
    /// \post <tt>norms(j) == (this->getVector[j])->dot (* (A.getVector[j]))</tt>
    ///
    /// The two-norm of a vector is the standard Euclidean norm, the
    /// square root of the sum of squares of the magnitudes of the
    /// vector's entries.  On exit, norms(k) is the two-norm of column
    /// k of this MultiVector.
    void
    norm2 (const Kokkos::View<mag_type*, Kokkos::HostSpace>& norms) const;

    template<class ViewType>
      typename std::enable_if<std::is_same<typename ViewType::value_type,mag_type>::value &&
                              std::is_same<typename ViewType::memory_space,typename device_type::memory_space>::value>::type
    norm2 (const ViewType& norms) const {
      // FIXME (mfh 11 Apr 2019) The enable_ifs make it useless for
      // this method to be templated.  (It only exists in case
      // HostSpace = device_type::memory_space.)
      using host_norms_view_type = Kokkos::View<mag_type*, Kokkos::HostSpace>;
      host_norms_view_type h_norms ("Tpetra::MV::h_norms", norms.extent (0));
      this->norm2 (h_norms);
      // DEEP_COPY REVIEW - NOT TESTED
      Kokkos::deep_copy (norms, h_norms);
    }

    /// \brief Compute the two-norm of each vector (column), storing
    ///   the result in a device view.
    ///
    /// See the above norm2() method for documentation.
    ///
    /// This method only exists if mag_type and T are different types.
    /// For example, if Teuchos::ScalarTraits<Scalar>::magnitudeType
    /// and mag_type differ, then this method ensures backwards
    /// compatibility with the previous interface (that returned norms
    /// as Teuchos::ScalarTraits<Scalar>::magnitudeType rather than as
    /// mag_type).  The complicated \c enable_if expression just
    /// ensures that the method only exists if mag_type and T are
    /// different types; the method still returns \c void, as above.
    template<typename T>
    typename std::enable_if< ! (std::is_same<mag_type, T>::value), void >::type
    norm2 (const Kokkos::View<T*, device_type>& norms) const
    {
      const size_t numNorms = norms.extent (0);
      Kokkos::View<mag_type*, device_type> theNorms ("MV::norm2 tmp", numNorms);
      // Call overload that takes a Kokkos::View<mag_type*, device_type>.
      this->norm2 (theNorms);
      // FIXME (mfh 14 Jul 2014) Does this actually work if mag_type
      // and T differ?  We would need a test for this, but only the
      // Sacado and Stokhos packages are likely to care about this use
      // case.  This could also come up with Kokkos::complex ->
      // std::complex conversion.
      // DEEP_COPY REVIEW - NOT TESTED
      Kokkos::deep_copy (norms, theNorms);
    }

    /// \brief Compute the two-norm of each vector (column).
    ///
    /// See the uppermost norm2() method above for documentation.
    void norm2 (const Teuchos::ArrayView<mag_type>& norms) const;

    /// \brief Compute the two-norm of each vector (column).
    /// \tparam T The output type of the norms.
    ///
    /// See the uppermost norm2() method above for documentation.
    ///
    /// This method only exists if mag_type and T are different types.
    /// For example, if Teuchos::ScalarTraits<Scalar>::magnitudeType
    /// and mag_type differ, then this method ensures backwards
    /// compatibility with the previous interface (that returned norms
    /// products as Teuchos::ScalarTraits<Scalar>::magnitudeType
    /// rather than as mag_type).  The complicated \c enable_if
    /// expression just ensures that the method only exists if
    /// mag_type and T are different types; the method still returns
    /// \c void, as above.
    template <typename T>
    typename std::enable_if< ! (std::is_same<mag_type,T>::value), void >::type
    norm2 (const Teuchos::ArrayView<T>& norms) const
    {
      typedef typename Teuchos::ArrayView<T>::size_type size_type;
      const size_type sz = norms.size ();
      Teuchos::Array<mag_type> theNorms (sz);
      this->norm2 (theNorms);
      for (size_type i = 0; i < sz; ++i) {
        // If T and mag_type differ, this does an implicit conversion.
        norms[i] = theNorms[i];
      }
    }

    /// \brief Compute the infinity-norm of each vector (column),
    ///   storing the result in a host View.
    ///
    /// The infinity-norm of a vector is the maximum of the magnitudes
    /// of the vector's entries.  On exit, norms(j) is the
    /// infinity-norm of column j of this MultiVector.
    void normInf (const Kokkos::View<mag_type*, Kokkos::HostSpace>& norms) const;

    template<class ViewType>
      typename std::enable_if<std::is_same<typename ViewType::value_type,mag_type>::value &&
                              std::is_same<typename ViewType::memory_space,typename device_type::memory_space>::value>::type
    normInf (const ViewType& norms) const {
      // FIXME (mfh 11 Apr 2019) The enable_ifs make it useless for
      // this method to be templated.  (It only exists in case
      // HostSpace = device_type::memory_space.)
      using host_norms_view_type = Kokkos::View<mag_type*, Kokkos::HostSpace>;
      host_norms_view_type h_norms ("Tpetra::MV::h_norms", norms.extent (0));
      this->normInf (h_norms);
      // DEEP_COPY REVIEW - HOST-TO-DEVICE
      Kokkos::deep_copy (execution_space(), norms, h_norms);
    }

    /// \brief Compute the infinity-norm of each vector (column),
    ///   storing the result in a device view.
    ///
    /// See the above normInf() method for documentation.
    ///
    /// This method only exists if mag_type and T are different types.
    /// For example, if Teuchos::ScalarTraits<Scalar>::magnitudeType
    /// and mag_type differ, then this method ensures backwards
    /// compatibility with the previous interface (that returned norms
    /// as Teuchos::ScalarTraits<Scalar>::magnitudeType rather than as
    /// mag_type).  The complicated \c enable_if expression just
    /// ensures that the method only exists if mag_type and T are
    /// different types; the method still returns \c void, as above.
    template<typename T>
    typename std::enable_if< ! (std::is_same<mag_type, T>::value), void >::type
    normInf (const Kokkos::View<T*, device_type>& norms) const
    {
      const size_t numNorms = norms.extent (0);
      Kokkos::View<mag_type*, device_type> theNorms ("MV::normInf tmp", numNorms);
      // Call overload that takes a Kokkos::View<mag_type*, device_type>.
      this->normInf (theNorms);
      // FIXME (mfh 15 Jul 2014) Does this actually work if mag_type
      // and T differ?  We would need a test for this, but only the
      // Sacado and Stokhos packages are likely to care about this use
      // case.  This could also come up with Kokkos::complex ->
      // std::complex conversion.
      // DEEP_COPY REVIEW - NOT TESTED
      Kokkos::deep_copy (norms, theNorms);
    }

    /// \brief Compute the infinity-norm of each vector (column),
    ///   storing the result in a Teuchos::ArrayView.
    ///
    /// See the uppermost normInf() method above for documentation.
    void normInf (const Teuchos::ArrayView<mag_type>& norms) const;

    /// \brief Compute the infinity-norm of each vector (column),
    ///   storing the result in a Teuchos::ArrayView.
    /// \tparam T The output type of the norms.
    ///
    /// See the uppermost normInf() method above for documentation.
    ///
    /// This method only exists if mag_type and T are different types.
    /// For example, if Teuchos::ScalarTraits<Scalar>::magnitudeType
    /// and mag_type differ, then this method ensures backwards
    /// compatibility with the previous interface (that returned norms
    /// products as Teuchos::ScalarTraits<Scalar>::magnitudeType
    /// rather than as mag_type).  The complicated \c enable_if
    /// expression just ensures that the method only exists if
    /// mag_type and T are different types; the method still returns
    /// \c void, as above.
    template <typename T>
    typename std::enable_if< ! (std::is_same<mag_type,T>::value), void >::type
    normInf (const Teuchos::ArrayView<T>& norms) const
    {
      typedef typename Teuchos::ArrayView<T>::size_type size_type;
      const size_type sz = norms.size ();
      Teuchos::Array<mag_type> theNorms (sz);
      this->norm2 (theNorms);
      for (size_type i = 0; i < sz; ++i) {
        // If T and mag_type differ, this does an implicit conversion.
        norms[i] = theNorms[i];
      }
    }


    /// \brief Compute mean (average) value of each column.
    ///
    /// The outcome of this routine is undefined for non-floating
    /// point scalar types (e.g., int).
    void meanValue (const Teuchos::ArrayView<impl_scalar_type>& means) const;

    template <typename T>
    typename std::enable_if<! std::is_same<impl_scalar_type, T>::value, void>::type
    meanValue (const Teuchos::ArrayView<T>& means) const
    {
      typedef typename Teuchos::Array<T>::size_type size_type;
      const size_type numMeans = means.size ();

      Teuchos::Array<impl_scalar_type> theMeans (numMeans);
      this->meanValue (theMeans ());
      for (size_type k = 0; k < numMeans; ++k) {
        means[k] = static_cast<T> (theMeans[k]);
      }
    }

    /// \brief Matrix-matrix multiplication: <tt>this = beta*this + alpha*op(A)*op(B)</tt>.
    ///
    /// If beta is zero, overwrite \c *this unconditionally, even if
    /// it contains NaN entries.  This imitates the semantics of
    /// analogous BLAS routines like DGEMM.
    void
    multiply (Teuchos::ETransp transA,
              Teuchos::ETransp transB,
              const Scalar& alpha,
              const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
              const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
              const Scalar& beta);

    /// \brief Multiply a Vector A elementwise by a MultiVector B.
    ///
    /// Compute <tt>this = scalarThis * this + scalarAB * B @ A</tt>
    /// where <tt>@</tt> denotes element-wise multiplication.  In
    /// pseudocode, if C denotes <tt>*this</tt> MultiVector:
    /// \code
    /// C(i,j) = scalarThis * C(i,j) + scalarAB * B(i,j) * A(i,1);
    /// \endcode
    /// for all rows i and columns j of C.
    ///
    /// B must have the same dimensions as <tt>*this</tt>, while A
    /// must have the same number of rows but a single column.
    ///
    /// We do not require that A, B, and <tt>*this</tt> have
    /// compatible Maps, as long as the number of rows in A, B, and
    /// <tt>*this</tt> on each process is the same.  For example, one
    /// or more of these vectors might have a locally replicated Map,
    /// or a Map with a local communicator (<tt>MPI_COMM_SELF</tt>).
    /// This case may occur in block relaxation algorithms when
    /// applying a diagonal scaling.
    void
    elementWiseMultiply (Scalar scalarAB,
                         const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                         const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
                         Scalar scalarThis);
    //@}
    //! @name Attribute access functions
    //@{

    //! Number of columns in the multivector.
    size_t getNumVectors() const;

    //! Local number of rows on the calling process.
    size_t getLocalLength() const;

    //! Global number of rows in the multivector.
    global_size_t getGlobalLength() const;

    /// \brief Stride between columns in the multivector.
    ///
    /// This is only meaningful if \c isConstantStride() returns true.
    ///
    /// \warning This may be different on different processes.
    size_t getStride() const;

    /// \brief Whether this multivector has constant stride between columns.
    ///
    /// \warning This may be different on different processes.
    bool isConstantStride() const;

    /// \brief Whether this multivector's memory might alias other. This is conservative: if either this or other
    ///     is not constant stride, then it simply checks whether the contiguous memory allocations overlap. It
    ///     doesn't check whether the sets of columns overlap. This is a symmetric relation: X.aliases(Y) == Y.aliases(X).
    bool aliases(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& other) const;

    //@}

    //! @name Overridden from Teuchos::Describable
    //@{

    //! A simple one-line description of this object.
    virtual std::string description() const override;

    /// \brief Print the object with the given verbosity level to a FancyOStream.
    ///
    /// \param out [out] Output stream to which to print.  For
    ///   verbosity levels VERB_LOW and lower, only the process with
    ///   rank 0 ("Proc 0") in the MultiVector's communicator prints.
    ///   For verbosity levels strictly higher than VERB_LOW, all
    ///   processes in the communicator need to be able to print to
    ///   the output stream.
    ///
    /// \param verbLevel [in] Verbosity level.  The default verbosity
    ///   (verbLevel=VERB_DEFAULT) is VERB_LOW.
    ///
    /// The amount and content of what this method prints depends on
    /// the verbosity level.  In the list below, each higher level
    /// includes all the content of the previous levels, as well as
    /// its own content.
    ///
    /// - VERB_LOW: Only Proc 0 prints; it prints the same thing as \c
    ///   description().
    /// - VERB_MEDIUM: Each process prints its local length (the
    ///   number of rows that it owns).
    /// - VERB_HIGH: Each process prints whether the multivector has
    ///   constant stride (see \c isConstantStride()), and if so, what
    ///   that stride is.  (Stride may differ on different processes.)
    /// - VERB_EXTREME: Each process prints the values in its local
    ///   part of the multivector.  This will print out as many rows
    ///   of data as the global number of rows in the multivector, so
    ///   beware.
    virtual void
    describe (Teuchos::FancyOStream& out,
              const Teuchos::EVerbosityLevel verbLevel =
              Teuchos::Describable::verbLevel_default) const override;
    //@}

    /// \brief Remove processes owning zero rows from the Map and their communicator.
    ///
    /// \warning This method is ONLY for use by experts.  We highly
    ///   recommend using the nonmember function of the same name
    ///   defined in Tpetra_DistObject_decl.hpp.
    ///
    /// \warning We make NO promises of backwards compatibility.
    ///   This method may change or disappear at any time.
    ///
    /// \param newMap [in] This <i>must</i> be the result of calling
    ///   the removeEmptyProcesses() method on the row Map.  If it
    ///   is not, this method's behavior is undefined.  This pointer
    ///   will be null on excluded processes.
    virtual void
    removeEmptyProcessesInPlace (const Teuchos::RCP<const map_type>& newMap) override;

    /// \brief Set whether this has copy (copyOrView = Teuchos::Copy)
    ///   or view (copyOrView = Teuchos::View) semantics.
    ///
    /// \warning The Kokkos refactor version of MultiVector
    ///   <i>only</i> implements view semantics.  If you attempt to
    ///   call this method with copyOrView == Teuchos::Copy, it will
    ///   throw std::invalid_argument.
    ///
    /// \warning This method is only for expert use.  It may change or
    ///   disappear at any time.
    void setCopyOrView (const Teuchos::DataAccess copyOrView) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        copyOrView == Teuchos::Copy, std::invalid_argument,
        "Tpetra::MultiVector::setCopyOrView: The Kokkos refactor version of "
        "MultiVector _only_ implements view semantics.  You may not call this "
        "method with copyOrView = Teuchos::Copy.  The only valid argument is "
        "Teuchos::View.");
    }

    /// \brief Get whether this has copy (copyOrView = Teuchos::Copy)
    ///   or view (copyOrView = Teuchos::View) semantics.
    ///
    // This method ONLY exists for the circa 2014 "Kokkos refactor"
    // effort.  It ALWAYS returns Teuchos::View.
    ///
    /// \warning This method is only for expert use.  It may change or
    ///   disappear at any time.
    Teuchos::DataAccess getCopyOrView () const {
      return Teuchos::View;
    }

    /// \brief Copy the contents of \c src into \c *this (deep copy).
    ///
    /// \param src [in] Source MultiVector (input of the deep copy).
    ///
    /// \pre <tt> ! src.getMap ().is_null () && ! this->getMap ().is_null () </tt>
    /// \pre <tt> src.getMap ()->isCompatible (* (this->getMap ()) </tt>
    ///
    /// \post Any outstanding views of \c src or \c *this remain valid.
    ///
   /// \note To implementers: The postcondition implies that the
    ///   implementation must not reallocate any memory of \c *this,
    ///   or otherwise change its dimensions.  This is <i>not</i> an
    ///   assignment operator; it does not change anything in \c *this
    ///   other than the contents of storage.
    void
    assign (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& src);

    /// \brief Return another MultiVector with the same entries, but
    ///   converted to a different Scalar type \c T.
    template <class T>
    Teuchos::RCP<MultiVector<T, LocalOrdinal, GlobalOrdinal, Node> >
    convert () const;


    // \brief Checks to see if the local length, number of vectors and size of Scalar type match
    /// \param src [in] MultiVector
    ///
    /// \pre <tt> ! vec.getMap ().is_null () && ! this->getMap ().is_null () </tt>
    /// \pre <tt> vec.getMap ()->isCompatible (* (this->getMap ()) </tt>
    ///
    /// \post Any outstanding views of \c src or \c *this remain valid.
    ///
    bool isSameSize(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & vec) const;

  private:
    //! The type of the base class of this class.
    using base_type = DistObject<scalar_type, local_ordinal_type,
                                 global_ordinal_type, node_type>;

  protected:
    template <class DS, class DL, class DG, class DN,
              class SS, class SL, class SG, class SN>
    friend void
    ::Tpetra::deep_copy (MultiVector<DS, DL, DG, DN>& dst,
                         const MultiVector<SS, SL, SG, SN>& src);

    /// \brief The Kokkos::DualView containing the MultiVector's data.
    ///
    /// This has to be declared \c mutable, so that get1dView() can
    /// retain its current \c const marking, even though it has always
    /// implied a device->host synchronization.  Lesson to the reader:
    /// Use \c const sparingly!
    mutable wrapped_dual_view_type view_;

    /// \brief Indices of columns this multivector is viewing.
    ///
    /// If this array has nonzero size, then this multivector is a
    /// view of another multivector (the "original" multivector).  In
    /// that case, whichVectors_ contains the indices of the columns
    /// of the original multivector.  Furthermore, isConstantStride()
    /// returns false in this case.
    ///
    /// If this array has zero size, then this multivector is not a
    /// view of any other multivector.  Furthermore, the stride
    /// between columns of this multivector is a constant: thus,
    /// isConstantStride() returns true.
    Teuchos::Array<size_t> whichVectors_;

    template<class SC, class LO, class GO, class NT>
    friend ::Teuchos::ArrayView<const size_t> getMultiVectorWhichVectors (const ::Tpetra::MultiVector<SC, LO, GO, NT>& X);

    //@}
    //! @name Misc. implementation details
    //@{

    /// \brief Implementation of description() for this class, and its
    ///   subclass Vector.
    ///
    /// \param className [in] Name of the class calling this method:
    ///   Either "Tpetra::MultiVector" or "Tpetra::Vector" (no quotes
    ///   in the string, in either case).
    std::string
    descriptionImpl (const std::string& className) const;

    /// \brief Print the calling process' verbose describe()
    ///   information to the returned string.
    ///
    /// This is an implementation detail of describe().
    ///
    /// \param vl [in] Verbosity level with which to print.
    std::string
    localDescribeToString (const Teuchos::EVerbosityLevel vl) const;

    /// \brief Implementation of describe() for this class, and its
    ///   subclass Vector.
    ///
    /// \param out [out] Output stream to which to write.  Only
    ///   Process 0 in this object's communicator may write to the
    ///   output stream.
    ///
    /// \param className [in] Name of the class calling this method.
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

    // Return true if and only if VectorIndex is a valid column index.
    bool vectorIndexOutOfRange (const size_t VectorIndex) const;

    /// \brief Persisting view of j-th column in the given ArrayRCP.
    ///
    /// This method considers isConstantStride().  The ArrayRCP may
    /// correspond either to a compute buffer or a host view.
    template <class T>
    Teuchos::ArrayRCP<T>
    getSubArrayRCP (Teuchos::ArrayRCP<T> arr, size_t j) const;

    //! "Original" number of rows in the (local) data.
    size_t getOrigNumLocalRows () const;

    //! "Original" number of columns in the (local) data.
    size_t getOrigNumLocalCols () const;

    //@}
    //! @name Implementation of Tpetra::DistObject
    //@{

    /// \typedef buffer_device_type
    /// \brief Kokkos::Device specialization for communication buffers.
    ///
    /// See #1088 for why this is not just <tt>device_type::device_type</tt>.
    using buffer_device_type =
      typename DistObject<scalar_type,
                          local_ordinal_type,
                          global_ordinal_type,
                          node_type>::buffer_device_type;

    /// \brief Whether data redistribution between \c sourceObj and this object is legal.
    ///
    /// This method is called in DistObject::doTransfer() to check
    /// whether data redistribution between the two objects is legal.
    virtual bool
    checkSizes (const SrcDistObject& sourceObj) override;

    //! Number of packets to send per LID
    virtual size_t constantNumberOfPackets () const override;

  // clang-format on
  virtual void copyAndPermute(
      const SrcDistObject &sourceObj, const size_t numSameIDs,
      const Kokkos::DualView<const local_ordinal_type *, buffer_device_type>
          &permuteToLIDs,
      const Kokkos::DualView<const local_ordinal_type *, buffer_device_type>
          &permuteFromLIDs,
      const CombineMode CM, const execution_space &space) override;

  virtual void copyAndPermute(
      const SrcDistObject &sourceObj, const size_t numSameIDs,
      const Kokkos::DualView<const local_ordinal_type *, buffer_device_type>
          &permuteToLIDs,
      const Kokkos::DualView<const local_ordinal_type *, buffer_device_type>
          &permuteFromLIDs,
      const CombineMode CM) override;
  // clang-format off

    virtual void
    packAndPrepare
    (const SrcDistObject& sourceObj,
     const Kokkos::DualView<
       const local_ordinal_type*,
       buffer_device_type>& exportLIDs,
     Kokkos::DualView<
       impl_scalar_type*,
       buffer_device_type>& exports,
     Kokkos::DualView<
       size_t*,
       buffer_device_type> /* numPacketsPerLID */,
     size_t& constantNumPackets,
     const execution_space &space) override;

    virtual void
    packAndPrepare
    (const SrcDistObject& sourceObj,
     const Kokkos::DualView<
       const local_ordinal_type*,
       buffer_device_type>& exportLIDs,
     Kokkos::DualView<
       impl_scalar_type*,
       buffer_device_type>& exports,
     Kokkos::DualView<
       size_t*,
       buffer_device_type> /* numPacketsPerLID */,
     size_t& constantNumPackets) override;

    virtual void
    unpackAndCombine
    (const Kokkos::DualView<
       const local_ordinal_type*,
       buffer_device_type>& importLIDs,
     Kokkos::DualView<
       impl_scalar_type*,
       buffer_device_type> imports,
     Kokkos::DualView<
       size_t*,
       buffer_device_type> /* numPacketsPerLID */,
     const size_t constantNumPackets,
     const CombineMode CM,
     const execution_space &space) override;

    virtual void
    unpackAndCombine
    (const Kokkos::DualView<
       const local_ordinal_type*,
       buffer_device_type>& importLIDs,
     Kokkos::DualView<
       impl_scalar_type*,
       buffer_device_type> imports,
     Kokkos::DualView<
       size_t*,
       buffer_device_type> /* numPacketsPerLID */,
     const size_t constantNumPackets,
     const CombineMode CM) override;

  private:

    // If comm buffers can be aliased to the data view, use this
    // implementation.
    template<class NO=Node>
    typename std::enable_if<std::is_same<typename Tpetra::Details::DefaultTypes::CommBufferMemorySpace<typename NO::execution_space>::type,
                                         typename NO::device_type::memory_space>::value, bool>::type
    reallocImportsIfNeededImpl (const size_t newSize,
                                 const bool verbose,
                                 const std::string* prefix,
                                 const bool areRemoteLIDsContiguous,
                                 const CombineMode CM);

    // If comm buffers cannot be aliased to the data view, use this
    // implementation. (Just calls DistObject::reallocImportsIfNeeded.)
    template<class NO=Node>
    typename std::enable_if<!std::is_same<typename Tpetra::Details::DefaultTypes::CommBufferMemorySpace<typename NO::execution_space>::type,
                                          typename NO::device_type::memory_space>::value, bool>::type
    reallocImportsIfNeededImpl (const size_t newSize,
                                 const bool verbose,
                                 const std::string* prefix,
                                 const bool areRemoteLIDsContiguous,
                                 const CombineMode CM);
  protected:

    virtual bool
    reallocImportsIfNeeded (const size_t newSize,
                                 const bool verbose,
                                 const std::string* prefix,
                                 const bool areRemoteLIDsContiguous=false,
                                 const CombineMode CM=INSERT) override;


  public:
    bool importsAreAliased();

  protected:
    Kokkos::DualView<impl_scalar_type*, buffer_device_type> unaliased_imports_;

    //@}
  }; // class MultiVector

  template<class SC, class LO, class GO, class NT>
  Teuchos::ArrayView<const size_t>
  getMultiVectorWhichVectors (const MultiVector<SC, LO, GO, NT>& X)
  {
    return X.whichVectors_ ();
  }


  /// \brief Specialization of deep_copy for MultiVector objects with
  ///   the same template parameters.
  template <class ST, class LO, class GO, class NT>
  void
  deep_copy (MultiVector<ST, LO, GO, NT>& dst,
             const MultiVector<ST, LO, GO, NT>& src)
  {
    // NOTE (mfh 11 Sep 2014) We can't implement deep_copy with
    // shallow-copy operator=, because that would invalidate existing
    // views of dst!
    dst.assign (src);
  }

  // Implementation of the most generic version of MultiVector deep_copy.
  template <class DS, class DL, class DG, class DN,
            class SS, class SL, class SG, class SN>
  void
  deep_copy (MultiVector<DS, DL, DG, DN>& dst,
             const MultiVector<SS, SL, SG, SN>& src)
  {
    using ::Tpetra::getMultiVectorWhichVectors;

    TEUCHOS_TEST_FOR_EXCEPTION(
      dst.getGlobalLength () != src.getGlobalLength () ||
      dst.getNumVectors () != src.getNumVectors (), std::invalid_argument,
      "Tpetra::deep_copy: Global dimensions of the two Tpetra::MultiVector "
      "objects do not match.  src has dimensions [" << src.getGlobalLength ()
      << "," << src.getNumVectors () << "], and dst has dimensions ["
      << dst.getGlobalLength () << "," << dst.getNumVectors () << "].");

    // FIXME (mfh 28 Jul 2014) Don't throw; just set a local error flag.
    TEUCHOS_TEST_FOR_EXCEPTION(
      dst.getLocalLength () != src.getLocalLength (), std::invalid_argument,
      "Tpetra::deep_copy: The local row counts of the two Tpetra::MultiVector "
      "objects do not match.  src has " << src.getLocalLength () << " row(s) "
      << " and dst has " << dst.getLocalLength () << " row(s).");

    const bool srcMostUpToDateOnDevice = ! src.need_sync_device ();

    if (src.isConstantStride () && dst.isConstantStride ()) {
      if (srcMostUpToDateOnDevice) {
        Details::localDeepCopyConstStride (
                 dst.getLocalViewDevice (Access::OverwriteAll),
                 src.getLocalViewDevice (Access::ReadOnly));
      }
      else {
        Details::localDeepCopyConstStride (
                 dst.getLocalViewDevice (Access::OverwriteAll),
                 src.getLocalViewHost (Access::ReadOnly));
      }
    }
    else {
      auto dstWhichVecs = getMultiVectorWhichVectors (dst);
      auto srcWhichVecs = getMultiVectorWhichVectors (src);

      if (srcMostUpToDateOnDevice) {
        Details::localDeepCopy (dst.getLocalViewDevice (Access::OverwriteAll),
                                src.getLocalViewDevice (Access::ReadOnly),
                                dst.isConstantStride (),
                                src.isConstantStride (),
                                dstWhichVecs,
                                srcWhichVecs);
      }
      else {
        Details::localDeepCopy (dst.getLocalViewDevice (Access::OverwriteAll),
                                src.getLocalViewHost (Access::ReadOnly),
                                dst.isConstantStride (),
                                src.isConstantStride (),
                                dstWhichVecs,
                                srcWhichVecs);
      }
    }
  }
} // namespace Tpetra


namespace Teuchos {

  // Give Teuchos::TypeNameTraits<Tpetra::MultiVector<...> > a
  // human-readable definition.
  template<class SC, class LO, class GO, class NT>
  class TypeNameTraits<Tpetra::MultiVector<SC, LO, GO, NT> > {
  public:
    static std::string name () {
      return std::string ("Tpetra::MultiVector<") +
        TypeNameTraits<SC>::name () + "," +
        TypeNameTraits<LO>::name () + "," +
        TypeNameTraits<GO>::name () + "," +
        TypeNameTraits<NT>::name () + ">";
    }

    static std::string
    concreteName (const Tpetra::MultiVector<SC, LO, GO, NT>&) {
      return name ();
    }
  };
} // namespace Teuchos

#endif // TPETRA_MULTIVECTOR_DECL_HPP
