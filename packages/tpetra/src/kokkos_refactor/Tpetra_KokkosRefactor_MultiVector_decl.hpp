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

#ifndef TPETRA_KOKKOS_REFACTOR_MULTIVECTOR_DECL_HPP
#define TPETRA_KOKKOS_REFACTOR_MULTIVECTOR_DECL_HPP

#include <Kokkos_DualView.hpp>
#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>
#include <Kokkos_ArithTraits.hpp>
#include <Tpetra_KokkosRefactor_DistObject.hpp>

#ifdef KOKKOS_HAVE_CXX11
#  include <type_traits>
#endif // KOKKOS_HAVE_CXX11

namespace KokkosClassic {
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration of DefaultArithmetic
  template<class KokkosMultiVectorType>
  class DefaultArithmetic;
#endif // DOXYGEN_SHOULD_SKIP_THIS
} // namespace KokkosClassic

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration of Vector, needed to prevent circular inclusions
  template<class S, class LO, class GO, class N> class Vector;

  // forward declaration of Map
  template<class LO, class GO, class N> class Map;
#endif // DOXYGEN_SHOULD_SKIP_THIS

  /// \class MultiVector
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
  ///  can use real-valued or complex-valued types here, unlike in
  ///  Epetra, where the scalar type is always \c double.)
  /// \tparam LocalOrdinal The type of local indices.  See the
  ///   documentation of Map for requirements.
  /// \tparam GlobalOrdinal The type of global indices.  See the
  ///   documentation of Map for requirements.
  /// \tparam Node The Kokkos Node type.  See the documentation of Map
  ///   for requirements.
  ///
  /// \section Kokkos_KR_MV_prereq Prerequisites
  ///
  /// Before reading the rest of this documentation, it helps to know
  /// something about the Teuchos memory management classes, in
  /// particular Teuchos::RCP, Teuchos::ArrayRCP, and
  /// Teuchos::ArrayView.  You may also want to know about the
  /// differences between BLAS 1, 2, and 3 operations, and learn a
  /// little bit about MPI (the Message Passing Interface for
  /// distributed-memory programming).  You won't have to use MPI
  /// directly to use MultiVector, but it helps to be familiar with
  /// the general idea of distributed storage of data over a
  /// communicator.
  ///
  /// \section Kokkos_KR_MV_layout Data layout and access
  ///
  /// \subsection Kokkos_KR_MV_noncontig Multivectors which view another multivector
  ///
  /// A multivector could be a <i>view</i> of some subset of another
  /// multivector's columns and rows.  A view is like a pointer; it
  /// provides access to the original multivector's data without
  /// copying the data.  There are no public constructors for creating
  /// a view, but any instance method with "view" in the name that
  /// returns an RCP<MultiVector> serves as a view constructor.
  ///
  /// The subset of columns in a view need not be contiguous.  For
  /// example, given a multivector X with 43 columns, it is possible
  /// to have a multivector Y which is a view of columns 1, 3, and 42
  /// (zero-based indices) of X.  We call such multivectors
  /// <i>noncontiguous</i>.  They have the the property that
  /// isConstantStride() returns false.
  ///
  /// Noncontiguous multivectors lose some performance advantages.
  /// For example, local computations may be slower, since Tpetra
  /// cannot use BLAS 3 routines (e.g., matrix-matrix multiply) on a
  /// noncontiguous multivectors without copying into temporary
  /// contiguous storage.  Noncontiguous multivectors also affect the
  /// ability to access the data in certain ways, which we will
  /// explain below.
  ///
  /// \subsection Kokkos_KR_MV_views Views of a multivector's data
  ///
  /// We have unfortunately overloaded the term "view."  In the
  /// section above, we explained the idea of a "multivector which is
  /// a view of another multivector."  This section is about "views of
  /// a multivector's data."  If you want to read or write the actual
  /// values in a multivector, this is what you want.  All the
  /// instance methods which return an ArrayRCP of Scalar data, or an
  /// ArrayRCP of ArrayRCP of Scalar data, return views to the data.
  /// These data are always <i>local</i> data, meaning that the
  /// corresponding rows of the multivector are owned by the calling
  /// process.  You can't use these methods to access remote data
  /// (rows that do not belong to the calling process).
  ///
  /// Data views may be either one-dimensional (1-D) or
  /// two-dimensional (2-D).  A 1-D view presents the data as a dense
  /// matrix in column-major order, returned as a single array.  On
  /// the calling process, the matrix has getLocalLength() rows,
  /// getNumVectors() columns, and column stride getStride().  You may
  /// not get a 1-D view of a noncontiguous multivector.  If you need
  /// the data of a noncontiguous multivector in a 1-D format, you may
  /// get a copy by calling get1dCopy().  A 2-D view presents the data
  /// as an array of arrays, one array per column (i.e., vector in the
  /// multivector).  The entries in each column are stored
  /// contiguously.  You may get a 2-D view of <i>any</i> multivector,
  /// whether or not it is noncontiguous.
  ///
  /// Views are not necessarily just encapsulated pointers.  The
  /// meaning of view depends in part on the Kokkos Node type (the
  /// Node template parameter).  This matters in particular if you are
  /// running on a Graphics Processing Unit (GPU) device.  You can
  /// tell at compile time whether you are running on a GPU by looking
  /// at the Kokkos Node type.  (Currently, the only GPU Node type we
  /// provide is KokkosClassic::ThrustGPUNode.  All other types are CPU
  /// Nodes.)  If the Kokkos Node is a GPU Node type, then views
  /// always reside in host (CPU) memory, rather than device (GPU)
  /// memory.  When you ask for a view, it copies data from the device
  /// to the host.
  ///
  /// What happens next to your view depends on whether the view is
  /// const (read-only) or nonconst (read and write).  Const views
  /// disappear (their host memory is deallocated) when the
  /// corresponding reference count (of the ArrayRCP) goes to zero.
  /// (Since the data were not changed, there is no need to update the
  /// original copy on the device.)  When a nonconst view's reference
  /// count goes to zero, the view's data are copied back to device
  /// memory, thus "pushing" your changes on the host back to the
  /// device.
  ///
  /// These device-host-device copy semantics on GPUs mean that we can
  /// only promise that a view is a snapshot of the multivector's data
  /// at the time the view was created.  If you create a const view,
  /// then a nonconst view, then modify the nonconst view, the
  /// contents of the const view are undefined.  For host-only (CPU
  /// only, no GPU) Kokkos Nodes, views may be just encapsulated
  /// pointers to the data, so modifying a nonconst view will change
  /// the original data.  For GPU Nodes, modifying a nonconst view
  /// will <i>not</i> change the original data until the view's
  /// reference count goes to zero.  Furthermore, if the nonconst
  /// view's reference count never goes to zero, the nonconst view
  /// will <i>never</i> be copied back to device memory, and thus the
  /// original data will never be changed.
  ///
  /// \subsection Kokkos_KR_MV_why_views Why won't you give me a raw pointer?
  ///
  /// Tpetra was designed to allow different data representations
  /// underneath the same interface.  This lets Tpetra run correctly
  /// and efficiently on many different kinds of hardware, including
  /// single-core CPUs, multicore CPUs with Non-Uniform Memory Access
  /// (NUMA), and even "discrete compute accelerators" like Graphics
  /// Processing Units (GPUs).  These different kinds of hardware all
  /// have in common the following:
  ///
  /// - Data layout matters a lot for performance
  /// - The right layout for your data depends on the hardware
  /// - Data may be distributed over different memory spaces in
  ///   hardware, and efficient code must respect this, whether or not
  ///   the programming model presents the different memories as a
  ///   single address space
  /// - Copying between different data layouts or memory spaces is
  ///   expensive and should be avoided whenever possible
  /// - Optimal data layout may require control over initialization of
  ///   storage
  ///
  /// These conclusions have practical consequences for the
  /// MultiVector interface.  In particular, we have deliberately made
  /// it difficult for you to access data directly by raw pointer.
  /// This is because the underlying layout may not be what you
  /// expect.  In some cases, you are not even allowed to dereference
  /// the raw pointer (for example, if it resides in GPU device
  /// memory, and you are working on the host CPU).  This is why we
  /// require accessing the data through views.
  ///
  /// \subsection Kokkos_KR_MV_why_no_square_brackets Why no operator[]?
  ///
  /// The above section also explains why we do not offer a <tt>Scalar&
  /// operator[]</tt> to access each entry of a vector directly.  Direct
  /// access on GPUs would require implicitly creating an internal
  /// host copy of the device data.  This would consume host memory,
  /// and it would not be clear when to write the resulting host data
  /// back to device memory.  The resulting operator would violate
  /// users' performance expectations, since it would be much slower
  /// than raw array access.  We have preferred in our design to
  /// expose what is expensive, by exposing data views and letting
  /// users control when to copy data between host and device.
  ///
  /// \subsection Kokkos_KR_MV_device_kernels How do I access the data directly?
  ///
  /// "Directly" here means without views, using a device kernel if
  /// the data reside on the GPU.
  ///
  /// There are two different options for direct access to the
  /// multivector's data.  One is to use the optional RTI (Reduction /
  /// Transformation Interface) subpackage of Tpetra.  You may enable
  /// this at Trilinos configure time by setting the CMake Boolean
  /// option \c Tpetra_ENABLE_RTI to \c ON.  Be aware that building
  /// and using RTI requires that your C++ compiler support the
  /// language features in the new C++11 standard.  RTI allows you to
  /// implement arbitrary element-wise operations over a vector,
  /// followed by arbitrary reductions over the elements of that
  /// vector.  We recommend RTI for most users.
  ///
  /// Another option is to access the local data through its Kokkos
  /// container data structure, KokkosClassic::MultiVector, and then use the
  /// \ref kokkos_node_api "Kokkos Node API" to implement arbitrary
  /// operations on the data.  We do not recommend this approach for
  /// most users.  In particular, the local data structures are likely
  /// to change over the next few releases.  If you find yourself
  /// wanting to try this option, please contact the Tpetra developers
  /// for recommendations.  We will be happy to work with you.
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
  /// these methods have the same blocking semantics as MPI_Allreduce.
  ///
  /// \warning Some computational methods, such as inner products and
  ///   norms, may return incorrect results if the MultiVector's Map
  ///   is overlapping (not one-to-one) but not locally replicated.
  ///   That is, if some but not all rows are shared by more than one
  ///   process in the communicator, then inner products and norms may
  ///   be wrong.  This behavior may change in future releases.
  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class DeviceType>
  class MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > :
    public DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >
  {
  private:
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> Node;
    typedef DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node> base_type;

  public:
    //! @name Typedefs to facilitate template metaprogramming.
    //@{

    //! The type of entries in the vector(s).
    typedef typename Kokkos::Details::ArithTraits<Scalar>::val_type scalar_type;
    //! The type of local indices.
    typedef LocalOrdinal local_ordinal_type;
    //! The type of global indices.
    typedef GlobalOrdinal global_ordinal_type;
    //! The Kokkos Node type.
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> node_type;

    /// \brief Type of an inner ("dot") product result.
    ///
    /// This is usually the same as <tt>scalar_type</tt>, but may
    /// differ if <tt>scalar_type</tt> is e.g., an uncertainty
    /// quantification type from the Stokhos package.
    typedef typename Kokkos::Details::InnerProductSpaceTraits<scalar_type>::dot_type dot_type;

    /// \brief Type of a norm result.
    ///
    /// This is usually the same as the type of the magnitude
    /// (absolute value) of <tt>scalar_type</tt>, but may differ if
    /// <tt>scalar_type</tt> is e.g., an uncertainty quantification
    /// type from the Stokhos package.
    typedef typename Kokkos::Details::ArithTraits<scalar_type>::mag_type mag_type;

    //! Type of the (new) Kokkos Device which implements parallel operations.
    typedef DeviceType device_type;

    //! Kokkos::DualView specialization used by this class.
    typedef Kokkos::DualView<scalar_type**, Kokkos::LayoutLeft, device_type> dual_view_type;

    //@}
    //! @name Constructors and destructor
    //@{

    //! Default constructor: makes a MultiVector with no rows or columns.
    MultiVector ();

    /// \brief Basic constuctor.
    ///
    /// \param map [in] Map describing the distribution of rows.
    /// \param NumVectors [in] Number of vectors (columns).
    /// \param zeroOut [in] Whether to initialize all the entries of
    ///   the MultiVector to zero.
    ///
    /// \note The Kokkos refactor version of MultiVector reserves the
    ///   right to initialize all entries of the MultiVector to zero,
    ///   regardless of the value of \c zeroOut.
    MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
                 size_t NumVectors,
                 bool zeroOut=true);

    //! Copy constructor (shallow copy!).
    MultiVector (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &source);

    /// \brief Copy constructor, with option to do deep or shallow copy.
    ///
    /// The Kokkos refactor version of Tpetra, unlike the "classic"
    /// version, always has view semantics.  Thus, copyOrView =
    /// Teuchos::View has no effect, and copyOrView = Teuchos::Copy
    /// does not mark this MultiVector as having copy semantics.
    /// However, copyOrView = Teuchos::Copy will make the resulting
    /// MultiVector a deep copy of the input MultiVector.
    MultiVector (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& source,
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
    MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
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
    MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
                 const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> >&ArrayOfPtrs,
                 const size_t NumVectors);

    /// \brief Expert mode constructor, that takes a Kokkos::DualView
    ///   of the MultiVector's data, and returns a MultiVector that
    ///   views those data.
    ///
    /// \warning This constructor is only for expert users.  We make
    ///   no promises about backwards compatibility for this interface.
    ///   It may change or go away at any time.
    ///
    /// \param map [in] Map describing the distribution of rows.
    /// \param view [in] Device view to the data (shallow copy).
    MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
                 const dual_view_type& view);

    /// \brief Expert mode constructor, that takes a Kokkos::DualView
    ///   of the MultiVector's data and the "original"
    ///   Kokkos::DualView of the data, and returns a MultiVector that
    ///   views those data.
    ///
    /// \warning This constructor is only for expert users.  We make
    ///   no promises about backwards compatibility for this interface.
    ///   It may change or go away at any time.
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
    MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
                 const dual_view_type& view,
                 const dual_view_type& origView);

    /// \brief Expert mode constructor for noncontiguous views.
    ///
    /// This constructor that takes a Kokkos::DualView of the
    /// MultiVector's data, and a list of the columns to view, and
    /// returns a MultiVector that views those data.  The resulting
    /// MultiVector does <i>not</i> have constant stride, that is,
    /// isConstantStride() returns false.
    ///
    /// \warning This constructor is only for expert users.  We make
    ///   no promises about backwards compatibility for this interface.
    ///   It may change or go away at any time.
    ///
    /// \param map [in] Map describing the distribution of rows.
    /// \param view [in] Device view to the data (shallow copy).
    /// \param whichVectors [in] Which columns (vectors) to view.
    MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
                 const dual_view_type& view,
                 const Teuchos::ArrayView<const size_t>& whichVectors);

    /// \brief Expert mode constructor for noncontiguous views, with
    ///   original view.
    ///
    /// This constructor that takes a Kokkos::DualView of the
    /// MultiVector's data, a view of the original data, and a list of
    /// the columns to view, and returns a MultiVector that views
    /// those data.  The resulting MultiVector does <i>not</i> have
    /// constant stride, that is, isConstantStride() returns false.
    ///
    /// \warning This constructor is only for expert users.  We make
    ///   no promises about backwards compatibility for this interface.
    ///   It may change or go away at any time.
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
    MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
                 const dual_view_type& view,
                 const dual_view_type& origView,
                 const Teuchos::ArrayView<const size_t>& whichVectors);

    //! Return a deep copy of <tt>*this</tt>, for a different Kokkos Node type.
    template <class Node2>
    Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node2> >
    clone (const Teuchos::RCP<Node2> &node2) const;

    //! Destructor (virtual for memory safety of derived classes).
    virtual ~MultiVector();

    //@}
    //! @name Post-construction modification routines
    //@{

    /// \brief Replace value, using global (row) index.
    ///
    /// Replace the current value at row \c globalRow (a global index)
    /// and column \c col with the given value.  The column index is
    /// zero based.
    ///
    /// \pre \c globalRow must be a valid global element on this
    ///   process, according to the row Map.
    void
    replaceGlobalValue (GlobalOrdinal globalRow,
                        size_t col,
                        const scalar_type& value);

    /// \brief Like the above replaceGlobalValue, but only enabled if
    ///   T differs from scalar_type.
    ///
    /// This method only exists if the template parameter T and
    /// scalar_type differ.  If C++11 is enabled, we further require
    /// that it be possible to convert T to scalar_type.
    ///
    /// This method is mainly useful for backwards compatibility, when
    /// Scalar differs from scalar_type.
    template<typename T>
#ifdef KOKKOS_HAVE_CXX11
    typename std::enable_if<! std::is_same<T, scalar_type>::value && std::is_convertible<T, scalar_type>::value, void>::type
#else
    typename Kokkos::Impl::enable_if<! Kokkos::Impl::is_same<T, scalar_type>::value, void>::type
#endif // KOKKOS_HAVE_CXX11
    replaceGlobalValue (GlobalOrdinal globalRow,
                        size_t col,
                        const T& value)
    {
      replaceGlobalValue (globalRow, col, static_cast<scalar_type> (value));
    }

    /// \brief Add value to existing value, using global (row) index.
    ///
    /// Add the given value to the existing value at row \c globalRow
    /// (a global index) and column \c col.  The column index is zero
    /// based.
    ///
    /// \pre \c globalRow must be a valid global element on this
    ///   process, according to the row Map.
    void
    sumIntoGlobalValue (GlobalOrdinal globalRow,
                        size_t col,
                        const scalar_type& value);

    /// \brief Like the above sumIntoGlobalValue, but only enabled if
    ///   T differs from scalar_type.
    ///
    /// This method only exists if the template parameter T and
    /// scalar_type differ.  If C++11 is enabled, we further require
    /// that it be possible to convert T to scalar_type.
    ///
    /// This method is mainly useful for backwards compatibility, when
    /// Scalar differs from scalar_type.
    template<typename T>
#ifdef KOKKOS_HAVE_CXX11
    typename std::enable_if<! std::is_same<T, scalar_type>::value && std::is_convertible<T, scalar_type>::value, void>::type
#else
    typename Kokkos::Impl::enable_if<! Kokkos::Impl::is_same<T, scalar_type>::value, void>::type
#endif // KOKKOS_HAVE_CXX11
    sumIntoGlobalValue (GlobalOrdinal globalRow,
                        size_t col,
                        const T& value)
    {
      sumIntoGlobalValue (globalRow, col, static_cast<scalar_type> (value));
    }

    /// \brief Replace value, using local (row) index.
    ///
    /// Replace the current value at row \c localRow (a local index)
    /// and column \c col with the given value.  The column index is
    /// zero based.
    ///
    /// \pre \c localRow must be a valid local element on this
    ///   process, according to the row Map.
    void
    replaceLocalValue (LocalOrdinal localRow,
                       size_t col,
                       const scalar_type& value);

    /// \brief Like the above replaceLocalValue, but only enabled if
    ///   T differs from scalar_type.
    ///
    /// This method only exists if the template parameter T and
    /// scalar_type differ.  If C++11 is enabled, we further require
    /// that it be possible to convert T to scalar_type.
    ///
    /// This method is mainly useful for backwards compatibility, when
    /// Scalar differs from scalar_type.
    template<typename T>
#ifdef KOKKOS_HAVE_CXX11
    typename std::enable_if<! std::is_same<T, scalar_type>::value && std::is_convertible<T, scalar_type>::value, void>::type
#else
    typename Kokkos::Impl::enable_if<! Kokkos::Impl::is_same<T, scalar_type>::value, void>::type
#endif // KOKKOS_HAVE_CXX11
    replaceLocalValue (LocalOrdinal localRow,
                       size_t col,
                       const T& value)
    {
      replaceLocalValue (localRow, col, static_cast<scalar_type> (value));
    }

    /// \brief Add value to existing value, using local (row) index.
    ///
    /// Add the given value to the existing value at row \c localRow
    /// (a local index) and column \c col.  The column index is zero
    /// based.
    ///
    /// \pre \c localRow must be a valid local element on this process,
    ///   according to the row Map.
    void
    sumIntoLocalValue (LocalOrdinal localRow,
                       size_t col,
                       const scalar_type& value);

    /// \brief Like the above sumIntoLocalValue, but only enabled if
    ///   T differs from scalar_type.
    ///
    /// This method only exists if the template parameter T and
    /// scalar_type differ.  If C++11 is enabled, we further require
    /// that it be possible to convert T to scalar_type.
    ///
    /// This method is mainly useful for backwards compatibility, when
    /// Scalar differs from scalar_type.
    template<typename T>
#ifdef KOKKOS_HAVE_CXX11
    typename std::enable_if<! std::is_same<T, scalar_type>::value && std::is_convertible<T, scalar_type>::value, void>::type
#else
    typename Kokkos::Impl::enable_if<! Kokkos::Impl::is_same<T, scalar_type>::value, void>::type
#endif // KOKKOS_HAVE_CXX11
    sumIntoLocalValue (LocalOrdinal localRow,
                       size_t col,
                       const T& value)
    {
      sumIntoLocalValue (localRow, col, static_cast<scalar_type> (value));
    }

    //! Set all values in the multivector with the given value.
    void putScalar (const Scalar &value);

    /// \brief Set all values in the multivector to pseudorandom numbers.
    ///
    /// \note The implementation of this method may depend on the
    ///   Kokkos Node type and on what third-party libraries you have
    ///   available.  Do not expect repeatable results.
    ///
    /// \note Behavior of this method may or may not depend on
    ///   external use of the C library routines \c srand() and \c
    ///   rand().
    ///
    /// \warning This method does <i>not</i> promise to use a
    ///   distributed-memory parallel pseudorandom number generator.
    ///   Corresponding values on different processes might be
    ///   correlated.  It also does not promise to use a high-quality
    ///   pseudorandom number generator within each process.
    void randomize();

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
    void replaceMap (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map);

    /// \brief Sum values of a locally replicated multivector across all processes.
    ///
    /// \warning This method may only be called for locally replicated
    ///   MultiVectors.
    ///
    /// \pre isDistributed() == false
    void reduce();

    /// \brief Assign the contents of \c source to this multivector (deep copy).
    ///
    /// \pre The two multivectors must have the same communicator.
    /// \pre The input multivector's Map must be compatible with this
    ///      multivector's Map.  That is, \code
    ///      this->getMap ()->isCompatible (source.getMap ());
    ///      \endcode
    /// \pre The two multivectors must have the same number of columns.
    ///
    /// \note This method must always be called as a collective
    ///   operation on all processes over which the multivector is
    ///   distributed.  This is because the method reserves the right
    ///   to check for compatibility of the two Maps, at least in
    ///   debug mode.
    MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&
    operator= (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& source);

    //@}

    //! @name Data Copy and View get methods
    /** These methods are used to get the data underlying the MultiVector. They return data in one of three forms:
      - a MultiVector with a subset of the columns of the target MultiVector
      - a raw C pointer or array of raw C pointers
      - one of the Teuchos memory management classes
      Not all of these methods are valid for a particular MultiVector. For instance, calling a method that accesses a
      view of the data in a 1-D format (i.e., get1dView) requires that the target MultiVector has constant stride.
     */
    //@{

    //! Return a MultiVector with copies of selected columns.
    Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
    subCopy (const Teuchos::Range1D &colRng) const;

    //! Return a MultiVector with copies of selected columns.
    Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
    subCopy (const Teuchos::ArrayView<const size_t> &cols) const;

    //! Return a const MultiVector with const views of selected columns.
    Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
    subView (const Teuchos::Range1D &colRng) const;

    //! Return a const MultiVector with const views of selected columns.
    Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
    subView (const Teuchos::ArrayView<const size_t> &cols) const;

    //! Return a MultiVector with views of selected columns.
    Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
    subViewNonConst (const Teuchos::Range1D &colRng);

    //! Return a MultiVector with views of selected columns.
    Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
    subViewNonConst (const Teuchos::ArrayView<const size_t> &cols);

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
    Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
    offsetView (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& subMap,
                size_t offset) const;

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
    Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
    offsetViewNonConst (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &subMap,
                        size_t offset);

    //! Return a Vector which is a const view of column j.
    Teuchos::RCP<const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
    getVector (size_t j) const;

    //! Return a Vector which is a nonconst view of column j.
    Teuchos::RCP<Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
    getVectorNonConst (size_t j);

    //! Const view of the local values in a particular vector of this multivector.
    Teuchos::ArrayRCP<const scalar_type> getData (size_t j) const;

    //! View of the local values in a particular vector of this multivector.
    Teuchos::ArrayRCP<scalar_type> getDataNonConst (size_t j);

    /// \brief Fill the given array with a copy of this multivector's local values.
    ///
    /// \param A [out] View of the array to fill.  We consider A as a
    ///   matrix with column-major storage.
    ///
    /// \param LDA [in] Leading dimension of the matrix A.
    void
    get1dCopy (const Teuchos::ArrayView<scalar_type>& A,
               const size_t LDA) const;

    template<typename T>
#ifdef KOKKOS_HAVE_CXX11
    typename std::enable_if<! std::is_same<T, scalar_type>::value &&
                            sizeof (T) == sizeof (scalar_type), void>::type
#else
    typename Kokkos::Impl::enable_if<! Kokkos::Impl::is_same<T, scalar_type>::value, void>::type
#endif // KOKKOS_HAVE_CXX11
    get1dCopy (const Teuchos::ArrayView<T>& A, const size_t LDA) const
    {
      Teuchos::ArrayView<scalar_type> A_r =
        Teuchos::av_reinterpret_cast<scalar_type> (A);
      this->get1dCopy (A_r, LDA);
    }

    /// \brief Fill the given array with a copy of this multivector's local values.
    ///
    /// \param ArrayOfPtrs [out] Array of arrays, one for each column
    ///   of the multivector.  On output, we fill ArrayOfPtrs[j] with
    ///   the data for column j of this multivector.
    void get2dCopy (Teuchos::ArrayView<const Teuchos::ArrayView<scalar_type> > ArrayOfPtrs) const;

    template<typename T>
#ifdef KOKKOS_HAVE_CXX11
    typename std::enable_if<! std::is_same<T, scalar_type>::value &&
                            sizeof (T) == sizeof (scalar_type), void>::type
#else
    typename Kokkos::Impl::enable_if<! Kokkos::Impl::is_same<T, scalar_type>::value, void>::type
#endif // KOKKOS_HAVE_CXX11
    get2dCopy (const Teuchos::ArrayView<const Teuchos::ArrayView<T> >& A,
               const size_t LDA) const
    {
      const size_t numVecs = this->getNumVectors ();
      const size_t numEnts = A.size ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        numVecs != numEnts, std::invalid_argument, "Tpetra::MultiVector::"
        "get2dCopy: this->getNumVectors() = " << numVecs << " != A.size() = "
        << numEnts << ".");
      for (size_t j = 0; j < numVecs; ++j) {
        Teuchos::ArrayView<T> A_j = Teuchos::av_reinterpret_cast<T> (A[j]);
        this->getVector (j)->get1dCopy (A_j);
      }
    }

    /// \brief Const persisting (1-D) view of this multivector's local values.
    ///
    /// This method assumes that the columns of the multivector are
    /// stored contiguously.  If not, this method throws
    /// std::runtime_error.
    Teuchos::ArrayRCP<const scalar_type> get1dView () const;

    //! Return const persisting pointers to values.
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<const scalar_type> > get2dView () const;

    /// \brief Nonconst persisting (1-D) view of this multivector's local values.
    ///
    /// This method assumes that the columns of the multivector are
    /// stored contiguously.  If not, this method throws
    /// std::runtime_error.
    Teuchos::ArrayRCP<scalar_type> get1dViewNonConst ();

    //! Return non-const persisting pointers to values.
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<scalar_type> > get2dViewNonConst ();

    /// \brief A view of the underlying KokkosClassic::MultiVector object.
    ///
    /// \brief This method is for expert users only.
    ///   It may change or be removed at any time.
    KokkosClassic::MultiVector<scalar_type, node_type> getLocalMV () const;

    /// \brief A nonconst reference to a view of the underlying
    ///   KokkosClassic::MultiVector object.
    ///
    /// \brief This method is for expert users only.
    ///   It may change or be removed at any time.
    ///
    /// \warning This method is DEPRECATED.  It may disappear at any
    ///   time.  Please call getLocalMV() instead.  There was never
    ///   actually a need for a getLocalMVNonConst() method, as far as
    ///   I can tell.
    TEUCHOS_DEPRECATED KokkosClassic::MultiVector<scalar_type, node_type>
    getLocalMVNonConst ();

    /// \brief Get the Kokkos::DualView which implements local storage.
    ///
    /// \warning This method is ONLY for expert developers.  Its
    ///   interface may change or it may disappear at any time.
    dual_view_type getDualView () const;

    /// \brief Update data on device or host only if data in the other
    ///   space has been marked as modified.
    ///
    /// If \c TargetDeviceType is the same as this MultiVector's
    /// device type, then copy data from host to device.  Otherwise,
    /// copy data from device to host.  In either case, only copy if
    /// the source of the copy has been modified.
    ///
    /// This is a one-way synchronization only.  If the target of the
    /// copy has been modified, this operation will discard those
    /// modifications.  It will also reset both device and host modified
    /// flags.
    ///
    /// \note This method doesn't know on its own whether you modified
    ///   the data in either memory space.  You must manually mark the
    ///   MultiVector as modified in the space in which you modified
    ///   it, by calling the modify() method with the appropriate
    ///   template parameter.
    template<class TargetDeviceType>
    void sync () {
      getDualView ().template sync<TargetDeviceType> ();
    }

    /// \brief Mark data as modified on the given device \c TargetDeviceType.
    ///
    /// If \c TargetDeviceType is the same as this MultiVector's
    /// device type, then mark the device's data as modified.
    /// Otherwise, mark the host's data as modified.
    template<class TargetDeviceType>
    void modify () {
      getDualView ().template modify<TargetDeviceType> ();
    }

    /// \brief Return a view of the local data on a specific device.
    /// \tparam TargetDeviceType The Kokkos Device type whose data to return.
    ///
    /// Please don't be afraid of the if_c expression in the return
    /// value's type.  That just tells the method what the return type
    /// should be: dual_view_type::t_dev if the \c TargetDeviceType
    /// template parameter matches this Tpetra object's device type,
    /// else dual_view_type::t_host.
    ///
    /// For example, suppose you create a Tpetra::MultiVector for the
    /// Kokkos::Cuda device, like this:
    /// \code
    /// typedef Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::Cuda> > node_type;
    /// typedef Tpetra::Map<int, int, node_type> map_type;
    /// typedef Tpetra::MultiVector<float, int, int, node_type> mv_type;
    ///
    /// RCP<const map_type> map = ...;
    /// mv_type DV (map, 3);
    /// \endcode
    /// If you want to get the CUDA device Kokkos::View, do this:
    /// \code
    /// typedef typename mv_type::dual_view_type dual_view_type;
    /// typedef typename dual_view_type::t_dev device_view_type;
    /// device_view_type cudaView = DV.getLocalView<Kokkos::Cuda> ();
    /// \endcode
    /// and if you want to get the host mirror of that View, do this:
    /// \code
    /// typedef typename dual_view_type::host_mirror_space host_device_type;
    /// typedef typename dual_view_type::t_host host_view_type;
    /// host_view_type hostView = DV.getLocalView<host_device_type> ();
    /// \endcode
    template<class TargetDeviceType>
    typename Kokkos::Impl::if_c<
      Kokkos::Impl::is_same<
        typename device_type::memory_space,
        typename TargetDeviceType::memory_space>::value,
      typename dual_view_type::t_dev,
      typename dual_view_type::t_host>::type
    getLocalView () const {
      return getDualView ().template view<TargetDeviceType> ();
    }

    //@}
    //! @name Mathematical methods
    //@{

    /// \brief Compute the dot product of each corresponding pair of
    ///   vectors (columns) in A and B.
    ///
    /// The "dot product" is the standard Euclidean inner product.  If
    /// the type of entries of the vectors (scalar_type) is complex,
    /// then A is transposed, not <tt>*this</tt>.  For example, if x
    /// and y each have one column, then <tt>x.dot (y, dots)</tt>
    /// computes \f$y^* x = \bar{y}^T x = \sum_i \bar{y}_i \cdot x_i\f$.
    ///
    /// \pre <tt>*this</tt> and A have the same number of columns (vectors).
    /// \pre \c dots has at least as many entries as the number of columns in A.
    ///
    /// \post <tt>dots[j] == (this->getVector[j])->dot (* (A.getVector[j]))</tt>
    void
    dot (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A,
         const Teuchos::ArrayView<dot_type>& dots) const;

    /// \brief Compute the dot product of each corresponding pair of
    ///   vectors (columns) in A and B.
    /// \tparam T The output type of the dot products.
    ///
    /// This method only exists if dot_type and T are different types.
    /// For example, if scalar_type and dot_type differ, then this
    /// method ensures backwards compatibility with the previous
    /// interface (that returned dot products as scalar_type rather
    /// than as dot_type).  The complicated \c enable_if expression
    /// just ensures that the method only exists if dot_type and T are
    /// different types; the method still returns \c void, as above.
    template <typename T>
    typename Kokkos::Impl::enable_if< !(Kokkos::Impl::is_same<dot_type, T>::value), void >::type
    dot (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A,
         const Teuchos::ArrayView<T> &dots) const
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
    /// the type of entries of the vectors (scalar_type) is complex,
    /// then A is transposed, not <tt>*this</tt>.  For example, if x
    /// and y each have one column, then <tt>x.dot (y, dots)</tt>
    /// computes \f$y^* x = \bar{y}^T x = \sum_i \bar{y}_i \cdot x_i\f$.
    ///
    /// \param A [in] MultiVector with which to dot \c *this.
    /// \param dots [out] Device View with getNumVectors() entries.
    ///
    /// \pre <tt>this->getNumVectors () == A.getNumVectors ()</tt>
    /// \pre <tt>dots.dimension_0 () == A.getNumVectors ()</tt>
    ///
    /// \post <tt>dots(j) == (this->getVector[j])->dot (* (A.getVector[j]))</tt>
    void
    dot (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A,
         const Kokkos::View<dot_type*, device_type>& dots) const;

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
    typename Kokkos::Impl::enable_if< !(Kokkos::Impl::is_same<dot_type, T>::value), void >::type
    dot (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A,
         const Kokkos::View<T*, device_type>& dots) const
    {
      const size_t numDots = dots.dimension_0 ();
      Kokkos::View<dot_type*, device_type> dts ("MV::dot tmp", numDots);
      // Call overload that takes a Kokkos::View<dot_type*, device_type>.
      this->dot (A, dts);
      // FIXME (mfh 14 Jul 2014) Does this actually work if dot_type
      // and T differ?  We would need a test for this, but only the
      // Sacado and Stokhos packages are likely to care about this use
      // case.  It could also come up for Kokkos::complex ->
      // std::complex conversions, but those two implementations
      // should generally be bitwise compatible.
      Kokkos::deep_copy (dots, dts);
    }

    //! Put element-wise absolute values of input Multi-vector in target: A = abs(this)
    void abs (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A);

    //! Put element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
    void reciprocal (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A);

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
    void scale (Teuchos::ArrayView<const Scalar> alpha);

    /// \brief Scale each column in place: <tt>this[j] = alpha[j]*this[j]</tt>.
    ///
    /// Replace each column j of this MultiVector with
    /// <tt>alpha[j]</tt> times the current column j of this
    /// MultiVector.  This method will always multiply, even if all
    /// the entries of alpha are zero.  That means, for example, that
    /// if \c *this contains NaN entries before calling this method,
    /// the NaN entries will remain after this method finishes.
    void scale (const Kokkos::View<const scalar_type*, device_type> alpha);

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
           const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A);

    /// \brief Update: <tt>this = beta*this + alpha*A</tt>.
    ///
    /// Update this MultiVector with scaled values of A.  If beta is
    /// zero, overwrite \c *this unconditionally, even if it contains
    /// NaN entries.  It is legal for the input A to alias this
    /// MultiVector.
    void
    update (const Scalar& alpha,
            const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A,
            const Scalar& beta);

    /// \brief Update: <tt>this = gamma*this + alpha*A + beta*B</tt>.
    ///
    /// Update this MultiVector with scaled values of A and B.  If
    /// gamma is zero, overwrite \c *this unconditionally, even if it
    /// contains NaN entries.  It is legal for the inputs A or B to
    /// alias this MultiVector.
    void
    update (const Scalar& alpha,
            const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A,
            const Scalar& beta,
            const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& B,
            const Scalar& gamma);

    /// \brief Compute the one-norm of each vector (column), storing
    ///   the result in a device view.
    ///
    /// The one-norm of a vector is the sum of squares of the
    /// magnitudes of the vector's entries.  On exit, norms(j) is the
    /// one-norm of column j of this MultiVector.
    ///
    /// \param norms [out] Device View with getNumVectors() entries.
    ///
    /// \pre <tt>norms.dimension_0 () == this->getNumVectors ()</tt>
    /// \post <tt>norms(j) == (this->getVector[j])->norm1 (* (A.getVector[j]))</tt>
    void norm1 (const Kokkos::View<mag_type*, device_type>& norms) const;

    /// \brief Compute the one-norm of each vector (column), storing
    ///   the result in a device view.
    /// \tparam T The output type of the dot products.
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
    typename Kokkos::Impl::enable_if< !(Kokkos::Impl::is_same<mag_type, T>::value), void >::type
    norm1 (const Kokkos::View<T*, device_type>& norms) const
    {
      const size_t numNorms = norms.dimension_0 ();
      Kokkos::View<mag_type*, device_type> tmpNorms ("MV::norm1 tmp", numNorms);
      // Call overload that takes a Kokkos::View<mag_type*, device_type>.
      this->norm1 (tmpNorms);
      // FIXME (mfh 15 Jul 2014) Does this actually work if mag_type
      // and T differ?  We would need a test for this, but only the
      // Sacado and Stokhos packages are likely to care about this use
      // case.  It could also come up with Kokkos::complex ->
      // std::complex conversion.
      Kokkos::deep_copy (norms, tmpNorms);
    }


    /// \brief Compute the one-norm of each vector (column).
    ///
    /// The one-norm of a vector is the sum of squares of the
    /// magnitudes of the vector's entries.  On exit, norms[j] is the
    /// one-norm of column j of this MultiVector.
    void norm1 (const Teuchos::ArrayView<mag_type>& norms) const;

    /// \brief Compute the one-norm of each vector (column).
    /// \tparam T The output type of the norms.
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
    typename Kokkos::Impl::enable_if< !(Kokkos::Impl::is_same<mag_type,T>::value), void >::type
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
    ///   the result in a device view.
    ///
    /// The two-norm of a vector is the standard Euclidean norm, the
    /// square root of the sum of squares of the magnitudes of the
    /// vector's entries.  On exit, norms(k) is the two-norm of column
    /// k of this MultiVector.
    ///
    /// \param norms [out] Device View with getNumVectors() entries.
    ///
    /// \pre <tt>norms.dimension_0 () == this->getNumVectors ()</tt>
    /// \post <tt>norms(j) == (this->getVector[j])->dot (* (A.getVector[j]))</tt>
    void norm2 (const Kokkos::View<mag_type*, device_type>& norms) const;

    /// \brief Compute the two-norm of each vector (column), storing
    ///   the result in a device view.
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
    typename Kokkos::Impl::enable_if< !(Kokkos::Impl::is_same<mag_type, T>::value), void >::type
    norm2 (const Kokkos::View<T*, device_type>& norms) const
    {
      const size_t numNorms = norms.dimension_0 ();
      Kokkos::View<mag_type*, device_type> theNorms ("MV::norm2 tmp", numNorms);
      // Call overload that takes a Kokkos::View<mag_type*, device_type>.
      this->norm2 (theNorms);
      // FIXME (mfh 14 Jul 2014) Does this actually work if mag_type
      // and T differ?  We would need a test for this, but only the
      // Sacado and Stokhos packages are likely to care about this use
      // case.  This could also come up with Kokkos::complex ->
      // std::complex conversion.
      Kokkos::deep_copy (norms, theNorms);
    }

    /// \brief Compute the two-norm of each vector (column).
    ///
    /// The two-norm of a vector is the standard Euclidean norm, the
    /// square root of the sum of squares of the magnitudes of the
    /// vector's entries.  On exit, norms[k] is the two-norm of column
    /// k of this MultiVector.
    void norm2 (const Teuchos::ArrayView<mag_type>& norms) const;

    /// \brief Compute the two-norm of each vector (column).
    /// \tparam T The output type of the norms.
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
    typename Kokkos::Impl::enable_if< !(Kokkos::Impl::is_same<mag_type,T>::value), void >::type
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
    ///   storing the result in a device view.
    ///
    /// The infinity-norm of a vector is the maximum of the magnitudes
    /// of the vector's entries.  On exit, norms(j) is the
    /// infinity-norm of column j of this MultiVector.
    void normInf (const Kokkos::View<mag_type*, device_type>& norms) const;

    /// \brief Compute the two-norm of each vector (column), storing
    ///   the result in a device view.
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
    typename Kokkos::Impl::enable_if< !(Kokkos::Impl::is_same<mag_type, T>::value), void >::type
    normInf (const Kokkos::View<T*, device_type>& norms) const
    {
      const size_t numNorms = norms.dimension_0 ();
      Kokkos::View<mag_type*, device_type> theNorms ("MV::normInf tmp", numNorms);
      // Call overload that takes a Kokkos::View<mag_type*, device_type>.
      this->normInf (theNorms);
      // FIXME (mfh 15 Jul 2014) Does this actually work if mag_type
      // and T differ?  We would need a test for this, but only the
      // Sacado and Stokhos packages are likely to care about this use
      // case.  This could also come up with Kokkos::complex ->
      // std::complex conversion.
      Kokkos::deep_copy (norms, theNorms);
    }

    /// \brief Compute the infinity-norm of each vector (column).
    ///
    /// The infinity-norm of a vector is the maximum of the magnitudes
    /// of the vector's entries.  On exit, norms[j] is the
    /// infinity-norm of column j of this MultiVector.
    void normInf (const Teuchos::ArrayView<mag_type>& norms) const;

    /// \brief Compute the infinity-norm of each vector (column).
    /// \tparam T The output type of the norms.
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
    typename Kokkos::Impl::enable_if< !(Kokkos::Impl::is_same<mag_type,T>::value), void >::type
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

    //! Compute Weighted 2-norm (RMS Norm) of each vector in multi-vector.
    //! The outcome of this routine is undefined for non-floating point scalar types (e.g., int).
    void
    normWeighted (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& weights,
                  const Teuchos::ArrayView<mag_type>& norms) const;

    /// \brief Compute the weighted 2-norm (RMS Norm) of each column.
    ///
    /// The outcome of this routine is undefined for non-floating
    /// point scalar types (e.g., int).
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
    typename Kokkos::Impl::enable_if< !(Kokkos::Impl::is_same<mag_type,T>::value), void >::type
    normWeighted (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& weights,
                  const Teuchos::ArrayView<T>& norms) const
    {
      typedef typename Teuchos::ArrayView<T>::size_type size_type;
      const size_type sz = norms.size ();
      Teuchos::Array<mag_type> theNorms (sz);
      this->normWeighted (weights, theNorms);
      for (size_type i = 0; i < sz; ++i) {
        // If T and mag_type differ, this does an implicit conversion.
        norms[i] = theNorms[i];
      }
    }

    /// \brief Compute mean (average) value of each column.
    ///
    /// The outcome of this routine is undefined for non-floating
    /// point scalar types (e.g., int).
    void meanValue (const Teuchos::ArrayView<scalar_type>& means) const;

    template <typename T>
    typename Kokkos::Impl::enable_if<! Kokkos::Impl::is_same<scalar_type, T>::value, void>::type
    meanValue (const Teuchos::ArrayView<T>& means) const
    {
      typedef typename Teuchos::Array<T>::size_type size_type;
      const size_type numMeans = means.size ();

      Teuchos::Array<scalar_type> theMeans (numMeans);
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
              const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A,
              const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& B,
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
                         const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A,
                         const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& B,
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

    //@}

    //! @name Overridden from Teuchos::Describable
    //@{

    //! A simple one-line description of this object.
    virtual std::string description() const;

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
              Teuchos::Describable::verbLevel_default) const;
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
    removeEmptyProcessesInPlace (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& newMap);

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
    /// \note The Kokkos refactor version of MultiVector <i>only</i>
    ///   implements view semantics.  This is not currently true for
    ///   the "classic" version of MultiVector, though that will
    ///   change in the near future.
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
    assign (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, node_type>& src);

  protected:
    template <class S, class LO, class GO, class D>
    friend MultiVector<S,LO,GO,Kokkos::Compat::KokkosDeviceWrapperNode<D> >
    createCopy (const MultiVector<S,LO,GO,Kokkos::Compat::KokkosDeviceWrapperNode<D> >& src);

    template <class DS, class DL, class DG, class DD, class SS, class SL, class SG, class SD>
    friend void
    deep_copy (MultiVector<DS,DL,DG,Kokkos::Compat::KokkosDeviceWrapperNode<DD> >& dst,
               const MultiVector<SS,SL,SG,Kokkos::Compat::KokkosDeviceWrapperNode<SD> >& src);

    typedef KokkosClassic::MultiVector<scalar_type, node_type> KMV;
    typedef KokkosClassic::DefaultArithmetic<KMV> MVT;

    /// \brief The KokkosClassic::MultiVector containing the MultiVector's data.
    ///
    /// This is actually just a view of the Kokkos::DualView \c view_.
    /// The latter owns the MultiVector's data.
    //KMV lclMV_;

    /// \brief The Kokkos::DualView containing the MultiVector's data.
    ///
    /// This has to be declared \c mutable, so that get1dView() can
    /// retain its current \c const marking, even though it has always
    /// implied a device->host synchronization.  Lesson to the reader:
    /// Use \c const sparingly!
    mutable dual_view_type view_;

    /// \brief The "original view" of the MultiVector's data.
    ///
    /// Methods like offsetView() return a view of a contiguous subset
    /// of rows.  At some point, we might like to get all of the rows
    /// back, by taking another view of a <i>super</i>set of rows.
    /// For example, we might like to get a column Map view of a
    /// (domain Map view of a (column Map MultiVector)).  Tpetra's
    /// implementation of Gauss-Seidel and SOR in CrsMatrix relies on
    /// this functionality.  However, Kokkos (rightfully) forbids us
    /// from taking a superset of rows of the current view.
    ///
    /// We deal with this at the Tpetra level by keeping around the
    /// original view of <i>all</i> the rows (and columns), which is
    /// \c origView_.  Methods like offsetView() then use origView_,
    /// not view_, to make the subview for the returned MultiVector.
    /// Furthermore, offsetView() can do error checking by getting the
    /// original number of rows from origView_.
    ///
    /// This may pose some problems for offsetView if it is given an
    /// offset other than zero, but that case is hardly exercised, so
    /// I am not going to worry about it for now.
    ///
    /// Note that the "original" view isn't always original.  It
    /// always has the original number of rows.  However, some special
    /// cases of constructors that take a whichVectors argument, when
    /// whichVectors.size() is 1, may point origView_ to the column to
    /// view.  Those constructors do this so that the resulting
    /// MultiVector has constant stride.  This special case does not
    /// affect correctness of offsetView and related methods.
    mutable dual_view_type origView_;

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

    //! @name View constructors, used only by nonmember constructors.
    //@{

    template <class S,class LO,class GO,class N>
    friend Teuchos::RCP<MultiVector<S,LO,GO,N> >
    createMultiVectorFromView (const Teuchos::RCP<const Map<LO, GO, N> >&,
                               const Teuchos::ArrayRCP<S>&,
                               const size_t, const size_t);

    bool vectorIndexOutOfRange (size_t VectorIndex) const;

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

    /// \brief Whether data redistribution between \c sourceObj and this object is legal.
    ///
    /// This method is called in DistObject::doTransfer() to check
    /// whether data redistribution between the two objects is legal.
    virtual bool
    checkSizes (const SrcDistObject& sourceObj);

    //! Number of packets to send per LID
    virtual size_t constantNumberOfPackets () const;

    /// \brief Whether lass implements old or new interface
    virtual bool useNewInterface () { return true; }

    virtual void
    copyAndPermuteNew (
      const SrcDistObject& sourceObj,
      size_t numSameIDs,
      const Kokkos::View<const local_ordinal_type*, device_type>& permuteToLIDs,
      const Kokkos::View<const local_ordinal_type*, device_type>& permuteFromLIDs);

    virtual void
    packAndPrepareNew (
      const SrcDistObject& sourceObj,
      const Kokkos::View<const local_ordinal_type*, device_type>& exportLIDs,
      Kokkos::View<scalar_type*, device_type>& exports,
      const Kokkos::View<size_t*, device_type>& numPacketsPerLID,
      size_t& constantNumPackets,
      Distributor &distor);

    virtual void
    unpackAndCombineNew (
      const Kokkos::View<const LocalOrdinal*, device_type>& importLIDs,
      const Kokkos::View<const scalar_type*, device_type>& imports,
      const Kokkos::View<size_t*, device_type>& numPacketsPerLID,
      size_t constantNumPackets,
      Distributor &distor,
      CombineMode CM);

    void createViews () const;
    void createViewsNonConst (KokkosClassic::ReadWriteOption rwo);
    void releaseViews () const;
    //@}

    typename dual_view_type::t_dev getKokkosView() const { return view_.d_view; }

  }; // class MultiVector

  namespace Details {

    // Partial specialization of MultiVectorCloner, for when the
    // source and destination MultiVector types are both Kokkos
    // refactor types, and when they both have the same Scalar type,
    // but all their other template parameters might be different.
    template<class ScalarType,
             class DstLocalOrdinalType, class DstGlobalOrdinalType, class DstDeviceType,
             class SrcLocalOrdinalType, class SrcGlobalOrdinalType, class SrcDeviceType>
    struct MultiVectorCloner< ::Tpetra::MultiVector<ScalarType,
                                                    DstLocalOrdinalType,
                                                    DstGlobalOrdinalType,
                                                    Kokkos::Compat::KokkosDeviceWrapperNode<DstDeviceType> >,
                              ::Tpetra::MultiVector<ScalarType,
                                                    SrcLocalOrdinalType,
                                                    SrcGlobalOrdinalType,
                                                    Kokkos::Compat::KokkosDeviceWrapperNode<SrcDeviceType> > >
    {
      typedef Kokkos::Compat::KokkosDeviceWrapperNode<DstDeviceType> dst_node_type;
      typedef Kokkos::Compat::KokkosDeviceWrapperNode<SrcDeviceType> src_node_type;
      typedef ::Tpetra::MultiVector<ScalarType, DstLocalOrdinalType,
                                    DstGlobalOrdinalType,
                                    dst_node_type> dst_mv_type;
      typedef ::Tpetra::MultiVector<ScalarType, SrcLocalOrdinalType,
                                    SrcGlobalOrdinalType,
                                    src_node_type> src_mv_type;

      static Teuchos::RCP<dst_mv_type>
      clone (const src_mv_type& X, const Teuchos::RCP<dst_node_type>& node2)
      {
        typedef typename src_mv_type::map_type src_map_type;
        typedef typename dst_mv_type::map_type dst_map_type;
        typedef typename dst_mv_type::node_type dst_node_type;
        typedef typename src_mv_type::dual_view_type src_dual_view_type;
        typedef typename dst_mv_type::dual_view_type dst_dual_view_type;

        // Clone X's Map to have the new Node type.
        RCP<const src_map_type> map1 = X.getMap ();
        RCP<const dst_map_type> map2 = map1.is_null () ?
          Teuchos::null : map1->template clone<dst_node_type> (node2);

        const size_t lclNumRows = X.getLocalLength ();
        const size_t numCols = X.getNumVectors ();
        src_dual_view_type X_view = X.getDualView ();
        dst_dual_view_type Y_view ("MV::dual_view", lclNumRows, numCols);

        RCP<dst_mv_type> Y = rcp (new dst_mv_type (map2, Y_view));
        // Let deep_copy do the work for us, to avoid code duplication.
        ::Tpetra::deep_copy (Y, X);
        return Y ;
      }
    };

    // Partial specialization for the Kokkos refactor specialization
    // of Tpetra::MultiVector.
    template<class S, class LO, class GO, class DeviceType>
    struct CreateMultiVectorFromView<MultiVector<S, LO, GO, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > > {
      typedef S scalar_type;
      typedef LO local_ordinal_type;
      typedef GO global_ordinal_type;
      typedef Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> node_type;
      typedef ::Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;
      typedef MultiVector<S, LO, GO, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > MultiVectorType;

      static Teuchos::RCP<MultiVectorType>
      create (const Teuchos::RCP<const map_type>& map,
              const Teuchos::ArrayRCP<scalar_type>& view,
              const size_t LDA,
              const size_t numVectors)
      {
        (void) map;
        (void) view;
        (void) LDA;
        (void) numVectors;

        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
      }
    };
  } // namespace Details

  template <class ST, class LO, class GO, class DT>
  MultiVector<ST, LO, GO, Kokkos::Compat::KokkosDeviceWrapperNode<DT> >
  createCopy (const MultiVector<ST, LO, GO, Kokkos::Compat::KokkosDeviceWrapperNode<DT> >& src);

  namespace { // (anonymous)
    template<class DstType, class SrcType, class IndexType, class DeviceType,
             const bool DstConstStride, const bool SrcConstStride>
    struct DeepCopySelectedVectors {
      typedef typename DeviceType::execution_space device_type;
      DstType dst_;
      SrcType src_;
      Kokkos::View<const IndexType*, DeviceType> whichVectorDst_;
      Kokkos::View<const IndexType*, DeviceType> whichVectorSrc_;
      const IndexType numVecs_;

      DeepCopySelectedVectors (DstType dst,
                               SrcType src,
                               const Kokkos::View<const IndexType*, DeviceType>& whichVectorDst,
                               const Kokkos::View<const IndexType*, DeviceType>& whichVectorSrc) :
        dst_ (dst),
        src_ (src),
        whichVectorDst_ (whichVectorDst),
        whichVectorSrc_ (whichVectorSrc),
        numVecs_ (whichVectorSrc_.dimension_0 ())
      {}

      DeepCopySelectedVectors (DstType dst, SrcType src) :
        dst_ (dst),
        src_ (src),
        numVecs_ (dst.dimension_1 ())
      {
        TEUCHOS_TEST_FOR_EXCEPTION(
          ! DstConstStride || ! SrcConstStride, std::logic_error,
          "Tpetra::DeepCopySelectedVectors: You may not use the constant-stride "
          "constructor if either of the Boolean template parameters is false.");
      }

      void KOKKOS_INLINE_FUNCTION operator () (const IndexType i) const {
        if (DstConstStride) {
          if (SrcConstStride) {
            for (IndexType j = 0; j < numVecs_; ++j) {
              dst_(i,j) = src_(i,j);
            }
          } else {
            for (IndexType j = 0; j < numVecs_; ++j) {
              dst_(i,j) = src_(i,whichVectorSrc_(j));
            }
          }
        } else {
          if (SrcConstStride) {
            for (IndexType j = 0; j < numVecs_; ++j) {
              dst_(i,whichVectorDst_(j)) = src_(i,j);
            }
          } else {
            for (IndexType j = 0; j < numVecs_; ++j) {
              dst_(i,whichVectorDst_(j)) = src_(i,whichVectorSrc_(j));
            }
          }
        }
      }
    };

    template<class DstViewType,
             class SrcViewType,
             class IndexType>
    struct LocalDeepCopy {
      typedef Kokkos::View<const IndexType*, typename DstViewType::memory_space> dst_wv_type;
      typedef Kokkos::View<const IndexType*, typename SrcViewType::memory_space> src_wv_type;

      static void
      run (const DstViewType& dst, const SrcViewType& src,
           const bool dstConstStride, const bool srcConstStride,
           const dst_wv_type& dstWhichVecs, const src_wv_type& srcWhichVecs)
      {
        typedef typename DstViewType::execution_space execution_space;

        // FIXME (mfh 22 Jul 2014, 10 Dec 2014) Currently, it doesn't
        // work to do a 2-D copy, even if both MultiVectors have
        // constant stride.  This is because Kokkos can't currently
        // tell the difference between padding (which permits a single
        // deep_copy for the whole 2-D View) and stride > numRows
        // (which does NOT permit a single deep_copy for the whole 2-D
        // View).  Carter is working on this, but for now, the
        // temporary fix is to copy one column at a time.
        const bool deepCopyWorks = false;

        if (dstConstStride) {
          if (srcConstStride) {
            if (deepCopyWorks) {
              Kokkos::deep_copy (dst, src);
            }
            else {
              typedef DeepCopySelectedVectors<DstViewType, SrcViewType, IndexType, execution_space, true, true> functor_type;
              functor_type f (dst, src, dstWhichVecs, srcWhichVecs);
              Kokkos::parallel_for (dst.dimension_0 (), f);
            }
          }
          else { // ! srcConstStride
            typedef DeepCopySelectedVectors<DstViewType, SrcViewType, IndexType, execution_space, true, false> functor_type;
            functor_type f (dst, src, dstWhichVecs, srcWhichVecs);
            Kokkos::parallel_for (dst.dimension_0 (), f);
          }
        }
        else { // ! dstConstStride
          if (srcConstStride) {
            typedef DeepCopySelectedVectors<DstViewType, SrcViewType, IndexType, execution_space, false, true> functor_type;
            functor_type f (dst, src, dstWhichVecs, srcWhichVecs);
            Kokkos::parallel_for (dst.dimension_0 (), f);
          }
          else { // ! srcConstStride
            typedef DeepCopySelectedVectors<DstViewType, SrcViewType, IndexType, execution_space, false, false> functor_type;
            functor_type f (dst, src, dstWhichVecs, srcWhichVecs);
            Kokkos::parallel_for (dst.dimension_0 (), f);
          }
        }
      }
    };

    template<class DstViewType,
             class SrcViewType,
             class IndexType>
    void
    localDeepCopy (const DstViewType& dst, const SrcViewType& src,
                   const bool dstConstStride, const bool srcConstStride,
                   const Kokkos::View<const IndexType*, typename DstViewType::memory_space>& dstWhichVecs,
                   const Kokkos::View<const IndexType*, typename SrcViewType::memory_space>& srcWhichVecs)
    {
      LocalDeepCopy<DstViewType, SrcViewType, IndexType>::run (dst, src, dstConstStride, srcConstStride, dstWhichVecs, srcWhichVecs);
    }
  } // namespace (anonymous)


  // NOTE (mfh 11 Sep 2014) Even though this partial specialization
  // looks redundant with the one in Tpetra_MultiVector_decl.hpp, it
  // needs to be here, else GCC 4.8.2 gives a compiler error saying
  // that calls to Tpetra::deep_copy are ambiguous.
  template <class ST, class LO, class GO, class DT>
  void
  deep_copy (MultiVector<ST, LO, GO, Kokkos::Compat::KokkosDeviceWrapperNode<DT> >& dst,
             const MultiVector<ST, LO, GO, Kokkos::Compat::KokkosDeviceWrapperNode<DT> >& src)
  {
    // NOTE (mfh 11 Sep 2014) We can't implement deep_copy with
    // shallow-copy operator=, because that would invalidate existing
    // views of dst!
    dst.assign (src);
  }


  template <class DS, class DL, class DG, class DD,
            class SS, class SL, class SG, class SD>
  void
  deep_copy (MultiVector<DS, DL, DG, Kokkos::Compat::KokkosDeviceWrapperNode<DD> >& dst,
             const MultiVector<SS, SL, SG, Kokkos::Compat::KokkosDeviceWrapperNode<SD> >& src)
  {
    using Kokkos::parallel_for;
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<SD> src_node_type;
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<DD> dst_node_type;
    typedef MultiVector<DS, DL, DG, dst_node_type> MVD;
    typedef const MultiVector<SS, SL, SG, src_node_type> MVS;

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

    if (src.isConstantStride () && dst.isConstantStride ()) {
      Kokkos::deep_copy (dst.getDualView (), src.getDualView ());
    }
    else {
      typedef Kokkos::DualView<SL*, DD> whichvecs_type;
      typedef typename whichvecs_type::host_mirror_space host_mirror_space ;

      if (dst.isConstantStride ()) {
        const SL numWhichVecs = static_cast<SL> (src.whichVectors_.size ());
        const std::string whichVecsLabel ("MV::deep_copy::whichVecs");

        // We can't sync src, since it is only an input argument.
        // Thus, we have to use the most recently modified version of
        // src, which coudl be either the device or host version.
        if (src.getDualView ().modified_device >= src.getDualView ().modified_host) {
          // Copy from the device version of src.
          //
          // whichVecs tells the kernel which vectors (columns) of src
          // to copy.  Fill whichVecs on the host, and sync to device.
          whichvecs_type whichVecs (whichVecsLabel, numWhichVecs);
          whichVecs.template modify<host_mirror_space> ();
          for (SL i = 0; i < numWhichVecs; ++i) {
            whichVecs.h_view(i) = static_cast<SL> (src.whichVectors_[i]);
          }
          // Sync the host version of whichVecs to the device.
          whichVecs.template sync<DD> ();

          // Mark the device version of dst's DualView as modified.
          dst.template modify<DD> ();
          // Copy from the selected vectors of src to dst, on the
          // device.  The functor ignores its 3rd arg in this case.
          typedef DeepCopySelectedVectors<typename MVD::dual_view_type::t_dev,
            typename MVS::dual_view_type::t_dev, SL, DD, true, false> functor_type;
          functor_type f (dst.getDualView ().template view<DD> (),
                          src.getDualView ().template view<DD> (),
                          whichVecs.d_view, whichVecs.d_view);
          Kokkos::parallel_for (src.getLocalLength (), f);
          // Sync dst's DualView to the host.  This is cheaper than
          // repeating the above copy from src to dst on the host.
          dst.template sync<host_mirror_space> ();
        }
        else { // host version of src was the most recently modified
          // Copy from the host version of src.
          //
          // whichVecs tells the kernel which vectors (columns) of src
          // to copy.  Fill whichVecs on the host, and use it there.
          typedef Kokkos::View<SL*, host_mirror_space> the_whichvecs_type;
          the_whichvecs_type whichVecs (whichVecsLabel, numWhichVecs);
          for (SL i = 0; i < numWhichVecs; ++i) {
            whichVecs(i) = static_cast<SL> (src.whichVectors_[i]);
          }
          // Copy from the selected vectors of src to dst, on the host.
          // The functor ignores its 3rd arg in this case.
          typedef DeepCopySelectedVectors<typename MVD::dual_view_type::t_host,
            typename MVS::dual_view_type::t_host, SL, host_mirror_space,
            true, false> functor_type;
          functor_type f (dst.getDualView ().template view<host_mirror_space> (),
                          src.getDualView ().template view<host_mirror_space> (),
                          whichVecs, whichVecs);
          Kokkos::parallel_for (src.getLocalLength (), f);
          // Sync dst back to the device, since we only copied on the host.
          dst.template sync<DD> ();
        }
      }
      else { // dst is NOT constant stride
        typedef typename Kokkos::ViewTraits<DL*,DD,void,void>::host_mirror_space
          host_mirror_space_type;

        if (src.isConstantStride ()) {
          if (src.getDualView ().modified_device >= src.getDualView ().modified_host) {
            // Copy from the device version of src.
            //
            // whichVecs tells the kernel which vectors (columns) of dst
            // to copy.  Fill whichVecs on the host, and sync to device.
            typedef Kokkos::DualView<DL*, DD> the_whichvecs_type;
            const std::string whichVecsLabel ("MV::deep_copy::whichVecs");
            const DL numWhichVecs = static_cast<DL> (dst.whichVectors_.size ());
            the_whichvecs_type whichVecs (whichVecsLabel, numWhichVecs);
            whichVecs.template modify<host_mirror_space_type> ();
            for (DL i = 0; i < numWhichVecs; ++i) {
              whichVecs.h_view(i) = dst.whichVectors_[i];
            }
            // Sync the host version of whichVecs to the device.
            whichVecs.template sync<DD> ();

            // Copy src to the selected vectors of dst, on the device.
            // The functor ignores its 4th arg in this case.
            typedef DeepCopySelectedVectors<typename MVD::dual_view_type::t_dev,
              typename MVS::dual_view_type::t_dev, DL, DD, false, true> functor_type;
            functor_type f (dst.getDualView ().template view<DD> (),
                            src.getDualView ().template view<DD> (),
                            whichVecs.d_view, whichVecs.d_view);
            Kokkos::parallel_for (src.getLocalLength (), f);
            // We can't sync src and repeat the above copy on the
            // host, so sync dst back to the host.
            //
            // FIXME (mfh 29 Jul 2014) This may overwrite columns that
            // don't actually belong to dst's view.
            dst.template sync<host_mirror_space_type> ();
          }
          else { // host version of src was the most recently modified
            // Copy from the host version of src.
            //
            // whichVecs tells the kernel which vectors (columns) of src
            // to copy.  Fill whichVecs on the host, and use it there.
            typedef Kokkos::View<DL*, host_mirror_space_type> the_whichvecs_type;
            const DL numWhichVecs = static_cast<DL> (dst.whichVectors_.size ());
            the_whichvecs_type whichVecs ("MV::deep_copy::whichVecs", numWhichVecs);
            for (DL i = 0; i < numWhichVecs; ++i) {
              whichVecs(i) = static_cast<DL> (dst.whichVectors_[i]);
            }
            // Copy from src to the selected vectors of dst, on the
            // host.  The functor ignores its 4th arg in this case.
            typedef DeepCopySelectedVectors<typename MVD::dual_view_type::t_host,
              typename MVS::dual_view_type::t_host, DL, host_mirror_space_type,
              false, true> functor_type;
            functor_type f (dst.getDualView ().template view<host_mirror_space_type> (),
                            src.getDualView ().template view<host_mirror_space_type> (),
                            whichVecs, whichVecs);
            Kokkos::parallel_for (src.getLocalLength (), f);
            // Sync dst back to the device, since we only copied on the host.
            //
            // FIXME (mfh 29 Jul 2014) This may overwrite columns that
            // don't actually belong to dst's view.
            dst.template sync<DD> ();
          }
        }
        else { // neither src nor dst have constant stride
          if (src.getDualView ().modified_device >= src.getDualView ().modified_host) {
            // Copy from the device version of src.
            //
            // whichVectorsDst tells the kernel which vectors
            // (columns) of dst to copy.  Fill it on the host, and
            // sync to device.
            const DL dstNumWhichVecs = static_cast<DL> (dst.whichVectors_.size ());
            Kokkos::DualView<DL*, DD> whichVecsDst ("MV::deep_copy::whichVecsDst",
                                                    dstNumWhichVecs);
            whichVecsDst.template modify<host_mirror_space> ();
            for (DL i = 0; i < dstNumWhichVecs; ++i) {
              whichVecsDst.h_view(i) = static_cast<DL> (dst.whichVectors_[i]);
            }
            // Sync the host version of whichVecsDst to the device.
            whichVecsDst.template sync<DD> ();

            // whichVectorsSrc tells the kernel which vectors
            // (columns) of src to copy.  Fill it on the host, and
            // sync to device.  Use the destination MultiVector's
            // LocalOrdinal type here.
            const DL srcNumWhichVecs = static_cast<DL> (src.whichVectors_.size ());
            Kokkos::DualView<DL*, DD> whichVecsSrc ("MV::deep_copy::whichVecsSrc",
                                                    srcNumWhichVecs);
            whichVecsSrc.template modify<host_mirror_space> ();
            for (DL i = 0; i < srcNumWhichVecs; ++i) {
              whichVecsSrc.h_view(i) = static_cast<DL> (src.whichVectors_[i]);
            }
            // Sync the host version of whichVecsSrc to the device.
            whichVecsSrc.template sync<DD> ();

            // Copy from the selected vectors of src to the selected
            // vectors of dst, on the device.
            typedef DeepCopySelectedVectors<typename MVD::dual_view_type::t_dev,
              typename MVS::dual_view_type::t_dev, DL, DD, false, false>
              functor_type;
            functor_type f (dst.getDualView ().template view<DD> (),
                            src.getDualView ().template view<DD> (),
                            whichVecsDst.d_view, whichVecsSrc.d_view);
            Kokkos::parallel_for (src.getLocalLength (), f);
          }
          else {
            const DL dstNumWhichVecs = static_cast<DL> (dst.whichVectors_.size ());
            Kokkos::View<DL*, host_mirror_space> whichVectorsDst ("dstWhichVecs", dstNumWhichVecs);
            for (DL i = 0; i < dstNumWhichVecs; ++i) {
              whichVectorsDst(i) = dst.whichVectors_[i];
            }

            // Use the destination MultiVector's LocalOrdinal type here.
            const DL srcNumWhichVecs = static_cast<DL> (src.whichVectors_.size ());
            Kokkos::View<DL*, host_mirror_space> whichVectorsSrc ("srcWhichVecs", srcNumWhichVecs);
            for (DL i = 0; i < srcNumWhichVecs; ++i) {
              whichVectorsSrc(i) = src.whichVectors_[i];
            }

            typedef DeepCopySelectedVectors<typename MVD::dual_view_type::t_host,
              typename MVS::dual_view_type::t_host,
              DL, host_mirror_space, false, false> functor_type;
            functor_type f (dst.getDualView ().template view<host_mirror_space> (),
                            src.getDualView ().template view<host_mirror_space> (),
                            whichVectorsDst, whichVectorsSrc);
            Kokkos::parallel_for (src.getLocalLength (), f);

            // We can't sync src and repeat the above copy on the
            // host, so sync dst back to the host.
            //
            // FIXME (mfh 29 Jul 2014) This may overwrite columns that
            // don't actually belong to dst's view.
            dst.template sync<host_mirror_space> ();
          }
        }
      }
    }
  }
} // namespace Tpetra

#endif // TPETRA_KOKKOS_REFACTOR_MULTIVECTOR_DECL_HPP
