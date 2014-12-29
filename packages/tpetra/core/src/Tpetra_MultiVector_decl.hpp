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

#ifndef TPETRA_MULTIVECTOR_DECL_HPP
#define TPETRA_MULTIVECTOR_DECL_HPP

#include <Teuchos_DataAccess.hpp>
#include <Teuchos_Range1D.hpp>
#include "Tpetra_ConfigDefs.hpp"
#include <Tpetra_Map_decl.hpp>
#if TPETRA_USE_KOKKOS_DISTOBJECT
#include "Tpetra_DistObjectKA.hpp"
#else
#include "Tpetra_DistObject.hpp"
#endif
#include "Tpetra_ViewAccepter.hpp"
#include <Kokkos_MultiVector.hpp>
#include <Teuchos_BLAS_types.hpp>

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
#endif // DOXYGEN_SHOULD_SKIP_THIS


  namespace Details {
    /// \brief Implementation of createMultiVectorFromView
    /// \tparam MultiVectorType A specialization of Tpetra::MultiVector.
    ///
    /// This struct lets us do partial specialization of the nonmember
    /// template function createMultiVectorFromView.  This is
    /// particularly useful so that we can partially specialize this
    /// function for the new Kokkos refactor specializations of
    /// MultiVector.
    template<class MultiVectorType>
    struct CreateMultiVectorFromView {
      typedef typename MultiVectorType::scalar_type scalar_type;
      typedef typename MultiVectorType::local_ordinal_type local_ordinal_type;
      typedef typename MultiVectorType::global_ordinal_type global_ordinal_type;
      typedef typename MultiVectorType::node_type node_type;
      typedef ::Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

      static Teuchos::RCP<MultiVectorType>
      create (const Teuchos::RCP<const map_type>& map,
              const Teuchos::ArrayRCP<scalar_type>& view,
              const size_t LDA,
              const size_t numVectors);
    };
  } // namespace Details


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
  ///   may use real-valued or complex-valued types here, unlike in
  ///   Epetra, where the scalar type is always \c double.)  The
  ///   default is \c double (real, double-precision floating-point
  ///   type).  You may use any type here that has a
  ///   Teuchos::ScalarTraits specialization.
  /// \tparam LocalOrdinal The type of local indices.  See the
  ///   documentation of Map for requirements.
  /// \tparam GlobalOrdinal The type of global indices.  See the
  ///   documentation of Map for requirements.
  /// \tparam Node The Kokkos Node type.  See the documentation of Map
  ///   for requirements.
  ///
  /// \section Kokkos_MV_prereq Prerequisites
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
  /// \section Kokkos_MV_layout Data layout and access
  ///
  /// \subsection Kokkos_MV_noncontig Multivectors which view another multivector
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
  /// \subsection Kokkos_MV_views Views of a multivector's data
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
  /// \subsection Kokkos_MV_why_views Why won't you give me a raw pointer?
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
  /// \subsection Kokkos_MV_why_no_square_brackets Why no operator[]?
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
  /// \subsection Kokkos_MV_device_kernels How do I access the data directly?
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
  /// \section Kokkos_MV_dist Parallel distribution of data
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
  template<class Scalar = double,
           class LocalOrdinal = Map<>::local_ordinal_type,
           class GlobalOrdinal = typename Map<LocalOrdinal>::global_ordinal_type,
           class Node = typename Map<LocalOrdinal, GlobalOrdinal>::node_type>
  class MultiVector :
#if TPETRA_USE_KOKKOS_DISTOBJECT
    public DistObjectKA<Scalar, LocalOrdinal, GlobalOrdinal, Node>
#else
    public DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node>
#endif
  {
  public:
    //! @name Typedefs to facilitate template metaprogramming.
    //@{

    //! The type of entries in the vector(s).
    typedef Scalar scalar_type;
    //! The type of local indices.
    typedef LocalOrdinal local_ordinal_type;
    //! The type of global indices.
    typedef GlobalOrdinal global_ordinal_type;
    //! The Kokkos Node type.
    typedef Node node_type;

    /// \brief Type of an inner ("dot") product result.
    ///
    /// This is usually the same as <tt>scalar_type</tt>, but may
    /// differ if <tt>scalar_type</tt> is e.g., an uncertainty
    /// quantification type from the Stokhos package.
    typedef scalar_type dot_type;

    /// \brief Type of a norm result.
    ///
    /// This is usually the same as the type of the magnitude
    /// (absolute value) of <tt>scalar_type</tt>, but may differ if
    /// <tt>scalar_type</tt> is e.g., an uncertainty quantification
    /// type from the Stokhos package.
    typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType mag_type;

#if TPETRA_USE_KOKKOS_DISTOBJECT
    typedef DistObjectKA<Scalar, LocalOrdinal, GlobalOrdinal, Node> DO;
    typedef typename DO::device_type device_type;
#else
    typedef DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node> DO;
#endif

    //@}
    //! \name Constructors and destructor
    //@{

    //! Default constructor: makes a MultiVector with no rows or columns.
    MultiVector ();

    /// \brief Basic constructor.
    ///
    /// \param map [in] The MultiVector's Map.  The Map describes the
    ///   distribution of rows over process(es) in the Map's
    ///   communicator.
    ///
    /// \param NumVectors [in] Number of columns in the MultiVector.
    ///
    /// \param zeroOut [in] If true (the default), require that all the
    ///   Vector's entries be zero on return.  If false, the Vector's
    ///   entries have undefined values on return, and must be set
    ///   explicitly.
    MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
                 size_t NumVectors,
                 bool zeroOut=true);

    /// \brief Copy constructor.
    ///
    /// Whether this does a deep copy or a shallow copy depends on
    /// whether \c source has "view semantics."  See discussion in the
    /// documentation of the two-argument copy constructor below.
    MultiVector (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &source);

    /// \brief Copy constructor, with option to do shallow copy and
    ///   mark the result as having "view semantics."
    ///
    /// If copyOrView is Teuchos::View, this constructor marks the
    /// result as having "view semantics."  This means that copy
    /// construction or assignment (operator=) with the resulting
    /// object will always do a shallow copy, and will transmit view
    /// semantics to the result of the shallow copy.  If copyOrView is
    /// Teuchos::Copy, this constructor <i>always</i> does a deep copy
    /// and marks the result as not having view semantics, whether or
    /// not \c source has view semantics.
    ///
    /// View semantics are a "forwards compatibility" measure for
    /// porting to the Kokkos refactor version of Tpetra.  The latter
    /// only ever has view semantics.  The "classic" version of Tpetra
    /// does not currently have view semantics by default, but this
    /// will change.
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
                 size_t LDA,
                 size_t NumVectors);

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
                 size_t NumVectors);

    /// \brief Expert mode constructor.
    ///
    /// \warning This constructor is only for expert users.  We make
    ///   no promises about backwards compatibility for this interface.
    ///   It may change or go away at any time.
    ///
    /// \param map [in] Map describing the distribution of rows.
    /// \param data [in] Device pointer to the data (column-major)
    /// \param LDA [in] Leading dimension (stride) of the data
    /// \param numVecs [in] Number of vectors (columns)
    MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
                 const Teuchos::ArrayRCP<Scalar>& data,
                 const size_t LDA,
                 const size_t numVectors);

    //! Create a cloned MultiVector for a different node type
    template <class Node2>
    Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node2> >
    clone(const Teuchos::RCP<Node2> &node2) const;

    //! Destructor (virtual for memory safety of derived classes).
    virtual ~MultiVector();

    //@}
    //! @name Post-construction modification routines
    //@{

    /// \brief Replace value, using global (row) index.
    ///
    /// Replace the current value at row globalRow (a global index)
    /// and column vectorIndex with the given value.  The column index
    /// is zero based.
    ///
    /// \pre \c globalRow must be a valid global element on this
    ///   process, according to the row Map.
    void
    replaceGlobalValue (GlobalOrdinal globalRow,
                        size_t vectorIndex,
                        const Scalar &value);

    /// \brief Add value to existing value, using global (row) index.
    ///
    /// Add the given value to the existing value at row globalRow (a
    /// global index) and column vectorIndex.  The column index is
    /// zero based.
    ///
    /// \pre \c globalRow must be a valid global element on this
    ///   process, according to the row Map.
    void
    sumIntoGlobalValue (GlobalOrdinal globalRow,
                        size_t vectorIndex,
                        const Scalar &value);

    /// \brief Replace value, using local (row) index.
    ///
    /// Replace the current value at row myRow (a local index) and
    /// column vectorIndex with the given value.  The column index is
    /// zero based.
    ///
    /// \pre \c myRow must be a valid local element on this process,
    ///   according to the row Map.
    void
    replaceLocalValue (LocalOrdinal myRow,
                       size_t vectorIndex,
                       const Scalar &value);

    /// \brief Add value to existing value, using local (row) index.
    ///
    /// Add the given value to the existing value at row myRow (a
    /// local index) and column vectorIndex.  The column index is
    /// zero based.
    ///
    /// \pre \c myRow must be a valid local element on this process,
    ///   according to the row Map.
    void
    sumIntoLocalValue (LocalOrdinal myRow,
                       size_t vectorIndex,
                       const Scalar &value);

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

    /// \brief Assignment operator.
    ///
    /// If this MultiVector (the left-hand side of the assignment) has
    /// view semantics (<tt>getCopyOrView() == Teuchos::View</tt>),
    /// then this does a shallow copy.  Otherwise, it does a deep
    /// copy.  The latter is the default behavior.
    ///
    /// A deep copy has the following prerequisites:
    ///
    /// \pre The input MultiVector's Map must be compatible with this
    ///      multivector's Map.  That is, \code
    ///      this->getMap ()->isCompatible (source.getMap ());
    ///      \endcode
    /// \pre Both MultiVectors must have the same number of columns.
    ///
    /// \note This method must always be called as a collective
    ///   operation on all processes over which the multivector is
    ///   distributed.  This is because the method reserves the right
    ///   to check for compatibility of the two Maps, at least in
    ///   debug mode.
    MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&
    operator= (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& source);

    //@}

    /// \name Data copy and view methods
    ///
    /// These methods are used to get the data underlying the
    /// MultiVector. They return data in one of two forms:
    /// <ul>
    /// <li> A MultiVector with a subset of the rows or columns of the
    ///      input MultiVector </li>
    /// <li> An array of data, wrapped in one of the Teuchos memory
    ///      management classes </li>
    /// </ul>
    /// Not all of these methods are valid for a particular
    /// MultiVector. For instance, calling a method that accesses a
    /// view of the data in a 1-D format (i.e., get1dView) requires
    /// that the target MultiVector has constant stride.
    //@{

    /// \brief Return a MultiVector with copies of selected columns.
    ///
    /// \param colRng [in] Inclusive, contiguous range of columns.
    ///   <tt>[colRng.lbound(), colRng.ubound()]</tt> defines the range.
    Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
    subCopy (const Teuchos::Range1D &colRng) const;

    //! Return a MultiVector with copies of selected columns.
    Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
    subCopy (const Teuchos::ArrayView<const size_t> &cols) const;

    /// \brief Return a MultiVector with const views of selected columns.
    ///
    /// \param colRng [in] Inclusive, contiguous range of columns.
    ///   <tt>[colRng.lbound(), colRng.ubound()]</tt> defines the range.
    Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
    subView (const Teuchos::Range1D &colRng) const;

    //! Return a const MultiVector with const views of selected columns.
    Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
    subView (const Teuchos::ArrayView<const size_t> &cols) const;

    /// \brief Return a MultiVector with views of selected columns.
    ///
    /// \param colRng [in] Inclusive, contiguous range of columns.
    ///   <tt>[colRng.lbound(), colRng.ubound()]</tt> defines the range.
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
    Teuchos::ArrayRCP<const Scalar> getData(size_t j) const;

    //! View of the local values in a particular vector of this multivector.
    Teuchos::ArrayRCP<Scalar> getDataNonConst(size_t j);

    /// \brief Fill the given array with a copy of this multivector's local values.
    ///
    /// \param A [out] View of the array to fill.  We consider A as a
    ///   matrix with column-major storage.
    ///
    /// \param LDA [in] Leading dimension of the matrix A.
    void get1dCopy (Teuchos::ArrayView<Scalar> A, size_t LDA) const;

    /// \brief Fill the given array with a copy of this multivector's local values.
    ///
    /// \param ArrayOfPtrs [out] Array of arrays, one for each column
    ///   of the multivector.  On output, we fill ArrayOfPtrs[j] with
    ///   the data for column j of this multivector.
    void get2dCopy (Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> > ArrayOfPtrs) const;

    /// \brief Const persisting (1-D) view of this multivector's local values.
    ///
    /// This method assumes that the columns of the multivector are
    /// stored contiguously.  If not, this method throws
    /// std::runtime_error.
    Teuchos::ArrayRCP<const Scalar> get1dView() const;

    //! Return const persisting pointers to values.
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > get2dView() const;

    /// \brief Nonconst persisting (1-D) view of this multivector's local values.
    ///
    /// This method assumes that the columns of the multivector are
    /// stored contiguously.  If not, this method throws
    /// std::runtime_error.
    Teuchos::ArrayRCP<Scalar> get1dViewNonConst();

    //! Return non-const persisting pointers to values.
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > get2dViewNonConst();

    /// \brief A view of the underlying KokkosClassic::MultiVector object.
    ///
    /// \brief This method is for expert users only.
    ///   It may change or be removed at any time.
    KokkosClassic::MultiVector<Scalar,Node> getLocalMV () const;

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
    TPETRA_DEPRECATED
    KokkosClassic::MultiVector<Scalar,Node>& getLocalMVNonConst ();

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
         const Teuchos::ArrayView<Scalar>& dots) const;

    //! Put element-wise absolute values of input Multi-vector in target: A = abs(this)
    void abs(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A);

    //! Put element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
    void reciprocal(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A);

    /// \brief Scale in place: <tt>this = alpha*this</tt>.
    ///
    /// Replace this MultiVector with alpha times this MultiVector.
    /// This method will always multiply, even if alpha is zero.  That
    /// means, for example, that if \c *this contains NaN entries
    /// before calling this method, the NaN entries will remain after
    /// this method finishes.
    void scale (const Scalar &alpha);

    /// \brief Scale each column in place: <tt>this[j] = alpha[j]*this[j]</tt>.
    ///
    /// Replace each column j of this MultiVector with
    /// <tt>alpha[j]</tt> times the current column j of this
    /// MultiVector.  This method will always multiply, even if all
    /// the entries of alpha are zero.  That means, for example, that
    /// if \c *this contains NaN entries before calling this method,
    /// the NaN entries will remain after this method finishes.
    void scale (Teuchos::ArrayView<const Scalar> alpha);

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

    /// \brief Compute the one-norm of each vector (column).
    ///
    /// The one-norm of a vector is the sum of squares of the
    /// magnitudes of the vector's entries.  On exit, norms[k] is the
    /// one-norm of column k of this MultiVector.
    void norm1(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const;

    /// \brief Compute the two-norm of each vector (column).
    ///
    /// The two-norm of a vector is the standard Euclidean norm, the
    /// square root of the sum of squares of the magnitudes of the
    /// vector's entries.  On exit, norms[k] is the two-norm of column
    /// k of this MultiVector.
    void norm2(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const;

    /// \brief Compute the infinity-norm of each vector (column).
    ///
    /// The infinity-norm of a vector is the maximum of the magnitudes
    /// of the vector's entries.  On exit, norms[k] is the
    /// infinity-norm of column k of this MultiVector.
    void normInf (const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& norms) const;

    //! Compute Weighted 2-norm (RMS Norm) of each vector in multi-vector.
    //! The outcome of this routine is undefined for non-floating point scalar types (e.g., int).
    void
    normWeighted (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& weights,
                  const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& norms) const;

    //! \brief Compute mean (average) value of each vector in multi-vector.
    //! The outcome of this routine is undefined for non-floating point scalar types (e.g., int).
    void meanValue (const Teuchos::ArrayView<Scalar>& means) const;

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
    /// \warning This method is only for expert use.  It may change or
    ///   disappear at any time.
    void setCopyOrView (const Teuchos::DataAccess copyOrView) {
      hasViewSemantics_ = (copyOrView == Teuchos::View);
    }

    /// \brief Get whether this has copy (copyOrView = Teuchos::Copy)
    ///   or view (copyOrView = Teuchos::View) semantics.
    ///
    /// \warning This method is only for expert use.  It may change or
    ///   disappear at any time.
    Teuchos::DataAccess getCopyOrView () const {
      return hasViewSemantics_ ? Teuchos::View : Teuchos::Copy;
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

  protected:
    typedef KokkosClassic::MultiVector<Scalar,Node> KMV;
    typedef KokkosClassic::DefaultArithmetic<KMV>   MVT;

    //! The KokkosClassic::MultiVector containing the compute buffer of data.
    KMV lclMV_;

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
    Array<size_t> whichVectors_;

    /// \brief Whether this MultiVector has view semantics.
    ///
    /// "View semantics" means that if this MultiVector is on the
    /// right side of an operator=, the left side gets a shallow copy,
    /// and acquires view semantics.  The Kokkos refactor version of
    /// MultiVector only ever has view semantics.  The "classic"
    /// version of MultiVector currently does not have view semantics
    /// by default, but this will change.
    ///
    /// You can set this for now by calling one of the constructors
    /// that accepts a Teuchos::DataAccess enum value.
    bool hasViewSemantics_;

    //! \name View constructors, used only by nonmember constructors.
    //@{

    // Implementation detail of the nonmember "constructor" function
    // createMultiVectorFromView.  Please consider this function
    // DEPRECATED.
    template <class MultiVectorType>
    friend struct Details::CreateMultiVectorFromView;

    /// \brief View constructor with user-allocated data.
    ///
    /// Please consider this constructor DEPRECATED.
    ///
    /// The tag says that views of the MultiVector are always host
    /// views, that is, they do not live on a separate device memory
    /// space (for example, on a GPU).
    ///
    /// This member constructor is meant to be called by its nonmember
    /// constructor friend; it is not meant to be called by users
    /// (hence it is protected).
    MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
                 const Teuchos::ArrayRCP<Scalar>& view,
                 size_t LDA,
                 size_t NumVectors,
                 EPrivateHostViewConstructor /* dummy */);

    bool vectorIndexOutOfRange (size_t VectorIndex) const;

    /// \fn getSubArrayRCP
    /// \brief Persisting view of j-th column in the given ArrayRCP.
    ///
    /// This method considers isConstantStride().  The ArrayRCP may
    /// correspond either to a compute buffer or a host view.
    template <class T>
    ArrayRCP<T> getSubArrayRCP(ArrayRCP<T> arr, size_t j) const;

    /// \brief Advanced constructorfor non-contiguous views.
    ///
    /// Please consider this constructor DEPRECATED.
    MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
                 Teuchos::ArrayRCP<Scalar> data,
                 size_t LDA,
                 Teuchos::ArrayView<const size_t> whichVectors,
                 EPrivateComputeViewConstructor /* dummy */);

    /// \brief Advanced constructor for contiguous views.
    ///
    /// Please consider this constructor DEPRECATED.
    MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
                 Teuchos::ArrayRCP<Scalar> data,
                 size_t LDA,
                 size_t NumVectors,
                 EPrivateComputeViewConstructor /* dummy */);

    /// \brief Advanced constructor for contiguous views.
    ///
    /// This version of the contiguous view constructor takes a
    /// previously constructed KokkosClassic::MultiVector, which views
    /// the local data.  The local multivector should have been made
    /// using the appropriate offsetView* method of
    /// KokkosClassic::MultiVector.
    MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
                 const KokkosClassic::MultiVector<Scalar, Node>& localMultiVector,
                 EPrivateComputeViewConstructor /* dummy */);

    /// \brief Advanced constructor for noncontiguous views.
    ///
    /// This version of the noncontiguous view constructor takes a
    /// previously constructed KokkosClassic::MultiVector, which is the
    /// correct view of the local data.  The local multivector should
    /// have been made using the appropriate offsetView* method of
    /// KokkosClassic::MultiVector.
    MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
                 const KokkosClassic::MultiVector<Scalar, Node>& localMultiVector,
                 Teuchos::ArrayView<const size_t> whichVectors,
                 EPrivateComputeViewConstructor /* dummy */);

    //@}
    //! \name Implementation of Tpetra::DistObject
    //@{

    /// \brief Whether data redistribution between \c sourceObj and this object is legal.
    ///
    /// This method is called in DistObject::doTransfer() to check
    /// whether data redistribution between the two objects is legal.
    virtual bool
    checkSizes (const SrcDistObject& sourceObj);

    //! Number of packets to send per LID
    virtual size_t constantNumberOfPackets () const;

#if TPETRA_USE_KOKKOS_DISTOBJECT
    virtual void
    copyAndPermute (
      const SrcDistObject& sourceObj,
      size_t numSameIDs,
      const Kokkos::View<const LocalOrdinal*, device_type> &permuteToLIDs,
      const Kokkos::View<const LocalOrdinal*, device_type> &permuteFromLIDs);

    virtual void
    packAndPrepare (
      const SrcDistObject& sourceObj,
      const Kokkos::View<const LocalOrdinal*, device_type> &exportLIDs,
      Kokkos::View<Scalar*, device_type> &exports,
      const Kokkos::View<size_t*, device_type> &numPacketsPerLID,
      size_t& constantNumPackets,
      Distributor &distor);

    virtual void
    unpackAndCombine (
      const Kokkos::View<const LocalOrdinal*, device_type> &importLIDs,
      const Kokkos::View<const Scalar*, device_type> &imports,
      const Kokkos::View<size_t*, device_type> &numPacketsPerLID,
      size_t constantNumPackets,
      Distributor &distor,
      CombineMode CM);
#else
    virtual void
    copyAndPermute (const SrcDistObject& sourceObj,
                    size_t numSameIDs,
                    const ArrayView<const LocalOrdinal>& permuteToLIDs,
                    const ArrayView<const LocalOrdinal>& permuteFromLIDs);

    virtual void
    packAndPrepare (const SrcDistObject& sourceObj,
                    const ArrayView<const LocalOrdinal>& exportLIDs,
                    Array<Scalar>& exports,
                    const ArrayView<size_t>& numExportPacketsPerLID,
                    size_t& constantNumPackets,
                    Distributor& distor);

    virtual void
    unpackAndCombine (const ArrayView<const LocalOrdinal>& importLIDs,
                      const ArrayView<const Scalar>& imports,
                      const ArrayView<size_t>& numPacketsPerLID,
                      size_t constantNumPackets,
                      Distributor& distor,
                      CombineMode CM);
#endif

    void createViews () const;
    void createViewsNonConst (KokkosClassic::ReadWriteOption rwo);
    void releaseViews () const;

    //! Nonconst host view created in createViewsNonConst().
    mutable ArrayRCP<Scalar> ncview_;

    //! Const host view created in createViews().
    mutable ArrayRCP<const Scalar> cview_;

    //@}

#if TPETRA_USE_KOKKOS_DISTOBJECT

    Kokkos::View<const Scalar*, device_type, Kokkos::MemoryUnmanaged>
    getKokkosView() const {
      Teuchos::ArrayRCP<const Scalar> buff = MVT::getValues (lclMV_);
      Kokkos::View<const Scalar*, device_type, Kokkos::MemoryUnmanaged> v(
        buff.getRawPtr(), buff.size());
      return v;
    }
    Kokkos::View<Scalar*, device_type, Kokkos::MemoryUnmanaged>
    getKokkosViewNonConst() {
      Teuchos::ArrayRCP<Scalar> buff = MVT::getValuesNonConst (lclMV_);
      Kokkos::View<Scalar*, device_type, Kokkos::MemoryUnmanaged> v(
        buff.getRawPtr(), buff.size());
      return v;
    }

#endif
  };

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
  deep_copy (MultiVector<DS,DL,DG,DN>& dst,
             const MultiVector<SS,SL,SG,SN>& src);
  // {
  //   TEUCHOS_TEST_FOR_EXCEPTION(
  //     true, std::logic_error, "The fully generic version of Tpetra::deep_copy "
  //     "is not implemented.");
  // }

  template <class SS, class SL, class SG, class SN>
  void
  deep_copy (MultiVector<SS,SL,SG,SN>& dst, const MultiVector<SS,SL,SG,SN>& src)
  {
    // NOTE (mfh 11 Sep 2014) We can't implement deep_copy with
    // shallow-copy operator=, because that would invalidate existing
    // views of dst!
    dst.assign (src);
  }

  /// \brief Return a deep copy of the MultiVector \c src.
  /// \relatesalso MultiVector
  ///
  /// Regarding Copy or View semantics: The returned MultiVector is
  /// always a deep copy of \c src, and <i>always</i> has view
  /// semantics.  This is because createCopy returns by value, and
  /// therefore it assumes that you want to pass MultiVector objects
  /// around by value, not by Teuchos::RCP.
  ///
  /// In the Kokkos refactor version of Tpetra, MultiVector always has
  /// View semantics.  However, the above remarks still apply.
  template <class ST, class LO, class GO, class NT>
  MultiVector<ST, LO, GO, NT>
  createCopy (const MultiVector<ST, LO, GO, NT>& src);

  namespace Details {
    /// \brief Implementation of ::Tpetra::MultiVector::clone().
    ///
    /// \tparam DstMultiVecType Specialization of
    ///   ::Tpetra::MultiVector, which is the result of (is returned
    ///   by) the clone() operation.
    ///
    /// \tparam SrcMultiVecType Specialization of
    ///   ::Tpetra::MultiVector, which is the source (input) of the
    ///   clone() operation.
    ///
    /// We provide partial specializations for the following cases:
    /// <ol>
    /// <li> Source and destination MultiVector types have the same
    ///      Scalar type, but all their other template parameters
    ///      might be different. </li>
    /// <li> Source and destination MultiVector types are the
    ///      same. </li>
    /// <li> Source and destination MultiVector types are both Kokkos
    ///      refactor types (we look at their Node types to determine
    ///      this), and have the same Scalar types, but all their
    ///      other template parameters might be different. </li>
    /// <li> Source and destination MultiVector types are both Kokkos
    ///      refactor types (we look at their Node types to determine
    ///      this), and both the same. </li>
    /// </ol>
    template<class DstMultiVecType, class SrcMultiVecType>
    struct MultiVectorCloner {
      typedef DstMultiVecType dst_mv_type;
      typedef SrcMultiVecType src_mv_type;

      static Teuchos::RCP<dst_mv_type>
      clone (const src_mv_type& X,
             const Teuchos::RCP<typename dst_mv_type::node_type>& node2);
    };

    // Partial specialization of MultiVectorCloner for when the source
    // and destination MultiVector types have the same Scalar type,
    // but all their other template parameters might be different.
    template<class ScalarType,
             class DstLocalOrdinalType, class DstGlobalOrdinalType, class DstNodeType,
             class SrcLocalOrdinalType, class SrcGlobalOrdinalType, class SrcNodeType>
    struct MultiVectorCloner< ::Tpetra::MultiVector<ScalarType, DstLocalOrdinalType, DstGlobalOrdinalType, DstNodeType>,
                              ::Tpetra::MultiVector<ScalarType, SrcLocalOrdinalType, SrcGlobalOrdinalType, SrcNodeType> >
    {
      typedef ::Tpetra::MultiVector<ScalarType, DstLocalOrdinalType, DstGlobalOrdinalType, DstNodeType> dst_mv_type;
      typedef ::Tpetra::MultiVector<ScalarType, SrcLocalOrdinalType, SrcGlobalOrdinalType, SrcNodeType> src_mv_type;

      static Teuchos::RCP<dst_mv_type>
      clone (const src_mv_type& X,
             const Teuchos::RCP<typename src_mv_type::node_type>& node2)
      {
        using Teuchos::ArrayRCP;
        using Teuchos::RCP;
        using Teuchos::rcp;
        typedef typename src_mv_type::map_type src_map_type;
        typedef typename dst_mv_type::map_type dst_map_type;
        typedef typename dst_mv_type::node_type dst_node_type;

        // Clone X's Map to have the new Node type.
        RCP<const src_map_type> map1 = X.getMap ();
        RCP<const dst_map_type> map2 = map1.is_null () ?
          Teuchos::null : map1->template clone<dst_node_type> (node2);

        const size_t lclNumRows = X.getLocalLength ();
        const size_t numCols = X.getNumVectors ();
        const size_t LDA = lclNumRows;

        // Get a host deep copy of X's data.
        ArrayRCP<typename src_mv_type::scalar_type> data1 (LDA * numCols);
        X.get1dCopy (data1 (), LDA);

        // Create the destination MultiVector.  This might do another
        // deep copy; I'm not really worried about that.  clone()
        // doesn't have to be super fast; it just can't be too slow.
        return rcp (new dst_mv_type (map2, data1 (), LDA, numCols));
      }
    };

    // Partial specialization of MultiVectorCloner for when the source
    // and destination MultiVector types are the same.
    template<class ScalarType, class LocalOrdinalType, class GlobalOrdinalType, class NodeType>
    struct MultiVectorCloner< ::Tpetra::MultiVector<ScalarType, LocalOrdinalType, GlobalOrdinalType, NodeType>,
                              ::Tpetra::MultiVector<ScalarType, LocalOrdinalType, GlobalOrdinalType, NodeType> >
    {
      typedef ::Tpetra::MultiVector<ScalarType, LocalOrdinalType, GlobalOrdinalType, NodeType> dst_mv_type;
      typedef dst_mv_type src_mv_type;

      static Teuchos::RCP<dst_mv_type>
      clone (const src_mv_type& X, const Teuchos::RCP<NodeType>& )
      {
        // Create a deep copy.
        RCP<dst_mv_type> X_clone = rcp (new dst_mv_type (X, Teuchos::Copy));
        // Set the cloned MultiVector to have the same copy-or-view
        // semantics as the input MultiVector X.
        X_clone->setCopyOrView (X.getCopyOrView ());

        return X_clone;
      }
    };

    template<class MultiVectorType>
    Teuchos::RCP<MultiVectorType>
    CreateMultiVectorFromView<MultiVectorType>::
    create (const Teuchos::RCP<const map_type>& map,
            const Teuchos::ArrayRCP<scalar_type>& view,
            const size_t LDA,
            const size_t numVectors)
    {
      using Teuchos::rcp;
      typedef Tpetra::details::ViewAccepter<node_type> VAN;

      // This uses a protected MultiVector constructor, but this
      // nonmember function was declared a friend of MultiVector.
      //
      // The ViewAccepter expression will fail to compile for
      // unsupported Kokkos Node types.
      return rcp (new MultiVectorType (map, VAN::template acceptView<scalar_type> (view),
                                       LDA, numVectors, HOST_VIEW_CONSTRUCTOR));
    }
  } // namespace Details


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  template <class Node2>
  Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node2> >
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  clone (const Teuchos::RCP<Node2>& node2) const
  {
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node2> dst_mv_type;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> src_mv_type;
    typedef Details::MultiVectorCloner<dst_mv_type, src_mv_type> cloner_type;
    return cloner_type::clone (*this, node2);
  }

  /// \brief Nonmember MultiVector constructor: make a MultiVector from a given Map.
  /// \relatesalso MultiVector
  /// \relatesalso Vector
  ///
  /// \param map [in] Map describing the distribution of rows of the
  ///   resulting MultiVector.
  /// \param numVectors [in] Number of columns of the resulting
  ///   MultiVector.
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  createMultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
                     const size_t numVectors)
  {
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
    return Teuchos::rcp (new MV (map, numVectors));
  }

} // namespace Tpetra

// Include KokkosRefactor partial specialization if enabled
#if defined(TPETRA_HAVE_KOKKOS_REFACTOR)
#include "Tpetra_KokkosRefactor_MultiVector_decl.hpp"
#endif

// Define createMultiVectorFromView after the above include, so that
// the function can pick up the partial specialization of
// CreateMultiVectorFromView.

namespace Tpetra {

  /// \brief Nonmember MultiVector constructor with view semantics using user-allocated data.
  /// \relatesalso MultiVector
  /// \relatesalso Vector
  ///
  /// \warning This function is not supported for all Kokkos Node
  ///   types.  Specifically, it is not typically supported for
  ///   GPU accelerator-based nodes like KokkosClassic::ThrustGPUNode.
  ///
  /// \param map [in] The Map describing the distribution of rows of
  ///   the multivector.
  /// \param view [in/out] A pointer to column-major dense matrix
  ///   data.  This will be the multivector's data on the calling
  ///   process.  The multivector will use the pointer directly,
  ///   without copying.
  /// \param LDA [in] The leading dimension (a.k.a. "stride") of the
  ///   column-major input data.
  /// \param numVectors [in] The number of columns in the input data.
  ///   This will be the number of vectors in the returned
  ///   multivector.
  ///
  /// \node To Kokkos and Tpetra developers: If you add a new Kokkos
  ///   Node type that is a host Node type (where memory lives in user
  ///   space, not in a different space as on a GPU), you will need to
  ///   add a specialization of Tpetra::details::ViewAccepter for your
  ///   new Node type.
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  TPETRA_DEPRECATED
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  createMultiVectorFromView (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
                             const Teuchos::ArrayRCP<Scalar>& view,
                             const size_t LDA,
                             const size_t numVectors)
  {
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> mv_type;
    typedef Details::CreateMultiVectorFromView<mv_type> impl_type;
    return impl_type::create (map, view, LDA, numVectors);
  }

} // namespace Tpetra

#endif // TPETRA_MULTIVECTOR_DECL_HPP
