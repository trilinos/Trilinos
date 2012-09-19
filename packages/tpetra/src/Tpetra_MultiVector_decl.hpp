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

#include <Teuchos_LabeledObject.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_BLAS_types.hpp>
#include <Teuchos_Range1D.hpp>

#include <Kokkos_MultiVector.hpp>
#include <Kokkos_DefaultArithmetic.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DistObject.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_ViewAccepter.hpp"

// TODO: add principal use case instructions for memory management interfaces (view/copy extraction)
// TODO: expand user-visible documentation

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration of Vector, needed to prevent circular inclusions
  template<class S, class LO, class GO, class N> class Vector;

  //template<class S, class LO, class GO, class N> class MultiVector;

  //template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  //RCP< MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  //createMultiVectorFromView(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map,
  //                          const ArrayRCP<Scalar> &view, size_t LDA, size_t numVectors);
#endif

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
  /// \tparam Scalar The type of the numerical entries of the vector(s).
  ///  (You can use real-valued or complex-valued types here, unlike in
  ///  Epetra, where the scalar type is always \c double.)
  ///
  /// \tparam LocalOrdinal The type of local indices.  Same as the \c
  ///   LocalOrdinal template parameter of \c Map objects used by this
  ///   matrix.  (In Epetra, this is just \c int.)  The default type is
  ///   \c int, which should suffice for most users.  This type must be
  ///   big enough to store the local (per process) number of rows.
  ///
  /// \tparam GlobalOrdinal The type of global indices.  Same as the
  ///   \c GlobalOrdinal template parameter of \c Map objects used by
  ///   this matrix.  (In Epetra, this is just \c int.  One advantage
  ///   of Tpetra over Epetra is that you can use a 64-bit integer
  ///   type here if you want to solve big problems.)  The default
  ///   type is <tt>LocalOrdinal</tt>.  This type must be big enough
  ///   to store the global (over all processes in the communicator)
  ///   number of rows or columns.
  ///
  /// \tparam Node A class implementing on-node shared-memory parallel
  ///   operations.  It must implement the
  ///   \ref kokkos_node_api "Kokkos Node API."
  ///   The default \c Node type should suffice for most users.
  ///   The actual default type depends on your Trilinos build options.
  ///
  /// \note If you use the default \c GlobalOrdinal type, which is
  ///   <tt>int</tt>, then the <i>global</i> number of rows or columns
  ///   in the matrix may be no more than \c INT_MAX, which for
  ///   typical 32-bit \c int is \f$2^{31} - 1\f$ (about two billion).
  ///   If you want to solve larger problems, you must use a 64-bit
  ///   integer type here.
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
  /// provide is Kokkos::ThrustGPUNode.  All other types are CPU
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
  /// container data structure, Kokkos::MultiVector, and then use the
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
  template<class Scalar,
           class LocalOrdinal=int,
           class GlobalOrdinal=LocalOrdinal,
           class Node=Kokkos::DefaultNode::DefaultNodeType>
  class MultiVector :
    public DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node>
  {
  public:
    //! @name Typedefs to facilitate template metaprogramming.
    //@{

    //! The type of entries in the vector(s).
    typedef Scalar        scalar_type;
    //! The type of local indices.
    typedef LocalOrdinal  local_ordinal_type;
    //! The type of global indices.
    typedef GlobalOrdinal global_ordinal_type;
    //! The Kokkos Node type.
    typedef Node          node_type;

    //@}
    //! @name Constructors and destructor
    //@{

    /// \brief Basic constuctor.
    ///
    /// \param map [in] Map describing the distribution of rows.
    /// \param NumVectors [in] Number of vectors (columns).
    /// \param zeroOut [in] Whether to initialize all the entries of
    ///   the MultiVector to zero.
    MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
                 size_t NumVectors,
                 bool zeroOut=true);

    //! Copy constructor (performs a deep copy).
    MultiVector (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &source);

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

    /// \brief Replace the underlying Map with a compatible one.
    ///
    /// This method relabels the rows of the multivector using the
    /// global IDs in the input Map.  Thus, it implicitly applies a
    /// permutation, without actually moving data.  This only works if
    /// the input Map is compatible (in the sense of \c
    /// Map::isCompatible()) with the multivector's current Map, so
    /// that the number of rows per process does not change.
    ///
    /// We only check for compatibility in debug mode (when Trilinos
    /// was built with the Trilinos_ENABLE_DEBUG option set).  In that
    /// case, if the input Map is <i>not</i> compatible, then this
    /// method throws \c std::invalid_argument.  We only check in
    /// debug mode because the check requires communication
    /// (\f$O(1)\f$ all-reduces).
    ///
    /// \note This method is <i>not</i> for arbitrary data
    ///   redistribution.  If you need to move data around, use \c
    ///   Import or \c Export.
    ///
    /// \note This method must always be called as a collective
    ///   operation on all processes over which the multivector is
    ///   distributed.  This is because the method reserves the right
    ///   to check for compatibility of the two Maps, at least in
    ///   debug mode.
    void replaceMap(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map);

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

    //! \brief Return a const MultiVector view of a subset of rows.
    /**
        Return a const view of this MultiVector consisting of a subset
        of the rows, as specified by an offset and a subset Map of
        this MultiVector's current row Map.

        \param In subMap - The row Map for the new MultiVector.
        \param In offset - The offset into the data of <tt>(*this)</tt>.

        \pre <tt>subMap->getNodeNumElements() + offset < this->getLocalLength()</tt>
        \pre The given Map must be a subset Map of this MultiVector's row Map.
     */
    Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
    offsetView (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& subMap,
                size_t offset) const;

    //! \brief Return a non-const MultiVector view of a subset of rows.
    /**
        Returns a non-const view of this MultiVector consisting of a
        subset of the rows, as specified by an offset and a subset Map
        of this MultiVector's current row Map.

        \param In subMap - The row Map for the new MultiVector.
        \param In offset - The offset into the data of <tt>(*this)</tt>.

        \pre <tt>subMap->getNodeNumElements() + offset < this->getLocalLength()</tt>
        \pre The given Map must be a subset Map of this MultiVector's row Map.
     */
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

    //! Return a const reference to the underlying Kokkos::MultiVector object (advanced use only)
    const Kokkos::MultiVector<Scalar,Node> & getLocalMV() const;

    //! Return a non-const reference to the underlying Kokkos::MultiVector object (advanced use only)
    Kokkos::MultiVector<Scalar,Node> & getLocalMVNonConst();

    //@}

    //! @name Mathematical methods
    //@{

    //! Compute dot product of each corresponding pair of vectors, dots[i] = this[i].dot(A[i])
    void dot(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const Teuchos::ArrayView<Scalar> &dots) const;

    //! Put element-wise absolute values of input Multi-vector in target: A = abs(this)
    void abs(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A);

    //! Put element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
    void reciprocal(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A);

    //! Scale the current values of a multi-vector, this = alpha*this.
    void scale(const Scalar &alpha);

    //! Scale the current values of a multi-vector, this[j] = alpha[j]*this[j].
    void scale(Teuchos::ArrayView<const Scalar> alpha);

    //! Replace multi-vector values with scaled values of A, this = alpha*A.
    void scale(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A);

    //! Update multi-vector values with scaled values of A, this = beta*this + alpha*A.
    void update(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const Scalar &beta);

    //! Update multi-vector with scaled values of A and B, this = gamma*this + alpha*A + beta*B.
    void update(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const Scalar &beta, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B, const Scalar &gamma);

    //! Compute 1-norm of each vector in multi-vector.
    void norm1(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const;

    //! Compute 2-norm of each vector in multi-vector.
    void norm2(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const;

    //! Compute Inf-norm of each vector in multi-vector.
    void normInf(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const;

    //! Compute Weighted 2-norm (RMS Norm) of each vector in multi-vector.
    void normWeighted(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &weights, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const;

    //! Compute mean (average) value of each vector in multi-vector.
    void meanValue(const Teuchos::ArrayView<Scalar> &means) const;

    //! Matrix-matrix multiplication: this = beta*this + alpha*op(A)*op(B).
    void
    multiply (Teuchos::ETransp transA,
              Teuchos::ETransp transB,
              const Scalar& alpha,
              const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A,
              const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& B,
              const Scalar& beta);

    //! Element-wise multiply of a Vector A with a MultiVector B.
    /** Forms this = scalarThis * this + scalarAB * B @ A
     *  where @ denotes element-wise multiplication.
     *  B must be the same shape (size and num-vectors) as this, while
     *  A is the same size but a single vector (column).
     *
     *  this = scalarThis * this(i,j) + scalarAB * B(i,j) * A(i,1) (A has only 1 column)
     *
     */
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
    std::string description() const;

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
    void
    describe (Teuchos::FancyOStream& out,
              const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

    //@}

  protected:

    typedef Kokkos::MultiVector<Scalar,Node> KMV;
    typedef Kokkos::DefaultArithmetic<KMV>   MVT;

    //! The Kokkos::MultiVector containing the compute buffer of data.
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

    //! @name View constructors, used only by nonmember constructors.
    //@{

    template <class S,class LO,class GO,class N>
    friend RCP<MultiVector<S,LO,GO,N> >
    createMultiVectorFromView (const Teuchos::RCP<const Map<LO,GO,N> >&, const Teuchos::ArrayRCP<S>&, size_t, size_t);

    /// \brief View constructor with user-allocated data, for CPU nodes only.
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

    inline bool vectorIndexOutOfRange(size_t VectorIndex) const {
      return (VectorIndex < 1 && VectorIndex != 0) || VectorIndex >= getNumVectors();
    }

    /// \fn getSubArrayRCP
    /// \brief Persisting view of j-th column in the given ArrayRCP.
    ///
    /// This method considers isConstantStride().  The ArrayRCP may
    /// correspond either to a compute buffer or a host view.
    template <class T>
    ArrayRCP<T> getSubArrayRCP(ArrayRCP<T> arr, size_t j) const;

    //! Advanced constructor for non-contiguous views.
    MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
                 Teuchos::ArrayRCP<Scalar> data,
                 size_t LDA,
                 Teuchos::ArrayView<const size_t> whichVectors,
                 EPrivateComputeViewConstructor /* dummy */);

    //! Advanced constructor for contiguous views.
    MultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
                 Teuchos::ArrayRCP<Scalar> data,
                 size_t LDA,
                 size_t NumVectors,
                 EPrivateComputeViewConstructor /* dummy */);

    //@}
    //! @name Implementation of Tpetra::DistObject
    //@{

    /// \brief Whether data redistribution between \c sourceObj and this object is legal.
    ///
    /// This method is called in \c DistObject::doTransfer() to check
    /// whether data redistribution between the two objects is legal.
    bool
    checkSizes (const DistObject<Scalar,LocalOrdinal,GlobalOrdinal,Node>& sourceObj);

    void
    copyAndPermute (const DistObject<Scalar,LocalOrdinal,GlobalOrdinal,Node>& sourceObj,
                    size_t numSameIDs,
                    const ArrayView<const LocalOrdinal>& permuteToLIDs,
                    const ArrayView<const LocalOrdinal>& permuteFromLIDs);

    void
    packAndPrepare (const DistObject<Scalar,LocalOrdinal,GlobalOrdinal,Node>& sourceObj,
                    const ArrayView<const LocalOrdinal>& exportLIDs,
                    Array<Scalar>& exports,
                    const ArrayView<size_t>& numExportPacketsPerLID,
                    size_t& constantNumPackets,
                    Distributor& distor);

    void
    unpackAndCombine (const ArrayView<const LocalOrdinal>& importLIDs,
                      const ArrayView<const Scalar>& imports,
                      const ArrayView<size_t>& numPacketsPerLID,
                      size_t constantNumPackets,
                      Distributor& distor,
                      CombineMode CM);

    void createViews () const;
    void createViewsNonConst (Kokkos::ReadWriteOption rwo);
    void releaseViews () const;

    //! Nonconst host view created in createViewsNonConst().
    mutable ArrayRCP<Scalar> ncview_;

    //! Const host view created in createViews().
    mutable ArrayRCP<const Scalar> cview_;

    //@}
    //! @name Instance data for tracking efficiency warnings
    //@{

    //! Whether releaseViews() has yet raised an efficiency warning.
    mutable bool releaseViewsRaisedEfficiencyWarning_;
    //! Whether createViews() has yet raised an efficiency warning.
    mutable bool createViewsRaisedEfficiencyWarning_;
    //! Whether createViewsNonConst() has yet raised an efficiency warning.
    mutable bool createViewsNonConstRaisedEfficiencyWarning_;

    //@}
  }; // class MultiVector

  /// \brief Nonmember MultiVector constructor: make a MultiVector from a given Map.
  /// \relatesalso MultiVector
  ///
  /// \param map [in] Map describing the distribution of rows of the
  ///   resulting MultiVector.
  /// \param numVectors [in] Number of columns of the resulting
  ///   MultiVector.
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  createMultiVector (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
                     size_t numVectors)
  {
    using Teuchos::rcp;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

    const bool initToZero = true;
    return rcp (new MV (map, numVectors, initToZero));
  }

  /// \brief Nonmember MultiVector constructor with view semantics using user-allocated data.
  /// \relatesalso MultiVector
  /// \relatesalso Vector
  ///
  /// \warning This function is not supported for all Kokkos Node
  ///   types.  Specifically, it is not typically supported for
  ///   GPU accelerator-based nodes like Kokkos::ThrustGPUNode.
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
  Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  createMultiVectorFromView (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
                             const Teuchos::ArrayRCP<Scalar>& view,
                             size_t LDA,
                             size_t numVectors)
  {
    using Teuchos::rcp;
    typedef Tpetra::details::ViewAccepter<Node> VAN;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

    // This uses a protected MultiVector constructor, but this
    // nonmember function was declared a friend of MultiVector.
    //
    // The ViewAccepter expression will fail to compile for
    // unsupported Kokkos Node types.
    return rcp (new MV (map, VAN::template acceptView<Scalar> (view),
                        LDA, numVectors, HOST_VIEW_CONSTRUCTOR));
  }

} // namespace Tpetra


#endif // TPETRA_MULTIVECTOR_DECL_HPP
