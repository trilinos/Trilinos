// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// clang-format off
#ifndef TPETRA_CRSMATRIX_DECL_HPP
#define TPETRA_CRSMATRIX_DECL_HPP

/// \file Tpetra_CrsMatrix_decl.hpp
/// \brief Declaration of the Tpetra::CrsMatrix class

#include "Tpetra_CrsMatrix_fwd.hpp"
#include "TpetraExt_MatrixMatrix_fwd.hpp"
#include "KokkosSparse_Utils.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#if KOKKOSKERNELS_VERSION >= 40299
#include "Tpetra_Details_MatrixApplyHelper.hpp"
#else
#include "Tpetra_LocalCrsMatrixOperator.hpp"
#endif
#include "Tpetra_RowMatrix_decl.hpp"
#include "Tpetra_Exceptions.hpp"
#include "Tpetra_DistObject.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Details_PackTraits.hpp" // unused here, could delete
#include "Tpetra_Details_ExecutionSpacesUser.hpp"
#include "Teuchos_DataAccess.hpp"


#include <memory> // std::shared_ptr

namespace Tpetra {

  // Forward declaration for CrsMatrix::swap() test
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node> class crsMatrix_Swap_Tester;

  /// \brief Nonmember CrsMatrix constructor that fuses Import and fillComplete().
  /// \relatesalso CrsMatrix
  /// \tparam CrsMatrixType A specialization of CrsMatrix.
  ///
  /// A common use case is to create an empty destination CrsMatrix,
  /// redistribute from a source CrsMatrix (by an Import or Export
  /// operation), then call fillComplete() on the destination
  /// CrsMatrix.  This constructor fuses these three cases, for an
  /// Import redistribution.
  ///
  /// Fusing redistribution and fillComplete() exposes potential
  /// optimizations.  For example, it may make constructing the column
  /// Map faster, and it may avoid intermediate unoptimized storage in
  /// the destination CrsMatrix.  These optimizations may improve
  /// performance for specialized kernels like sparse matrix-matrix
  /// multiply, as well as for redistributing data after doing load
  /// balancing.
  ///
  /// The resulting matrix is fill complete (in the sense of
  /// isFillComplete()) and has optimized storage (in the sense of
  /// isStorageOptimized()).  By default, its domain Map is the domain
  /// Map of the source matrix, and its range Map is the range Map of
  /// the source matrix.
  ///
  /// \warning If the target Map of the Import is a subset of the
  ///   source Map of the Import, then you cannot use the default
  ///   range Map.  You should instead construct a nonoverlapping
  ///   version of the target Map and supply that as the nondefault
  ///   value of the range Map.
  ///
  /// \param sourceMatrix [in] The source matrix from which to
  ///   import.  The source of an Import must have a nonoverlapping
  ///   distribution.
  ///
  /// \param importer [in] The Import instance containing a
  ///   precomputed redistribution plan.  The source Map of the
  ///   Import must be the same as the rowMap of sourceMatrix unless
  ///   the "Reverse Mode" option on the params list, in which case
  ///   the targetMap of Import must match the rowMap of the sourceMatrix
  ///
  /// \param domainMap [in] Domain Map of the returned matrix.  If
  ///   null, we use the default, which is the domain Map of the
  ///   source matrix.
  ///
  /// \param rangeMap [in] Range Map of the returned matrix.  If
  ///   null, we use the default, which is the range Map of the
  ///   source matrix.
  ///
  /// \param params [in/out] Optional list of parameters.  If not
  ///   null, any missing parameters will be filled in with their
  ///   default values.
  template<class CrsMatrixType>
  Teuchos::RCP<CrsMatrixType>
  importAndFillCompleteCrsMatrix (const Teuchos::RCP<const CrsMatrixType>& sourceMatrix,
                                  const Import<typename CrsMatrixType::local_ordinal_type,
                                               typename CrsMatrixType::global_ordinal_type,
                                               typename CrsMatrixType::node_type>& importer,
                                  const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                               typename CrsMatrixType::global_ordinal_type,
                                                               typename CrsMatrixType::node_type> >& domainMap = Teuchos::null,
                                  const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                               typename CrsMatrixType::global_ordinal_type,
                                                               typename CrsMatrixType::node_type> >& rangeMap = Teuchos::null,
                                  const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

  /// \brief Nonmember CrsMatrix constructor that fuses Import and fillComplete().
  /// \relatesalso CrsMatrix
  /// \tparam CrsMatrixType A specialization of CrsMatrix.
  ///
  /// A common use case is to create an empty destination CrsMatrix,
  /// redistribute from a source CrsMatrix (by an Import or Export
  /// operation), then call fillComplete() on the destination
  /// CrsMatrix.  This constructor fuses these three cases, for an
  /// Import redistribution.
  ///
  /// Fusing redistribution and fillComplete() exposes potential
  /// optimizations.  For example, it may make constructing the column
  /// Map faster, and it may avoid intermediate unoptimized storage in
  /// the destination CrsMatrix.  These optimizations may improve
  /// performance for specialized kernels like sparse matrix-matrix
  /// multiply, as well as for redistributing data after doing load
  /// balancing.
  ///
  /// The resulting matrix is fill complete (in the sense of
  /// isFillComplete()) and has optimized storage (in the sense of
  /// isStorageOptimized()).  By default, its domain Map is the domain
  /// Map of the source matrix, and its range Map is the range Map of
  /// the source matrix.
  ///
  /// \warning If the target Map of the Import is a subset of the
  ///   source Map of the Import, then you cannot use the default
  ///   range Map.  You should instead construct a nonoverlapping
  ///   version of the target Map and supply that as the nondefault
  ///   value of the range Map.
  ///
  /// \param sourceMatrix [in] The source matrix from which to
  ///   import.  The source of an Import must have a nonoverlapping
  ///   distribution.
  ///
  /// \param rowImporter [in] The Import instance containing a
  ///   precomputed redistribution plan.  The source Map of the
  ///   Import must be the same as the rowMap of sourceMatrix unless
  ///   the "Reverse Mode" option on the params list, in which case
  ///   the targetMap of Import must match the rowMap of the sourceMatrix
  ///
  /// \param domainImporter [in] The Import instance containing a
  ///   precomputed redistribution plan.  The source Map of the
  ///   Import must be the same as the domainMap of sourceMatrix unless
  ///   the "Reverse Mode" option on the params list, in which case
  ///   the targetMap of Import must match the domainMap of the sourceMatrix
  ///
  /// \param domainMap [in] Domain Map of the returned matrix.
  ///
  /// \param rangeMap [in] Range Map of the returned matrix.
  ///
  /// \param params [in/out] Optional list of parameters.  If not
  ///   null, any missing parameters will be filled in with their
  ///   default values.
  template<class CrsMatrixType>
  Teuchos::RCP<CrsMatrixType>
  importAndFillCompleteCrsMatrix (const Teuchos::RCP<const CrsMatrixType>& sourceMatrix,
                                  const Import<typename CrsMatrixType::local_ordinal_type,
                                               typename CrsMatrixType::global_ordinal_type,
                                               typename CrsMatrixType::node_type>& rowImporter,
                                  const Import<typename CrsMatrixType::local_ordinal_type,
                                              typename CrsMatrixType::global_ordinal_type,
                                              typename CrsMatrixType::node_type>& domainImporter,
                                  const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                               typename CrsMatrixType::global_ordinal_type,
                                                               typename CrsMatrixType::node_type> >& domainMap,
                                  const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                               typename CrsMatrixType::global_ordinal_type,
                                                               typename CrsMatrixType::node_type> >& rangeMap,
                                  const Teuchos::RCP<Teuchos::ParameterList>& params);

  /// \brief Nonmember CrsMatrix constructor that fuses Export and fillComplete().
  /// \relatesalso CrsMatrix
  /// \tparam CrsMatrixType A specialization of CrsMatrix.
  ///
  /// For justification, see the documentation of
  /// importAndFillCompleteCrsMatrix() (which is the Import analog of
  /// this function).
  ///
  /// The resulting matrix is fill complete (in the sense of
  /// isFillComplete()) and has optimized storage (in the sense of
  /// isStorageOptimized()).  By default, its domain Map is the domain
  /// Map of the source matrix, and its range Map is the range Map of
  /// the source matrix.
  ///
  /// \param sourceMatrix [in] The source matrix from which to
  ///   export.  Its row Map may be overlapping, since the source of
  ///   an Export may be overlapping.
  ///
  /// \param exporter [in] The Export instance containing a
  ///   precomputed redistribution plan.  The source Map of the
  ///   Export must be the same as the row Map of sourceMatrix.
  ///
  /// \param domainMap [in] Domain Map of the returned matrix.  If
  ///   null, we use the default, which is the domain Map of the
  ///   source matrix.
  ///
  /// \param rangeMap [in] Range Map of the returned matrix.  If
  ///   null, we use the default, which is the range Map of the
  ///   source matrix.
  ///
  /// \param params [in/out] Optional list of parameters.  If not
  ///   null, any missing parameters will be filled in with their
  ///   default values.
  template<class CrsMatrixType>
  Teuchos::RCP<CrsMatrixType>
  exportAndFillCompleteCrsMatrix (const Teuchos::RCP<const CrsMatrixType>& sourceMatrix,
                                  const Export<typename CrsMatrixType::local_ordinal_type,
                                               typename CrsMatrixType::global_ordinal_type,
                                               typename CrsMatrixType::node_type>& exporter,
                                  const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                               typename CrsMatrixType::global_ordinal_type,
                                                               typename CrsMatrixType::node_type> >& domainMap = Teuchos::null,
                                  const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                               typename CrsMatrixType::global_ordinal_type,
                                                               typename CrsMatrixType::node_type> >& rangeMap = Teuchos::null,
                                  const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

  /// \brief Nonmember CrsMatrix constructor that fuses Export and fillComplete().
  /// \relatesalso CrsMatrix
  /// \tparam CrsMatrixType A specialization of CrsMatrix.
  ///
  /// For justification, see the documentation of
  /// importAndFillCompleteCrsMatrix() (which is the Import analog of
  /// this function).
  ///
  /// The resulting matrix is fill complete (in the sense of
  /// isFillComplete()) and has optimized storage (in the sense of
  /// isStorageOptimized()).  By default, its domain Map is the domain
  /// Map of the source matrix, and its range Map is the range Map of
  /// the source matrix.
  ///
  /// \param sourceMatrix [in] The source matrix from which to
  ///   export.  Its row Map may be overlapping, since the source of
  ///   an Export may be overlapping.
  ///
  /// \param rowExporter [in] The Export instance containing a
  ///   precomputed redistribution plan.  The source Map of the
  ///   Export must be the same as the row Map of sourceMatrix.
  ///
  /// \param domainExporter [in] The Export instance containing a
  ///   precomputed redistribution plan.  The source Map of the
  ///   Export must be the same as the domain Map of sourceMatrix.
  ///
  /// \param domainMap [in] Domain Map of the returned matrix.
  ///
  /// \param rangeMap [in] Range Map of the returned matrix.
  ///
  /// \param params [in/out] Optional list of parameters.  If not
  ///   null, any missing parameters will be filled in with their
  ///   default values.
  template<class CrsMatrixType>
  Teuchos::RCP<CrsMatrixType>
  exportAndFillCompleteCrsMatrix (const Teuchos::RCP<const CrsMatrixType>& sourceMatrix,
                                  const Export<typename CrsMatrixType::local_ordinal_type,
                                               typename CrsMatrixType::global_ordinal_type,
                                               typename CrsMatrixType::node_type>& rowExporter,
                                  const Export<typename CrsMatrixType::local_ordinal_type,
                                               typename CrsMatrixType::global_ordinal_type,
                                               typename CrsMatrixType::node_type>& domainExporter,
                                  const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                               typename CrsMatrixType::global_ordinal_type,
                                                               typename CrsMatrixType::node_type> >& domainMap,
                                  const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                               typename CrsMatrixType::global_ordinal_type,
                                                               typename CrsMatrixType::node_type> >& rangeMap,
                                  const Teuchos::RCP<Teuchos::ParameterList>& params);

  /// \brief Nonmember function that computes a residual
  /// Computes R = B - A * X
  namespace Details {
    template<class SC, class LO, class GO, class NO>
    void residual(const Operator<SC,LO,GO,NO> &   A,
                  const MultiVector<SC,LO,GO,NO> & X,
                  const MultiVector<SC,LO,GO,NO> & B,
                  MultiVector<SC,LO,GO,NO> & R);
  }

  /// \class CrsMatrix
  /// \brief Sparse matrix that presents a row-oriented interface that
  ///   lets users read or modify entries.
  ///
  /// \tparam Scalar The type of the numerical entries of the matrix.
  ///   (You can use real-valued or complex-valued types here, unlike
  ///   in Epetra, where the scalar type is always \c double.)
  /// \tparam LocalOrdinal The type of local indices.  See the
  ///   documentation of Map for requirements.
  /// \tparam GlobalOrdinal The type of global indices.  See the
  ///   documentation of Map for requirements.
  /// \tparam Node The Kokkos Node type.  See the documentation of Map
  ///   for requirements.
  ///
  /// This class implements a distributed-memory parallel sparse matrix,
  /// and provides sparse matrix-vector multiply (including transpose)
  /// and sparse triangular solve operations.  It provides access by rows
  /// to the elements of the matrix, as if the local data were stored in
  /// compressed sparse row format.  (Implementations are <i>not</i>
  /// required to store the data in this way internally.)  This class has
  /// an interface like that of Epetra_CrsMatrix, but also allows
  /// insertion of data into nonowned rows, much like Epetra_FECrsMatrix.
  ///
  /// \section Tpetra_CrsMatrix_prereq Prerequisites
  ///
  /// Before reading the rest of this documentation, it helps to know
  /// something about the Teuchos memory management classes, in
  /// particular Teuchos::RCP, Teuchos::ArrayRCP, and Teuchos::ArrayView.
  /// You should also know a little bit about MPI (the Message Passing
  /// Interface for distributed-memory programming).  You won't have to
  /// use MPI directly to use CrsMatrix, but it helps to be familiar with
  /// the general idea of distributed storage of data over a
  /// communicator.  Finally, you should read the documentation of Map
  /// and MultiVector.
  ///
  /// \section Tpetra_CrsMatrix_local_vs_global Local and global indices
  ///
  /// The distinction between local and global indices might confuse new
  /// Tpetra users.  Please refer to the documentation of Map for a
  /// detailed explanation.  This is important because many of
  /// CrsMatrix's methods for adding, modifying, or accessing entries
  /// come in versions that take either local or global indices.  The
  /// matrix itself may store indices either as local or global, and the
  /// same matrix may use global indices or local indices at different
  /// points in its life.  You should only use the method version
  /// corresponding to the current state of the matrix.  For example,
  /// getGlobalRowView() returns a view to the indices represented as
  /// global; it is incorrect to call this method if the matrix is
  /// storing indices as local.  Call isGloballyIndexed() or
  /// isLocallyIndexed() to find out whether the matrix currently stores
  /// indices as local or global.
  ///
  /// It may also help to read CrsGraph's documentation.
  ///
  /// \section Tpetra_CrsMatrix_insertion_into_nonowned_rows Insertion into nonowned rows
  ///
  /// All methods (except for insertGlobalValues() and
  /// sumIntoGlobalValues(); see below) that work with global indices
  /// only allow operations on indices owned by the calling process.  For
  /// example, methods that take a global row index expect that row to be
  /// owned by the calling process.  Access to <i>nonowned rows</i>, that
  /// is, rows <i>not</i> owned by the calling process, requires
  /// performing an explicit communication via the Import / Export
  /// capabilities of the CrsMatrix object.  See the documentation of
  /// DistObject for more details.
  ///
  /// The methods insertGlobalValues() and sumIntoGlobalValues() are
  /// exceptions to this rule.  They both allows you to add data to
  /// nonowned rows.  These data are stored locally and communicated to
  /// the appropriate process on the next call to globalAssemble() or
  /// fillComplete().  This means that CrsMatrix provides the same
  /// nonowned insertion functionality that Epetra provides via
  /// Epetra_FECrsMatrix.
  ///
  /// \section Tpetra_DistObject_MultDist Note for developers on DistObject
  ///
  /// DistObject only takes a single Map as input to its constructor.
  /// MultiVector is an example of a subclass for which a single Map
  /// suffices to describe its data distribution.  In that case,
  /// DistObject's getMap() method obviously must return that Map.
  /// CrsMatrix is an example of a subclass that requires two Map
  /// objects: a row Map and a column Map.  For CrsMatrix, getMap()
  /// returns the row Map.  This means that doTransfer() (which
  /// CrsMatrix does not override) uses the row Map objects of the
  /// source and target CrsMatrix objects.  CrsMatrix in turn uses its
  /// column Map (if it has one) to "filter" incoming sparse matrix
  /// entries whose column indices are not in that process' column
  /// Map.  This means that CrsMatrix may perform extra communication,
  /// though the Import and Export operations are still correct.
  ///
  /// This is necessary if the CrsMatrix does not yet have a column
  /// Map.  Other processes might have added new entries to the
  /// matrix; the calling process has to see them in order to accept
  /// them.  However, the CrsMatrix may already have a column Map, for
  /// example, if it was created with the constructor that takes both
  /// a row and a column Map, or if it is fill complete (which creates
  /// the column Map if the matrix does not yet have one).  In this
  /// case, it could be possible to "filter" on the sender (instead of
  /// on the receiver, as CrsMatrix currently does) and avoid sending
  /// data corresponding to columns that the receiver does not own.
  /// Doing this would require revising the Import or Export object
  /// (instead of the incoming data) using the column Map, to remove
  /// global indices and their target process ranks from the send
  /// lists if the target process does not own those columns, and to
  /// remove global indices and their source process ranks from the
  /// receive lists if the calling process does not own those columns.
  /// (Abstractly, this is a kind of set difference between an Import
  /// or Export object for the row Maps, and the Import resp. Export
  /// object for the column Maps.)  This could be done separate from
  /// DistObject, by creating a new "filtered" Import or Export
  /// object, that keeps the same source and target Map objects but
  /// has a different communication plan.  We have not yet implemented
  /// this optimization.
  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  class CrsMatrix :
    public RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>,
    public DistObject<char, LocalOrdinal, GlobalOrdinal, Node>,
    public Details::Spaces::User
  {
  // clang-format on
private:
  using dist_object_type =
      DistObject<char, LocalOrdinal, GlobalOrdinal,
                 Node>; ///< Type of the DistObject specialization from which
                        ///< this class inherits.
  // clang-format off

  public:
    //! @name Typedefs
    //@{

    //! The type of each entry in the matrix.
    using scalar_type = Scalar;
    //! The type of each local index in the matrix.
    using local_ordinal_type = LocalOrdinal;
    //! The type of each global index in the matrix.
    using global_ordinal_type = GlobalOrdinal;
    //! The Kokkos device type.
    using device_type = typename Node::device_type;
    //! The Kokkos execution space.
    using execution_space = typename device_type::execution_space;
    //! The Kokkos memory space.
    using memory_space = typename device_type::memory_space;

    /// \brief This class' Kokkos Node type.
    ///
    /// This is a leftover that will be deprecated and removed.
    /// See e.g., GitHub Issue #57.
    using node_type = Node;

    //! The Map specialization suitable for this CrsMatrix specialization.
    using map_type = Map<LocalOrdinal, GlobalOrdinal, Node>;

    //! The Import specialization suitable for this CrsMatrix specialization.
    using import_type = Import<LocalOrdinal, GlobalOrdinal, Node>;

    //! The Export specialization suitable for this CrsMatrix specialization.
    using export_type = Export<LocalOrdinal, GlobalOrdinal, Node>;

    //! The RowMatrix representing the base class of CrsMatrix
    using row_matrix_type = RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

    /// \brief The type used internally in place of \c Scalar.
    ///
    /// Some \c Scalar types might not work with Kokkos on all
    /// execution spaces, due to missing CUDA device macros or
    /// volatile overloads.  The C++ standard type std::complex<T> has
    /// this problem.  To fix this, we replace std::complex<T> values
    /// internally with the (usually) bitwise identical type
    /// Kokkos::complex<T>.  The latter is the \c impl_scalar_type
    /// corresponding to \c Scalar = std::complex.
    using impl_scalar_type = typename row_matrix_type::impl_scalar_type;
    /// \brief Type of a norm result.
    ///
    /// This is usually the same as the type of the magnitude
    /// (absolute value) of <tt>Scalar</tt>, but may differ for
    /// certain <tt>Scalar</tt> types.
    using mag_type = typename Kokkos::ArithTraits<impl_scalar_type>::mag_type;

    //! The CrsGraph specialization suitable for this CrsMatrix specialization.
    using crs_graph_type = CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;

    //! The part of the sparse matrix's graph on each MPI process.
    using local_graph_device_type = typename crs_graph_type::local_graph_device_type;
    using local_graph_host_type = typename crs_graph_type::local_graph_host_type;

    /// \brief The specialization of Kokkos::CrsMatrix that represents
    ///   the part of the sparse matrix on each MPI process.
    using local_matrix_device_type =
      KokkosSparse::CrsMatrix<impl_scalar_type,
                              local_ordinal_type,
                              device_type,
                              void,
                              typename local_graph_device_type::size_type>;
    using local_matrix_host_type = 
          typename local_matrix_device_type::HostMirror;

#if KOKKOSKERNELS_VERSION < 40299
    /// \brief The type of the local matrix-vector operator (a wrapper of \c KokkosSparse::CrsMatrix )
    using local_multiply_op_type =
      LocalCrsMatrixOperator<scalar_type,
                             scalar_type,
                             device_type>;
#endif

    using row_ptrs_device_view_type = 
          typename row_matrix_type::row_ptrs_device_view_type;
    using row_ptrs_host_view_type = 
          typename row_matrix_type::row_ptrs_host_view_type;


    using local_inds_device_view_type = 
          typename row_matrix_type::local_inds_device_view_type;
    using local_inds_host_view_type = 
          typename row_matrix_type::local_inds_host_view_type;
    using nonconst_local_inds_host_view_type = 
          typename row_matrix_type::nonconst_local_inds_host_view_type;

    using global_inds_device_view_type = 
          typename row_matrix_type::global_inds_device_view_type;
    using global_inds_host_view_type = 
          typename row_matrix_type::global_inds_host_view_type;
    using nonconst_global_inds_host_view_type = 
          typename row_matrix_type::nonconst_global_inds_host_view_type;

    using values_device_view_type = 
          typename row_matrix_type::values_device_view_type;
    using values_host_view_type = 
          typename row_matrix_type::values_host_view_type;
    using nonconst_values_host_view_type = 
          typename row_matrix_type::nonconst_values_host_view_type;

    //@}
    //! @name Constructors and destructor
    //@{

    //! Copy constructor.
    CrsMatrix (const CrsMatrix<Scalar, LocalOrdinal,
                               GlobalOrdinal, Node>&) = default;

    //! Move constructor.
    CrsMatrix (CrsMatrix<Scalar, LocalOrdinal,
                         GlobalOrdinal, Node>&&) = default;

    //! Copy assignment.
    CrsMatrix&
    operator= (const CrsMatrix<Scalar, LocalOrdinal,
                               GlobalOrdinal, Node>&) = default;

    //! Move assignment.
    CrsMatrix&
    operator= (CrsMatrix<Scalar, LocalOrdinal,
                         GlobalOrdinal, Node>&&) = default;

    /// \brief Constructor specifying the maximum number of entries
    ///   that any row on the process can take.
    ///
    /// \param rowMap [in] Distribution of rows of the matrix.
    ///
    /// \param maxNumEntriesPerRow [in] Maximum number of matrix
    ///   entries per row.  This is a strict upper bound.  It may
    ///   differ on different processes.
    ///
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
               const size_t maxNumEntriesPerRow,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying (possibly different) number of entries in each row.
    ///
    /// \param rowMap [in] Distribution of rows of the matrix.
    ///
    /// \param numEntPerRowToAlloc [in] Maximum number of matrix
    ///   entries to allocate for each row.  This is a strict upper
    ///   bound.  It may differ on different processes.
    ///
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
               const Teuchos::ArrayView<const size_t>& numEntPerRowToAlloc,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);


    /// \brief Constructor specifying column Map and fixed number of entries for each row.
    ///
    /// The column Map will be used to filter any matrix entries
    /// inserted using insertLocalValues() or insertGlobalValues().
    ///
    /// \param rowMap [in] Distribution of rows of the matrix.
    ///
    /// \param colMap [in] Distribution of columns of the matrix.
    ///   See replaceColMap() for the requirements.
    ///
    /// \param maxNumEntPerRow [in] Maximum number of matrix entries
    ///   per row.  This is a strict upper bound.  It may differ on
    ///   different processes.
    ///
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
               const Teuchos::RCP<const map_type>& colMap,
               const size_t maxNumEntPerRow,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying column Map and number of entries in each row.
    ///
    /// The column Map will be used to filter any matrix indices
    /// inserted using insertLocalValues() or insertGlobalValues().
    ///
    /// \param rowMap [in] Distribution of rows of the matrix.
    ///
    /// \param colMap [in] Distribution of columns of the matrix.
    ///   See replaceColMap() for the requirements.
    ///
    /// \param numEntPerRowToAlloc [in] Maximum number of matrix
    ///   entries to allocate for each row.  This is a strict upper
    ///   bound.  It may differ on different processes.
    ///
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
               const Teuchos::RCP<const map_type>& colMap,
               const Teuchos::ArrayView<const size_t>& numEntPerRowToAlloc,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying a matrix and a previously
    ///   constructed graph, presumably a subset of the matrix's graph.
    ///   This matrix will alias the first N values of the passed-in
    ///   matrix, where N is the number of entries in the graph.
    ///
    /// Calling this constructor fixes the graph structure of the
    /// sparse matrix.  We say in this case that the matrix has a
    /// "static graph."  If you create a CrsMatrix with this
    /// constructor, you are not allowed to insert new entries into
    /// the matrix, but you are allowed to change values in the
    /// matrix.
    ///
    /// The given graph must be fill complete.  Note that calling
    /// resumeFill() on the graph makes it not fill complete, even if
    /// you had previously called fillComplete() on the graph.  In
    /// that case, you must call fillComplete() on the graph again
    /// before invoking this CrsMatrix constructor.
    ///
    /// This constructor is marked \c explicit so that you can't
    /// create a CrsMatrix by accident when passing a CrsGraph into a
    /// function that takes a CrsMatrix.
    ///
    /// \param matrix [in] The existing matrix whose values this one will alias.
    /// \param graph [in] The graph structure of the sparse matrix.
    ///   The graph <i>must</i> be fill complete.
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    explicit CrsMatrix (CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& matrix,
                        const Teuchos::RCP<const crs_graph_type>& graph,
                        const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying a previously constructed graph.
    ///
    /// Calling this constructor fixes the graph structure of the
    /// sparse matrix.  We say in this case that the matrix has a
    /// "static graph."  If you create a CrsMatrix with this
    /// constructor, you are not allowed to insert new entries into
    /// the matrix, but you are allowed to change values in the
    /// matrix.
    ///
    /// The given graph must be fill complete.  Note that calling
    /// resumeFill() on the graph makes it not fill complete, even if
    /// you had previously called fillComplete() on the graph.  In
    /// that case, you must call fillComplete() on the graph again
    /// before invoking this CrsMatrix constructor.
    ///
    /// This constructor is marked \c explicit so that you can't
    /// create a CrsMatrix by accident when passing a CrsGraph into a
    /// function that takes a CrsMatrix.
    ///
    /// \param graph [in] The graph structure of the sparse matrix.
    ///   The graph <i>must</i> be fill complete.
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    explicit CrsMatrix (const Teuchos::RCP<const crs_graph_type>& graph,
                        const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying a previously constructed graph
    ///   and entries array
    ///
    /// Calling this constructor fixes the graph structure of the
    /// sparse matrix.  We say in this case that the matrix has a
    /// "static graph."  If you create a CrsMatrix with this
    /// constructor, you are not allowed to insert new entries into
    /// the matrix, but you are allowed to change values in the
    /// matrix.
    ///
    /// The given graph must be fill complete.  Note that calling
    /// resumeFill() on the graph makes it not fill complete, even if
    /// you had previously called fillComplete() on the graph.  In
    /// that case, you must call fillComplete() on the graph again
    /// before invoking this CrsMatrix constructor.
    ///
    /// This constructor is marked \c explicit so that you can't
    /// create a CrsMatrix by accident when passing a CrsGraph into a
    /// function that takes a CrsMatrix.
    ///
    /// \param graph [in] The graph structure of the sparse matrix.
    ///   The graph <i>must</i> be fill complete.
    /// \param values [in] The local entries in the matrix,
    ///   as in a CSR "vals" array.  The length of this vector
    ///   should be equal to the number of unknowns in the matrix.
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    explicit CrsMatrix (const Teuchos::RCP<const crs_graph_type>& graph,
                        const typename local_matrix_device_type::values_type& values,
                        const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying column Map and arrays containing
    ///   the matrix in local indices.  In almost all cases the indices
    ///   must be sorted on input, but if they aren't sorted,
    ///   "sorted" must be set to false in params.
    ///
    /// \param rowMap [in] Distribution of rows of the matrix.
    ///
    /// \param colMap [in] Distribution of columns of the matrix.
    ///   See replaceColMap() for the requirements.
    ///
    /// \param rowPointers [in] The beginning of each row in the matrix,
    ///   as in a CSR "rowptr" array.  The length of this vector should be
    ///   equal to the number of rows in the graph, plus one.  This last
    ///   entry should store the nunber of nonzeros in the matrix.
    ///
    /// \param columnIndices [in] The local indices of the columns,
    ///   as in a CSR "colind" array.  The length of this vector
    ///   should be equal to the number of unknowns in the matrix.
    ///
    /// \param values [in] The local entries in the matrix,
    ///   as in a CSR "vals" array.  The length of this vector
    ///   should be equal to the number of unknowns in the matrix.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
               const Teuchos::RCP<const map_type>& colMap,
               const typename local_graph_device_type::row_map_type& rowPointers,
               const typename local_graph_device_type::entries_type::non_const_type& columnIndices,
               const typename local_matrix_device_type::values_type& values,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying column Map and arrays containing
    ///   the local matrix.  In almost all cases the local matrix
    ///   must be sorted on input, but if it isn't sorted,
    ///   "sorted" must be set to false in params.
    ///
    /// \param rowMap [in] Distribution of rows of the matrix.
    ///
    /// \param colMap [in] Distribution of columns of the matrix.
    ///   See replaceColMap() for the requirements.
    ///
    /// \param rowPointers [in] The beginning of each row in the matrix,
    ///   as in a CSR "rowptr" array.  The length of this vector should be
    ///   equal to the number of rows in the graph, plus one.  This last
    ///   entry should store the nunber of nonzeros in the matrix.
    ///
    /// \param columnIndices [in] The local indices of the columns,
    ///   as in a CSR "colind" array.  The length of this vector
    ///   should be equal to the number of unknowns in the matrix.
    ///
    /// \param values [in] The local entries in the matrix,
    ///   as in a CSR "vals" array.  The length of this vector
    ///   should be equal to the number of unknowns in the matrix.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
               const Teuchos::RCP<const map_type>& colMap,
               const Teuchos::ArrayRCP<size_t>& rowPointers,
               const Teuchos::ArrayRCP<LocalOrdinal>& columnIndices,
               const Teuchos::ArrayRCP<Scalar>& values,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying column Map and a local matrix,
    ///   which the resulting CrsMatrix views.
    ///
    /// Unlike most other CrsMatrix constructors, successful
    /// completion of this constructor will result in a fill-complete
    /// matrix.
    ///
    /// \param rowMap [in] Distribution of rows of the matrix.
    ///
    /// \param colMap [in] Distribution of columns of the matrix.
    ///   See replaceColMap() for the requirements.
    ///
    /// \param lclMatrix [in] A local CrsMatrix containing all local
    ///    matrix values as well as a local graph.  The graph's local
    ///    row indices must come from the specified row Map, and its
    ///    local column indices must come from the specified column
    ///    Map.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
               const Teuchos::RCP<const map_type>& colMap,
               const local_matrix_device_type& lclMatrix,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying column, domain and range Maps,
    ///   and a local matrix, which the resulting CrsMatrix views.
    ///   In almost all cases the local matrix
    ///   must be sorted on input, but if it isn't sorted,
    ///   "sorted" must be set to false in params.
    ///
    /// Unlike most other CrsMatrix constructors, successful
    /// completion of this constructor will result in a fill-complete
    /// matrix.
    ///
    /// \param rowMap [in] Distribution of rows of the matrix.
    ///
    /// \param colMap [in] Distribution of columns of the matrix.
    ///   See replaceColMap() for the requirements.
    ///
    /// \param domainMap [in] The matrix's domain Map.  MUST be one to
    ///   one!
    ///
    /// \param rangeMap [in] The matrix's range Map.  MUST be one to
    ///   one!  May be, but need not be, the same as the domain Map.
    ///
    /// \param lclMatrix [in] A local CrsMatrix containing all local
    ///    matrix values as well as a local graph.  The graph's local
    ///    row indices must come from the specified row Map, and its
    ///    local column indices must come from the specified column
    ///    Map.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsMatrix (const local_matrix_device_type& lclMatrix,
               const Teuchos::RCP<const map_type>& rowMap,
               const Teuchos::RCP<const map_type>& colMap,
               const Teuchos::RCP<const map_type>& domainMap = Teuchos::null,
               const Teuchos::RCP<const map_type>& rangeMap = Teuchos::null,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    //! Create a fill-complete CrsMatrix from all the things it needs.
    ///   \param lclMatrix [in] In almost all cases the local matrix
    ///     must be sorted on input, but if it isn't sorted,
    ///     "sorted" must be set to false in params.
    CrsMatrix (const local_matrix_device_type& lclMatrix,
               const Teuchos::RCP<const map_type>& rowMap,
               const Teuchos::RCP<const map_type>& colMap,
               const Teuchos::RCP<const map_type>& domainMap,
               const Teuchos::RCP<const map_type>& rangeMap,
               const Teuchos::RCP<const import_type>& importer,
               const Teuchos::RCP<const export_type>& exporter,
               const Teuchos::RCP<Teuchos::ParameterList>& params =
                 Teuchos::null);

    /// \brief Copy constructor, with option to do deep or shallow copy.
    // This function in 'Copy' mode is only guaranteed to work correctly for matrices
    // which are fillComplete.
    CrsMatrix (const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& source,
               const Teuchos::DataAccess copyOrView);

    /// \brief Destructor (virtual for memory safety of derived classes).
    ///
    /// \note To Tpetra developers: See the C++ Core Guidelines C.21
    ///   ("If you define or <tt>=delete</tt> any default operation,
    ///   define or <tt>=delete</tt> them all"), in particular the
    ///   AbstractBase example, for why this destructor declaration
    ///   implies that we need the above four <tt>=default</tt>
    ///   declarations for copy construction, move construction, copy
    ///   assignment, and move assignment.
    virtual ~CrsMatrix () = default;

    // This friend declaration makes the clone() method work.
    template <class S2, class LO2, class GO2, class N2>
    friend class CrsMatrix;

    // This friend declaration allows for fused residual calculation
    template <class S2, class LO2, class GO2, class N2>
    friend void Details::residual(const Operator<S2,LO2,GO2,N2> &   A,
                                  const MultiVector<S2,LO2,GO2,N2> & X,
                                  const MultiVector<S2,LO2,GO2,N2> & B,
                                  MultiVector<S2,LO2,GO2,N2> & R);

    // This friend declaration allows for batching of apply calls
    template <class MatrixArray, class MultiVectorArray>
    friend void batchedApply(const MatrixArray &Matrices,
                             const typename std::remove_pointer<typename MultiVectorArray::value_type>::type & X,
                             MultiVectorArray &Y,
                             typename std::remove_pointer<typename MatrixArray::value_type>::type::scalar_type alpha,
                             typename std::remove_pointer<typename MatrixArray::value_type>::type::scalar_type beta,
                             Teuchos::RCP<Teuchos::ParameterList> params);

  public:
    //@}
    //! @name Methods for inserting, modifying, or removing entries
    //@{

    /// \brief Insert one or more entries into the matrix, using
    ///   global column indices.
    ///
    /// \param globalRow [in] Global index of the row into which to
    ///   insert the entries.
    /// \param cols [in] Global indices of the columns into which
    ///   to insert the entries.
    /// \param vals [in] Values to insert into the above columns.
    ///
    /// For all k in 0, ..., <tt>col.size()-1</tt>, insert the value
    /// <tt>values[k]</tt> into entry <tt>(globalRow, cols[k])</tt> of
    /// the matrix.  If that entry already exists, add the new value
    /// to the old value.
    ///
    /// This is a local operation.  It does not communicate (using
    /// MPI).  If row \c globalRow is owned by the calling process,
    /// the entries will be inserted immediately.  Otherwise, if that
    /// row is <i>not</i> owned by the calling process, then the
    /// entries will be stored locally for now, and only communicated
    /// to the process that owns the row when either fillComplete() or
    /// globalAssemble() is called.  If that process already has an
    /// entry, the incoming value will be added to the old value, just
    /// as if it were inserted on the owning process.
    //
    /// If the matrix has a column Map (<tt>hasColMap() == true</tt>),
    /// and if globalRow is owned by process p, then it is forbidden
    /// to insert column indices that are not in the column Map on
    /// process p.  Tpetra will test the input column indices to
    /// ensure this is the case, but if \c globalRow is not owned by
    /// the calling process, the test will be deferred until the next
    /// call to globalAssemble() or fillComplete().
    ///
    /// \warning The behavior described in the above paragraph differs
    ///   from that of Epetra.  If the matrix has a column Map,
    ///   Epetra_CrsMatrix "filters" column indices not in the column
    ///   Map.  Many users found this confusing, so we changed it so
    ///   that nonowned column indices are forbidden.
    ///
    /// It is legal to call this method whether the matrix's column
    /// indices are globally or locally indexed.  If the matrix's
    /// column indices are locally indexed (<tt>isLocallyIndexed() ==
    /// true</tt>), then this method will convert the input global
    /// column indices to local column indices.
    ///
    /// For better performance when filling entries into a sparse
    /// matrix, consider the following tips:
    /// <ol>
    /// <li>Use local indices (e.g., insertLocalValues()) if you know
    ///   the column Map in advance.  Converting global indices to
    ///   local indices is expensive.  Of course, if you don't know
    ///   the column Map in advance, you must use global indices.</li>
    /// <li>When invoking the CrsMatrix constructor, give the best
    ///   possible upper bounds on the number of entries in each row
    ///   of the matrix.  This will avoid expensive reallocation if
    ///   your bound was not large enough.</li>
    /// <li>If you plan to reuse a matrix's graph structure, but
    ///   change its values, in repeated fillComplete() / resumeFill()
    ///   cycles, you can get the best performance by creating the
    ///   matrix with a const CrsGraph.  Do this by using the
    ///   CrsMatrix constructor that accepts an RCP of a const
    ///   CrsGraph.  If you do this, you must use the "replace" or
    ///   "sumInto" methods to change the values of the matrix; you
    ///   may not use insertGlobalValues() or
    ///   insertLocalValues().</li>
    /// </ol>
    void
    insertGlobalValues (const GlobalOrdinal globalRow,
                        const Teuchos::ArrayView<const GlobalOrdinal>& cols,
                        const Teuchos::ArrayView<const Scalar>& vals);

    /// \brief Epetra compatibility version of insertGlobalValues (see
    ///   above) that takes arguments as raw pointers, rather than
    ///   Teuchos::ArrayView.
    ///
    /// Arguments are the same and in the same order as
    /// Epetra_CrsMatrix::InsertGlobalValues.
    ///
    /// \param globalRow [in] Global index of the row into which to
    ///   insert the entries.
    /// \param numEnt [in] Number of entries to insert; number of
    ///   valid entries in \c vals and \c inds.
    /// \param vals [in] Values to insert.
    /// \param inds [in] Global indices of the columns into which
    ///   to insert the entries.
    void
    insertGlobalValues (const GlobalOrdinal globalRow,
                        const LocalOrdinal numEnt,
                        const Scalar vals[],
                        const GlobalOrdinal inds[]);

    /// \brief Insert one or more entries into the matrix, using local
    ///   column indices.
    ///
    /// \param localRow [in] Local index of the row into which to
    ///   insert the entries.  It must be owned by the row Map on the
    ///   calling process.
    /// \param cols [in] Local indices of the columns into which to
    ///   insert the entries.  All of the column indices must be owned
    ///   by the column Map on the calling process.
    /// \param vals [in] Values to insert into the above columns.
    /// \param CM [in] How values should be inserted. Valid options
    ///   are: ADD (default) inserts values that are not yet in the
    ///   matrix graph, and sums values that are already present. INSERT
    ///   inserts values that are not yet in the matrix graph, and
    ///   replaces values that are already present.
    ///
    /// For all k in 0, ..., <tt>cols.size()-1</tt>, insert the value
    /// <tt>values[k]</tt> into entry <tt>(globalRow, cols[k])</tt> of
    /// the matrix.  If that entry already exists, add the new value
    /// to the old value, if <tt>CM=ADD</tt>, otherwise replace
    /// the old value.
    ///
    /// In order to call this method, the matrix must be locally
    /// indexed, and it must have a column Map.
    ///
    /// For better performance when filling entries into a sparse
    /// matrix, consider the following tips:
    /// <ol>
    /// <li>When invoking the CrsMatrix constructor, give the best
    ///   possible upper bounds on the number of entries in each row
    ///   of the matrix.  This will avoid expensive reallocation if
    ///   your bound was not large enough.</li>
    /// <li>If you plan to reuse a matrix's graph structure, but
    ///   change its values, in repeated fillComplete() / resumeFill()
    ///   cycles, you can get the best performance by creating the
    ///   matrix with a const CrsGraph.  Do this by using the
    ///   CrsMatrix constructor that accepts an RCP of a const
    ///   CrsGraph.  If you do this, you must use the "replace" or
    ///   "sumInto" methods to change the values of the matrix; you
    ///   may not use insertGlobalValues() or
    ///   insertLocalValues().</li>
    /// </ol>
    void
    insertLocalValues (const LocalOrdinal localRow,
                       const Teuchos::ArrayView<const LocalOrdinal> &cols,
                       const Teuchos::ArrayView<const Scalar> &vals,
                       const CombineMode CM=ADD);

    /// \brief Epetra compatibility version of insertLocalValues (see
    ///   above) that takes arguments as raw pointers, rather than
    ///   Teuchos::ArrayView.
    ///
    /// Arguments are the same and in the same order as
    /// Epetra_CrsMatrix::InsertMyValues.
    ///
    /// \param localRow [in] Local index of the row into which to
    ///   insert the entries.
    /// \param numEnt [in] Number of entries to insert; number of
    ///   valid entries in \c vals and \c cols.
    /// \param vals [in] Values to insert.
    /// \param cols [in] Global indices of the columns into which
    ///   to insert the entries.
    /// \param CM [in] How values should be inserted. Valid options
    ///   are: ADD (default) inserts values that are not yet in the
    ///   matrix graph, and sums values that are already present. INSERT
    ///   inserts values that are not yet in the matrix graph, and
    ///   replaces values that are already present.
    void
    insertLocalValues (const LocalOrdinal localRow,
                       const LocalOrdinal numEnt,
                       const Scalar vals[],
                       const LocalOrdinal cols[],
                       const CombineMode CM=ADD);

  protected:
    /// \brief Implementation detail of replaceGlobalValues.
    ///
    /// \param rowVals [in/out] On input: Values of the row of the
    ///   sparse matrix to modify.  On output: The modified values.
    /// \param graph [in] The matrix's graph.
    /// \param rowInfo [in] Result of graph.getRowInfo on the index of
    ///   the local row of the matrix to modify.
    /// \param inds [in] Global column indices of that row to modify.
    /// \param newVals [in] For each k, replace the value in rowVals
    ///   corresponding to local column index inds[k] with newVals[k].
    virtual LocalOrdinal
    replaceGlobalValuesImpl (impl_scalar_type rowVals[],
                             const crs_graph_type& graph,
                             const RowInfo& rowInfo,
                             const GlobalOrdinal inds[],
                             const impl_scalar_type newVals[],
                             const LocalOrdinal numElts);

  public:
    /// \brief Replace one or more entries' values, using global indices.
    ///
    /// \param globalRow [in] Global index of the row in which to
    ///   replace the entries.  This row <i>must</i> be owned by the
    ///   calling process.
    /// \param inputInds [in] Kokkos::View of the global indices of
    ///   the columns in which to replace the entries.
    /// \param inputVals [in] Kokkos::View of the values to use for
    ///   replacing the entries.
    ///
    /// For all k in 0, ..., <tt>inputInds.extent(0)-1</tt>,
    /// replace the value at entry <tt>(globalRow, inputInds(k))</tt>
    /// of the matrix with <tt>inputVals(k)</tt>.  That entry must
    /// exist in the matrix already.
    ///
    /// If <tt>(globalRow, inputInds(k))</tt> corresponds to an entry
    /// that is duplicated in this matrix row (likely because it was
    /// inserted more than once and fillComplete() has not been called
    /// in the interim), the behavior of this method is not defined.
    ///
    /// \return The number of indices for which values were actually
    ///   replaced; the number of "correct" indices.
    ///
    /// If the returned value N satisfies
    ///
    /// <tt>0 <= N < inputInds.extent(0)</tt>,
    ///
    /// then <tt>inputInds.extent(0) - N</tt> of the entries of
    /// <tt>cols</tt> are not valid global column indices.  If the
    /// returned value is
    /// <tt>Teuchos::OrdinalTraits<LocalOrdinal>::invalid()</tt>, then
    /// at least one of the following is true:
    /// <ul>
    /// <li> <tt>! isFillActive ()</tt> </li>
    /// <li> <tt> inputInds.extent (0) != inputVals.extent (0)</tt> </li>
    /// </ul>
    local_ordinal_type
    replaceGlobalValues(
      const global_ordinal_type globalRow,
      const Kokkos::View<const global_ordinal_type*, Kokkos::AnonymousSpace>& inputInds,
      const Kokkos::View<const impl_scalar_type*, Kokkos::AnonymousSpace>& inputVals);

    /// \brief Overload of replaceGlobalValues (see above), that takes
    ///   Teuchos::ArrayView (host pointers) instead of Kokkos::View.
    LocalOrdinal
    replaceGlobalValues (const GlobalOrdinal globalRow,
                         const Teuchos::ArrayView<const GlobalOrdinal>& cols,
                         const Teuchos::ArrayView<const Scalar>& vals);

    /// \brief Overload of replaceGlobalValues (see above), that takes
    ///   raw pointers instead of Kokkos::View.
    ///
    /// This version of the method takes the same arguments in the
    /// same order as Epetra_CrsMatrix::ReplaceGlobalValues.
    ///
    /// \param globalRow [in] Global index of the row in which to
    ///   replace the entries.  This row <i>must</i> be owned by the
    ///   calling process.
    /// \param numEnt [in] Number of entries to replace; number of
    ///   valid entries in \c vals and \c cols.
    /// \param vals [in] Values to use for replacing the entries.
    /// \param cols [in] Global indices of the columns in which to
    ///   replace the entries.
    LocalOrdinal
    replaceGlobalValues (const GlobalOrdinal globalRow,
                         const LocalOrdinal numEnt,
                         const Scalar vals[],
                         const GlobalOrdinal cols[]);

  protected:
    /// \brief Implementation detail of replaceLocalValues.
    ///
    /// \param rowVals [in/out] On input: Values of the row of the
    ///   sparse matrix to modify.  On output: The modified values.
    /// \param graph [in] The matrix's graph.
    /// \param rowInfo [in] Result of graph.getRowInfo on the index of
    ///   the local row of the matrix to modify.
    /// \param inds [in] Local column indices of that row to modify.
    /// \param newVals [in] For each k, replace the value in rowVals
    ///   corresponding to local column index inds[k] with newVals[k].
    virtual LocalOrdinal
    replaceLocalValuesImpl (impl_scalar_type rowVals[],
                            const crs_graph_type& graph,
                            const RowInfo& rowInfo,
                            const LocalOrdinal inds[],
                            const impl_scalar_type newVals[],
                            const LocalOrdinal numElts);

  public:
    /// \brief Replace one or more entries' values, using local
    ///   row and column indices.
    ///
    /// \param localRow [in] local index of the row in which to
    ///   replace the entries.  This row <i>must</i> be owned by the
    ///   calling process.
    /// \param cols [in] Local indices of the columns in which to
    ///   replace the entries.
    /// \param vals [in] Values to use for replacing the entries.
    ///
    /// For local row index \c localRow and local column indices
    /// <tt>cols</tt>, do <tt>A(localRow, cols(k)) = vals(k)</tt>.
    /// The row index and column indices must be valid on the calling
    /// process, and all matrix entries <tt>A(localRow, cols(k))</tt>
    /// must already exist.  (This method does <i>not</i> change the
    /// matrix's structure.)  If the row index is valid, any invalid
    /// column indices are ignored, but counted in the return value.
    ///
    /// \return The number of indices for which values were actually
    ///   replaced; the number of "correct" indices.
    ///
    /// If the returned value N satisfies
    ///
    /// <tt>0 <= N < cols.extent(0)</tt>,
    ///
    /// then <tt>cols.extent(0) - N</tt> of the entries of
    /// <tt>cols</tt> are not valid local column indices.  If the
    /// returned value is
    /// <tt>Teuchos::OrdinalTraits<LocalOrdinal>::invalid()</tt>,
    /// then at least one of the following is true:
    ///   <ul>
    ///   <li> <tt>! isFillActive ()</tt> </li>
    ///   <li> <tt>! hasColMap ()</tt> </li>
    ///   <li> <tt> cols.extent (0) != vals.extent (0)</tt> </li>
    ///   </ul>
    local_ordinal_type
    replaceLocalValues(
      const local_ordinal_type localRow,
      const Kokkos::View<const local_ordinal_type*, Kokkos::AnonymousSpace>& inputInds,
      const Kokkos::View<const impl_scalar_type*, Kokkos::AnonymousSpace>& inputVals);

    /// \brief Backwards compatibility version of replaceLocalValues
    ///   (see above), that takes Teuchos::ArrayView (host pointers)
    ///   instead of Kokkos::View.
    LocalOrdinal
    replaceLocalValues (const LocalOrdinal localRow,
                        const Teuchos::ArrayView<const LocalOrdinal>& cols,
                        const Teuchos::ArrayView<const Scalar>& vals);

    /// \brief Epetra compatibility version of replaceLocalValues,
    ///   that takes raw pointers instead of Kokkos::View.
    ///
    /// This version of the method takes the same arguments in the
    /// same order as Epetra_CrsMatrix::ReplaceMyValues.
    ///
    /// \param localRow [in] local index of the row in which to
    ///   replace the entries.  This row <i>must</i> be owned by the
    ///   calling process.
    /// \param numEnt [in] Number of entries to replace; number of
    ///   valid entries in \c inputVals and \c inputCols.
    /// \param inputVals [in] Values to use for replacing the entries.
    /// \param inputCols [in] Local indices of the columns in which to
    ///   replace the entries.
    ///
    /// \return The number of indices for which values were actually
    ///   replaced; the number of "correct" indices.
    LocalOrdinal
    replaceLocalValues (const LocalOrdinal localRow,
                        const LocalOrdinal numEnt,
                        const Scalar inputVals[],
                        const LocalOrdinal inputCols[]);

  private:
    /// \brief Whether sumIntoLocalValues and sumIntoGlobalValues
    ///   should use atomic updates by default.
    ///
    /// \warning This is an implementation detail.
    static const bool useAtomicUpdatesByDefault =
#ifdef KOKKOS_ENABLE_SERIAL
      ! std::is_same<execution_space, Kokkos::Serial>::value;
#else
      true;
#endif // KOKKOS_ENABLE_SERIAL

    /// \brief Implementation detail of sumIntoGlobalValues.
    ///
    /// \tparam InputMemorySpace Kokkos memory space / device in which
    ///   the input data live.  This may differ from the memory space
    ///   in which the current matrix values (rowVals) live.
    /// \tparam ValsMemorySpace Kokkos memory space / device in which
    ///   the matrix's current values live.  This may differ from the
    ///   memory space in which the input data (inds and newVals)
    ///   live.
    ///
    /// \param rowVals [in/out] On input: Values of the row of the
    ///   sparse matrix to modify.  On output: The modified values.
    /// \param graph [in] The matrix's graph.
    /// \param rowInfo [in] Result of getRowInfo on the index of the
    ///   local row of the matrix to modify.
    /// \param inds [in] Global column indices of that row to modify.
    /// \param newVals [in] For each k, increment the value in rowVals
    ///   corresponding to global column index inds[k] by newVals[k].
    ///
    /// \return The number of valid input column indices.  In case of
    ///   error other than one or more invalid column indices, this
    ///   method returns
    ///   Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
  protected:
    virtual LocalOrdinal
    sumIntoGlobalValuesImpl (impl_scalar_type rowVals[],
                             const crs_graph_type& graph,
                             const RowInfo& rowInfo,
                             const GlobalOrdinal inds[],
                             const impl_scalar_type newVals[],
                             const LocalOrdinal numElts,
                             const bool atomic = useAtomicUpdatesByDefault);

  public:
    /// \brief Sum into one or more sparse matrix entries, using
    ///   global indices.
    ///
    /// This is a local operation; it does not involve communication.
    /// However, if you sum into rows not owned by the calling
    /// process, it may result in future communication in
    /// globalAssemble() (which is called by fillComplete()).
    ///
    /// If \c globalRow is owned by the calling process, then this
    /// method performs the sum-into operation right away.  Otherwise,
    /// if the row is <i>not</i> owned by the calling process, this
    /// method defers the sum-into operation until globalAssemble().
    /// That method communicates data for nonowned rows to the
    /// processes that own those rows.  Then, globalAssemble() does
    /// one of the following:
    /// <ul>
    /// <li> It calls insertGlobalValues() for that data if the matrix
    ///      has a dynamic graph. </li>
    /// <li> It calls sumIntoGlobalValues() for that data if the matrix
    ///      has a static graph.  The matrix silently ignores
    ///      (row,column) pairs that do not exist in the graph.
    /// </ul>
    ///
    /// \param globalRow [in] The global index of the row in which to
    ///   sum into the matrix entries.
    /// \param cols [in] One or more column indices.
    /// \param vals [in] One or more values corresponding to those
    ///   column indices.  <tt>vals[k]</tt> corresponds to
    ///   <tt>cols[k]</tt>.
    /// \param atomic [in] Whether to use atomic updates.
    ///
    /// \return The number of indices for which values were actually
    ///   modified; the number of "correct" indices.
    ///
    /// This method has the same preconditions and return value
    /// meaning as replaceGlobalValues() (which see).
    LocalOrdinal
    sumIntoGlobalValues (const GlobalOrdinal globalRow,
                         const Teuchos::ArrayView<const GlobalOrdinal>& cols,
                         const Teuchos::ArrayView<const Scalar>& vals,
                         const bool atomic = useAtomicUpdatesByDefault);

    /// \brief Epetra compatibility version of sumIntoGlobalValues
    ///   (see above), that takes input as raw pointers instead of
    ///   Kokkos::View.
    ///
    /// Arguments are the same and in the same order as those of
    /// Epetra_CrsMatrix::SumIntoGlobalValues, except for \c atomic,
    /// which is as above.
    ///
    /// \param globalRow [in] The global index of the row in which to
    ///   sum into the matrix entries.
    /// \param numEnt [in] Number of valid entries in \c vals and
    ///   \c cols.  This has type \c LocalOrdinal because we assume
    ///   that users will never want to insert more column indices
    ///   in one call than the matrix has columns.
    /// \param vals [in] \c numEnt values corresponding to the column
    ///   indices in \c cols.  That is, \c vals[k] is the value
    ///   corresponding to \c cols[k].
    /// \param cols [in] \c numEnt global column indices.
    /// \param atomic [in] Whether to use atomic updates.
    ///
    /// \return The number of indices for which values were actually
    ///   modified; the number of "correct" indices.
    LocalOrdinal
    sumIntoGlobalValues (const GlobalOrdinal globalRow,
                         const LocalOrdinal numEnt,
                         const Scalar vals[],
                         const GlobalOrdinal cols[],
                         const bool atomic = useAtomicUpdatesByDefault);

  protected:
    /// \brief Implementation detail of sumIntoLocalValues.
    ///
    /// \param rowVals [in/out] On input: Values of the row of the
    ///   sparse matrix to modify.  On output: The modified values.
    /// \param graph [in] The matrix's graph.
    /// \param rowInfo [in] Result of graph.getRowInfo on the index of
    ///   the local row of the matrix to modify.
    /// \param inds [in] Local column indices of that row to modify.
    /// \param newVals [in] For each k, increment the value in rowVals
    ///   corresponding to local column index inds[k] by newVals[k].
    /// \param atomic [in] Whether to use atomic updates (+=) when
    ///   incrementing values.
    virtual LocalOrdinal
    sumIntoLocalValuesImpl (impl_scalar_type rowVals[],
                            const crs_graph_type& graph,
                            const RowInfo& rowInfo,
                            const LocalOrdinal inds[],
                            const impl_scalar_type newVals[],
                            const LocalOrdinal numElts,
                            const bool atomic = useAtomicUpdatesByDefault);

  public:
    /// \brief Sum into one or more sparse matrix entries, using local
    ///   row and column indices.
    ///
    /// For local row index \c localRow and local column indices
    /// <tt>cols</tt>, perform the update <tt>A(localRow, cols[k]) +=
    /// vals[k]</tt>.  The row index and column indices must be valid
    /// on the calling process, and all matrix entries <tt>A(localRow,
    /// cols[k])</tt> must already exist.  (This method does
    /// <i>not</i> change the matrix's structure.)  If the row index
    /// is valid, any invalid column indices are ignored, but counted
    /// in the return value.
    ///
    /// This overload of the method takes the column indices and
    /// values as Kokkos::View.  See below for an overload that takes
    /// Teuchos::ArrayView instead.
    ///
    /// \tparam LocalIndicesViewType Kokkos::View specialization that
    ///   is a 1-D array of LocalOrdinal.
    /// \tparam ImplScalarViewType Kokkos::View specialization that is
    ///   a 1-D array of impl_scalar_type (usually the same as Scalar,
    ///   unless Scalar is std::complex<T> for some T, in which case
    ///   it is Kokkos::complex<T>).
    ///
    /// \param localRow [in] Local index of a row.  This row
    ///   <i>must</i> be owned by the calling process.
    /// \param cols [in] Local indices of the columns whose entries we
    ///   want to modify.
    /// \param vals [in] Values corresponding to the above column
    ///   indices.  <tt>vals(k)</tt> corresponds to <tt>cols(k)</tt>.
    /// \param atomic [in] Whether to use atomic updates.
    ///
    /// \return The number of indices for which values were actually
    ///   modified; the number of "correct" indices.
    ///
    /// This method has the same preconditions and return value
    /// meaning as replaceLocalValues() (which see).
    local_ordinal_type
    sumIntoLocalValues(
      const local_ordinal_type localRow,
      const Kokkos::View<const local_ordinal_type*, Kokkos::AnonymousSpace>& inputInds,
      const Kokkos::View<const impl_scalar_type*, Kokkos::AnonymousSpace>& inputVals,
      const bool atomic = useAtomicUpdatesByDefault);

    /// \brief Sum into one or more sparse matrix entries, using local
    ///   row and column indices.
    ///
    /// For local row index \c localRow and local column indices
    /// <tt>cols</tt>, perform the update <tt>A(localRow, cols[k]) +=
    /// vals[k]</tt>.  The row index and column indices must be valid
    /// on the calling process, and all matrix entries <tt>A(localRow,
    /// cols[k])</tt> must already exist.  (This method does
    /// <i>not</i> change the matrix's structure.)  If the row index
    /// is valid, any invalid column indices are ignored, but counted
    /// in the return value.
    ///
    /// This overload of the method takes the column indices and
    /// values as Teuchos::ArrayView.  See above for an overload that
    /// takes Kokkos::View instead.
    ///
    /// \param localRow [in] Local index of a row.  This row
    ///   <i>must</i> be owned by the calling process.
    /// \param cols [in] Local indices of the columns whose entries we
    ///   want to modify.
    /// \param vals [in] Values corresponding to the above column
    ///   indices.  <tt>vals[k]</tt> corresponds to <tt>cols[k]</tt>.
    /// \param atomic [in] Whether to use atomic updates.
    ///
    /// \return The number of indices for which values were actually
    ///   modified; the number of "correct" indices.
    ///
    /// This method has the same preconditions and return value
    /// meaning as replaceLocalValues() (which see).
    LocalOrdinal
    sumIntoLocalValues (const LocalOrdinal localRow,
                        const Teuchos::ArrayView<const LocalOrdinal>& cols,
                        const Teuchos::ArrayView<const Scalar>& vals,
                        const bool atomic = useAtomicUpdatesByDefault);

    /// \brief Epetra compatibility version of sumIntoLocalValues (see
    ///   above) that takes raw pointers instead of Kokkos::View.
    ///
    /// Arguments are the same and in the same order as
    /// Epetra_CrsMatrix::SumIntoMyValues, except for the \c atomic
    /// last argument, which is as above.
    ///
    /// \param localRow [in] The local index of the row in which to
    ///   sum into the matrix entries.
    /// \param numEnt [in] Number of valid entries in \c vals and
    ///   \c cols.  This has type \c LocalOrdinal because we assume
    ///   that users will never want to insert more column indices
    ///   in one call than the matrix has columns.
    /// \param vals [in] \c numEnt values corresponding to the column
    ///   indices in \c cols.  That is, \c vals[k] is the value
    ///   corresponding to \c cols[k].
    /// \param cols [in] \c numEnt local column indices.
    /// \param atomic [in] Whether to use atomic updates.
    ///
    /// \return The number of indices for which values were actually
    ///   modified; the number of "correct" indices.
    LocalOrdinal
    sumIntoLocalValues (const LocalOrdinal localRow,
                        const LocalOrdinal numEnt,
                        const Scalar vals[],
                        const LocalOrdinal cols[],
                        const bool atomic = useAtomicUpdatesByDefault);

  private:
    /// \brief Transform the given values using local indices.
    ///
    /// \param rowVals [in/out] The values to be transformed.  They
    ///   correspond to the row indicated by rowInfo.
    /// \param graph [in] The matrix's graph; <tt>*(this->staticGraph_)</tt>.
    /// \param rowInfo [in] Result of graph.getRowInfo(lclRow), where lclRow
    ///   is the local index of the row in which to transform values.
    ///   For <tt>rowInfo = getRowInfo(lclRow)</tt>,
    ///   <tt>rowInfo.localRow == lclRow</tt>.
    /// \param inds [in] (Local) column indices, for which to
    ///   transform the corresponding values in rowVals.
    /// \param newVals [in] Values to use for transforming rowVals.
    ///   It's probably OK for these to alias rowVals.
    /// \param numElts [in] Number of entries in inds and newVals.
    /// \param f [in] A binary function used to transform rowVals.  It
    ///   takes two <tt>const impl_scalar_type&</tt> arguments, and
    ///   returns impl_scalar_type.
    ///
    /// This method transforms the values using the expression
    /// \code
    /// newVals[k] = f (rowVals[k], newVals[j]);
    /// \endcode
    /// where k is the local index corresponding to <tt>inds[j]</tt>.
    /// It ignores invalid local column indices, but they are counted
    /// in the return value.
    ///
    /// \return The number of valid local column indices.  In case of
    ///   error other than one or more invalid column indices, this
    ///   method returns
    ///   Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
    LocalOrdinal
    transformLocalValues (impl_scalar_type rowVals[],
                          const crs_graph_type& graph,
                          const RowInfo& rowInfo,
                          const LocalOrdinal inds[],
                          const impl_scalar_type newVals[],
                          const LocalOrdinal numElts,
                          std::function<impl_scalar_type (const impl_scalar_type&, const impl_scalar_type&) > f,
                          const bool atomic = useAtomicUpdatesByDefault);

    /// \brief Transform the given values using global indices.
    ///
    /// \param rowVals [in/out] The values to be transformed.  They
    ///   correspond to the row indicated by rowInfo.
    /// \param graph [in] The matrix's graph; <tt>*(this->staticGraph_)</tt>.
    /// \param rowInfo [in] Result of graph.getRowInfo(lclRow), where lclRow
    ///   is the local index of the row in which to transform values.
    ///   For <tt>rowInfo = getRowInfo(lclRow)</tt>,
    ///   <tt>rowInfo.localRow == lclRow</tt>.
    /// \param inds [in] (Global) column indices, for which to
    ///   transform the corresponding values in rowVals.
    /// \param newVals [in] Values to use for transforming rowVals.
    ///   It's probably OK for these to alias rowVals.
    /// \param numElts [in] Number of entries in inds and newVals.
    /// \param f [in] A binary function used to transform rowVals.  It
    ///   takes two <tt>const impl_scalar_type&</tt> arguments, and
    ///   returns impl_scalar_type.
    ///
    /// This method transforms the values using the expression
    /// \code
    /// newVals[k] = f (rowVals[k], newVals[j]);
    /// \endcode
    /// where k is the local index corresponding to <tt>inds[j]</tt>.
    /// It ignores invalid input column indices, but they are counted
    /// in the return value.
    ///
    /// \return The number of valid input column indices.  In case of
    ///   error other than one or more invalid column indices, this
    ///   method returns
    ///   Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
    LocalOrdinal
    transformGlobalValues (impl_scalar_type rowVals[],
                           const crs_graph_type& graph,
                           const RowInfo& rowInfo,
                           const GlobalOrdinal inds[],
                           const impl_scalar_type newVals[],
                           const LocalOrdinal numElts,
                           std::function<impl_scalar_type (const impl_scalar_type&, const impl_scalar_type&) > f,
                           const bool atomic = useAtomicUpdatesByDefault);

    /// \brief Transform the given values using local indices.
    ///
    /// \param lclRow [in] Local index of the row in which to transform.
    /// \param numInputEnt [in] Number of entries in inputVals and
    ///   inputCols.
    /// \param inputVals [in] Values to use for transforming the
    ///   values.
    /// \param inputCols [in] (Local) column indices, for which to
    ///   transform the corresponding values.
    /// \param f [in] A binary function used to transform rowVals.  It
    ///   takes two <tt>const impl_scalar_type&</tt> arguments, and
    ///   returns impl_scalar_type.
    /// \param atomic [in] Whether to use atomic updates.
    ///
    /// This method transforms the values using the expression
    /// \code
    /// newVals[k] = f (rowVals[k], newVals[j]);
    /// \endcode
    /// where k is the local index corresponding to
    /// <tt>inputInds[j]</tt>.  It ignores invalid input column
    /// indices, but they are counted in the return value.
    ///
    /// \return The number of valid input column indices.  In case of
    ///   error other than one or more invalid column indices, this
    ///   method returns
    ///   Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
    LocalOrdinal
    transformLocalValues (const LocalOrdinal lclRow,
                          const LocalOrdinal numInputEnt,
                          const impl_scalar_type inputVals[],
                          const LocalOrdinal inputCols[],
                          std::function<impl_scalar_type (const impl_scalar_type&, const impl_scalar_type&) > f,
                          const bool atomic = useAtomicUpdatesByDefault);

    /// \brief Transform the given values using global indices.
    ///
    /// \param gblRow [in] Global index of the row in which to transform.
    /// \param numInputEnt [in] Number of entries in inputVals and
    ///   inputCols.
    /// \param inputVals [in] Values to use for transforming the
    ///   values.
    /// \param inputCols [in] (Global) column indices, for which to
    ///   transform the corresponding values.
    /// \param f [in] A binary function used to transform rowVals.  It
    ///   takes two <tt>const impl_scalar_type&</tt> arguments, and
    ///   returns impl_scalar_type.
    /// \param atomic [in] Whether to use atomic updates.
    ///
    /// This method transforms the values using the expression
    /// \code
    /// newVals[k] = f (rowVals[k], newVals[j]);
    /// \endcode
    /// where k is the local index corresponding to
    /// <tt>inputInds[j]</tt>.  It ignores invalid input column
    /// indices, but they are counted in the return value.
    ///
    /// \return The number of valid local column indices.  In case of
    ///   error other than one or more invalid column indices, this
    ///   method returns
    ///   Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
    LocalOrdinal
    transformGlobalValues (const GlobalOrdinal gblRow,
                           const LocalOrdinal numInputEnt,
                           const impl_scalar_type inputVals[],
                           const GlobalOrdinal inputCols[],
                           std::function<impl_scalar_type (const impl_scalar_type&, const impl_scalar_type&) > f,
                           const bool atomic = useAtomicUpdatesByDefault);

  public:
    /// \brief Transform CrsMatrix entries in place, using local
    ///   indices to select the entries in the row to transform.
    ///
    /// For every entry \f$A(i,j)\f$ to transform, if \f$v_{ij}\f$ is
    /// the corresponding entry of the \c inputVals array, then we
    /// apply the binary function f to \f$A(i,j)\f$ as follows:
    /// \f[
    ///   A(i,j) := f(A(i,j), v_{ij}).
    /// \f]
    /// For example, BinaryFunction = std::plus<impl_scalar_type> does
    /// the same thing as sumIntoLocalValues, and BinaryFunction =
    /// project2nd<impl_scalar_type,impl_scalar_type> does the same
    /// thing as replaceLocalValues.  (It is generally more efficient
    /// to call sumIntoLocalValues resp. replaceLocalValues than to do
    /// this.)
    ///
    /// This overload of the method takes the column indices and
    /// values as Kokkos::View.  See below for an overload that takes
    /// Teuchos::ArrayView instead.
    ///
    /// \tparam LocalIndicesViewType Kokkos::View specialization that
    ///   is a 1-D array of LocalOrdinal.
    /// \tparam ImplScalarViewType Kokkos::View specialization that is
    ///   a 1-D array of impl_scalar_type (usually the same as Scalar,
    ///   unless Scalar is std::complex<T> for some T, in which case
    ///   it is Kokkos::complex<T>).
    /// \tparam BinaryFunction The type of the binary function f to
    ///   use for updating the sparse matrix's value(s).  This should
    ///   be convertible to
    ///   std::function<impl_scalar_type (const impl_scalar_type&,
    ///                                   const impl_scalar_type&)>.
    ///
    /// \param lclRow [in] (Local) index of the row to modify.
    ///   This row <i>must</t> be owned by the calling process.  (This
    ///   is a stricter requirement than for sumIntoGlobalValues.)
    /// \param inputInds [in] (Local) indices in the row to modify.
    ///   Indices not in the row on the calling process, and their
    ///   corresponding values, will be ignored.
    /// \param inputVals [in] Values to use for modification.
    /// \param f [in] The binary function to use for updating the
    ///   sparse matrix's value.  It takes two \c impl_scalar_type
    ///   values and returns \c impl_scalar_type.
    /// \param atomic [in] Whether to use atomic updates.
    template<class LocalIndicesViewType,
             class ImplScalarViewType,
             class BinaryFunction>
    LocalOrdinal
    transformLocalValues (const LocalOrdinal lclRow,
                          const typename UnmanagedView<LocalIndicesViewType>::type& inputInds,
                          const typename UnmanagedView<ImplScalarViewType>::type& inputVals,
                          BinaryFunction f,
                          const bool atomic = useAtomicUpdatesByDefault)
    {
      // We use static_assert here to check the template parameters,
      // rather than std::enable_if (e.g., on the return value, to
      // enable compilation only if the template parameters match the
      // desired attributes).  This turns obscure link errors into
      // clear compilation errors.  It also makes the return value a
      // lot easier to see.
      static_assert (Kokkos::is_view<LocalIndicesViewType>::value,
                     "First template parameter LocalIndicesViewType must be "
                     "a Kokkos::View.");
      static_assert (Kokkos::is_view<ImplScalarViewType>::value,
                     "Second template parameter ImplScalarViewType must be a "
                     "Kokkos::View.");
      static_assert (static_cast<int> (LocalIndicesViewType::rank) == 1,
                     "First template parameter LocalIndicesViewType must "
                     "have rank 1.");
      static_assert (static_cast<int> (ImplScalarViewType::rank) == 1,
                     "Second template parameter ImplScalarViewType must have "
                     "rank 1.");
      static_assert (std::is_same<
                       typename LocalIndicesViewType::non_const_value_type,
                       local_ordinal_type>::value,
                     "First template parameter LocalIndicesViewType must "
                     "contain values of type local_ordinal_type.");
      static_assert (std::is_same<
                       typename ImplScalarViewType::non_const_value_type,
                       impl_scalar_type>::value,
                     "Second template parameter ImplScalarViewType must "
                     "contain values of type impl_scalar_type.");
      typedef LocalOrdinal LO;
      const LO numInputEnt = inputInds.extent (0);
      if (static_cast<LO> (inputVals.extent (0)) != numInputEnt) {
        return Teuchos::OrdinalTraits<LO>::invalid ();
      }
      return this->transformLocalValues (lclRow,
                                         numInputEnt,
                                         inputVals.data (),
                                         inputInds.data (),
                                         f,
                                         atomic);
    }

    /// \brief Transform CrsMatrix entries in place, using global
    ///   indices to select the entries in the row to transform.
    ///
    /// For every entry \f$A(i,j)\f$ to transform, if \f$v_{ij}\f$ is
    /// the corresponding entry of the \c inputVals array, then we
    /// apply the binary function f to \f$A(i,j)\f$ as follows:
    /// \f[
    ///   A(i,j) := f(A(i,j), v_{ij}).
    /// \f]
    /// For example, BinaryFunction = std::plus<impl_scalar_type> does
    /// the same thing as sumIntoLocalValues, and BinaryFunction =
    /// project2nd<impl_scalar_type,impl_scalar_type> does the same
    /// thing as replaceLocalValues.  (It is generally more efficient
    /// to call sumIntoLocalValues resp. replaceLocalValues than to do
    /// this.)
    ///
    /// \tparam BinaryFunction The type of the binary function f to
    ///   use for updating the sparse matrix's value(s).  This should
    ///   be convertible to
    ///   std::function<impl_scalar_type (const impl_scalar_type&,
    ///                                   const impl_scalar_type&)>.
    /// \tparam InputMemorySpace Kokkos memory space / device in which
    ///   the input data live.  This may differ from the memory space
    ///   in which the current matrix's row's values live.
    ///
    /// \param gblRow [in] (Global) index of the row to modify.  This
    ///   row <i>must</t> be owned by the calling process.  (This is a
    ///   stricter requirement than for sumIntoGlobalValues.)
    /// \param inputInds [in] (Global) indices in the row to modify.
    ///   Indices not in the row on the calling process, and their
    ///   corresponding values, will be ignored.
    /// \param inputVals [in] Values to use for modification.
    /// \param f [in] The binary function to use for updating the
    ///   sparse matrix's value.  It takes two \c impl_scalar_type
    ///   values and returns \c impl_scalar_type.
    /// \param atomic [in] Whether to use atomic updates.
    ///
    /// This method works whether indices are local or global.
    /// However, it will cost more if indices are local, since it will
    /// have to convert the input global indices to local indices in
    /// that case.
    template<class BinaryFunction, class InputMemorySpace>
    LocalOrdinal
    transformGlobalValues (const GlobalOrdinal gblRow,
                           const Kokkos::View<const GlobalOrdinal*,
                             InputMemorySpace,
                             Kokkos::MemoryUnmanaged>& inputInds,
                           const Kokkos::View<const impl_scalar_type*,
                             InputMemorySpace,
                             Kokkos::MemoryUnmanaged>& inputVals,
                           BinaryFunction f,
                           const bool atomic = useAtomicUpdatesByDefault)
    {
      typedef LocalOrdinal LO;
      const LO numInputEnt = inputInds.extent (0);
      if (static_cast<LO> (inputVals.extent (0)) != numInputEnt) {
        return Teuchos::OrdinalTraits<LO>::invalid ();
      }
      return this->transformGlobalValues (gblRow,
                                          numInputEnt,
                                          inputVals.data (),
                                          inputInds.data (),
                                          f,
                                          atomic);
    }

    //! Set all matrix entries equal to \c alpha.
    void setAllToScalar (const Scalar& alpha);

    //! Scale the matrix's values: <tt>this := alpha*this</tt>.
    void scale (const Scalar& alpha);

    /// \brief Set the local matrix using three (compressed sparse row) arrays.
    ///
    /// \pre ind is sorted within each row
    /// \pre <tt>hasColMap() == true</tt>
    /// \pre <tt>getGraph() != Teuchos::null</tt>
    /// \pre No insert/sum routines have been called
    ///
    /// \warning This is for EXPERT USE ONLY.  We make NO PROMISES of
    ///   backwards compatibility.
    ///
    /// This method behaves like the CrsMatrix constructor that takes
    /// a const CrsGraph.  It fixes the matrix's graph, but does not
    /// call fillComplete on the matrix.  The graph might not
    /// necessarily be fill complete, but it must have a local graph.
    ///
    /// The input arguments might be used directly (shallow copy), or
    /// they might be (deep) copied.
    ///
    /// \param ptr [in] Array of row offsets.
    /// \param ind [in] Array of (local) column indices.
    /// \param val [in/out] Array of values.  This is in/out because
    ///   the matrix reserves the right to take this argument by
    ///   shallow copy.  Any method that changes the matrix's values
    ///   may then change this.
    void
    setAllValues (const typename local_graph_device_type::row_map_type& ptr,
                  const typename local_graph_device_type::entries_type::non_const_type& ind,
                  const typename local_matrix_device_type::values_type& val);

    /// \brief Set the local matrix using an existing local matrix.
    ///
    /// \pre column indices are sorted within each row
    /// \pre <tt>hasColMap() == true</tt>
    /// \pre <tt>getGraph() != Teuchos::null</tt>
    /// \pre No insert/sum routines have been called
    ///
    /// \warning This is for EXPERT USE ONLY.  We make NO PROMISES of
    ///   backwards compatibility.
    ///
    /// This method simply calls the method setAllValues that accepts three
    /// compressed sparse row arrays.
    ///
    /// The input argument might be used directly (shallow copy), or
    /// it might be (deep) copied.
    ///
    /// \param ptr [in/out] Kokkos sparse local matrix
    ///   This is in/out because the matrix reserves the right to take this argument by
    ///   shallow copy.  Any method that changes the matrix's values
    ///   may then change this.
    void
    setAllValues (const local_matrix_device_type& localMatrix);

    /// \brief Set the local matrix using three (compressed sparse row) arrays.
    ///
    /// \pre ind is sorted within each row
    /// \pre <tt>hasColMap() == true</tt>
    /// \pre <tt>getGraph() != Teuchos::null</tt>
    /// \pre No insert/sum routines have been called
    ///
    /// \warning This is for EXPERT USE ONLY.  We make NO PROMISES of
    ///   backwards compatibility.
    ///
    /// This method behaves like the CrsMatrix constructor that takes
    /// a const CrsGraph.  It fixes the matrix's graph, but does not
    /// call fillComplete on the matrix.  The graph might not
    /// necessarily be fill complete, but it must have a local graph.
    ///
    /// The input arguments might be used directly (shallow copy), or
    /// they might be (deep) copied.
    ///
    /// \param ptr [in] Array of row offsets.
    /// \param ind [in] Array of (local) column indices.
    /// \param val [in/out] Array of values.  This is in/out because
    ///   the matrix reserves the right to take this argument by
    ///   shallow copy.  Any method that changes the matrix's values
    ///   may then change this.
    void
    setAllValues (const Teuchos::ArrayRCP<size_t>& ptr,
                  const Teuchos::ArrayRCP<LocalOrdinal>& ind,
                  const Teuchos::ArrayRCP<Scalar>& val);

    /// \brief Get a host view of the CRS packed row pointers
    row_ptrs_host_view_type getLocalRowPtrsHost () const
    { return getCrsGraph()->getLocalRowPtrsHost(); }

    /// \brief Get a device view of the CRS packed row pointers
    row_ptrs_device_view_type getLocalRowPtrsDevice () const
    { return getCrsGraph()->getLocalRowPtrsDevice(); }

    /// \brief Get a host view of the CRS packed column indicies
    local_inds_host_view_type getLocalIndicesHost () const
    { return getCrsGraph()->getLocalIndicesHost(); }

    /// \brief Get a device_view of the CRS packed column indicies
    local_inds_device_view_type getLocalIndicesDevice () const
    { return getCrsGraph()->getLocalIndicesDevice(); }

    //@}
    //! @name Transformational methods
    //@{

    /// \brief Communicate nonlocal contributions to other processes.
    ///
    /// Users do not normally need to call this method.  fillComplete
    /// always calls this method, unless you specifically tell
    /// fillComplete to do otherwise by setting its "No Nonlocal
    /// Changes" parameter to \c true.  Thus, it suffices to call
    /// fillComplete.
    ///
    /// Methods like insertGlobalValues and sumIntoGlobalValues let
    /// you add or modify entries in rows that are not owned by the
    /// calling process.  These entries are called "nonlocal
    /// contributions."  The methods that allow nonlocal contributions
    /// store the entries on the calling process, until globalAssemble
    /// is called.  globalAssemble sends these nonlocal contributions
    /// to the process(es) that own them, where they then become part
    /// of the matrix.
    ///
    /// This method only does global assembly if there are nonlocal
    /// entries on at least one process.  It does an all-reduce to
    /// find that out.  If not, it returns early, without doing any
    /// more communication or work.
    ///
    /// If you previously inserted into a row which is not owned by
    /// <i>any</i> process in the row Map, the behavior of this method
    /// is undefined.  It may detect the invalid row indices and throw
    /// an exception, or it may silently drop the entries inserted
    /// into invalid rows.  Behavior may vary, depending on whether
    /// Tpetra was built with debug checking enabled.
    void globalAssemble();

    /// \brief Resume operations that may change the values or
    ///   structure of the matrix.
    ///
    /// This method must be called as a collective operation.
    ///
    /// Calling fillComplete "freezes" both the values and the
    /// structure of the matrix.  If you want to modify the matrix
    /// again, you must first call resumeFill.  You then may not call
    /// resumeFill again on that matrix until you first call
    /// fillComplete.  You may make sequences of fillComplete,
    /// resumeFill calls as many times as you wish.
    ///
    /// \post <tt>isFillActive() && ! isFillComplete()</tt>
    void resumeFill (const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Tell the matrix that you are done changing its
    ///   structure or values, and that you are ready to do
    ///   computational kernels (e.g., sparse matrix-vector multiply)
    ///   with it.
    ///
    /// This tells the graph to optimize its data structures for
    /// computational kernels, and to prepare (MPI) communication
    /// patterns.
    ///
    /// Off-process indices are distributed (via globalAssemble()),
    /// indices are sorted, redundant indices are fused, and global
    /// indices are transformed to local indices.
    ///
    /// \warning The domain Map and row Map arguments to this method
    ///   MUST be one to one!  If you have Maps that are not one to
    ///   one, and you do not know how to make a Map that covers the
    ///   same global indices but <i>is</i> one to one, then you may
    ///   call Tpetra::createOneToOne() (see Map's header file) to
    ///   make a one-to-one version of your Map.
    ///
    /// \pre  <tt>   isFillActive() && ! isFillComplete() </tt>
    /// \post <tt> ! isFillActive() &&   isFillComplete() </tt>
    ///
    /// \param domainMap [in] The matrix's domain Map.  MUST be one to
    ///   one!
    /// \param rangeMap [in] The matrix's range Map.  MUST be one to
    ///   one!  May be, but need not be, the same as the domain Map.
    /// \param params [in/out] List of parameters controlling this
    ///   method's behavior.  See below for valid parameters.
    ///
    /// List of valid parameters in <tt>params</tt>:
    /// <ul>
    /// <li> "No Nonlocal Changes" (\c bool): Default is false.  If
    ///      true, the caller promises that no modifications to
    ///      nonowned rows have happened on any process since the last
    ///      call to fillComplete.  This saves a global all-reduce to
    ///      check whether any process did a nonlocal insert.
    ///      Nonlocal changes include any sumIntoGlobalValues or
    ///      insertGlobalValues call with a row index that is not in
    ///      the row Map of the calling process.
    /// </li>
    ///
    /// <li> "Sort column Map ghost GIDs" (\c bool): Default is true.
    ///      makeColMap() (which fillComplete may call) always groups
    ///      remote GIDs by process rank, so that all remote GIDs with
    ///      the same owning rank occur contiguously.  By default, it
    ///      always sorts remote GIDs in increasing order within those
    ///      groups.  This behavior differs from Epetra, which does
    ///      not sort remote GIDs with the same owning process.  If
    ///      you don't want to sort (for compatibility with Epetra),
    ///      set this parameter to \c false.  This parameter only
    ///      takes effect if the matrix owns the graph.  This is an
    ///      expert mode parameter ONLY.  We make no promises about
    ///      backwards compatibility of this parameter.  It may change
    ///      or disappear at any time.
    /// </li>
    /// </ul>
    void
    fillComplete (const Teuchos::RCP<const map_type>& domainMap,
                  const Teuchos::RCP<const map_type>& rangeMap,
                  const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Tell the matrix that you are done changing its
    ///   structure or values, and that you are ready to do
    ///   computational kernels (e.g., sparse matrix-vector multiply)
    ///   with it.  Set default domain and range Maps.
    ///
    /// See above three-argument version of fillComplete for full
    /// documentation.  If the matrix does not yet have domain and
    /// range Maps (i.e., if fillComplete has not yet been called on
    /// this matrix at least once), then this method uses the matrix's
    /// row Map (result of this->getRowMap()) as both the domain Map
    /// and the range Map.  Otherwise, this method uses the matrix's
    /// existing domain and range Maps.
    ///
    /// \warning It is only valid to call this overload of
    ///   fillComplete if the row Map is one to one!  If the row Map
    ///   is NOT one to one, you must call the above three-argument
    ///   version of fillComplete, and supply one-to-one domain and
    ///   range Maps.  If you have Maps that are not one to one, and
    ///   you do not know how to make a Map that covers the same
    ///   global indices but <i>is</i> one to one, then you may call
    ///   Tpetra::createOneToOne() (see Map's header file) to make a
    ///   one-to-one version of your Map.
    ///
    /// \param params [in/out] List of parameters controlling this
    ///   method's behavior.  See documentation of the three-argument
    ///   version of fillComplete (above) for valid parameters.
    void
    fillComplete (const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Perform a fillComplete on a matrix that already has data.
    ///
    /// The matrix must already have filled local 1-D storage
    /// (k_clInds1D_ and k_rowPtrs_ for the graph, and k_values1D_ in
    /// the matrix).  If the matrix has been constructed in any other
    /// way, this method will throw an exception.  This routine is
    /// needed to support other Trilinos packages and should not be
    /// called by ordinary users.
    ///
    /// \warning This method is intended for expert developer use
    ///   only, and should never be called by user code.
    ///
    /// \param domainMap [in] The matrix's domain Map.  MUST be one to
    ///   one!
    /// \param rangeMap [in] The matrix's range Map.  MUST be one to
    ///   one!  May be, but need not be, the same as the domain Map.
    /// \param importer [in] Import from the matrix's domain Map to
    ///   its column Map.  If no Import is necessary (i.e., if the
    ///   domain and column Maps are the same, in the sense of
    ///   Tpetra::Map::isSameAs), then this may be Teuchos::null.
    /// \param exporter [in] Export from the matrix's row Map to its
    ///   range Map.  If no Export is necessary (i.e., if the row and
    ///   range Maps are the same, in the sense of
    ///   Tpetra::Map::isSameAs), then this may be Teuchos::null.
    /// \param params [in/out] List of parameters controlling this
    ///   method's behavior.
    void
    expertStaticFillComplete (const Teuchos::RCP<const map_type>& domainMap,
                              const Teuchos::RCP<const map_type>& rangeMap,
                              const Teuchos::RCP<const import_type>& importer = Teuchos::null,
                              const Teuchos::RCP<const export_type>& exporter = Teuchos::null,
                              const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Replace the matrix's column Map with the given Map.
    ///
    /// \param newColMap [in] New column Map.  Must be nonnull.
    ///   Within Tpetra, there are no particular restrictions on the column map.
    ///   However, if this graph will be used in Xpetra, Ifpack2, or MueLu,
    ///   the column map's list of global indices must follow "Aztec ordering":
    ///   locally owned GIDs (same order as in domain map), followed by remote GIDs
    ///   (in order of owning proc, and sorted within each proc).
    ///
    ///   It is strongly recommended to use \c Tpetra::Details::makeColMap()
    ///   to create the column map. makeColMap() follows Aztec ordering by default.
    ///
    /// \pre The matrix must have no entries inserted yet, on any
    ///   process in the row Map's communicator.
    ///
    /// \pre The matrix must not have been created with a constant
    ///   (a.k.a. "static") CrsGraph.
    void
    replaceColMap (const Teuchos::RCP<const map_type>& newColMap);

    /// \brief Reindex the column indices in place, and replace the
    ///   column Map.  Optionally, replace the Import object as well.
    ///
    /// \pre The matrix is <i>not</i> fill complete:
    ///   <tt>! this->isFillComplete() </tt>.
    /// \pre Either the input graph is \c NULL, or it is <i>not</i>
    ///   fill complete:
    ///   <tt>graph == NULL || ! graph->isFillComplete()</tt>.
    /// \pre On every calling process, every index owned by the
    ///   current column Map must also be owned by the new column Map.
    /// \pre If the new Import object is provided, the new Import
    ///   object's source Map must be the same as the current domain
    ///   Map, and the new Import's target Map must be the same as the
    ///   new column Map.
    ///
    /// \param graph [in] The matrix's graph.  If you don't provide
    ///   this (i.e., if <tt>graph == NULL</tt>), then the matrix must
    ///   own its graph, which will be modified in place.  (That is,
    ///   you must <i>not</i> have created the matrix with a constant
    ///   graph.)  If you <i>do</i> provide this, then the method will
    ///   assume that it is the same graph as the matrix's graph, and
    ///   the provided graph will be modified in place.
    /// \param newColMap [in] New column Map.  Must be nonnull.
    /// \param newImport [in] New Import object.  Optional; computed
    ///   if not provided or if null.  Computing an Import is
    ///   expensive, so it is worth providing this if you can.
    /// \param sortEachRow [in] If true, sort the indices (and their
    ///   corresponding values) in each row after reindexing.
    ///
    /// Why would you want to use this method?  Well, for example, you
    /// might need to use an Ifpack2 preconditioner that only accepts
    /// a matrix with a certain kind of column Map.  Your matrix has
    /// the wrong kind of column Map, but you know how to compute the
    /// right kind of column Map.  You might also know an efficient
    /// way to compute an Import object from the current domain Map to
    /// the new column Map.  (For an instance of the latter, see the
    /// Details::makeOptimizedColMapAndImport function in
    /// Tpetra_Details_makeOptimizedColMap.hpp.)
    ///
    /// Suppose that you created this CrsMatrix with a constant graph;
    /// that is, that you called the CrsMatrix constructor that takes
    /// a CrsGraph as input:
    ///
    /// \code
    /// RCP<CrsGraph<> > G (new CrsGraph<> (rowMap, origColMap, ...));
    /// // ... fill G ...
    /// G->fillComplete (domMap, ranMap);
    /// CrsMatrix<> A (G);
    /// // ... fill A ...
    /// \endcode
    ///
    /// Now suppose that you want to give A to a preconditioner that
    /// can't handle a matrix with an arbitrary column Map (in the
    /// example above, <tt>origColMap</tt>).  You first must create a
    /// new suitable column Map <tt>newColMap</tt>, and optionally a
    /// new Import object <tt>newImport</tt> from the matrix's current
    /// domain Map to the new column Map.  Then, call this method,
    /// passing in G (which must <i>not</i> be fill complete) while
    /// the matrix is <i>not</i> fill complete.  Be sure to save the
    /// graph's <i>original</i> Import object; you'll need that later.
    ///
    /// \code
    /// RCP<const CrsGraph<>::import_type> origImport = G->getImporter ();
    /// G->resumeFill ();
    /// A.reindexColumns (G.getRawPtr (), newColMap, newImport);
    /// G.fillComplete (domMap, ranMap);
    /// A.fillComplete (domMap, ranMap);
    /// \endcode
    ///
    /// Now you may give the matrix A to the preconditioner in
    /// question.  After doing so, and after you solve the linear
    /// system using the preconditioner, you might want to put the
    /// matrix back like it originally was.  You can do that, too!
    ///
    /// \code
    /// A.resumeFill ();
    /// G->resumeFill ();
    /// A.reindexColumns (G.getRawPtr (), origColMap, origImport);
    /// G->fillComplete (domMap, ranMap);
    /// A->fillComplete (domMap, ranMap);
    /// \endcode
    void
    reindexColumns (crs_graph_type* const graph,
                    const Teuchos::RCP<const map_type>& newColMap,
                    const Teuchos::RCP<const import_type>& newImport = Teuchos::null,
                    const bool sortEachRow = true);

    /// \brief Replace the current domain Map with the given objects.
    ///
    /// The matrix's Import object will be recomputed if needed.
    ///
    /// \param newDomainMap [in] New domain Map.  Must be nonnull.
    ///
    /// \pre The matrix must be fill complete:
    ///   <tt>isFillComplete() == true</tt>.
    /// 
    void
    replaceDomainMap (const Teuchos::RCP<const map_type>& newDomainMap);

    /// \brief Replace the current domain Map and Import with the given objects.
    ///
    /// \param newDomainMap [in] New domain Map.  Must be nonnull.
    /// \param newImporter [in] Optional Import object.  If null, the new Domain Map must equal the matrix's Column Map
    ///
    /// \pre The matrix must be fill complete:
    ///   <tt>isFillComplete() == true</tt>.
    /// \pre If the Import is provided, its target Map must be the
    ///   same as the column Map of the matrix.
    /// \pre If the Import is provided, its source Map must be the
    ///   same as the provided new domain Map.
    /// \pre If the Import is not provided, the new Domain Map must be the
    ///   same as the matrix's Column Map.
    void
    replaceDomainMapAndImporter (const Teuchos::RCP<const map_type>& newDomainMap,
                                 Teuchos::RCP<const import_type>& newImporter);

    /// \brief Replace the current range Map with the given objects.
    ///
    /// The matrix's Export object will be recomputed if needed.
    ///
    /// \param newRangeMap [in] New Range Map.  Must be nonnull.
    ///
    /// \pre The matrix must be fill complete:
    ///   <tt>isFillComplete() == true</tt>.
    /// 
    void
    replaceRangeMap (const Teuchos::RCP<const map_type>& newRangeMap);

    /// \brief Replace the current Range Map and Export with the given objects.
    ///
    /// \param newRangeMap [in] New domain Map.  Must be nonnull.
    /// \param newExporter [in] Optional Export object.  If null, the new Range Map must equal the matrix's Row Map
    ///
    /// \pre The matrix must be fill complete:
    ///   <tt>isFillComplete() == true</tt>.
    /// \pre If the Export is provided, its target Map must be the
    ///   same as the new Range Map of the matrix.
    /// \pre If the Export is provided, its source Map must be the
    ///   same as the Row Map of the matrix
    /// \pre If the Export is not provided, the new Range Map must be the
    ///   same as the matrix's Row Map.
    void
    replaceRangeMapAndExporter (const Teuchos::RCP<const map_type>& newRangeMap,
                                Teuchos::RCP<const export_type>& newExporter);

    /// \brief Remove processes owning zero rows from the Maps and their communicator.
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

    //@}
    //! @name Methods implementing RowMatrix
    //@{

    //! The communicator over which the matrix is distributed.
    Teuchos::RCP<const Teuchos::Comm<int> > getComm() const override;


    //! The Map that describes the row distribution in this matrix.
    Teuchos::RCP<const map_type> getRowMap () const override;

    //! The Map that describes the column distribution in this matrix.
    Teuchos::RCP<const map_type> getColMap () const override;

    //! This matrix's graph, as a RowGraph.
    Teuchos::RCP<const RowGraph<LocalOrdinal, GlobalOrdinal, Node> >
    getGraph () const override;

    //! This matrix's graph, as a CrsGraph.
    Teuchos::RCP<const crs_graph_type> getCrsGraph () const;

  private:
    /// \brief Const reference to this matrix's graph, as a CrsGraph.
    ///
    /// This is a thread-safe version of getCrsGraph() (see above).
    /// Teuchos::RCP's copy constructor, assignment operator
    /// (operator=), and destructor are not currently thread safe (as
    /// of 17 May 2017).  Thus, if we want to write
    /// host-thread-parallel code, it's important to avoid creating or
    /// destroying Teuchos::RCP instances.  This method lets CrsMatrix
    /// access its graph, without creating an Teuchos::RCP instance
    /// (as the return value of getCrsGraph() does do).
    const crs_graph_type& getCrsGraphRef () const;

  public:
#if __armclang_major__ == 22 && __armclang_minor__ == 1
   // On Stria, PR 13052 caused a 25% performance regression in the
   // CGSolve performance test that is fixed by forcing
   // getLocalMatrixDevice to always be inlined. Restrict the fix
   // to the specific toolchain where the problem was observed
#define TPETRA_DETAILS_ALWAYS_INLINE __attribute__((always_inline))
#else
#define TPETRA_DETAILS_ALWAYS_INLINE
#endif
     /// \brief The local sparse matrix.
     ///
     /// \warning It is only valid to call this method under certain
    /// \brief The local sparse matrix.
    ///
    /// \warning It is only valid to call this method under certain
    ///   circumstances.  In particular, either the CrsMatrix must
    ///   have been created with a \c local_matrix_type object, or
    ///   fillComplete must have been called on this CrsMatrix at
    ///   least once.  This method will do no error checking, so you
    ///   are responsible for knowing when it is safe to call this
    ///   method.
    TPETRA_DETAILS_ALWAYS_INLINE local_matrix_device_type
    getLocalMatrixDevice () const;
    local_matrix_host_type getLocalMatrixHost () const;
#undef TPETRA_DETAILS_ALWAYS_INLINE

#if KOKKOSKERNELS_VERSION < 40299
    /// \brief The local sparse matrix operator 
    ///   (a wrapper of \c getLocalMatrixDevice()
    ///   that supports local matrix-vector multiply)
    ///
    /// \warning It is only valid to call this method if this->isFillComplete().
    std::shared_ptr<local_multiply_op_type> getLocalMultiplyOperator () const;
#endif

    /// \brief Number of global elements in the row map of this matrix.
    ///
    /// This is <it>not</it> the number of rows in the matrix as a
    /// mathematical object.  This method returns the global sum of
    /// the number of local elements in the row map on each processor,
    /// which is the row map's getGlobalNumElements().  Since the row
    /// map is not one-to-one in general, that global sum could be
    /// different than the number of rows in the matrix.  If you want
    /// the number of rows in the matrix, ask the range map for its
    /// global number of elements, using the following code:
    /// <code>
    /// global_size_t globalNumRows = getRangeMap()->getGlobalNumElements();
    /// </code>
    /// This method retains the behavior of Epetra, which also asks
    /// the row map for the global number of rows, rather than asking
    /// the range map.
    ///
    /// \warning Undefined if isFillActive().
    ///
    global_size_t getGlobalNumRows() const override;

    /// \brief The number of global columns in the matrix.
    ///
    /// This equals the number of entries in the matrix's domain Map.
    ///
    /// \warning Undefined if isFillActive().
    global_size_t getGlobalNumCols() const override;

    /// \brief The number of matrix rows owned by the calling process.
    ///
    /// Note that the sum of all the return values over all processes
    /// in the row Map's communicator does not necessarily equal the
    /// global number of rows in the matrix, if the row Map is
    /// overlapping.
    size_t getLocalNumRows() const override;

    /// \brief The number of columns connected to the locally owned rows of this matrix.
    ///
    /// Throws std::runtime_error if <tt>! hasColMap ()</tt>.
    size_t getLocalNumCols() const override;

    //! The index base for global indices for this matrix.
    GlobalOrdinal getIndexBase() const override;

    //! The global number of entries in this matrix.
    global_size_t getGlobalNumEntries() const override;

    //! The local number of entries in this matrix.
    size_t getLocalNumEntries() const override;

    /// \brief Number of entries in the sparse matrix in the given
    ///   global row, on the calling (MPI) process.
    ///
    /// \return <tt>OrdinalTraits<size_t>::invalid()</tt>if the
    ///   specified global row index is invalid on the calling
    ///   process, else the number of entries in the given row.
    size_t getNumEntriesInGlobalRow (GlobalOrdinal globalRow) const override;

    /// \brief Number of entries in the sparse matrix in the given
    ///   local row, on the calling (MPI) process.
    ///
    /// \return <tt>OrdinalTraits<size_t>::invalid()</tt>if the
    ///   specified local row index is invalid on the calling process,
    ///   else the number of entries in the given row.
    size_t getNumEntriesInLocalRow (local_ordinal_type localRow) const override;

    /// \brief Maximum number of entries in any row of the matrix,
    ///   over all processes in the matrix's communicator.
    ///
    /// \pre <tt>! isFillActive()</tt>
    ///
    /// This method only uses the matrix's graph.  Explicitly stored
    /// zeros count as "entries."
    size_t getGlobalMaxNumRowEntries () const override;

    /// \brief Maximum number of entries in any row of the matrix,
    ///   on this process.
    ///
    /// \pre <tt>! isFillActive()</tt>
    ///
    /// This method only uses the matrix's graph.  Explicitly stored
    /// zeros count as "entries."
    size_t getLocalMaxNumRowEntries () const override;

    //! The number of degrees of freedom per mesh point.
    virtual LocalOrdinal getBlockSize () const override { return 1; }

    //! Whether the matrix has a well-defined column Map.
    bool hasColMap () const override;


    /// \brief Whether the matrix is locally indexed on the calling
    ///   process.
    ///
    /// The matrix is locally indexed on the calling process if and
    /// only if all of the following hold:
    /// <ol>
    /// <li> The matrix is not empty on the calling process </li>
    /// <li> The matrix has a column Map </li>
    /// </ol>
    ///
    /// The following is always true:
    /// \code
    /// (! locallyIndexed() && ! globallyIndexed()) || (locallyIndexed() || globallyIndexed());
    /// \endcode
    /// That is, a matrix may be neither locally nor globally indexed,
    /// but it can never be both.  Furthermore a matrix that is not
    /// fill complete, might have some processes that are neither
    /// locally nor globally indexed, and some processes that are
    /// globally indexed.  The processes that are neither do not have
    /// any entries.
    bool isLocallyIndexed() const override;

    /// \brief Whether the matrix is globally indexed on the calling process.
    ///
    /// The matrix is globally indexed on the calling process if and
    /// only if all of the following hold:
    /// <ol>
    /// <li> The matrix is not empty on the calling process </li>
    /// <li> The matrix does not yet have a column Map </li>
    /// </ol>
    ///
    /// The following is always true:
    /// \code
    /// (! locallyIndexed() && ! globallyIndexed()) ||
    ///   (locallyIndexed() || globallyIndexed());
    /// \endcode
    /// That is, a matrix may be neither locally nor globally indexed,
    /// but it can never be both.  Furthermore a matrix that is not
    /// fill complete, might have some processes that are neither
    /// locally nor globally indexed, and some processes that are
    /// globally indexed.  The processes that are neither do not have
    /// any entries.
    bool isGloballyIndexed() const override;

    /// \brief Whether the matrix is fill complete.
    ///
    /// A matrix is <i>fill complete</i> (or "in compute mode") when
    /// fillComplete() has been called without an intervening call to
    /// resumeFill().  A matrix must be fill complete in order to call
    /// computational kernels like sparse matrix-vector multiply and
    /// sparse triangular solve.  A matrix must be <i>not</i> fill
    /// complete ("in edit mode") in order to call methods that
    /// insert, modify, or remove entries.
    ///
    /// The following are always true:
    /// <ul>
    /// <li> <tt> isFillActive() == ! isFillComplete() </tt>
    /// <li> <tt> isFillActive() || isFillComplete() </tt>
    /// </ul>
    ///
    /// A matrix starts out (after its constructor returns) as not
    /// fill complete.  It becomes fill complete after fillComplete()
    /// returns, and becomes not fill complete again if resumeFill()
    /// is called.  Some methods like clone() and some of the
    /// "nonmember constructors" (like importAndFillComplete() and
    /// exportAndFillComplete()) may return a fill-complete matrix.
    bool isFillComplete() const override;

    /// \brief Whether the matrix is not fill complete.
    ///
    /// A matrix is <i>fill complete</i> (or "in compute mode") when
    /// fillComplete() has been called without an intervening call to
    /// resumeFill().  A matrix must be fill complete in order to call
    /// computational kernels like sparse matrix-vector multiply and
    /// sparse triangular solve.  A matrix must be <i>not</i> fill
    /// complete ("in edit mode") in order to call methods that
    /// insert, modify, or remove entries.
    ///
    /// The following are always true:
    /// <ul>
    /// <li> <tt> isFillActive() == ! isFillComplete() </tt>
    /// <li> <tt> isFillActive() || isFillComplete() </tt>
    /// </ul>
    ///
    /// A matrix starts out (after its constructor returns) as not
    /// fill complete.  It becomes fill complete after fillComplete()
    /// returns, and becomes not fill complete again if resumeFill()
    /// is called.  Some methods like clone() and some of the
    /// "nonmember constructors" (like importAndFillComplete() and
    /// exportAndFillComplete()) may return a fill-complete matrix.
    bool isFillActive() const;

    //! \brief Returns \c true if storage has been optimized.
    /**
       Optimized storage means that the allocation of each row is equal to the
       number of entries. The effect is that a pass through the matrix, i.e.,
       during a mat-vec, requires minimal memory traffic. One limitation of
       optimized storage is that no new indices can be added to the matrix.
    */
    bool isStorageOptimized () const;

    //! Indicates that the graph is static, so that new entries cannot be added to this matrix.
    bool isStaticGraph () const;

    /// \brief Compute and return the Frobenius norm of the matrix.
    ///
    /// The Frobenius norm of the matrix is defined as
    /// \f\[
    ///   \|A\|_F = \sqrt{\sum_{i,j} \|A(i,j)\|^2}.
    /// \f\].
    ///
    mag_type getFrobeniusNorm () const override;

    /// \brief Return \c true if getLocalRowView() and
    ///   getGlobalRowView() are valid for this object.
    virtual bool supportsRowViews () const override;

protected:
    using values_dualv_type =
          Kokkos::DualView<impl_scalar_type*, device_type>;
    using values_wdv_type = 
          Details::WrappedDualView<values_dualv_type>;
    values_wdv_type valuesUnpacked_wdv;
    mutable values_wdv_type valuesPacked_wdv;

#if KOKKOSKERNELS_VERSION < 40299
    using ordinal_rowptrs_type = typename local_multiply_op_type::ordinal_view_type;
    /// \brief local_ordinal typed version of local matrix's rowptrs.
    ///   This allows the LocalCrsMatrixOperator to have rowptrs and entries be the same type,
    ///   so cuSPARSE SpMV (including merge-path) can be used for apply.
    ///   This is allocated and populated lazily in getLocalMultiplyOperator(), only if all 4 conditions are met:
    ///     - node_type is KokkosCudaWrapperNode
    ///     - the cuSPARSE TPL is enabled
    ///     - local_ordinal_type can represent getLocalNumEntries()
    mutable ordinal_rowptrs_type ordinalRowptrs;
#endif

public:

    /// \brief Fill given arrays with a deep copy of the locally owned
    ///   entries of the matrix in a given row, using global column
    ///   indices.
    ///
    /// \param GlobalRow [in] Global index of the row for which to
    ///   return entries.
    /// \param Indices [out] Global column indices corresponding to
    ///   values.
    /// \param Values [out] Matrix values.
    /// \param NumEntries [out] Number of entries.
    ///
    /// \note To Tpetra developers: Discussion of whether to use
    ///   <tt>Scalar</tt> or <tt>impl_scalar_type</tt> for output
    ///   array of matrix values.
    ///
    /// If \c Scalar differs from <tt>impl_scalar_type</tt>, as for
    /// example with std::complex<T> and Kokkos::complex<T>, we must
    /// choose which type to use.  We must make the same choice as
    /// RowMatrix does, else CrsMatrix won't compile, because it won't
    /// implement a pure virtual method.  We choose <tt>Scalar</tt>,
    /// for the following reasons.  First, <tt>Scalar</tt> is the
    /// user's preferred type, and <tt>impl_scalar_type</tt> an
    /// implementation detail that makes Tpetra work with Kokkos.
    /// Second, Tpetra's public interface provides a host-only
    /// interface, which eliminates some reasons for requiring
    /// implementation-specific types like Kokkos::complex.
    ///
    /// We do eventually want to put Tpetra methods in Kokkos kernels,
    /// but we only <i>need</i> to put them in host kernels, since
    /// Tpetra is a host-only interface.  Users can still manually
    /// handle conversion from <tt>Scalar</tt> to
    /// <tt>impl_scalar_type</tt> for reductions.
    ///
    /// The right thing to do would be to rewrite RowMatrix so that
    /// getGlobalRowCopy is NOT inherited, but is implemented by a
    /// pure virtual "hook" getGlobalRowCopyImpl.  The latter takes
    /// raw pointers.  That would give us the freedom to overload
    /// getGlobalRowCopy, which one normally can't do with virtual
    /// methods.  It would make sense for one getGlobalRowCopyImpl
    /// method to implement both Teuchos::ArrayView and Kokos::View
    /// versions of getGlobalRowCopy.
    ///
    /// Note: A std::runtime_error exception is thrown if either
    /// <tt>Indices</tt> or <tt>Values</tt> is not large enough to
    /// hold the data associated with row \c GlobalRow. If row
    /// <tt>GlobalRow</tt> is not owned by the calling process, then
    /// \c Indices and \c Values are unchanged and \c NumIndices is
    /// returned as Teuchos::OrdinalTraits<size_t>::invalid().
    void
    getGlobalRowCopy (GlobalOrdinal GlobalRow,
                      nonconst_global_inds_host_view_type &Indices,
                      nonconst_values_host_view_type &Values,
                      size_t& NumEntries) const override;
    /// \brief Fill given arrays with a deep copy of the locally owned
    ///   entries of the matrix in a given row, using local column
    ///   indices.
    ///
    /// \param LocalRow   [in]  Local index of the row for which to return entries.
    /// \param Indices    [out] Local column indices corresponding to values.
    /// \param Values       [out] Matrix values.
    /// \param NumEntries [out] Number of entries returned.
    ///
    /// Note: A std::runtime_error exception is thrown if either
    /// <tt>colInds</tt> or \c vals is not large enough to hold the
    /// data associated with row \c localRow. If row \c localRow is
    /// not owned by the calling process, then <tt>colInds</tt> and
    /// <tt>vals</tt> are unchanged and <tt>numEntries</tt> is
    /// returned as Teuchos::OrdinalTraits<size_t>::invalid().
    void
    getLocalRowCopy (LocalOrdinal LocalRow,
                     nonconst_local_inds_host_view_type &Indices,
                     nonconst_values_host_view_type &Values,
                     size_t& NumEntries) const override;

    /// \brief Get a constant, nonpersisting view of a row of this
    ///   matrix, using global row and column indices.
    ///
    /// \param GlobalRow [in]  Global index of the row to view.
    /// \param indices   [out] On output: view of the global column indices in the row.
    /// \param values    [out] On output: view of the values in the row.
    ///
    /// \pre <tt>isLocallyIndexed () == false</tt>
    /// \post <tt>indices.size () == this->getNumEntriesInGlobalRow (GlobalRow)</tt>
    ///
    /// If \c GlobalRow is not a valid global row index on the calling
    /// process, then \c indices is set to null.

    void
    getGlobalRowView (GlobalOrdinal GlobalRow,
                      global_inds_host_view_type &indices,
                      values_host_view_type &values) const override;

    /// \brief Get a constant view of a row of this
    ///   matrix, using local row and column indices.
    ///
    /// \param LocalRow [in]  Local index of the row to view.
    /// \param indices  [out] On output: view of the local column indices in the row.
    /// \param values   [out] On output: view of the values in the row.
    ///
    /// \pre <tt>isGloballyIndexed () == false</tt>
    /// \post <tt>indices.size () == this->getNumEntriesInLocalRow (LocalRow)</tt>
    ///
    /// If \c LocalRow is not a valid local row index on the calling
    /// process, then \c indices is set to null.
    void
    getLocalRowView(LocalOrdinal LocalRow,
                    local_inds_host_view_type &indices,
                    values_host_view_type &values) const override;

    /// \brief Get a constant, nonpersisting view of a row of this
    ///   matrix, using local row and column indices, with raw
    ///   pointers.
    ///
    /// This overload exists only if Scalar differs from
    /// impl_scalar_type.  In that case, this overload takes a Scalar
    /// pointer.

    /// \brief Get a copy of the diagonal entries of the matrix.
    ///
    /// This method returns a Vector with the same Map as this
    /// matrix's row Map.  On each process, it contains the diagonal
    /// entries owned by the calling process. If the matrix has an empty
    /// row, the diagonal entry contains a zero.
    void
    getLocalDiagCopy (Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& diag) const override;

    /// \brief Get offsets of the diagonal entries in the matrix.
    ///
    /// \warning This method is DEPRECATED.  Call
    ///   CrsGraph::getLocalDiagOffsets, in particular the overload
    ///   that returns the offsets as a Kokkos::View.
    ///
    /// \warning This method is only for expert users.
    /// \warning We make no promises about backwards compatibility
    ///   for this method.  It may disappear or change at any time.
    /// \warning This method must be called collectively.  We reserve
    ///   the right to do extra checking in a debug build that will
    ///   require collectives.
    ///
    /// \pre The matrix must be locally indexed (which means that it
    ///   has a column Map).
    /// \pre All diagonal entries of the matrix's graph must be
    ///   populated on this process.  Results are undefined otherwise.
    /// \post <tt>offsets.size() == getLocalNumRows()</tt>
    ///
    /// This method creates an array of offsets of the local diagonal
    /// entries in the matrix.  This array is suitable for use in the
    /// two-argument version of getLocalDiagCopy().  However, its
    /// contents are not defined in any other context.  For example,
    /// you should not rely on offsets[i] being the index of the
    /// diagonal entry in the views returned by getLocalRowView().
    /// This may be the case, but it need not be.  (For example, we
    /// may choose to optimize the lookups down to the optimized
    /// storage level, in which case the offsets will be computed with
    /// respect to the underlying storage format, rather than with
    /// respect to the views.)
    ///
    /// Calling any of the following invalidates the output array:
    /// <ul>
    /// <li> insertGlobalValues() </li>
    /// <li> insertLocalValues() </li>
    /// <li> fillComplete() (with a dynamic graph) </li>
    /// <li> resumeFill() (with a dynamic graph) </li>
    /// </ul>
    ///
    /// If the matrix has a const ("static") graph, and if that graph
    /// is fill complete, then the offsets array remains valid through
    /// calls to fillComplete() and resumeFill().  "Invalidates" means
    /// that you must call this method again to recompute the offsets.
    void getLocalDiagOffsets (Teuchos::ArrayRCP<size_t>& offsets) const;

    /// \brief Variant of getLocalDiagCopy() that uses precomputed offsets.
    ///
    /// \warning This method is only for expert users.
    /// \warning We make no promises about backwards compatibility
    ///   for this method.  It may disappear or change at any time.
    ///
    /// This method uses the offsets of the diagonal entries, as
    /// precomputed by the Kokkos::View overload of
    /// getLocalDiagOffsets(), to speed up copying the diagonal of the
    /// matrix.  The offsets must be recomputed if any of the
    /// following methods are called:
    /// <ul>
    /// <li> insertGlobalValues() </li>
    /// <li> insertLocalValues() </li>
    /// <li> fillComplete() (with a dynamic graph) </li>
    /// <li> resumeFill() (with a dynamic graph) </li>
    /// </ul>
    ///
    /// If the matrix has a const ("static") graph, and if that graph
    /// is fill complete, then the offsets array remains valid through
    /// calls to fillComplete() and resumeFill().
    void
    getLocalDiagCopy (Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& diag,
                      const Kokkos::View<const size_t*, device_type,
                        Kokkos::MemoryUnmanaged>& offsets) const;

    /// \brief Variant of getLocalDiagCopy() that uses precomputed offsets.
    ///
    /// \warning This overload of the method is DEPRECATED.  Call the
    ///   overload above that returns the offsets as a Kokkos::View.
    /// \warning This method is only for expert users.
    /// \warning We make no promises about backwards compatibility
    ///   for this method.  It may disappear or change at any time.
    ///
    /// This method uses the offsets of the diagonal entries, as
    /// precomputed by getLocalDiagOffsets(), to speed up copying the
    /// diagonal of the matrix.  The offsets must be recomputed if any
    /// of the following methods are called:
    /// <ul>
    /// <li> insertGlobalValues() </li>
    /// <li> insertLocalValues() </li>
    /// <li> fillComplete() (with a dynamic graph) </li>
    /// <li> resumeFill() (with a dynamic graph) </li>
    /// </ul>
    ///
    /// If the matrix has a const ("static") graph, and if that graph
    /// is fill complete, then the offsets array remains valid through
    /// calls to fillComplete() and resumeFill().
    void
    getLocalDiagCopy (Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& diag,
                      const Teuchos::ArrayView<const size_t>& offsets) const;

    /// \brief Scale the matrix on the left with the given Vector.
    ///
    /// On return, for all entries i,j in the matrix,
    /// \f$A(i,j) = x(i)*A(i,j)\f$.
    void
    leftScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) override;

    /// \brief Scale the matrix on the right with the given Vector.
    ///
    /// On return, for all entries i,j in the matrix,
    /// \f$A(i,j) = x(j)*A(i,j)\f$.
    void
    rightScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) override;

    //@}
    //! @name Local apply
    //@{

    /// \brief Compute the local part of a sparse matrix-(Multi)Vector
    ///   multiply.
    ///
    /// Most Tpetra users want the apply() method (which see), not
    /// this method.
    ///
    /// This method computes <tt>Y := beta*Y + alpha*Op(A)*X</tt>,
    /// where <tt>Op(A)</tt> is either \f$A\f$, \f$A^T\f$ (the
    /// transpose), or \f$A^H\f$ (the conjugate transpose).
    ///
    /// "The local part" means that this method does no communication
    /// between processes, even if this is necessary for correctness
    /// of the matrix-vector multiply.  Use the apply() method if you
    /// want to compute mathematical sparse matrix-vector multiply.
    ///
    /// This method is mainly of use to Tpetra developers, though some
    /// users may find it helpful if they plan to reuse the result of
    /// doing an Import on the input MultiVector for several sparse
    /// matrix-vector multiplies with matrices that have the same
    /// column Map.
    ///
    /// If you want to do global (not just local) sparse matrix-vector
    /// multiplies with MultiVectors that have different Scalar type
    /// than the CrsMatrix, use CrsMatrixMultiplyOp.  If you want to
    /// do local sparse matrix-vector multiplies with MultiVectors
    /// that have different Scalar type than the CrsMatrix (i.e., if
    /// you're wondering where the templated localMultiply method
    /// went), use LocalCrsMatrixOperator.
    ///
    /// The Map of X and \c mode must satisfy the following:
    /// \code
    /// mode == Teuchos::NO_TRANS &&
    ///   X.getMap ()->isSameAs(* (this->getColMap ())) ||
    /// mode != Teuchos::NO_TRANS &&
    ///   X.getMap ()->isSameAs(* (this->getRowMap ()));
    /// \endcode
    ///
    /// The Map of Y and \c mode must satisfy the following:
    /// \code
    /// mode == Teuchos::NO_TRANS &&
    ///   Y.getMap ()->isSameAs(* (this->getRowMap ())) ||
    /// mode != Teuchos::NO_TRANS &&
    ///   Y.getMap ()->isSameAs(* (this->getColMap ()));
    /// \endcode
    ///
    /// We say that X is "post-Imported," and that Y is
    /// "pre-Exported."
    ///
    /// Both X and Y must have constant stride, and they may not alias
    /// one another (that is, occupy overlapping space in memory).  We
    /// may not necessarily check for aliasing, and if we do, we will
    /// only do this in a debug build.  Aliasing X and Y may cause
    /// nondeterministically incorrect results.
    ///
    /// If <tt>beta == 0</tt>, this operation will enjoy overwrite
    /// semantics: Y's entries will be ignored, and Y will be
    /// overwritten with the result of the multiplication, even if it
    /// contains <tt>Inf</tt> or <tt>NaN</tt> floating-point entries.
    /// Likewise, if <tt>alpha == 0</tt>, this operation will ignore A
    /// and X, even if they contain <tt>Inf</tt> or <tt>NaN</tt>
    /// floating-point entries.
    void
    localApply (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
                MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&Y,
                const Teuchos::ETransp mode = Teuchos::NO_TRANS,
                const Scalar& alpha = Teuchos::ScalarTraits<Scalar>::one (),
                const Scalar& beta = Teuchos::ScalarTraits<Scalar>::zero ()) const;

    /// \brief Return another CrsMatrix with the same entries, but
    ///   converted to a different Scalar type \c T.
    template <class T>
    Teuchos::RCP<CrsMatrix<T, LocalOrdinal, GlobalOrdinal, Node> >
    convert () const;

    //@}
    //! @name Methods implementing Operator
    //@{

    /// \brief Compute a sparse matrix-MultiVector multiply.
    ///
    /// This method computes <tt>Y := beta*Y + alpha*Op(A)*X</tt>,
    /// where <tt>Op(A)</tt> is either \f$A\f$, \f$A^T\f$ (the
    /// transpose), or \f$A^H\f$ (the conjugate transpose).
    ///
    /// If <tt>beta == 0</tt>, this operation will enjoy overwrite
    /// semantics: Y's entries will be ignored, and Y will be
    /// overwritten with the result of the multiplication, even if it
    /// contains <tt>NaN</tt> (not-a-number) floating-point entries.
    void
    apply (const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
           MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&Y,
           Teuchos::ETransp mode = Teuchos::NO_TRANS,
           Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
           Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const override;

    /// \brief Whether apply() allows applying the transpose or
    ///   conjugate transpose.
    bool hasTransposeApply () const override;

    /// \brief The domain Map of this matrix.
    ///
    /// This method implements Tpetra::Operator.  If fillComplete()
    /// has not yet been called at least once on this matrix, or if
    /// the matrix was not constructed with a domain Map, then this
    /// method returns Teuchos::null.
    Teuchos::RCP<const map_type> getDomainMap () const override;

    /// \brief The range Map of this matrix.
    ///
    /// This method implements Tpetra::Operator.  If fillComplete()
    /// has not yet been called at least once on this matrix, or if
    /// the matrix was not constructed with a domain Map, then this
    /// method returns Teuchos::null.
    Teuchos::RCP<const map_type> getRangeMap () const override;

    //@}
    //! @name Other "apply"-like methods
    //@{

    /// \brief Implementation of RowMatrix::add: return <tt>alpha*A + beta*this</tt>.
    ///
    /// This override of the default implementation ensures that, when
    /// called on a CrsMatrix, this method always returns a CrsMatrix
    /// of exactly the same type as <tt>*this</tt>.  "Exactly the same
    /// type" means that all the template parameters match, including
    /// the fifth template parameter.  The input matrix A need not
    /// necessarily be a CrsMatrix or a CrsMatrix of the same type as
    /// <tt>*this</tt>, though this method may be able to optimize
    /// further in that case.
    virtual Teuchos::RCP<RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    add (const Scalar& alpha,
         const RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
         const Scalar& beta,
         const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& domainMap,
         const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& rangeMap,
         const Teuchos::RCP<Teuchos::ParameterList>& params) const override;

    //@}
    //! @name Implementation of Teuchos::Describable interface
    //@{

    //! A one-line description of this object.
    std::string description () const override;

    /// \brief Print this object with the given verbosity level to the
    ///   given output stream.
    void
    describe (Teuchos::FancyOStream& out,
              const Teuchos::EVerbosityLevel verbLevel =
              Teuchos::Describable::verbLevel_default) const override;

    //@}
    //! @name Implementation of DistObject interface
    //@{

    /// \typedef buffer_device_type
    /// \brief Kokkos::Device specialization for communication buffers.
    ///
    /// See #1088 for why this is not just <tt>device_type::device_type</tt>.
    typedef typename DistObject<Scalar, LocalOrdinal, GlobalOrdinal,
                                Node>::buffer_device_type buffer_device_type;

    virtual bool
    checkSizes (const SrcDistObject& source) override;

    void
    applyCrsPadding(
      const typename crs_graph_type::padding_type& padding,
      const bool verbose);

  private:
    void
    copyAndPermuteStaticGraph(
      const RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& source,
      const size_t numSameIDs,
      const LocalOrdinal permuteToLIDs[],
      const LocalOrdinal permuteFromLIDs[],
      const size_t numPermutes);

    void
    copyAndPermuteNonStaticGraph(
      const RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& source,
      const size_t numSameIDs,
      const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& permuteToLIDs_dv,
      const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& permuteFromLIDs_dv,
      const size_t numPermutes);

  protected:

  // clang-format on
  using dist_object_type::
      copyAndPermute; /// DistObject copyAndPermute has multiple overloads --
                      /// use copyAndPermutes for anything we don't override
  // clang-format off

    virtual void
    copyAndPermute
    (const SrcDistObject& source,
     const size_t numSameIDs,
     const Kokkos::DualView<
       const local_ordinal_type*,
       buffer_device_type>& permuteToLIDs,
     const Kokkos::DualView<
       const local_ordinal_type*,
       buffer_device_type>& permuteFromLIDs,
     const CombineMode CM) override;

    virtual void
    packAndPrepare
    (const SrcDistObject& source,
     const Kokkos::DualView<
       const local_ordinal_type*,
       buffer_device_type>& exportLIDs,
     Kokkos::DualView<char*, buffer_device_type>& exports,
     Kokkos::DualView<size_t*, buffer_device_type> numPacketsPerLID,
     size_t& constantNumPackets) override;

  // clang-format on
  using dist_object_type::packAndPrepare; ///< DistObject overloads
                                          ///< packAndPrepare. Explicitly use
                                          ///< DistObject's packAndPrepare for
                                          ///< anything we don't override
                                          // clang-format off

  private:
    /// \brief Unpack the imported column indices and values, and
    ///   combine into matrix.
    void
    unpackAndCombineImpl(
      const Kokkos::DualView<const local_ordinal_type*,
        buffer_device_type>& importLIDs,
      Kokkos::DualView<char*, buffer_device_type> imports,
      Kokkos::DualView<size_t*, buffer_device_type> numPacketsPerLID,
      const size_t constantNumPackets,
      const CombineMode combineMode,
      const bool verbose);

    /// \brief Implementation of unpackAndCombineImpl for when the
    ///   target matrix's structure may change.
    void
    unpackAndCombineImplNonStatic(
      const Kokkos::DualView<const local_ordinal_type*,
        buffer_device_type>& importLIDs,
      Kokkos::DualView<char*, buffer_device_type> imports,
      Kokkos::DualView<size_t*, buffer_device_type> numPacketsPerLID,
      const size_t constantNumPackets,
      const CombineMode combineMode);

  public:
    /// \brief Unpack the imported column indices and values, and
    ///   combine into matrix.
    ///
    /// \warning The allowed CombineMode depends on whether the
    ///   matrix's graph is static or dynamic.  ADD, REPLACE, and
    ///   ABSMAX are valid for a static graph, but INSERT is not.  ADD
    ///   and INSERT are valid for a dynamic graph; ABSMAX and REPLACE
    ///   have not yet been implemented (and would require serious
    ///   changes to matrix assembly in order to implement sensibly).
    void
    unpackAndCombine
    (const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& importLIDs,
     Kokkos::DualView<char*, buffer_device_type> imports,
     Kokkos::DualView<size_t*, buffer_device_type> numPacketsPerLID,
     const size_t constantNumPackets,
     const CombineMode CM) override;

  // clang-format on
  using dist_object_type::unpackAndCombine; ///< DistObject has overloaded
                                            ///< unpackAndCombine, use the
                                            ///< DistObject's implementation for
                                            ///< anything we don't override.
                                            // clang-format off

    /// \brief Pack this object's data for an Import or Export.
    ///
    /// \warning To be called only by the packAndPrepare method of
    ///   appropriate classes of DistObject.
    ///
    /// \param exportLIDs [in] Local indices of the rows to pack.
    /// \param exports [out] On output: array of packed matrix
    ///   entries; allocated by method.
    /// \param numPacketsPerLID [out] On output: numPacketsPerLID[i]
    ///   is the number of bytes of the \c exports array used for
    ///   storing packed local row \c exportLIDs[i].
    /// \param constantNumPackets [out] If zero on output, the packed
    ///   rows may have different numbers of entries.  If nonzero on
    ///   output, then that number gives the constant number of
    ///   entries for all packed rows <i>on all processes in the
    ///   matrix's communicator</i>.
    ///
    /// \subsection Tpetra_CrsMatrix_packNew_summary Packing scheme
    ///
    /// The number of "packets" per row is the number of bytes per
    /// row.  Each row has the following storage format:
    ///
    /// <tt>[numEnt, vals, inds]</tt>,
    ///
    /// where:
    /// <ul>
    /// <li> \c numEnt (\c LocalOrdinal): number of entries in the
    ///      row. </li>
    /// <li> \c vals: array of \c Scalar.  For the k-th entry in the
    ///      row, \c vals[k] is its value and \c inds[k] its global
    ///      column index. </li>
    /// <li> \c inds: array of \c GlobalOrdinal.  For the k-th entry
    ///      in the row, \c vals[k] is its value and \c inds[k] its
    ///      global column index. </li>
    /// </ul>
    ///
    /// We reserve the right to pad for alignment in the future.  In
    /// that case, the number of bytes reported by \c numPacketsPerLID
    /// will reflect padding to align each datum to its size, and the
    /// row will have final padding as well to ensure that the
    /// <i>next</i> row is aligned.  Rows with zero entries will still
    /// take zero bytes, however.
    ///
    /// RowMatrix::pack will always use the same packing scheme as
    /// this method.  This ensures correct Import / Export from a
    /// RowMatrix to a CrsMatrix.
    ///
    /// We do <i>not</i> recommend relying on the details of this
    /// packing scheme.  We describe it here more for Tpetra
    /// developers and less for users.
    ///
    /// \subsection Tpetra_CrsMatrix_packNew_disc Discussion
    ///
    /// DistObject requires packing an object's entries as type
    /// <tt>Packet</tt>, which is the first template parameter of
    /// DistObject.  Since sparse matrices have both values and
    /// indices, we use <tt>Packet=char</tt> and pack them into
    /// buffers of <tt>char</tt> (really "byte").  Indices are stored
    /// as global indices, in case the source and target matrices have
    /// different column Maps (or don't have a column Map yet).
    ///
    /// Currently, we only pack values and column indices.  Row
    /// indices are stored implicitly as the local indices (LIDs) to
    /// pack (see \c exportLIDs).  This is because a DistObject
    /// instance only has one Map, and currently we use the row Map
    /// for CrsMatrix (and RowMatrix).  This makes redistribution of
    /// matrices with 2-D distributions less efficient, but it works
    /// for now.  This may change in the future.
    ///
    /// On output, \c numPacketsPerLID[i] gives the number of bytes
    /// used to pack local row \c exportLIDs[i] of \c this object (the
    /// source object of an Import or Export).  If \c offset is the
    /// exclusive prefix sum-scan of \c numPacketsPerLID, then on
    /// output, <tt>exports[offset[i] .. offset[i+1]]</tt>
    /// (half-exclusive range) contains the packed entries for local
    /// row \c exportLIDs[i].
    ///
    /// Entries for each row use a "struct of arrays" pattern to match
    /// how sparse matrices actually store their data.  The number of
    /// entries in the row goes first, all values go next, and all
    /// column indices (stored as global indices) go last.  Values and
    /// column indices occur in the same order.  Rows with zero
    /// entries always take zero bytes (we do not store their number
    /// of entries explicitly).  This ensures sparsity of storage and
    /// communication in case most rows are empty.
    ///
    /// \subsection Tpetra_CrsMatrix_packNew_why Justification
    ///
    /// GCC >= 4.9 and recent-future versions of the Intel compiler
    /// implement stricter aliasing rules that forbid unaligned type
    /// punning.  If we were to pack as an "array of structs" -- in
    /// this case, an array of <tt>(Scalar, GlobalOrdinal)</tt> pairs
    /// -- then we would either have to pad each matrix entry for
    /// alignment, or call memcpy twice per matrix entry to pack and
    /// unpack.  The "struct of arrays" storage scheme reduces the
    /// padding requirement to a constant per row, or reduces the
    /// number of memcpy calls to two per row.
    ///
    /// We include the number of entries in each row in that row's
    /// packed data, to make unpacking easier.  This saves us from an
    /// error-prone computation to find the number of entries from the
    /// number of bytes.  That computation gets even more difficult if
    /// we have to introduce padding for alignment in the future.
    /// Knowing the number of entries for each row also makes
    /// parallelizing packing and unpacking easier.
    void
    packNew (const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& exportLIDs,
             Kokkos::DualView<char*, buffer_device_type>& exports,
             const Kokkos::DualView<size_t*, buffer_device_type>& numPacketsPerLID,
             size_t& constantNumPackets) const;

  private:
    /// \brief Pack this matrix (part of implementation of packAndPrepare).
    ///
    /// This method helps implement packAndPrepare.
    ///
    /// Call this only when this matrix (which is the source matrix to
    /// pack) does not yet have a KokkosSparse::CrsMatrix.
    void
    packNonStaticNew (const Kokkos::DualView<const local_ordinal_type*, buffer_device_type>& exportLIDs,
                      Kokkos::DualView<char*, buffer_device_type>& exports,
                      const Kokkos::DualView<size_t*, buffer_device_type>& numPacketsPerLID,
                      size_t& constantNumPackets) const;

    /// \brief Pack data for the current row to send.
    ///
    /// \param exports [out] The entire array of packed data to
    ///   "export," that is, to send out from this process.
    /// \param offset [in] Offset into \c exports (which see), at
    ///   which to begin writing this row's packed data.
    /// \param numEnt [in] Number of entries in the row.
    /// \param gidsIn [in] Array of global column indices in the row,
    ///   to pack into the \c exports output buffer.
    /// \param valsIn [in] Array of values in the row, to pack into
    ///   the \c exports output buffer.
    /// \param numBytesPerValue [in] Number of bytes to use for the
    ///   packed representation of a single \c impl_scalar_type matrix
    ///   value.
    ///
    /// This method, like the rest of Tpetra, assumes that all values
    /// of the same type have the same number of bytes in their packed
    /// representation.
    ///
    /// This method does not allocate temporary storage.  We intend
    /// for this to be safe to call in a thread-parallel way at some
    /// point, though it is currently not, due to thread safety issues
    /// with Teuchos::RCP (always) and Teuchos::ArrayView (in a debug
    /// build).
    ///
    /// \return The number of bytes used in the row's packed
    ///   representation.  If numEnt is zero, then the row always uses
    ///   zero bytes (we don't even pack the number of entries in that
    ///   case).
    size_t
    packRow (char exports[],
             const size_t offset,
             const size_t numEnt,
             const GlobalOrdinal gidsIn[],
             const impl_scalar_type valsIn[],
             const size_t numBytesPerValue) const;

    /// \brief Pack data for the current row to send, if the matrix's
    ///   graph is known to be static (and therefore fill complete,
    ///   and locally indexed).
    ///
    /// \param numEntOut [out] Where to write the number of entries in
    ///   the row.
    /// \param valOut [out] Output (packed) array of matrix values.
    /// \param indOut [out] Output (packed) array of matrix column
    ///   indices (as global indices).
    /// \param numEnt [in] Number of entries in the row.
    /// \param lclRow [in] Local index of the row.
    ///
    /// This method does not allocate temporary storage.  We intend
    /// for this to be safe to call in a thread-parallel way on host
    /// (not in CUDA).
    ///
    /// \return \c true if the method succeeded, else \c false.
    ///
    /// \warning (mfh 24 Mar 2017) The current implementation of this
    ///   kernel assumes CUDA UVM.  If we want to fix that, we need to
    ///   write a pack kernel for the whole matrix, that runs on
    ///   device.  As a work-around, consider a fence before and after
    ///   packing.
    bool
    packRowStatic (char* const numEntOut,
                   char* const valOut,
                   char* const indOut,
                   const size_t numEnt,
                   const LocalOrdinal lclRow) const;

    /// \brief Unpack and combine received data for the current row.
    ///
    /// \param gidsOut [out] On output: The row's global column indices.
    /// \param valsOut [out] On output: The row's values.  valsOut[k]
    ///   is the value corresponding to global column index gidsOut[k]
    ///   in this row.
    /// \param imports [in] The entire array of "imported" packed
    ///   data; that is, all the data received from other processes.
    /// \param offset [in] Offset into \c imports (which see), at
    ///   which to begin reading this row's packed data.
    /// \param numBytes [in] Number of bytes of data available to
    ///   unpack for this row.
    /// \param numEnt [in] Number of entries in the row.
    /// \param numBytesPerValue [in] Number of bytes to use for the
    ///   packed representation of a single \c impl_scalar_type matrix
    ///   value.
    ///
    /// This method, like the rest of Tpetra, assumes that all values
    /// of the same type have the same number of bytes in their packed
    /// representation.
    ///
    /// \return The number of bytes used in the row's packed
    ///   representation.  If numEnt is zero, then the row always uses
    ///   zero bytes (we don't even pack the number of entries in that
    ///   case).
    size_t
    unpackRow (GlobalOrdinal gidsOut[],
               impl_scalar_type valsOut[],
               const char imports[],
               const size_t offset,
               const size_t numBytes,
               const size_t numEnt,
               const size_t numBytesPerValue);

    /// \brief Allocate space for packNew() to pack entries to send.
    ///
    /// This is part of the implementation of packAndPrepare, which
    /// helps implement the "new" DistObject interface.
    ///
    /// \param exports [in/out] Pack buffer to (re)allocate.
    /// \param totalNumEntries [out] Total number of entries to send.
    /// \param exportLIDs [in] Local indices of the rows to send.
    void
    allocatePackSpaceNew (Kokkos::DualView<char*, buffer_device_type>& exports,
                          size_t& totalNumEntries,
                          const Kokkos::DualView<const local_ordinal_type*,
                            buffer_device_type>& exportLIDs) const;
    //@}

  public:
    //! Get the Kokkos local values on host, read only
    typename local_matrix_host_type::values_type::const_type
    getLocalValuesHost (Access::ReadOnlyStruct s) const
    {
      return valuesPacked_wdv.getHostView(s);
    }

    //! Get the Kokkos local values on host, read write
    typename local_matrix_host_type::values_type
    getLocalValuesHost (Access::ReadWriteStruct s)
    {
      return valuesPacked_wdv.getHostView(s);
    }

    //! Get the Kokkos local values on host, overwrite all
    typename local_matrix_host_type::values_type
    getLocalValuesHost (Access::OverwriteAllStruct s)
    {
      return valuesPacked_wdv.getHostView(s);
    }

    //! Get the Kokkos local values on device, read only
    typename local_matrix_device_type::values_type::const_type
    getLocalValuesDevice (Access::ReadOnlyStruct s) const
    {
      return valuesPacked_wdv.getDeviceView(s);
    }

    //! Get the Kokkos local values on device, read write
    typename local_matrix_device_type::values_type
    getLocalValuesDevice (Access::ReadWriteStruct s)
    {
      return valuesPacked_wdv.getDeviceView(s);
    }

    //! Get the Kokkos local values on device, overwrite all
    typename local_matrix_device_type::values_type
    getLocalValuesDevice (Access::OverwriteAllStruct s)
    {
      return valuesPacked_wdv.getDeviceView(s);
    }

  private:
    // Friend declaration for nonmember function.
    template<class CrsMatrixType>
    friend Teuchos::RCP<CrsMatrixType>
    Tpetra::importAndFillCompleteCrsMatrix (const Teuchos::RCP<const CrsMatrixType>& sourceMatrix,
                                    const Import<typename CrsMatrixType::local_ordinal_type,
                                                 typename CrsMatrixType::global_ordinal_type,
                                                 typename CrsMatrixType::node_type>& importer,
                                    const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                                 typename CrsMatrixType::global_ordinal_type,
                                                                 typename CrsMatrixType::node_type> >& domainMap,
                                    const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                                 typename CrsMatrixType::global_ordinal_type,
                                                                 typename CrsMatrixType::node_type> >& rangeMap,
                                    const Teuchos::RCP<Teuchos::ParameterList>& params);

    // Friend declaration for nonmember function.
    template<class CrsMatrixType>
    friend Teuchos::RCP<CrsMatrixType>
    Tpetra::importAndFillCompleteCrsMatrix (const Teuchos::RCP<const CrsMatrixType>& sourceMatrix,
                                    const Import<typename CrsMatrixType::local_ordinal_type,
                                                 typename CrsMatrixType::global_ordinal_type,
                                                 typename CrsMatrixType::node_type>& rowImporter,
                                   const Import<typename CrsMatrixType::local_ordinal_type,
                                                typename CrsMatrixType::global_ordinal_type,
                                                typename CrsMatrixType::node_type>& domainImporter,
                                    const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                                 typename CrsMatrixType::global_ordinal_type,
                                                                 typename CrsMatrixType::node_type> >& domainMap,
                                    const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                                 typename CrsMatrixType::global_ordinal_type,
                                                                 typename CrsMatrixType::node_type> >& rangeMap,
                                    const Teuchos::RCP<Teuchos::ParameterList>& params);


    // Friend declaration for nonmember function.
    template<class CrsMatrixType>
    friend Teuchos::RCP<CrsMatrixType>
    Tpetra::exportAndFillCompleteCrsMatrix (const Teuchos::RCP<const CrsMatrixType>& sourceMatrix,
                                    const Export<typename CrsMatrixType::local_ordinal_type,
                                                 typename CrsMatrixType::global_ordinal_type,
                                                 typename CrsMatrixType::node_type>& exporter,
                                    const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                                 typename CrsMatrixType::global_ordinal_type,
                                                                 typename CrsMatrixType::node_type> >& domainMap,
                                    const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                                 typename CrsMatrixType::global_ordinal_type,
                                                                 typename CrsMatrixType::node_type> >& rangeMap,
                                    const Teuchos::RCP<Teuchos::ParameterList>& params);

    // Friend declaration for nonmember function.
    template<class CrsMatrixType>
    friend Teuchos::RCP<CrsMatrixType>
    Tpetra::exportAndFillCompleteCrsMatrix (const Teuchos::RCP<const CrsMatrixType>& sourceMatrix,
                                    const Export<typename CrsMatrixType::local_ordinal_type,
                                                 typename CrsMatrixType::global_ordinal_type,
                                                 typename CrsMatrixType::node_type>& rowExporter,
                                    const Export<typename CrsMatrixType::local_ordinal_type,
                                                 typename CrsMatrixType::global_ordinal_type,
                                                 typename CrsMatrixType::node_type>& domainExporter,
                                    const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                                 typename CrsMatrixType::global_ordinal_type,
                                                                 typename CrsMatrixType::node_type> >& domainMap,
                                    const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                                 typename CrsMatrixType::global_ordinal_type,
                                                                 typename CrsMatrixType::node_type> >& rangeMap,
                                    const Teuchos::RCP<Teuchos::ParameterList>& params);

  public:
    /// \brief Import from <tt>this</tt> to the given destination
    ///   matrix, and make the result fill complete.
    ///
    /// If destMatrix.is_null(), this creates a new matrix as the
    /// destination.  (This is why destMatrix is passed in by nonconst
    /// reference to RCP.)  Otherwise it checks for "pristine" status
    /// and throws if that is not the case.  "Pristine" means that the
    /// matrix has no entries and is not fill complete.
    ///
    /// Use of the "non-member constructor" version of this method,
    /// exportAndFillCompleteCrsMatrix, is preferred for user
    /// applications.
    ///
    /// \warning This method is intended for expert developer use
    ///   only, and should never be called by user code.
    void
    importAndFillComplete (Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& destMatrix,
                           const import_type& importer,
                           const Teuchos::RCP<const map_type>& domainMap,
                           const Teuchos::RCP<const map_type>& rangeMap,
                           const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) const;

    /// \brief Import from <tt>this</tt> to the given destination
    ///   matrix, and make the result fill complete.
    ///
    /// If destMatrix.is_null(), this creates a new matrix as the
    /// destination.  (This is why destMatrix is passed in by nonconst
    /// reference to RCP.)  Otherwise it checks for "pristine" status
    /// and throws if that is not the case.  "Pristine" means that the
    /// matrix has no entries and is not fill complete.
    ///
    /// Use of the "non-member constructor" version of this method,
    /// exportAndFillCompleteCrsMatrix, is preferred for user
    /// applications.
    ///
    /// \warning This method is intended for expert developer use
    ///   only, and should never be called by user code.
    void
    importAndFillComplete (Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& destMatrix,
                           const import_type& rowImporter,
                           const import_type& domainImporter,
                           const Teuchos::RCP<const map_type>& domainMap,
                           const Teuchos::RCP<const map_type>& rangeMap,
                           const Teuchos::RCP<Teuchos::ParameterList>& params) const;


    /// \brief Export from <tt>this</tt> to the given destination
    ///   matrix, and make the result fill complete.
    ///
    /// If destMatrix.is_null(), this creates a new matrix as the
    /// destination.  (This is why destMatrix is passed in by nonconst
    /// reference to RCP.)  Otherwise it checks for "pristine" status
    /// and throws if that is not the case.  "Pristine" means that the
    /// matrix has no entries and is not fill complete.
    ///
    /// Use of the "non-member constructor" version of this method,
    /// exportAndFillCompleteCrsMatrix, is preferred for user
    /// applications.
    ///
    /// \warning This method is intended for expert developer use
    ///   only, and should never be called by user code.
    void
    exportAndFillComplete (Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& destMatrix,
                           const export_type& exporter,
                           const Teuchos::RCP<const map_type>& domainMap = Teuchos::null,
                           const Teuchos::RCP<const map_type>& rangeMap = Teuchos::null,
                           const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) const;

    /// \brief Export from <tt>this</tt> to the given destination
    ///   matrix, and make the result fill complete.
    ///
    /// If destMatrix.is_null(), this creates a new matrix as the
    /// destination.  (This is why destMatrix is passed in by nonconst
    /// reference to RCP.)  Otherwise it checks for "pristine" status
    /// and throws if that is not the case.  "Pristine" means that the
    /// matrix has no entries and is not fill complete.
    ///
    /// Use of the "non-member constructor" version of this method,
    /// exportAndFillCompleteCrsMatrix, is preferred for user
    /// applications.
    ///
    /// \warning This method is intended for expert developer use
    ///   only, and should never be called by user code.
    void
    exportAndFillComplete (Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& destMatrix,
                           const export_type& rowExporter,
                           const export_type& domainExporter,
                           const Teuchos::RCP<const map_type>& domainMap,
                           const Teuchos::RCP<const map_type>& rangeMap,
                           const Teuchos::RCP<Teuchos::ParameterList>& params) const;


  private:
    /// \brief Transfer (e.g. Import/Export) from <tt>this</tt> to the
    ///   given destination matrix, and make the result fill complete.
    ///
    /// If destMat.is_null(), this creates a new matrix, otherwise it
    /// checks for "pristine" status and throws if that is not the
    /// case.  This method implements importAndFillComplete and
    /// exportAndFillComplete, which in turn implemment the nonmember
    /// "constructors" importAndFillCompleteCrsMatrix and
    /// exportAndFillCompleteCrsMatrix.  It's convenient to put those
    /// nonmember constructors' implementations inside the CrsMatrix
    /// class, so that we don't have to put much code in the _decl
    /// header file.
    ///
    /// The point of this method is to fuse three tasks:
    ///
    ///   1. Create a destination matrix (CrsMatrix constructor)
    ///   2. Import or Export this matrix to the destination matrix
    ///   3. Call fillComplete on the destination matrix
    ///
    /// Fusing these tasks can avoid some communication and work.
    void
    transferAndFillComplete (Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& destMatrix,
                             const ::Tpetra::Details::Transfer<LocalOrdinal, GlobalOrdinal, Node>& rowTransfer,
                             const Teuchos::RCP<const ::Tpetra::Details::Transfer<LocalOrdinal, GlobalOrdinal, Node> > & domainTransfer,
                             const Teuchos::RCP<const map_type>& domainMap = Teuchos::null,
                             const Teuchos::RCP<const map_type>& rangeMap = Teuchos::null,
                             const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) const;

    /// \brief Common implementation detail of insertGlobalValues and
    ///   insertGlobalValuesFiltered.
    ///
    /// \pre <tt>&graph == this->getCrsGraph ().getRawPtr ()</tt>
    /// \pre <tt>rowInfo == graph.getRowInfo (rowInfo.localRow)</tt>
    /// \pre <tt>graph.getRowMap ()->isNodeLocalElement (rowInfo.localRow)</tt></li>
    /// \pre <tt>! this->isStaticGraph ()</tt>
    /// \pre If graph has a column Map, then all entries of gblColInds
    ///      are in the column Map on the calling process.  That is, the
    ///      entries of gblColInds (and their corresponding vals entries)
    ///      are "prefiltered," if we needed to filter them.
  protected:
    virtual void
    insertGlobalValuesImpl (crs_graph_type& graph,
                            RowInfo& rowInfo,
                            const GlobalOrdinal gblColInds[],
                            const impl_scalar_type vals[],
                            const size_t numInputEnt);

  private:
    /// \brief Like insertGlobalValues(), but with column filtering.
    ///
    /// "Column filtering" means that if the matrix has a column Map,
    /// then this method ignores entries in columns that are not in
    /// the column Map.
    ///
    /// See discussion in the documentation of getGlobalRowCopy()
    /// about why we use \c Scalar and not \c impl_scalar_type here
    /// for the input array type.
    void
    insertGlobalValuesFiltered(
      const GlobalOrdinal globalRow,
      const Teuchos::ArrayView<const GlobalOrdinal>& indices,
      const Teuchos::ArrayView<const Scalar>& values,
      const bool debug);

    /// \brief Wrapper for insertGlobalValuesFiltered that prints
    ///   helpful error messages if insertGlobalValuesFiltered throws.
    void
    insertGlobalValuesFilteredChecked(
      const GlobalOrdinal globalRow,
      const Teuchos::ArrayView<const GlobalOrdinal>& indices,
      const Teuchos::ArrayView<const Scalar>& values,
      const char* const prefix,
      const bool debug,
      const bool verbose);

    /// \brief Combine in the data using the given combine mode.
    ///
    /// The copyAndPermute() and unpackAndCombine() methods use this
    /// function to combine incoming entries from the source matrix
    /// with the target matrix's current data.  This method's behavior
    /// depends on whether the target matrix (that is, this matrix)
    /// has a static graph.
    ///
    /// See discussion in the documentation of getGlobalRowCopy()
    /// about why we use \c Scalar and not \c impl_scalar_type here
    /// for the input array type.
    void
    combineGlobalValues(
      const GlobalOrdinal globalRowIndex,
      const Teuchos::ArrayView<const GlobalOrdinal>& columnIndices,
      const Teuchos::ArrayView<const Scalar>& values,
      const Tpetra::CombineMode combineMode,
      const char* const prefix,
      const bool debug,
      const bool verbose);

    /// \brief Combine in the data using the given combine mode.
    ///
    /// The copyAndPermute() and unpackAndCombine() methods may this
    /// function to combine incoming entries from the source matrix
    /// with the target matrix's current data.  This method's behavior
    /// depends on whether the target matrix (that is, this matrix)
    /// has a static graph.
    ///
    /// \param lclRow [in] <i>Local</i> row index of the row to modify.
    /// \param numEnt [in] Number of entries in the input data.
    /// \param vals [in] Input values to combine.
    /// \param cols [in] Input (global) column indices corresponding
    ///   to the above values.
    /// \param combineMode [in] The CombineMode to use.
    /// \param prefix [in] Prefix for verbose debugging output; must
    ///   be nonnull if verbose is true.
    /// \param debug [in] Whether to do debug checking.
    /// \param verbose [in] Whether to print verbose debugging output.
    ///
    /// \return The number of modified entries.  No error if and only
    ///   if equal to numEnt.
    LocalOrdinal
    combineGlobalValuesRaw(const LocalOrdinal lclRow,
                           const LocalOrdinal numEnt,
                           const impl_scalar_type vals[],
                           const GlobalOrdinal cols[],
                           const Tpetra::CombineMode combineMode,
                           const char* const prefix,
                           const bool debug,
                           const bool verbose);

    /// \brief Transform CrsMatrix entries, using global indices;
    ///   backwards compatibility version that takes
    ///   Teuchos::ArrayView instead of Kokkos::View.
    ///
    /// See above overload of transformGlobalValues for full documentation.
    ///
    /// \tparam BinaryFunction The type of binary function to apply.
    ///
    /// \param globalRow [in] (Global) index of the row to modify.
    /// \param indices [in] (Global) indices in the row to modify.
    /// \param values [in] Values to use for modification.
    template<class BinaryFunction>
    LocalOrdinal
    transformGlobalValues (const GlobalOrdinal globalRow,
                           const Teuchos::ArrayView<const GlobalOrdinal>& indices,
                           const Teuchos::ArrayView<const Scalar>& values,
                           BinaryFunction f,
                           const bool atomic = useAtomicUpdatesByDefault)
    {
      typedef impl_scalar_type IST;
      typedef LocalOrdinal LO;
      typedef GlobalOrdinal GO;

      const LO numInputEnt = static_cast<LO> (indices.size ());
      if (static_cast<LO> (values.size ()) != numInputEnt) {
        return Teuchos::OrdinalTraits<LO>::invalid ();
      }

      const GO* const inputCols = indices.getRawPtr ();
      const IST* const inputVals =
        reinterpret_cast<const IST*> (values.getRawPtr ());
      return this->transformGlobalValues (globalRow, numInputEnt, inputVals,
                                          inputCols, f, atomic);
    }

    /// \brief Special case of insertGlobalValues for when globalRow
    ///   is <i>not<i> owned by the calling process.
    ///
    /// See discussion in the documentation of getGlobalRowCopy()
    /// about why we use \c Scalar and not \c impl_scalar_type here
    /// for the input array type.
    void
    insertNonownedGlobalValues (const GlobalOrdinal globalRow,
                                const Teuchos::ArrayView<const GlobalOrdinal>& indices,
                                const Teuchos::ArrayView<const Scalar>& values);

    /// \brief Insert indices and their values into the given row.
    ///
    /// \tparam Scalar The type of a single value.  When this method
    ///   is called by CrsMatrix, \c Scalar corresponds to the first
    ///   template parameter of CrsMatrix.
    ///
    /// \pre <tt>! (lg == LocalIndices && I == GlobalIndices)</tt>.
    ///   It does not make sense to give this method local column
    ///   indices (meaning that the graph has a column Map), yet to
    ///   ask it to store global indices.
    ///
    /// \param graph [in/out] The graph into which to insert indices.
    ///
    /// \param rowInfo [in/out] On input: Result of the graph's
    ///   getRowInfo() or updateAllocAndValues() methods, for the
    ///   locally owned row (whose local index is
    ///   <tt>rowInfo.localRow</tt>) for which you want to insert
    ///   indices.  On output: RowInfo with updated newEntries field.
    ///
    /// \param newInds [in] View of the column indices to insert.  If
    ///   <tt>lg == GlobalIndices</tt>, then newInds.ginds, a
    ///   <tt>Teuchos::ArrayView<const GlobalOrdinal></tt>, contains
    ///   the (global) column indices to insert.  Otherwise, if <tt>lg
    ///   == LocalIndices</tt>, then newInds.linds, a
    ///   <tt>Teuchos::ArrayView<const LocalOrdinal></tt>, contains
    ///   the (local) column indices to insert.
    ///
    /// \param oldRowVals [out] View of the current values.  They will
    ///   be overwritten with the new values.
    ///
    /// \param newRowVals [in] View of the new values.  They will be
    ///   copied over the old values.
    ///
    /// \param lg If <tt>lg == GlobalIndices</tt>, then the input
    ///   indices (in \c newInds) are global indices.  Otherwise, if
    ///   <tt>lg == LocalIndices</tt>, the input indices are local
    ///   indices.
    ///
    /// \param I If <tt>lg == GlobalIndices</tt>, then this method
    ///   will store the input indices as global indices.  Otherwise,
    ///   if <tt>I == LocalIndices</tt>, this method will store the
    ///   input indices as local indices.
    void
    insertIndicesAndValues (crs_graph_type& graph,
                            RowInfo& rowInfo,
                            const typename crs_graph_type::SLocalGlobalViews& newInds,
                            const Teuchos::ArrayView<impl_scalar_type>& oldRowVals,
                            const Teuchos::ArrayView<const impl_scalar_type>& newRowVals,
                            const ELocalGlobal lg,
                            const ELocalGlobal I);

  protected:
    // useful typedefs
    typedef Teuchos::OrdinalTraits<LocalOrdinal> OTL;
    typedef Kokkos::ArithTraits<impl_scalar_type> STS;
    typedef Kokkos::ArithTraits<mag_type> STM;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>      V;
    typedef crs_graph_type Graph;

    // Enums
    enum GraphAllocationStatus {
      GraphAlreadyAllocated,
      GraphNotYetAllocated
    };

  protected:
    /// \brief Allocate values (and optionally indices) using the Node.
    ///
    /// \param gas [in] If GraphNotYetAllocated, allocate the
    ///   indices of \c myGraph_ via \c allocateIndices(lg) before
    ///   allocating values.
    ///
    /// \param lg [in] Argument passed into \c
    ///   myGraph_->allocateIndices(), if applicable.
    ///
    /// \param verbose [in] Whether to print verbose debugging output.
    ///
    /// \pre If the graph (that is, staticGraph_) indices are
    ///   already allocated, then gas must be GraphAlreadyAllocated.
    ///   Otherwise, gas must be GraphNotYetAllocated.  We only
    ///   check for this precondition in debug mode.
    ///
    /// \pre If the graph indices are not already allocated, then
    ///   the graph must be owned by the matrix.
    void allocateValues (ELocalGlobal lg, GraphAllocationStatus gas,
                         const bool verbose);

    /// \brief Merge duplicate row indices in the given row, along
    ///   with their corresponding values.
    ///
    /// This method is only called by sortAndMergeIndicesAndValues(),
    /// and only when the matrix owns the graph, not when the matrix
    /// was constructed with a const graph.
    ///
    /// \pre The graph is not already storage optimized:
    ///   <tt>isStorageOptimized() == false</tt>
    /// \return The new row length, after merging.
    static size_t
    mergeRowIndicesAndValues (size_t rowLen, local_ordinal_type* cols, impl_scalar_type* vals);

    /// \brief Sort and merge duplicate local column indices in all
    ///   rows on the calling process, along with their corresponding
    ///   values.
    ///
    /// \pre The matrix is locally indexed (more precisely, not
    ///   globally indexed).
    /// \pre The matrix owns its graph.
    /// \pre The matrix's graph is not already storage optimized:
    ///   <tt>isStorageOptimized() == false</tt>.
    ///
    /// \param sorted [in] If true, the column indices in each row on
    ///   the calling process are already sorted.
    /// \param merged [in] If true, the column indices in each row on
    ///   the calling process are already merged.
    void
    sortAndMergeIndicesAndValues (const bool sorted,
                                  const bool merged);

  public:

    //! Returns true if globalConstants have been computed; false otherwise
    bool haveGlobalConstants() const;

  protected:
    /// \brief Column Map MultiVector used in apply().
    ///
    /// This is a column Map MultiVector.  It is used as the target of
    /// the forward mode Import operation (if necessary) in apply(),
    /// and the source of the reverse mode Export
    /// operation (if necessary) in these methods.  Both of these
    /// methods create this MultiVector on demand if needed, and reuse
    /// it (if possible) for subsequent calls.
    ///
    /// This is declared <tt>mutable</tt> because the methods in
    /// question are const, yet want to cache the MultiVector for
    /// later use.
    mutable Teuchos::RCP<MV> importMV_;

    /// \brief Row Map MultiVector used in apply().
    ///
    /// This is a row Map MultiVector.  It is uses as the source of
    /// the forward mode Export operation (if necessary) in apply(),
    /// and the target of the reverse mode Import
    /// operation (if necessary) in these methods.  Both of these
    /// methods create this MultiVector on demand if needed, and reuse
    /// it (if possible) for subsequent calls.
    ///
    /// This is declared <tt>mutable</tt> because the methods in
    /// question are const, yet want to cache the MultiVector for
    /// later use.
    mutable Teuchos::RCP<MV> exportMV_;

    /// \brief Create a (or fetch a cached) column Map MultiVector.
    ///
    /// \param X_domainMap [in] A domain Map Multivector.  The
    ///   returned MultiVector, if nonnull, will have the same number
    ///   of columns as Y_domainMap.
    ///
    /// \param force [in] Force creating the MultiVector if it hasn't
    ///   been created already.
    ///
    /// The \c force parameter is helpful when the domain Map and the
    /// column Map are the same (so that normally we wouldn't need the
    /// column Map MultiVector), but the following (for example)
    /// holds:
    ///
    /// 1. The kernel needs a constant stride input MultiVector, but
    ///    the given input MultiVector is not constant stride.
    ///
    /// We don't test for the above in this method, because it depends
    /// on the specific kernel.
    Teuchos::RCP<MV>
    getColumnMapMultiVector (const MV& X_domainMap,
                             const bool force = false) const;

    /// \brief Create a (or fetch a cached) row Map MultiVector.
    ///
    /// \param Y_rangeMap [in] A range Map Multivector.  The returned
    ///   MultiVector, if nonnull, will have the same number of
    ///   columns as Y_rangeMap.
    ///
    /// \param force [in] Force creating the MultiVector if it hasn't
    ///   been created already.
    ///
    /// The \c force parameter is helpful when the range Map and the
    /// row Map are the same (so that normally we wouldn't need the
    /// row Map MultiVector), but one of the following holds:
    ///
    /// 1. The kernel needs a constant stride output MultiVector,
    ///    but the given output MultiVector is not constant stride.
    ///
    /// 2. The kernel does not permit aliasing of its input and output
    ///    MultiVector arguments, but they do alias each other.
    ///
    /// We don't test for the above in this method, because it depends
    /// on the specific kernel.
    Teuchos::RCP<MV>
    getRowMapMultiVector (const MV& Y_rangeMap,
                          const bool force = false) const;

    //! Special case of apply() for <tt>mode == Teuchos::NO_TRANS</tt>.
    void
    applyNonTranspose (const MV& X_in,
                       MV& Y_in,
                       Scalar alpha,
                       Scalar beta) const;

    //! Special case of apply() for <tt>mode != Teuchos::NO_TRANS</tt>.
    void
    applyTranspose (const MV& X_in,
                    MV& Y_in,
                    const Teuchos::ETransp mode,
                    Scalar alpha,
                    Scalar beta) const;

    // matrix data accessors

    /// \brief Get a const Host view of the locally owned values
    ///  row myRow, such that rowinfo = getRowInfo(myRow).
    typename values_dualv_type::t_host::const_type
    getValuesViewHost (const RowInfo& rowinfo) const;

    /// \brief Get a const Device view of the locally owned values
    ///  row myRow, such that rowinfo = getRowInfo(myRow).
    typename values_dualv_type::t_dev::const_type
    getValuesViewDevice (const RowInfo& rowinfo) const;

    /// \brief Get a non-const Host view of the locally owned values
    ///  row myRow, such that rowinfo = getRowInfo(myRow).
    typename values_dualv_type::t_host
    getValuesViewHostNonConst (const RowInfo& rowinfo);

    /// \brief Get a non-const Device view of the locally owned values
    ///  row myRow, such that rowinfo = getRowInfo(myRow).
    typename values_dualv_type::t_dev
    getValuesViewDeviceNonConst (const RowInfo& rowinfo);

#if KOKKOSKERNELS_VERSION >= 40299
private:
    // TODO: When KokkosKernels 4.4 is released, local_matrix_device_type can be permanently modified to use the default_size_type
    // of KK. This is always a type that is enabled by KK's ETI (preferring int if both or neither int and size_t are enabled).
    //
    // At that point the ApplyHelper can be replaced with just a SPMVHandle.
    using local_matrix_int_rowptrs_device_type =
      KokkosSparse::CrsMatrix<impl_scalar_type,
                              local_ordinal_type,
                              device_type,
                              void,
                              int>;

    /// The specialization of Details::MatrixApplyHelper used by this class in apply().
    using ApplyHelper = Details::MatrixApplyHelper<
      local_matrix_device_type,
      local_matrix_int_rowptrs_device_type, 
      typename MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::device_view_type>;


    std::shared_ptr<ApplyHelper> getApplyHelper() const {
      if (!applyHelper) {
        auto A_lcl = getLocalMatrixDevice();
        applyHelper = std::make_shared<ApplyHelper>(A_lcl.nnz(), A_lcl.graph.row_map);
      }
      return applyHelper;
    }
#endif

  protected:

    // Friend the tester for CrsMatrix::swap
    friend class Tpetra::crsMatrix_Swap_Tester<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

    // Friend the matrix multiply kernels so they can access internally-cached integer
    // row pointers without making them part of the CrsMatrix interface
    template<typename S, typename LO, typename GO, typename NODE, typename LOV> friend struct Tpetra::MMdetails::KernelWrappers;
    template<typename S, typename LO, typename GO, typename NODE, typename LOV> friend struct Tpetra::MMdetails::KernelWrappers2;


    // friend Matrix Matrix utility function that needs to access integer-typed rowptrs
    friend 
    void Tpetra::MMdetails::import_and_extract_views<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
    const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>&   A,
    Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >   targetMap,
    CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>&   Aview,
    Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> > prototypeImporter,
    bool                                                          userAssertsThereAreNoRemotes,
    const std::string&                                            label,
    const Teuchos::RCP<Teuchos::ParameterList>&                   params);

    /// \brief Swaps the data from *this with the data and maps from crsMatrix
    ///
    /// \param matrix [in/out] a crsMatrix
    void swap(CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & matrix);


  protected:

    /// \brief Fill data into the local matrix.
    ///
    /// This method is only called in fillComplete(), and it is only
    /// called if the graph's structure is already fixed (that is, if
    /// the matrix does not own the graph).
    void fillLocalMatrix (const Teuchos::RCP<Teuchos::ParameterList>& params);

    /// \brief Fill data into the local graph and matrix.
    ///
    /// This method is only called in fillComplete(), and it is only
    /// called if the graph's structure is <i>not</i> already fixed
    /// (that is, if the matrix <i>does</i> own the graph).
    void fillLocalGraphAndMatrix (const Teuchos::RCP<Teuchos::ParameterList>& params);

    //! Check that this object's state is sane; throw if it's not.
    void checkInternalState () const;

    /// \name (Global) graph pointers
    ///
    /// We keep two graph pointers in order to maintain const
    /// correctness.  myGraph_ is a graph which we create internally.
    /// Operations that change the sparsity structure also modify
    /// myGraph_.  If myGraph_ != null, then staticGraph_ == myGraph_
    /// pointerwise (we set the pointers equal to each other when we
    /// create myGraph_).  myGraph_ is only null if this CrsMatrix was
    /// created using the constructor with a const CrsGraph input
    /// argument.  In this case, staticGraph_ is set to the input
    /// CrsGraph.
    //@{
    Teuchos::RCP<const Graph> staticGraph_;
    Teuchos::RCP<      Graph>     myGraph_;
    //@}

protected:
    /// \brief Status of the matrix's storage, when not in a
    ///   fill-complete state.
    ///
    /// The phrase "When not in a fill-complete state" is important.
    /// When the matrix is fill complete, it <i>always</i> uses 1-D
    /// "packed" storage.  However, if the "Optimize Storage"
    /// parameter to fillComplete was false, the matrix may keep
    /// unpacked 1-D storage around and resume it on the next
    /// resumeFill call.
    Details::EStorageStatus storageStatus_ =
      Details::STORAGE_1D_UNPACKED;

    //! Whether the matrix is fill complete.
    bool fillComplete_ = false;

    /// \brief Nonlocal data added using insertGlobalValues().
    ///
    /// These data are cleared by globalAssemble(), once it finishes
    /// redistributing them to their owning processes.
    ///
    /// For a given nonowned global row gRow which was given to
    /// insertGlobalValues() or sumIntoGlobalValues(),
    /// <tt>nonlocals_[gRow].first[k]</tt> is the column index of an
    /// inserted entry, and <tt>nonlocals_[gRow].second[k]</tt> is its
    /// value.  Duplicate column indices for the same row index are
    /// allowed and will be summed during globalAssemble().
    ///
    /// This used to be a map from GlobalOrdinal to (GlobalOrdinal,
    /// Scalar) pairs.  This makes gcc issue a "note" about the ABI of
    /// structs containing std::complex members changing.  CDash
    /// reports this as a warning, even though it's a "note," not a
    /// warning.  However, I don't want it to show up, so I rearranged
    /// the map's value type to a pair of arrays, rather than an array
    /// of pairs.
    ///
    /// \note For Epetra developers: Tpetra::CrsMatrix corresponds
    ///   more to Epetra_FECrsMatrix than to Epetra_CrsMatrix.  The
    ///   insertGlobalValues() method in Tpetra::CrsMatrix, unlike
    ///   its corresponding method in Epetra_CrsMatrix, allows
    ///   insertion into rows which are not owned by the calling
    ///   process.  The globalAssemble() method redistributes these
    ///   to their owning processes.
    std::map<GlobalOrdinal, std::pair<Teuchos::Array<GlobalOrdinal>,
                                      Teuchos::Array<Scalar> > > nonlocals_;

  private:
#if KOKKOSKERNELS_VERSION >= 40299
    /// The apply helper is lazily created in apply(), and reset when resumeFill is called.
    /// It performs 3 functions:
    /// - Decides whether a version of the local matrix with int-typed rowptrs can and should be used to enable spmv TPLs
    /// - Keeps SPMVHandles for both the regular local matrix, and the int-typed version
    /// - Stores the int-typed rowptrs (if they can all be represented by int)
    mutable std::shared_ptr<ApplyHelper> applyHelper;
#endif

  public:
    // FIXME (mfh 24 Feb 2014) Is it _really_ necessary to make this a
    // public inner class of CrsMatrix?  It looks like it doesn't
    // depend on any implementation details of CrsMatrix at all.  It
    // should really be declared and defined outside of CrsMatrix.
    template<class DestViewType, class SrcViewType,
             class DestOffsetViewType, class SrcOffsetViewType>
    struct pack_functor {
      typedef typename DestViewType::execution_space execution_space;
      SrcViewType src_;
      DestViewType dst_;
      SrcOffsetViewType src_offset_;
      DestOffsetViewType dst_offset_;
      typedef typename DestOffsetViewType::non_const_value_type scalar_index_type;

      pack_functor (DestViewType dst, 
                    const SrcViewType src,
                    DestOffsetViewType dst_offset,
                    const SrcOffsetViewType src_offset) :
        src_ (src),
        dst_ (dst),
        src_offset_ (src_offset),
        dst_offset_ (dst_offset)
      {}

      KOKKOS_INLINE_FUNCTION
      void operator () (const LocalOrdinal row) const {
        scalar_index_type srcPos = src_offset_(row);
        const scalar_index_type dstEnd = dst_offset_(row+1);
        scalar_index_type dstPos = dst_offset_(row);
        for ( ; dstPos < dstEnd; ++dstPos, ++srcPos) {
          dst_(dstPos) = src_(srcPos);
        }
      }
    };
  }; // class CrsMatrix

  /// \brief Create an empty CrsMatrix given a row map and a single
  ///   integer upper bound on the number of stored entries per row.
  ///
  /// \relatesalso CrsMatrix
  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class Node>
  Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  createCrsMatrix(
    const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& map,
    const size_t maxNumEntriesPerRow = 0,
    const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
  {
    using matrix_type =
      CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    return Teuchos::rcp(new matrix_type(map, maxNumEntriesPerRow,
                                        params));
  }

  template<class CrsMatrixType>
  Teuchos::RCP<CrsMatrixType>
  importAndFillCompleteCrsMatrix (const Teuchos::RCP<const CrsMatrixType>& sourceMatrix,
                                  const Import<typename CrsMatrixType::local_ordinal_type,
                                               typename CrsMatrixType::global_ordinal_type,
                                               typename CrsMatrixType::node_type>& importer,
                                  const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                               typename CrsMatrixType::global_ordinal_type,
                                                               typename CrsMatrixType::node_type> >& domainMap,
                                  const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                               typename CrsMatrixType::global_ordinal_type,
                                                               typename CrsMatrixType::node_type> >& rangeMap,
                                  const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    Teuchos::RCP<CrsMatrixType> destMatrix;
    sourceMatrix->importAndFillComplete (destMatrix, importer, domainMap, rangeMap, params);
    return destMatrix;
  }

  template<class CrsMatrixType>
  Teuchos::RCP<CrsMatrixType>
  importAndFillCompleteCrsMatrix (const Teuchos::RCP<const CrsMatrixType>& sourceMatrix,
                                  const Import<typename CrsMatrixType::local_ordinal_type,
                                               typename CrsMatrixType::global_ordinal_type,
                                               typename CrsMatrixType::node_type>& rowImporter,
                                  const Import<typename CrsMatrixType::local_ordinal_type,
                                              typename CrsMatrixType::global_ordinal_type,
                                              typename CrsMatrixType::node_type>& domainImporter,
                                  const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                               typename CrsMatrixType::global_ordinal_type,
                                                               typename CrsMatrixType::node_type> >& domainMap,
                                  const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                               typename CrsMatrixType::global_ordinal_type,
                                                               typename CrsMatrixType::node_type> >& rangeMap,
                                  const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    Teuchos::RCP<CrsMatrixType> destMatrix;
    sourceMatrix->importAndFillComplete (destMatrix, rowImporter, domainImporter, domainMap, rangeMap, params);
    return destMatrix;
  }

  template<class CrsMatrixType>
  Teuchos::RCP<CrsMatrixType>
  exportAndFillCompleteCrsMatrix (const Teuchos::RCP<const CrsMatrixType>& sourceMatrix,
                                  const Export<typename CrsMatrixType::local_ordinal_type,
                                               typename CrsMatrixType::global_ordinal_type,
                                               typename CrsMatrixType::node_type>& exporter,
                                  const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                               typename CrsMatrixType::global_ordinal_type,
                                                               typename CrsMatrixType::node_type> >& domainMap,
                                  const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                               typename CrsMatrixType::global_ordinal_type,
                                                               typename CrsMatrixType::node_type> >& rangeMap,
                                  const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    Teuchos::RCP<CrsMatrixType> destMatrix;
    sourceMatrix->exportAndFillComplete (destMatrix, exporter, domainMap, rangeMap, params);
    return destMatrix;
  }

  template<class CrsMatrixType>
  Teuchos::RCP<CrsMatrixType>
  exportAndFillCompleteCrsMatrix (const Teuchos::RCP<const CrsMatrixType>& sourceMatrix,
                                  const Export<typename CrsMatrixType::local_ordinal_type,
                                               typename CrsMatrixType::global_ordinal_type,
                                               typename CrsMatrixType::node_type>& rowExporter,
                                  const Export<typename CrsMatrixType::local_ordinal_type,
                                               typename CrsMatrixType::global_ordinal_type,
                                               typename CrsMatrixType::node_type>& domainExporter,
                                  const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                               typename CrsMatrixType::global_ordinal_type,
                                                               typename CrsMatrixType::node_type> >& domainMap,
                                  const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                               typename CrsMatrixType::global_ordinal_type,
                                                               typename CrsMatrixType::node_type> >& rangeMap,
                                  const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    Teuchos::RCP<CrsMatrixType> destMatrix;
    sourceMatrix->exportAndFillComplete (destMatrix, rowExporter, domainExporter, domainMap, rangeMap, params);
    return destMatrix;
  }

  /// \brief Remove zero entries from a matrix.
  ///
  /// \param matrix    [in/out] CrsMatrix
  /// \param threshold [in]     magnitude threshold below which an entry is deemed to be zero
  ///
  /// \relatesalso CrsMatrix
  template<class CrsMatrixType>
  void
  removeCrsMatrixZeros(CrsMatrixType& matrix,
                       typename Teuchos::ScalarTraits<typename CrsMatrixType::scalar_type>::magnitudeType const & threshold =
                       Teuchos::ScalarTraits<typename CrsMatrixType::scalar_type>::magnitude( Teuchos::ScalarTraits<typename CrsMatrixType::scalar_type>::zero() ))
  {
    auto localMatrix = matrix.getLocalMatrixDevice();
    size_t nnzBefore = localMatrix.nnz();
    localMatrix = KokkosSparse::removeCrsMatrixZeros(localMatrix,threshold);
    size_t localNNZRemoved = nnzBefore - localMatrix.nnz();
    //Skip the expertStaticFillComplete if no entries were removed on any process.
    //The fill complete can perform MPI collectives, so it can only be skipped on all processes or none.
    size_t globalNNZRemoved = 0;
    Teuchos::reduceAll<int, size_t> (*(matrix.getComm()), Teuchos::REDUCE_SUM, 1, &localNNZRemoved, &globalNNZRemoved);
    if(globalNNZRemoved != size_t(0)) {
      matrix.resumeFill();
      matrix.setAllValues(localMatrix);
      matrix.expertStaticFillComplete(matrix.getDomainMap(),matrix.getRangeMap());
    }
  }

} // namespace Tpetra

/**
  \example CrsMatrix_NonlocalAfterResume.hpp
  \brief An example for inserting non-local entries into a
    Tpetra::CrsMatrix using Tpetra::CrsMatrix::insertGlobalValues(),
    with multiple calls to Tpetra::CrsMatrix::fillComplete().
 */

#endif // TPETRA_CRSMATRIX_DECL_HPP
