// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_CONTAINER_DECL_HPP
#define IFPACK2_CONTAINER_DECL_HPP

/// \file Ifpack2_Container.hpp
/// \brief Ifpack2::Container class declaration

#include "Ifpack2_ConfigDefs.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Teuchos_Describable.hpp"
#include <Tpetra_Map.hpp>
#include <Tpetra_BlockCrsMatrix.hpp>
#include <Teuchos_ParameterList.hpp>
#include <iostream>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Teuchos {
  // Forward declaration to avoid include.
  class ParameterList;
} // namespace Teuchos
#endif // DOXYGEN_SHOULD_SKIP_THIS

namespace Ifpack2 {

/// \class Container
/// \brief Interface for creating and solving a set of local linear problems.
/// \tparam MatrixType A specialization of Tpetra::RowMatrix.
///
/// This class is mainly useful for the implementation of
/// BlockRelaxation, and other preconditioners that need to solve
/// linear systems with diagonal blocks of a sparse matrix.
///
/// Users of BlockRelaxation (and any analogous preconditioners) do
/// not normally need to interact with the Container interface.
/// However, they do need to specify a specific Container subclass to
/// use, for example as the second template parameter
/// (<tt>ContainerType</tt>) of BlockRelaxation.  Implementations of
/// Container specify
/// - the kind of data structure used to store the local matrix, and
/// - how to solve a linear system with the local matrix.
///
/// For example, the SparseContainer subclass uses sparse matrices
/// (Tpetra::CrsMatrix) to store each diagonal block, and
/// can use any given Ifpack2 Preconditioner subclass to solve linear
/// systems.
///
/// A Container can create, populate, and solve local linear
/// systems. A local linear system matrix, B, is a submatrix of the
/// local components of the distributed matrix, A. The idea of Container
/// is to specify the rows of A that are contained in each B, then extract
/// B from A, and compute all it is necessary to solve a linear system
/// in B.  Then, set the initial guess (if necessary) and right-hand
/// sides for B, and solve the linear systems in B.
///
/// If you are writing a class (comparable to BlockRelaxation) that
/// uses Container, you should use it in the following way:
/// <ol>
/// <li> Create a Container object, specifying the global matrix A and
///      the indices of the local rows of A that are contained in each B.
///      The latter indices come from a Partitioner object.</li>
/// <li> Optionally, set linear solve parameters using setParameters().</li>
/// <li> Initialize the container by calling initialize().</li>
/// <li> Prepare the linear system solver using compute().</li>
/// <li> Solve a linear system using apply() for each
///      block index (from 0 to numBlocks_) as needed.</li>
/// </ol>
/// For an example of Steps 1-5 above, see
/// ExtractSubmatrices() and ApplyInverseGS() in
/// Ifpack2_BlockRelaxation_def.hpp.
template<class MatrixType>
class Container : public Teuchos::Describable
{
protected:
  using scalar_type = typename MatrixType::scalar_type;
  using local_ordinal_type = typename MatrixType::local_ordinal_type;
  using global_ordinal_type = typename MatrixType::global_ordinal_type;
  using node_type = typename MatrixType::node_type;
  using SC = scalar_type;
  using LO = local_ordinal_type;
  using GO = global_ordinal_type;
  using NO = node_type;
  using import_type = Tpetra::Import<LO, GO, NO>;
  using mv_type = Tpetra::MultiVector<SC, LO, GO, NO>;
  using vector_type = Tpetra::Vector<SC, LO, GO, NO>;
  using map_type = Tpetra::Map<LO, GO, NO>;
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NO>;
  using block_crs_matrix_type = Tpetra::BlockCrsMatrix<SC, LO, GO, NO>;
  using row_matrix_type = Tpetra::RowMatrix<SC, LO, GO, NO>;

  static_assert(std::is_same<MatrixType, row_matrix_type>::value,
                "Ifpack2::Container: Please use MatrixType = Tpetra::RowMatrix.");

  //! Internal representation of Scalar in Kokkos::View
  using ISC = typename Kokkos::ArithTraits<SC>::val_type;

  //! HostView (the host-space internal representation for Tpetra::Multivector) is the
  //! type of the vector arguments of DoJacobi, DoGaussSeidel, and DoSGS.
  using HostView = typename mv_type::dual_view_type::t_host;
  using ConstHostView = typename HostView::const_type;

public:
  /// \brief Constructor.
  ///
  /// \param matrix [in] The original input matrix.  This Container
  ///   will construct local diagonal blocks from its rows according to
  ///   <tt>partitions</tt>.
  /// \param partitioner [in] The Partitioner object that assigns
  ///   local rows of the input matrix to blocks.
  /// \param pointIndexed [in] If the input matrix is a \c Tpetra::BlockCrsMatrix,
  ///    whether elements of \c partitions[k] identify rows within blocks (true) or
  ///    whole blocks (false).
  Container (const Teuchos::RCP<const row_matrix_type>& matrix,
             const Teuchos::Array<Teuchos::Array<LO> >& partitions,
             bool pointIndexed);

  //! Destructor.
  virtual ~Container();

  /// \brief Local indices of the rows of the input matrix that belong to this block.
  ///
  /// The set of (local) rows assigned to this Container is defined by
  /// passing in a set of indices <tt>blockRows[i] = j</tt> to the
  /// constructor, where i (from 0 to <tt>getNumRows() - 1</tt>)
  /// indicates the Container's row, and j indicates the local row in
  /// the calling process.  Subclasses must always pass along these
  /// indices to the base class.
  ///
  /// The indices are usually used to reorder the local row index (on
  /// the calling process) of the i-th row in the Container.
  ///
  /// For an example of how to use these indices, see the
  /// implementation of BlockRelaxation::ExtractSubmatrices() in
  /// Ifpack2_BlockRelaxation_def.hpp.
  Teuchos::ArrayView<const LO> getBlockRows(int blockIndex) const;

  /// \brief Do all set-up operations that only require matrix structure.
  ///
  /// If the input matrix's structure changes, you must call this
  /// method before you may call compute().  You must then call
  /// compute() before you may call apply() or weightedApply().
  ///
  /// "Structure" refers to the graph of the matrix: the local and
  /// global dimensions, and the populated entries in each row.
  virtual void initialize () = 0;

  //! Initialize arrays with information about block sizes.
  void setBlockSizes(const Teuchos::Array<Teuchos::Array<LO> >& partitions);

  void getMatDiag() const;

  //! Whether the container has been successfully initialized.
  bool isInitialized () const;

  //! Whether the container has been successfully computed.
  bool isComputed () const;

  /// \brief Extract the local diagonal blocks and prepare the solver.
  ///
  /// If any entries' values in the input matrix have changed, you
  /// must call this method before you may call apply() or
  /// weightedApply().
  ///
  /// If DOF decoupling is to be used, it must be enabled with enableDecoupling() 
  /// before calling compute().
  virtual void compute () = 0;

  void DoJacobi(ConstHostView X, HostView Y, SC dampingFactor) const;
  void DoOverlappingJacobi(ConstHostView X, HostView Y, ConstHostView W, SC dampingFactor, bool nonsymScaling) const;
  void DoGaussSeidel(ConstHostView X, HostView Y, HostView Y2, SC dampingFactor) const;
  void DoSGS(ConstHostView X, HostView Y, HostView Y2, SC dampingFactor) const;

  //! Set parameters, if any.
  virtual void setParameters (const Teuchos::ParameterList& List) = 0;

  /// \brief Compute <tt>Y := alpha * M^{-1} X + beta*Y</tt>.
  ///
  /// X is in the domain Map of the original matrix (the argument to
  /// compute()), and Y is in the range Map of the original matrix.
  /// This method only reads resp. modifies the permuted subset of
  /// entries of X resp. Y related to the diagonal block M.  That
  /// permuted subset is defined by the indices passed into the
  /// constructor.
  ///
  /// This method is marked \c const for compatibility with
  /// Tpetra::Operator's method of the same name.  This might require
  /// subclasses to mark some of their instance data as \c mutable.
  virtual void
  apply(ConstHostView X,
        HostView Y,
        int blockIndex,
        Teuchos::ETransp mode = Teuchos::NO_TRANS,
        SC alpha = Teuchos::ScalarTraits<SC>::one(),
        SC beta = Teuchos::ScalarTraits<SC>::zero()) const = 0;

  //! Compute <tt>Y := alpha * diag(D) * M^{-1} (diag(D) * X) + beta*Y</tt>.
  virtual void
  weightedApply(ConstHostView X,
                HostView Y,
                ConstHostView D,
                int blockIndex,
                Teuchos::ETransp mode = Teuchos::NO_TRANS,
                SC alpha = Teuchos::ScalarTraits<SC>::one(),
                SC beta = Teuchos::ScalarTraits<SC>::zero()) const = 0;

  /// \brief Compute <tt>Y := (1 - a) Y + a D^{-1} (X - R*Y)</tt>.
  ///
  // The underlying container implements the splitting <tt>A = D + R</tt>. Only
  // it can have efficient access to D and R, as these are constructed in the
  // symbolic and numeric phases.
  //
  // This is the first performance-portable implementation of a block
  // relaxation, and it is supported currently only by BlockTriDiContainer.
  virtual void applyInverseJacobi (const mv_type& /* X */, mv_type& /* Y */,
                                   SC dampingFactor,
                                   bool /* zeroStartingSolution = false */,
                                   int /* numSweeps = 1 */) const = 0;

  //! Wrapper for apply with MultiVector
  virtual void applyMV (const mv_type& X, mv_type& Y) const;

  //! Wrapper for weightedApply with MultiVector
  virtual void weightedApplyMV (const mv_type& X,
                        mv_type& Y,
                        vector_type& W) const;

  virtual void clearBlocks();

  //! Print basic information about the container to \c os.
  virtual std::ostream& print (std::ostream& os) const = 0;

  //! Returns string describing the container.
  //! See <tt>Details::ContainerFactory</tt>.
  static std::string getName();

protected:

  //! Do one step of Gauss-Seidel on block i (used by DoGaussSeidel and DoSGS)
  virtual void DoGSBlock(ConstHostView X, HostView Y, HostView Y2, HostView Resid,
      SC dampingFactor, LO i) const;

  //! The input matrix to the constructor.
  Teuchos::RCP<const row_matrix_type> inputMatrix_;

  //! The input matrix, dynamic cast to CrsMatrix. May be null.
  Teuchos::RCP<const crs_matrix_type> inputCrsMatrix_;

  //! The input matrix, dynamic cast to BlockCrsMatrix. May be null.
  Teuchos::RCP<const block_crs_matrix_type> inputBlockMatrix_;

  //! The number of blocks (partitions) in the container.
  int numBlocks_;
  //! Local indices of the rows of the input matrix that belong to this block.
  Teuchos::Array<LO> blockRows_;      //size: total # of local rows (in all local blocks)
  //! Number of rows in each block.
  Teuchos::Array<LO> blockSizes_;     //size: # of blocks
  //! Starting index in blockRows_ of local row indices for each block.
  Teuchos::Array<LO> blockOffsets_;
  //! Diagonal elements.
  mutable Teuchos::RCP<vector_type> Diag_;
  //! Whether the problem is distributed across multiple MPI processes.
  bool IsParallel_;
  //! Number of local rows in input matrix.
  LO NumLocalRows_;
  //! Number of global rows in input matrix.
  GO NumGlobalRows_;
  //! Number of nonzeros in input matrix.
  GO NumGlobalNonzeros_;
  //! Whether the input matrix is a BlockCRS matrix.
  bool hasBlockCrs_;
  //! If hasBlockCrs_, the number of DOFs per vertex. Otherwise 1.
  int bcrsBlockSize_;
  //! (If hasBlockCrs_) Whether the blocks are described using sub-block row indices instead of full block rows.
  bool pointIndexed_;
  //! Number of scalars corresponding to one element of blockRows_
  //! If pointIndexed_, always 1. Otherwise if hasBlockCrs_, then bcrsBlockSize_. Otherwise 1.
  LO scalarsPerRow_;

  //! If \c true, the container has been successfully initialized.
  bool IsInitialized_;
  //! If \c true, the container has been successfully computed.
  bool IsComputed_;

  LO maxBlockSize_;
};

namespace Details
{
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  struct StridedRowView;
}

/// \class ContainerImpl
/// \brief The implementation of the numerical features of
/// Container (Jacobi, Gauss-Seidel, SGS). This class 
/// allows a custom scalar type (LocalScalarType) to be 
/// used for storing blocks and solving the block systems.
/// Hiding this template parameter from the Container
/// interface simplifies the 
/// BlockRelaxation and ContainerFactory classes.

template<class MatrixType, class LocalScalarType>
class ContainerImpl : public Container<MatrixType>
{
  //! @name Internal typedefs (protected)
  //@{
protected:
  using local_scalar_type = LocalScalarType;
  using SC = typename Container<MatrixType>::scalar_type;
  using LO = typename Container<MatrixType>::local_ordinal_type;
  using GO = typename Container<MatrixType>::global_ordinal_type;
  using NO = typename Container<MatrixType>::node_type;
  using StridedRowView = Details::StridedRowView<SC, LO, GO, NO>;
  using typename Container<MatrixType>::import_type;
  using typename Container<MatrixType>::row_matrix_type;
  using typename Container<MatrixType>::crs_matrix_type;
  using typename Container<MatrixType>::block_crs_matrix_type;
  using typename Container<MatrixType>::mv_type;
  using typename Container<MatrixType>::vector_type;
  using typename Container<MatrixType>::map_type;
  using typename Container<MatrixType>::ISC;
  //! The internal representation of LocalScalarType in Kokkos::View
  using LSC = LocalScalarType;
  using LISC = typename Kokkos::ArithTraits<LSC>::val_type;

  using local_mv_type = Tpetra::MultiVector<LSC, LO, GO, NO>;

  using typename Container<MatrixType>::HostView;
  using typename Container<MatrixType>::ConstHostView;
  using HostViewLocal = typename local_mv_type::dual_view_type::t_host;
  using HostSubviewLocal = Kokkos::View<LISC**, Kokkos::LayoutStride, typename HostViewLocal::memory_space>;
  using ConstHostSubviewLocal = Kokkos::View<const LISC**, Kokkos::LayoutStride, typename HostViewLocal::memory_space>;

  static_assert(std::is_same<MatrixType, row_matrix_type>::value,
                "Ifpack2::Container: Please use MatrixType = Tpetra::RowMatrix.");
  //@}
  //

  //! Translate local row (if `pointIndexed_`, then a row within a block row/node)
  //! to the corresponding column, so that the intersection of the
  //! row and column is on the global diagonal. Translates to global ID first
  //! using the row map, and then back to a local column ID using the col map.
  //! If the corresponding column is off-process, then returns
  //! OrdinalTraits<LO>::invalid() and does not throw an exception.
  LO translateRowToCol(LO row);

  //! View a row of the input matrix.
  Details::StridedRowView<SC, LO, GO, NO> getInputRowView(LO row) const;

  /// \brief Extract the submatrices identified by the local indices set by the constructor. BlockMatrix may be any type that
  /// supports direct entry access: `Scalar& operator()(size_t row, size_t col)`.
  /// \param diagBlocks [in/out] The diagonal block matrices. Its size must be `this->numBlocks_`.
  ///   Each BlockMatrix must be ready for entries to be assigned.
  virtual void extract() = 0;

public:

  ContainerImpl (const Teuchos::RCP<const row_matrix_type>& matrix,
                 const Teuchos::Array<Teuchos::Array<LO> >& partitions,
                 bool pointIndexed);

  //! Destructor.
  virtual ~ContainerImpl();

  /// \brief Do all set-up operations that only require matrix structure.
  ///
  /// If the input matrix's structure changes, you must call this
  /// method before you may call compute().  You must then call
  /// compute() before you may call apply() or weightedApply().
  ///
  /// "Structure" refers to the graph of the matrix: the local and
  /// global dimensions, and the populated entries in each row.
  virtual void initialize () = 0;

  /// \brief Extract the local diagonal blocks and prepare the solver.
  ///
  /// If any entries' values in the input matrix have changed, you
  /// must call this method before you may call apply() or
  /// weightedApply().
  ///
  /// If DOF decoupling is to be used, it must be enabled with enableDecoupling() 
  /// before calling compute().
  virtual void compute () = 0;

  //! Set parameters, if any.
  virtual void setParameters (const Teuchos::ParameterList& List);

  /// \brief Compute <tt>Y := alpha * M^{-1} X + beta*Y</tt>.
  ///
  /// X is in the domain Map of the original matrix (the argument to
  /// compute()), and Y is in the range Map of the original matrix.
  /// This method only reads resp. modifies the permuted subset of
  /// entries of X resp. Y related to the diagonal block M.  That
  /// permuted subset is defined by the indices passed into the
  /// constructor.
  ///
  /// This method is marked \c const for compatibility with
  /// Tpetra::Operator's method of the same name.  This might require
  /// subclasses to mark some of their instance data as \c mutable.
  virtual void
  apply(ConstHostView X,
        HostView Y,
        int blockIndex,
        Teuchos::ETransp mode = Teuchos::NO_TRANS,
        SC alpha = Teuchos::ScalarTraits<SC>::one(),
        SC beta = Teuchos::ScalarTraits<SC>::zero()) const;

  //! Compute <tt>Y := alpha * diag(D) * M^{-1} (diag(D) * X) + beta*Y</tt>.
  virtual void
  weightedApply(ConstHostView X,
                HostView Y,
                ConstHostView D,
                int blockIndex,
                Teuchos::ETransp mode = Teuchos::NO_TRANS,
                SC alpha = Teuchos::ScalarTraits<SC>::one(),
                SC beta = Teuchos::ScalarTraits<SC>::zero()) const;

  /// \brief Compute <tt>Y := (1 - a) Y + a D^{-1} (X - R*Y)</tt>.
  ///
  // The underlying container implements the splitting <tt>A = D + R</tt>. Only
  // it can have efficient access to D and R, as these are constructed in the
  // symbolic and numeric phases.
  //
  // This is the first performance-portable implementation of a block
  // relaxation, and it is supported currently only by BlockTriDiContainer.
  virtual void applyInverseJacobi (const mv_type& /* X */, mv_type& /* Y */,
                                   SC dampingFactor,
                                   bool /* zeroStartingSolution = false */,
                                   int /* numSweeps = 1 */) const;

  //! Wrapper for apply with MVs, used in unit tests (never called by BlockRelaxation)
  void applyMV (const mv_type& X, mv_type& Y) const;

  //! Wrapper for weightedApply with MVs, used in unit tests (never called by BlockRelaxation)
  void weightedApplyMV (const mv_type& X,
                        mv_type& Y,
                        vector_type& W) const;

  virtual void clearBlocks();

  //! Print basic information about the container to \c os.
  virtual std::ostream& print (std::ostream& os) const = 0;

  //! Returns string describing the container.
  //! See <tt>Details::ContainerFactory</tt>.
  static std::string getName();

protected:
  //Do Gauss-Seidel on only block i (this is used by DoGaussSeidel and DoSGS)
  void DoGSBlock(ConstHostView X, HostView Y, HostView Y2, HostView Resid,
      SC dampingFactor, LO i) const;

  //! Exactly solves the linear system By = x, where B is a diagonal block matrix
  //! (blockIndex), and X, Y are multivector subviews with the same length as B's dimensions.
  //!
  //! The Dense, Banded and TriDi containers all implement this and it is used in ContainerImpl::apply().
  //! The Sparse and BlockTriDi containers have their own implementation of apply() and do not use solveBlock.
  virtual void
  solveBlock(ConstHostSubviewLocal X,
             HostSubviewLocal Y,
             int blockIndex,
             Teuchos::ETransp mode,
             const LSC alpha,
             const LSC beta) const;

  //! Scratch vectors used in apply().
  mutable HostViewLocal X_local_;  //length: blockRows_.size()
  mutable HostViewLocal Y_local_;  //length: blockRows_.size()

  //! \c applyScratch provides space for temporary block-sized vectors
  //! in \c weightedApply(), so that full Kokkos::Views don't
  //! need to be created for every block apply (slow).
  //!
  //! Layout of applyScratch_:
  //! | Name     | Range                                   |
  //! |----------|-----------------------------------------|
  //! | D_local  | 0 to maxBlockSize_                      |
  //! | X_scaled | maxBlockSize_  to 2 * maxBlockSize_     |
  //! | Y_temp   | 2 * maxBlockSize_  to 3 * maxBlockSize_ |
  mutable HostViewLocal weightedApplyScratch_;

  //! Views for holding pieces of X corresponding to each block
  mutable std::vector<HostSubviewLocal> X_localBlocks_;

  //! Views for holding pieces of Y corresponding to each block
  mutable std::vector<HostSubviewLocal> Y_localBlocks_;
};

namespace Details {
  /// \brief Structure for read-only views of general matrix rows
  ///
  /// Supports rows within the nodes of a BlockCrsMatrix (point indexing).
  /// Use of getLocalRowCopy
  /// This is required for extracting diagonal blocks, and decoupling DOFs.
  template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
  struct StridedRowView
  {
    using SC = Scalar;
    using LO = LocalOrdinal;

    using block_crs_matrix_type = Tpetra::BlockCrsMatrix<SC, LO, GlobalOrdinal, Node>;

    using h_inds_type = typename block_crs_matrix_type::local_inds_host_view_type;
    using h_vals_type = typename block_crs_matrix_type::values_host_view_type;
    //! Constructor for row views (preferred)
    StridedRowView(h_vals_type vals_, h_inds_type inds_, int blockSize_, size_t nnz_);

    //! Constructor for row views 
    //    StridedRowView(const SC* vals_, const LO* inds_, int blockSize_, size_t nnz_);

    //! Constructor for deep copy (fallback, if matrix doesn't support row views)
    StridedRowView(Teuchos::Array<SC>& vals_, Teuchos::Array<LO>& inds_);
        
    SC val(size_t i) const;
    LO ind(size_t i) const;

    size_t size() const;

    private:
    h_vals_type vals;
    h_inds_type inds;
    int blockSize;
    size_t nnz;
    //These arrays are only used if the inputMatrix_ doesn't support row views.
    Teuchos::Array<SC> valsCopy;
    Teuchos::Array<LO> indsCopy;
  };
} // namespace Details

} // namespace Ifpack2

//! Print information about the given Container to the output stream \c os.
template <class MatrixType>
std::ostream& operator<<(std::ostream& os, const Ifpack2::Container<MatrixType>& obj);

namespace Teuchos {

/// \brief Partial specialization of TypeNameTraits for Ifpack2::Container.
///
/// \tparam MatrixType The template parameter of Ifpack2::Container.
///   Must be a Tpetra::RowMatrix specialization.
template<class MatrixType>
class TEUCHOSCORE_LIB_DLL_EXPORT TypeNameTraits< ::Ifpack2::Container<MatrixType> >
{
 public:
  static std::string name () {
    return std::string ("Ifpack2::Container<") +
      TypeNameTraits<MatrixType>::name () + ">";
  }

  static std::string concreteName (const ::Ifpack2::Container<MatrixType>&) {
    return name ();
  }
};

} // namespace Teuchos

#endif // IFPACK2_CONTAINER_HPP
