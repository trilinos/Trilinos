/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_CONTAINER_HPP
#define IFPACK2_CONTAINER_HPP

/// \file Ifpack2_Container.hpp
/// \brief Ifpack2::Container class declaration

#include "Ifpack2_ConfigDefs.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Teuchos_Describable.hpp"
#include <Ifpack2_Partitioner.hpp>
#include <Ifpack2_Details_MultiVectorLocalGatherScatter.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_BlockCrsMatrix.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Time.hpp>
#include <iostream>
#include <type_traits>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Teuchos {
  // Forward declaration to avoid include.
  class ParameterList;
} // namespace Teuchos
#endif // DOXYGEN_SHOULD_SKIP_THIS

namespace Ifpack2 {

/// \class Container
/// \brief Interface for creating and solving a local linear problem.
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
/// For example, the SparseContainer subclass uses a sparse matrix (in
/// particular, Tpetra::CrsMatrix) to store each diagonal block, and
/// can use any given Ifpack2 Preconditioner subclass to solve linear
/// systems.
///
/// A Container can create, populate, and solve a local linear
/// system. The local linear system matrix, B, is a submatrix of the
/// local components of a distributed matrix, A. The idea of Container
/// is to specify the rows of A that are contained in B, then extract
/// B from A, and compute all it is necessary to solve a linear system
/// in B.  Then, set the initial guess (if necessary) and right-hand
/// side for B, and solve the linear system in B.
///
/// If you are writing a class (comparable to BlockRelaxation) that
/// uses Container, you should use it in the following way:
/// <ol>
/// <li> Create a Container object, specifying the global matrix A and
///      the indices of the local rows of A that are contained in B.
///      The latter indices come from a Partitioner object.</li>
/// <li> Optionally, set linear solve parameters using setParameters().</li>
/// <li> Initialize the container by calling initialize().</li>
/// <li> Prepare the linear system solver using compute().</li>
/// <li> Solve the linear system using apply().</li>
/// </ol>
/// For an example of Steps 1-5 above, see the implementation of
/// BlockRelaxation::ExtractSubmatrices() in
/// Ifpack2_BlockRelaxation_def.hpp.
template<class MatrixType>
class Container : public Teuchos::Describable
{
public:
  typedef typename MatrixType::scalar_type scalar_type;
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;
  typedef typename MatrixType::node_type node_type;
  typedef Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type> import_type;
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> mv_type;
  typedef Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> vector_type;
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> crs_matrix_type;
  typedef Tpetra::BlockCrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> block_crs_matrix_type;
  typedef Tpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> row_matrix_type;

  static_assert(std::is_same<MatrixType, row_matrix_type>::value,
                "Ifpack2::Container: Please use MatrixType = Tpetra::RowMatrix.");

  //! Internal representation of Scalar in Kokkos::View
  typedef typename Kokkos::Details::ArithTraits<scalar_type>::val_type impl_scalar_type;

  //! HostView (the host-space internal representation for Tpetra::Multivector) is the
  //! type of the vector arguments of DoJacobi, DoGaussSeidel, and DoSGS.
  typedef typename mv_type::dual_view_type::t_host HostView;

  /// \brief Constructor.
  ///
  /// \brief matrix [in] The original input matrix.  This Container
  ///   will construct local diagonal blocks from the rows given by
  ///   <tt>partitioner</tt>.
  ///
  /// \param partitioner [in] The Partitioner object that assigns
  ///   local rows of the input matrix to blocks.
  Container (const Teuchos::RCP<const row_matrix_type>& matrix,
             const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
             bool pointIndexed) :
    inputMatrix_ (matrix),
    inputCrsMatrix_ (Teuchos::rcp_dynamic_cast<const crs_matrix_type>(inputMatrix_)),
    inputBlockMatrix_ (Teuchos::rcp_dynamic_cast<const block_crs_matrix_type>(inputMatrix_)),
    pointIndexed_(pointIndexed),
    IsInitialized_(false),
    IsComputed_(false)
  {
    using Teuchos::Ptr;
    using Teuchos::RCP;
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::Comm;
    NumLocalRows_ = inputMatrix_->getNodeNumRows();
    NumGlobalRows_ = inputMatrix_->getGlobalNumRows();
    NumGlobalNonzeros_ = inputMatrix_->getGlobalNumEntries();
    IsParallel_ = inputMatrix_->getRowMap()->getComm()->getSize() != 1;
    hasBlockCrs_ = !inputBlockMatrix_.is_null();
    if(hasBlockCrs_)
      bcrsBlockSize_ = inputBlockMatrix_->getBlockSize();
    else
      bcrsBlockSize_ = 1;
    if(hasBlockCrs_ && !pointIndexed_)
      scalarsPerRow_ = bcrsBlockSize_;
    else
      scalarsPerRow_ = 1;
    setBlockSizes(partitions);
    //Sanity check the partitions
    #ifdef HAVE_IFPACK2_DEBUG
    // Check whether the input set of local row indices is correct.
    const map_type& rowMap = *inputMatrix_->getRowMap();
    for(int i = 0; i < numBlocks_; i++)
    {
      Teuchos::ArrayView<const local_ordinal_type> blockRows = getBlockRows(i);
      for(local_ordinal_type j = 0; j < blockSizes_[i]; j++)
      {
        local_ordinal_type row = blockRows[j];
        if(pointIndexed)
        {
          //convert the point row to the corresponding block row
          row /= bcrsBlockSize_;
        }
        TEUCHOS_TEST_FOR_EXCEPTION(
          !rowMap.isNodeLocalElement(row),
          std::invalid_argument, "Ifpack2::Container: "
          "On process " << rowMap.getComm()->getRank() << " of "
          << rowMap.getComm()->getSize() << ", in the given set of local row "
          "indices blockRows = " << Teuchos::toString(blockRows) << ", the following "
          "entries is not valid local row index on the calling process: "
          << row << ".");
      }
    }
    #endif
  }

  /// \brief Constructor for single block (used in unit tests)
  ///
  /// \brief matrix [in] The original input matrix.  This Container
  ///   will construct a local diagonal block from the rows given by
  ///   <tt>blockRows</tt>.
  ///
  /// \param blockRows [in] The set of (local) rows assigned to this
  ///   container.  <tt>blockRows[i] == j</tt>, where i (from 0 to
  ///   <tt>getNumRows() - 1</tt>) indicates the Container's row, and
  ///   j indicates the local row in the calling process.  Subclasses
  ///   must always pass along these indices to the base class.
  Container (const Teuchos::RCP<const row_matrix_type>& matrix,
             const Teuchos::Array<local_ordinal_type>& blockRows,
             bool pointIndexed) :
    inputMatrix_ (matrix),
    inputCrsMatrix_ (Teuchos::rcp_dynamic_cast<const crs_matrix_type>(inputMatrix_)),
    inputBlockMatrix_ (Teuchos::rcp_dynamic_cast<const block_crs_matrix_type>(inputMatrix_)),
    numBlocks_ (1),
    blockRows_ (blockRows),
    blockSizes_ (1, blockRows.size()),
    blockOffsets_ (1, 0),
    pointIndexed_(pointIndexed),
    IsInitialized_(false),
    IsComputed_(false)
  {
    NumLocalRows_ = inputMatrix_->getNodeNumRows();
    NumGlobalRows_ = inputMatrix_->getGlobalNumRows();
    NumGlobalNonzeros_ = inputMatrix_->getGlobalNumEntries();
    IsParallel_ = inputMatrix_->getRowMap()->getComm()->getSize() > 1;
    hasBlockCrs_ = !inputBlockMatrix_.is_null();
    if(hasBlockCrs_)
      bcrsBlockSize_ = inputBlockMatrix_->getBlockSize();
    else
      bcrsBlockSize_ = 1;
    if(hasBlockCrs_ && !pointIndexed_)
      scalarsPerRow_ = bcrsBlockSize_;
    else
      scalarsPerRow_ = 1;
    maxBlockSize_ = blockSizes_[0] * scalarsPerRow_;
    //Sanity check the rows in the block
    #ifdef HAVE_IFPACK2_DEBUG
    // Check whether the input set of local row indices is correct.
    const map_type& rowMap = *inputMatrix_->getRowMap();
    for(size_t i = 0; i < blockRows_.size(); i++)
    {
      local_ordinal_type row = blockRows_[i];
      if(pointIndexed_)
      {
        row /= bcrsBlockSize_;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(
        !rowMap.isNodeLocalElement(row),
        std::invalid_argument, "Ifpack2::Container: "
        "On process " << rowMap.getComm()->getRank() << " of "
        << rowMap.getComm()->getSize() << ", in the given set of local row "
        "indices blockRows = " << Teuchos::toString(blockRows_) << ", the following "
        "entries is not valid local row index on the calling process: "
        << row << ".");
    }
    #endif
  }

  //! Destructor.
  virtual ~Container() {};

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
  Teuchos::ArrayView<const local_ordinal_type> getBlockRows(int blockIndex) const
  {
    return Teuchos::ArrayView<const local_ordinal_type>
      (&blockRows_[blockOffsets_[blockIndex]], blockSizes_[blockIndex]);
  }

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
  void setBlockSizes(const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions)
  {
    //First, create a grand total of all the rows in all the blocks
    //Note: If partitioner allowed overlap, this could be greater than the # of local rows
    local_ordinal_type totalBlockRows = 0;
    numBlocks_ = partitions.size();
    blockSizes_.resize(numBlocks_);
    blockOffsets_.resize(numBlocks_);
    maxBlockSize_ = 0;
    for(int i = 0; i < numBlocks_; i++)
    {
      local_ordinal_type rowsInBlock = partitions[i].size();
      blockSizes_[i] = rowsInBlock;
      blockOffsets_[i] = totalBlockRows;
      totalBlockRows += rowsInBlock;
      maxBlockSize_ = std::max(maxBlockSize_, rowsInBlock * scalarsPerRow_);
    }
    blockRows_.resize(totalBlockRows);
    //set blockRows_: each entry is the partition/block of the row
    local_ordinal_type iter = 0;
    for(int i = 0; i < numBlocks_; i++)
    {
      for(int j = 0; j < blockSizes_[i]; j++)
      {
        blockRows_[iter++] = partitions[i][j];
      }
    }
  }

  void getMatDiag() const
  {
    if(Diag_.is_null())
    {
      Diag_ = rcp(new vector_type(inputMatrix_->getDomainMap()));
      inputMatrix_->getLocalDiagCopy(*Diag_);
    }
  }

  //! Whether the container has been successfully initialized.
  bool isInitialized () const {
    return IsInitialized_;
  }

  //! Whether the container has been successfully computed.
  bool isComputed () const {
    return IsComputed_;
  }

  /// \brief Extract the local diagonal blocks and prepare the solver.
  ///
  /// If any entries' values in the input matrix have changed, you
  /// must call this method before you may call apply() or
  /// weightedApply().
  ///
  /// If DOF decoupling is to be used, it must be enabled with enableDecoupling() 
  /// before calling compute().
  virtual void compute () = 0;

  void DoJacobi(HostView& X, HostView& Y, scalar_type dampingFactor) const;
  void DoOverlappingJacobi(HostView& X, HostView& Y, HostView& W, scalar_type dampingFactor) const;
  void DoGaussSeidel(HostView& X, HostView& Y, HostView& Y2, scalar_type dampingFactor) const;
  void DoSGS(HostView& X, HostView& Y, HostView& Y2, scalar_type dampingFactor) const;

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
  apply(HostView& X,
        HostView& Y,
        int blockIndex,
        Teuchos::ETransp mode = Teuchos::NO_TRANS,
        scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
        scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const = 0;

  //! Compute <tt>Y := alpha * diag(D) * M^{-1} (diag(D) * X) + beta*Y</tt>.
  virtual void
  weightedApply(HostView& X,
                HostView& Y,
                HostView& D,
                int blockIndex,
                Teuchos::ETransp mode = Teuchos::NO_TRANS,
                scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
                scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const = 0;

  /// \brief Compute <tt>Y := (1 - a) Y + a D^{-1} (X - R*Y)</tt>.
  ///
  // The underlying container implements the splitting <tt>A = D + R</tt>. Only
  // it can have efficient access to D and R, as these are constructed in the
  // symbolic and numeric phases.
  //
  // This is the first performance-portable implementation of a block
  // relaxation, and it is supported currently only by BlockTriDiContainer.
  virtual void applyInverseJacobi (const mv_type& /* X */, mv_type& /* Y */,
                                   scalar_type dampingFactor,
                                   bool /* zeroStartingSolution = false */,
                                   int /* numSweeps = 1 */) const = 0;

  //! Wrapper for apply with MVs, used in unit tests (never called by BlockRelaxation)
  virtual void applyMV (mv_type& X, mv_type& Y) const
  {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Not implemented.");
  }

  //! Wrapper for weightedApply with MVs, used in unit tests (never called by BlockRelaxation)
  virtual void weightedApplyMV (mv_type& X,
                        mv_type& Y,
                        vector_type& W)
  {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Not implemented.");
  }

  virtual void clearBlocks();

  //! Print basic information about the container to \c os.
  virtual std::ostream& print (std::ostream& os) const = 0;

  //! Returns string describing the container.
  //! See <tt>Details::ContainerFactory</tt>.
  static std::string getName()
  {
    return "Generic";
  }

protected:

  //Do Gauss-Seidel on only block i (this is used by DoGaussSeidel and DoSGS)
  virtual void DoGSBlock(HostView& X, HostView& Y, HostView& Y2, HostView& Resid,
      scalar_type dampingFactor, local_ordinal_type i) const
  {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Not implemented.");
  }

  //! The input matrix to the constructor.
  Teuchos::RCP<const row_matrix_type> inputMatrix_;

  //! The input matrix, dynamic cast to CrsMatrix. May be null.
  Teuchos::RCP<const crs_matrix_type> inputCrsMatrix_;

  //! The input matrix, dynamic cast to BlockCrsMatrix. May be null.
  Teuchos::RCP<const block_crs_matrix_type> inputBlockMatrix_;

  //! The number of blocks (partitions) in the container.
  int numBlocks_;
  //! Local indices of the rows of the input matrix that belong to this block.
  Teuchos::Array<local_ordinal_type> blockRows_;      //size: total # of local rows (in all local blocks)
  //! Number of rows in each block.
  Teuchos::Array<local_ordinal_type> blockSizes_;     //size: # of blocks
  //! Starting index in blockRows_ of local row indices for each block.
  Teuchos::Array<local_ordinal_type> blockOffsets_;
  //! Diagonal elements.
  mutable Teuchos::RCP<vector_type> Diag_;
  //! Whether the problem is distributed across multiple MPI processes.
  bool IsParallel_;
  //! Number of local rows in input matrix.
  local_ordinal_type NumLocalRows_;
  //! Number of global rows in input matrix.
  global_ordinal_type NumGlobalRows_;
  //! Number of nonzeros in input matrix.
  global_ordinal_type NumGlobalNonzeros_;
  //! Whether the input matrix is a BlockCRS matrix.
  bool hasBlockCrs_;
  //! If hasBlockCrs_, the number of DOFs per vertex. Otherwise 1.
  int bcrsBlockSize_;
  //! (If hasBlockCrs_) Whether the blocks are described using sub-block row indices instead of full block rows.
  bool pointIndexed_;
  //! Number of scalars corresponding to one element of blockRows_
  //! If pointIndexed_, always 1. Otherwise if hasBlockCrs_, then bcrsBlockSize_. Otherwise 1.
  local_ordinal_type scalarsPerRow_;

  //! If \c true, the container has been successfully initialized.
  bool IsInitialized_;
  //! If \c true, the container has been successfully computed.
  bool IsComputed_;

  local_ordinal_type maxBlockSize_;
};

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
  typedef LocalScalarType local_scalar_type;
  typedef typename Container<MatrixType>::scalar_type scalar_type;
  typedef typename Container<MatrixType>::local_ordinal_type local_ordinal_type;
  typedef typename Container<MatrixType>::global_ordinal_type global_ordinal_type;
  typedef typename Container<MatrixType>::node_type node_type;
  typedef typename Container<MatrixType>::import_type import_type;
  typedef typename Container<MatrixType>::row_matrix_type row_matrix_type;
  typedef typename Container<MatrixType>::crs_matrix_type crs_matrix_type;
  typedef typename Container<MatrixType>::block_crs_matrix_type block_crs_matrix_type;
  typedef typename Container<MatrixType>::mv_type mv_type;
  typedef typename Container<MatrixType>::vector_type vector_type;
  typedef typename Container<MatrixType>::map_type map_type;
  typedef typename Container<MatrixType>::STS STS;
  typedef typename Container<MatrixType>::impl_scalar_type impl_scalar_type;
  typedef Tpetra::MultiVector<local_scalar_type, local_ordinal_type, global_ordinal_type, node_type> local_mv_type;

  //! The internal representation of LocalScalarType in Kokkos::View
  typedef typename Kokkos::Details::ArithTraits<local_scalar_type>::val_type local_impl_scalar_type;

  typedef typename mv_type::dual_view_type::t_host HostView;
  typedef typename local_mv_type::dual_view_type::t_host HostViewLocal;
  typedef typename Kokkos::View<local_impl_scalar_type**, Kokkos::LayoutStride, Kokkos::HostSpace> HostSubview;

  static_assert(std::is_same<MatrixType, row_matrix_type>::value,
                "Ifpack2::Container: Please use MatrixType = Tpetra::RowMatrix.");
  //@}
  //

  /// \brief Internal structure for read-only views of general matrix rows
  ///
  /// Supports rows within the nodes of a BlockCrsMatrix
  /// (required for extracting diagonal blocks, and decoupling DOFs)
  /// In that case, values is viewed as a 1-D array and indices refer to columns within blocks.
  struct StridedRowView
  {
    //! Constructor for row views (preferred)
    StridedRowView(const scalar_type* vals_, const local_ordinal_type* inds_, int blockSize_, size_t nnz_)
      : vals(vals_), inds(inds_), blockSize(blockSize_), nnz(nnz_)
    {}

    //! Constructor for deep copy (fallback, if matrix doesn't support row views)
    StridedRowView(Teuchos::Array<scalar_type>& vals_, Teuchos::Array<local_ordinal_type>& inds_)
      : vals(nullptr), inds(nullptr), blockSize(1), nnz(vals_.size())
    {
      valsCopy.swap(vals_);
      indsCopy.swap(inds_);
    }
    
    scalar_type val(size_t i) const
    {
      #ifdef HAVE_IFPACK2_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION(i >= nnz, std::runtime_error,
            "Out-of-bounds access into Ifpack2::Container::StridedRowView");
      #endif
      if(vals)
      {
        if(blockSize == 1)
          return vals[i];
        else
          return vals[i * blockSize];
      }
      else
        return valsCopy[i];
    }
    local_ordinal_type ind(size_t i) const
    {
      #ifdef HAVE_IFPACK2_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION(i >= nnz, std::runtime_error,
            "Out-of-bounds access into Ifpack2::Container::StridedRowView");
      #endif
      //inds is smaller than vals by a factor of the block size (dofs/node)
      if(inds)
      {
        if(blockSize == 1)
          return inds[i];
        else
          return inds[i / blockSize] * blockSize + i % blockSize;
      }
      else
        return indsCopy[i];
    }
    size_t size()
    {
      return nnz;
    }
    private:
      const scalar_type* vals;
      const local_ordinal_type* inds;
      int blockSize;
      size_t nnz;
      //These arrays are only used if the inputMatrix_ doesn't support row views.
      Teuchos::Array<scalar_type> valsCopy;
      Teuchos::Array<local_ordinal_type> indsCopy;
  };

  //! View a row of the input matrix.
  StridedRowView getInputRowView(local_ordinal_type row) const
  {
    if(this->hasBlockCrs_)
    {
      const local_ordinal_type* colinds;
      scalar_type* values;
      local_ordinal_type numEntries;
      this->inputBlockMatrix_->getLocalRowView(row / this->bcrsBlockSize_, colinds, values, numEntries);
      return StridedRowView(values + row % this->bcrsBlockSize_, colinds, this->bcrsBlockSize_, numEntries * this->bcrsBlockSize_);
    }
    else if(!this->inputMatrix_->supportsRowViews())
    {
      size_t maxEntries = this->inputMatrix_->getNodeMaxNumRowEntries();
      Teuchos::Array<local_ordinal_type> indsCopy(maxEntries);
      Teuchos::Array<scalar_type> valsCopy(maxEntries);
      size_t numEntries;
      this->inputMatrix_->getLocalRowCopy(row, indsCopy, valsCopy, numEntries);
      indsCopy.resize(numEntries);
      valsCopy.resize(numEntries);
      return StridedRowView(valsCopy, indsCopy);
    }
    else
    {
      const local_ordinal_type* colinds;
      const scalar_type* values;
      local_ordinal_type numEntries;
      this->inputMatrix_->getLocalRowViewRaw(row, numEntries, colinds, values);
      return StridedRowView(values, colinds, 1, numEntries);
    }
  }
  
  //! Translate local row (if `pointIndexed_`, then a row within a block row/node)
  //! to the corresponding column, so that the intersection of the
  //! row and column is on the global diagonal. Translates to global ID first
  //! using the row map, and then back to a local column ID using the col map.
  //! If the corresponding column is off-process, then returns
  //! OrdinalTraits<local_ordinal_type>::invalid() and does not throw an exception.
  local_ordinal_type translateRowToCol(local_ordinal_type row)
  {
    auto LO_INVALID = Teuchos::OrdinalTraits<local_ordinal_type>::invalid();
    auto GO_INVALID = Teuchos::OrdinalTraits<global_ordinal_type>::invalid();
    const map_type& globalRowMap = *(this->inputMatrix_->getRowMap());
    const map_type& globalColMap = *(this->inputMatrix_->getColMap());
    local_ordinal_type rowLID = row;
    local_ordinal_type dofOffset = 0;
    if(this->pointIndexed_)
    {
      rowLID = row / this->bcrsBlockSize_;
      dofOffset = row % this->bcrsBlockSize_;
    }
    global_ordinal_type diagGID = globalRowMap.getGlobalElement(rowLID);
    TEUCHOS_TEST_FOR_EXCEPTION(
      diagGID == GO_INVALID,
      std::runtime_error, "Ifpack2::Container::translateRowToCol: "
      "On process " << this->inputMatrix_->getRowMap()->getComm()->getRank() <<
      ", at least one row index in the set of local "
      "row indices given to the constructor is not a valid local row index in "
      "the input matrix's row Map on this process.  This should be impossible "
      "because the constructor checks for this case.  Here is the complete set "
      "of invalid local row indices: " << rowLID << ".  "
      "Please report this bug to the Ifpack2 developers.");
    //now, can translate diagGID (both a global row AND global col ID) to local column
    local_ordinal_type colLID = globalColMap.getLocalElement(diagGID);
    TEUCHOS_TEST_FOR_EXCEPTION(
      colLID == LO_INVALID,
      std::runtime_error, "Ifpack2::Container::translateRowToCol: "
      "On process " << this->inputMatrix_->getRowMap()->getComm()->getRank() << ", "
      "at least one row index in the set of row indices given to the constructor "
      "does not have a corresponding column index in the input matrix's column "
      "Map.  This probably means that the column(s) in question is/are empty on "
      "this process, which would make the submatrix to extract structurally "
      "singular. The invalid global column index is " << diagGID << ".");
    //colLID could identify a block column - translate to split column if needed
    if(this->pointIndexed_)
      return colLID * this->bcrsBlockSize_ + dofOffset;
    return colLID;
  }

  //! Extract the submatrices identified by the local indices set by the constructor. BlockMatrix may be any type that
  //! supports direct entry access: `Scalar& operator()(size_t row, size_t col)`.
  //! \pre `diagBlocks.size() == this->numBlocks_`
  //! \pre If diagBlocks elements do not manage their own
  //!     memory, there must be space allocated for all entries
  template<typename BlockMatrix>
  void extract(std::vector<BlockMatrix>& diagBlocks)
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    const auto INVALID = Teuchos::OrdinalTraits<local_ordinal_type>::invalid();
    //To extract diagonal blocks, need to translate local rows to local columns.
    //Strategy: make a lookup table that translates local cols in the matrix to offsets in blockRows_:
    //blockOffsets_[b] <= offset < blockOffsets_[b+1]: tests whether the column is in block b.
    //offset - blockOffsets_[b]: gives the column within block b.
    //
    //This provides the block and col within a block in O(1).
    if(this->scalarsPerRow_ > 1)
    {
      Array<local_ordinal_type> colToBlockOffset(this->inputBlockMatrix_->getNodeNumCols(), INVALID);
      for(int i = 0; i < this->numBlocks_; i++)
      {
        //Get the interval where block i is defined in blockRows_
        local_ordinal_type blockStart = this->blockOffsets_[i];
        local_ordinal_type blockEnd = blockStart + this->blockSizes_[i];
        ArrayView<const local_ordinal_type> localRows = this->getBlockRows(i);
        //Set the lookup table entries for the columns appearing in block i.
        //If OverlapLevel_ > 0, then this may overwrite values for previous blocks, but
        //this is OK. The values updated here are only needed to process block i's entries.
        for(size_t j = 0; j < (size_t) localRows.size(); j++)
        {
          local_ordinal_type localCol = translateRowToCol(localRows[j]);
          colToBlockOffset[localCol] = blockStart + j;
        }
        for(local_ordinal_type blockRow = 0; blockRow < (local_ordinal_type) localRows.size(); blockRow++)
        {
          //get a raw view of the whole block row
          const local_ordinal_type* indices;
          scalar_type* values;
          local_ordinal_type numEntries;
          local_ordinal_type inputRow = this->blockRows_[blockStart + blockRow];
          this->inputBlockMatrix_->getLocalRowView(inputRow, indices, values, numEntries);
          for(local_ordinal_type k = 0; k < numEntries; k++)
          {
            local_ordinal_type colOffset = colToBlockOffset[indices[k]];
            if(blockStart <= colOffset && colOffset < blockEnd)
            {
              //This entry does appear in the diagonal block.
              //(br, bc) identifies the scalar's position in the BlockCrs block.
              //Convert this to (r, c) which is its position in the container block.
              local_ordinal_type blockCol = colOffset - blockStart;
              for(local_ordinal_type bc = 0; bc < this->bcrsBlockSize_; bc++)
              {
                for(local_ordinal_type br = 0; br < this->bcrsBlockSize_; br++)
                {
                  local_ordinal_type r = this->bcrsBlockSize_ * blockRow + br;
                  local_ordinal_type c = this->bcrsBlockSize_ * blockCol + bc;
                  auto val = values[k * (this->bcrsBlockSize_ * this->bcrsBlockSize_) + (br + this->bcrsBlockSize_ * bc)];
                  if(val != 0)
                    diagBlocks[i](r, c) = val;
                }
              }
            }
          }
        }
      }
    }
    else
    {
      //get the mapping from point-indexed matrix columns to offsets in blockRows_
      //(this includes regular CrsMatrix columns, in which case bcrsBlockSize_ == 1)
      Array<local_ordinal_type> colToBlockOffset(this->inputMatrix_->getNodeNumCols() * this->bcrsBlockSize_, INVALID);
      for(int i = 0; i < this->numBlocks_; i++)
      {
        //Get the interval where block i is defined in blockRows_
        local_ordinal_type blockStart = this->blockOffsets_[i];
        local_ordinal_type blockEnd = blockStart + this->blockSizes_[i];
        ArrayView<const local_ordinal_type> localRows = this->getBlockRows(i);
        //Set the lookup table entries for the columns appearing in block i.
        //If OverlapLevel_ > 0, then this may overwrite values for previous blocks, but
        //this is OK. The values updated here are only needed to process block i's entries.
        for(size_t j = 0; j < (size_t) localRows.size(); j++)
        {
          //translateRowToCol will return the corresponding split column
          local_ordinal_type localCol = translateRowToCol(localRows[j]);
          colToBlockOffset[localCol] = blockStart + j;
        }
        for(size_t blockRow = 0; blockRow < (size_t) localRows.size(); blockRow++)
        {
          //get a view of the split row
          local_ordinal_type inputPointRow = this->blockRows_[blockStart + blockRow];
          auto rowView = getInputRowView(inputPointRow);
          for(size_t k = 0; k < rowView.size(); k++)
          {
            local_ordinal_type colOffset = colToBlockOffset[rowView.ind(k)];
            if(blockStart <= colOffset && colOffset < blockEnd)
            {
              local_ordinal_type blockCol = colOffset - blockStart;
              auto val = rowView.val(k);
              if(val != 0)
                diagBlocks[i](blockRow, blockCol) = rowView.val(k);
            }
          }
        }
      }
    }
  }

public:

  ContainerImpl (const Teuchos::RCP<const row_matrix_type>& matrix,
                 const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
                 bool pointIndexed)
    : Container<MatrixType>(matrix, partitions, pointIndexed)
  {}

  ContainerImpl (const Teuchos::RCP<const row_matrix_type>& matrix,
                 const Teuchos::Array<local_ordinal_type>& blockRows,
                 bool pointIndexed)
    : Container<MatrixType>(matrix, blockRows, pointIndexed)
  {}

  //! Destructor.
  virtual ~ContainerImpl() {};

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
  virtual void setParameters (const Teuchos::ParameterList& List) {}

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
  apply(HostView& X,
        HostView& Y,
        int blockIndex,
        Teuchos::ETransp mode = Teuchos::NO_TRANS,
        scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
        scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //! Compute <tt>Y := alpha * diag(D) * M^{-1} (diag(D) * X) + beta*Y</tt>.
  virtual void
  weightedApply(HostView& X,
                HostView& Y,
                HostView& D,
                int blockIndex,
                Teuchos::ETransp mode = Teuchos::NO_TRANS,
                scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
                scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  /// \brief Compute <tt>Y := (1 - a) Y + a D^{-1} (X - R*Y)</tt>.
  ///
  // The underlying container implements the splitting <tt>A = D + R</tt>. Only
  // it can have efficient access to D and R, as these are constructed in the
  // symbolic and numeric phases.
  //
  // This is the first performance-portable implementation of a block
  // relaxation, and it is supported currently only by BlockTriDiContainer.
  virtual void applyInverseJacobi (const mv_type& /* X */, mv_type& /* Y */,
                                   scalar_type dampingFactor,
                                   bool /* zeroStartingSolution = false */,
                                   int /* numSweeps = 1 */) const
  { TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Not implemented."); }

  //! Wrapper for apply with MVs, used in unit tests (never called by BlockRelaxation)
  void applyMV (mv_type& X, mv_type& Y) const
  {
    HostView XView = X.getLocalViewHost();
    HostView YView = Y.getLocalViewHost();
    this->apply (XView, YView, 0);
  }

  //! Wrapper for weightedApply with MVs, used in unit tests (never called by BlockRelaxation)
  void weightedApplyMV (mv_type& X,
                        mv_type& Y,
                        vector_type& W)
  {
    HostView XView = X.getLocalViewHost();
    HostView YView = Y.getLocalViewHost();
    HostView WView = W.getLocalViewHost();
    weightedApply (XView, YView, WView, 0);
  }

  virtual void clearBlocks();

  //! Print basic information about the container to \c os.
  virtual std::ostream& print (std::ostream& os) const = 0;

  //! Returns string describing the container.
  //! See <tt>Details::ContainerFactory</tt>.
  static std::string getName()
  {
    return "Generic";
  }

protected:
  //Do Gauss-Seidel on only block i (this is used by DoGaussSeidel and DoSGS)
  void DoGSBlock(HostView& X, HostView& Y, HostView& Y2, HostView& Resid,
      scalar_type dampingFactor, local_ordinal_type i) const;

  //! Exactly solves the linear system By = x, where B is a diagonal block matrix
  //! (blockIndex), and X, Y are multivector subviews with the same length as B's dimensions.
  //!
  //! The Dense, Banded and TriDi containers all implement this and it is used in ContainerImpl::apply().
  //! The Sparse and BlockTriDi containers have their own implementation of apply() and do not use solveBlock.
  virtual void
  solveBlock(HostSubview& X,
             HostSubview& Y,
             int blockIndex,
             Teuchos::ETransp mode,
             const local_scalar_type alpha,
             const local_scalar_type beta) const
  {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Not implemented.");
  }

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
  mutable std::vector<HostSubview> X_localBlocks_;

  //! Views for holding pieces of Y corresponding to each block
  mutable std::vector<HostSubview> Y_localBlocks_;
};

//! Print information about the given Container to the output stream \c os.
template <class MatrixType>
inline std::ostream&
operator<< (std::ostream& os, const Ifpack2::Container<MatrixType>& obj)
{
  return obj.print (os);
}

template <class MatrixType>
void Container<MatrixType>::DoJacobi(HostView& X, HostView& Y, scalar_type dampingFactor) const
{
  const scalar_type one = STS::one();
  // Note: Flop counts copied naively from Ifpack.
  // use blockRows_ and blockSizes_
  size_t numVecs = X.extent(1);
  // Non-overlapping Jacobi
  for (local_ordinal_type i = 0; i < numBlocks_; i++)
  {
    // may happen that a partition is empty
    if(blockSizes_[i] != 1 || hasBlockCrs_)
    {
      if(blockSizes_[i] == 0 )
        continue;
      apply(X, Y, i, Teuchos::NO_TRANS, dampingFactor, one);
    }
    else    // singleton, can't access Containers_[i] as it was never filled and may be null.
    {
      local_ordinal_type LRID = blockRows_[blockOffsets_[i]];
      getMatDiag();
      HostView diagView = Diag_->getLocalViewHost();
      impl_scalar_type d = (impl_scalar_type) one / diagView(LRID, 0);
      for(size_t nv = 0; nv < numVecs; nv++)
      {
        impl_scalar_type x = X(LRID, nv);
        Y(LRID, nv) = x * d;
      }
    }
  }
}

template <class MatrixType>
void Container<MatrixType>::DoOverlappingJacobi(HostView& X, HostView& Y, HostView& W, scalar_type dampingFactor) const
{
  // Overlapping Jacobi
  for(local_ordinal_type i = 0; i < numBlocks_; i++)
  {
    // may happen that a partition is empty
    if(blockSizes_[i] == 0)
      continue;
    if(blockSizes_[i] != 1)
      weightedApply(X, Y, W, i, Teuchos::NO_TRANS, dampingFactor, STS::one());
  }
}

//Do Gauss-Seidel with just block i
//This is used 3 times: once in DoGaussSeidel and twice in DoSGS
template<class MatrixType, typename LocalScalarType>
void ContainerImpl<MatrixType, LocalScalarType>::DoGSBlock(
    HostView& X, HostView& Y, HostView& Y2, HostView& Resid,
    scalar_type dampingFactor, local_ordinal_type i) const
{
  using Teuchos::ArrayView;
  size_t numVecs = X.extent(1);
  const scalar_type one = STS::one();
  if(this->blockSizes_[i] == 0)
    return; // Skip empty partitions
  if(this->hasBlockCrs_ && !this->pointIndexed_)
  {
    //Use efficient blocked version
    ArrayView<const local_ordinal_type> blockRows = this->getBlockRows(i);
    const size_t localNumRows = this->blockSizes_[i];
    for(size_t j = 0; j < localNumRows; j++)
    {
      local_ordinal_type row = blockRows[j]; // Containers_[i]->ID (j);
      local_ordinal_type numEntries;
      scalar_type* values;
      const local_ordinal_type* colinds;
      this->inputBlockMatrix_->getLocalRowView(row, colinds, values, numEntries);
      for(size_t m = 0; m < numVecs; m++)
      {
        for (int localR = 0; localR < this->bcrsBlockSize_; localR++)
          Resid(row * this->bcrsBlockSize_ + localR, m) = X(row * this->bcrsBlockSize_ + localR, m);
        for (local_ordinal_type k = 0; k < numEntries; ++k)
        {
          const local_ordinal_type col = colinds[k];
          for(int localR = 0; localR < this->bcrsBlockSize_; localR++)
          {
            for(int localC = 0; localC < this->bcrsBlockSize_; localC++)
            {
              Resid(row * this->bcrsBlockSize_ + localR, m) -=
                values[k * this->bcrsBlockSize_ * this->bcrsBlockSize_ + localR + localC * this->bcrsBlockSize_]
                * Y2(col * this->bcrsBlockSize_ + localC, m); }
          }
        }
      }
    }
    // solve with this block
    //
    // Note: I'm abusing the ordering information, knowing that X/Y
    // and Y2 have the same ordering for on-proc unknowns.
    //
    // Note: Add flop counts for inverse apply
    apply(Resid, Y2, i, Teuchos::NO_TRANS, dampingFactor, one);
  }
  else if(!this->hasBlockCrs_ && this->blockSizes_[i] == 1)
  {
    // singleton, can't access Containers_[i] as it was never filled and may be null.
    // a singleton calculation (just using matrix diagonal) is exact, all residuals should be zero.
    local_ordinal_type LRID = this->blockOffsets_[i];  // by definition, a singleton 1 row in block.
    HostView diagView = this->Diag_->getLocalViewHost();
    impl_scalar_type d = (impl_scalar_type) one / diagView(LRID, 0);
    for(size_t m = 0; m < numVecs; m++)
    {
      impl_scalar_type x = X(LRID, m);
      impl_scalar_type newy = x * d;
      Y2(LRID, m) = newy;
    }
  }
  else if(!this->inputCrsMatrix_.is_null() &&
      std::is_same<typename crs_matrix_type::device_type::memory_space, Kokkos::HostSpace>::value)
  {
    //Use the KokkosSparse internal matrix for low-overhead values/indices access
    //But, can only do this if the matrix is accessible directly from host, since it's not a DualView
    auto localA = this->inputCrsMatrix_->getLocalMatrix();
    typedef typename decltype(localA)::size_type size_type;
    const auto& rowmap = localA.graph.row_map;
    const auto& entries = localA.graph.entries;
    const auto& values = localA.values;
    ArrayView<const local_ordinal_type> blockRows = this->getBlockRows(i);
    for(size_t j = 0; j < (size_t) blockRows.size(); j++)
    {
      const local_ordinal_type row = blockRows[j];
      for(size_t m = 0; m < numVecs; m++)
      {
        scalar_type r = X(row, m);
        for(size_type k = rowmap(row); k < rowmap(row + 1); k++)
        {
          const local_ordinal_type col = entries(k);
          r -= values(k) * Y2(col, m);
        }
        Resid(row, m) = r;
      }
    }
    // solve with this block
    //
    // Note: I'm abusing the ordering information, knowing that X/Y
    // and Y2 have the same ordering for on-proc unknowns.
    //
    // Note: Add flop counts for inverse apply
    apply(Resid, Y2, i, Teuchos::NO_TRANS, dampingFactor, one);
  }
  else
  {
    //Either a point-indexed block matrix, or a normal row matrix
    //that doesn't support getLocalMatrix
    ArrayView<const local_ordinal_type> blockRows = this->getBlockRows(i);
    for(size_t j = 0; j < (size_t) blockRows.size(); j++)
    {
      const local_ordinal_type row = blockRows[j];
      auto rowView = getInputRowView(row);
      for(size_t m = 0; m < numVecs; m++)
      {
        Resid(row, m) = X(row, m);
        for (size_t k = 0; k < rowView.size(); ++k)
        {
          const local_ordinal_type col = rowView.ind(k);
          Resid(row, m) -= rowView.val(k) * Y2(col, m);
        }
      }
    }
    // solve with this block
    //
    // Note: I'm abusing the ordering information, knowing that X/Y
    // and Y2 have the same ordering for on-proc unknowns.
    //
    // Note: Add flop counts for inverse apply
    apply(Resid, Y2, i, Teuchos::NO_TRANS, dampingFactor, one);
  }
}

template<class MatrixType>
void Container<MatrixType>::
DoGaussSeidel(HostView& X, HostView& Y, HostView& Y2, scalar_type dampingFactor) const
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::Ptr;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  //This function just extracts the diagonal if it hasn't already.
  getMatDiag();
  // Note: Flop counts copied naively from Ifpack.
  auto numVecs = X.extent(1);
  // X = RHS, Y = initial guess
  HostView Resid("", X.extent(0), X.extent(1));
  for(local_ordinal_type i = 0; i < numBlocks_; i++)
  {
    DoGSBlock(X, Y, Y2, Resid, dampingFactor, i);
  }
  if(IsParallel_)
  {
    auto numMyRows = inputMatrix_->getNodeNumRows();
    for (size_t m = 0; m < numVecs; ++m)
    {
      for (size_t i = 0; i < numMyRows * bcrsBlockSize_; ++i)
      {
        Y(i, m) = Y2(i, m);
      }
    }
  }
}

template<class MatrixType>
void Container<MatrixType>::
DoSGS(HostView& X, HostView& Y, HostView& Y2, scalar_type dampingFactor) const
{
  // X = RHS, Y = initial guess
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::Ptr;
  using Teuchos::ptr;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  auto numVecs = X.extent(1);
  HostView Resid("", X.extent(0), X.extent(1));
  // Forward Sweep
  for(local_ordinal_type i = 0; i < numBlocks_; i++)
  {
    DoGSBlock(X, Y, Y2, Resid, dampingFactor, i);
  }
  static_assert(std::is_signed<local_ordinal_type>::value,
      "Local ordinal must be signed (unsigned breaks reverse iteration to 0)");
  // Reverse Sweep
  for(local_ordinal_type i = numBlocks_ - 1; i >= 0; --i)
  {
    DoGSBlock(X, Y, Y2, Resid, dampingFactor, i);
  }
  if(IsParallel_)
  {
    auto numMyRows = inputMatrix_->getNodeNumRows();
    for (size_t m = 0; m < numVecs; ++m)
    {
      for (size_t i = 0; i < numMyRows * bcrsBlockSize_; ++i)
      {
        Y(i, m) = Y2(i, m);
      }
    }
  }
}

template<class MatrixType, class LocalScalarType>
void ContainerImpl<MatrixType, LocalScalarType>::
apply (HostView& X,
       HostView& Y,
       int blockIndex,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // The local operator might have a different Scalar type than
  // MatrixType.  This means that we might have to convert X and Y to
  // the Tpetra::MultiVector specialization that the local operator
  // wants.  This class' X_ and Y_ internal fields are of the right
  // type for the local operator, so we can use those as targets.

  Details::MultiVectorLocalGatherScatter<mv_type, local_mv_type> mvgs;

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! this->IsComputed_, std::runtime_error, "Ifpack2::Container::apply: "
    "You must have called the compute() method before you may call apply().  "
    "You may call the apply() method as many times as you want after calling "
    "compute() once, but you must have called compute() at least once.");

  const size_t numVecs = X.extent(1);

  if(numVecs == 0) {
    return; // done! nothing to do
  }

  // The local operator works on a permuted subset of the local parts
  // of X and Y.  The subset and permutation are defined by the index
  // array returned by getBlockRows().  If the permutation is trivial
  // and the subset is exactly equal to the local indices, then we
  // could use the local parts of X and Y exactly, without needing to
  // permute.  Otherwise, we have to use temporary storage to permute
  // X and Y.  For now, we always use temporary storage.
  //
  // Create temporary permuted versions of the input and output.
  // (Re)allocate X_ and/or Y_ only if necessary.  We'll use them to
  // store the permuted versions of X resp. Y.  Note that X_local has
  // the domain Map of the operator, which may be a permuted subset of
  // the local Map corresponding to X.getMap().  Similarly, Y_local
  // has the range Map of the operator, which may be a permuted subset
  // of the local Map corresponding to Y.getMap().  numRows_ here
  // gives the number of rows in the row Map of the local Inverse_
  // operator.
  //
  // FIXME (mfh 20 Aug 2013) There might be an implicit assumption
  // here that the row Map and the range Map of that operator are
  // the same.
  //
  // FIXME (mfh 20 Aug 2013) This "local permutation" functionality
  // really belongs in Tpetra.

  if(X_localBlocks_.size() == 0 || X.extent(1) != X_local_.extent(1))
  {
    //need to resize (or create for the first time) the three scratch arrays
    X_localBlocks_.clear();
    Y_localBlocks_.clear();
    X_localBlocks_.reserve(this->numBlocks_);
    Y_localBlocks_.reserve(this->numBlocks_);

    X_local_ = HostViewLocal("X_local", this->blockRows_.size() * this->scalarsPerRow_, numVecs);
    Y_local_ = HostViewLocal("Y_local", this->blockRows_.size() * this->scalarsPerRow_, numVecs);

    //create all X_local and Y_local managed Views at once, are
    //reused in subsequent apply() calls
    for(int i = 0; i < this->numBlocks_; i++)
    {
      auto blockBounds = std::make_pair(this->blockOffsets_[i] * this->scalarsPerRow_,
          (this->blockOffsets_[i] + this->blockSizes_[i]) * this->scalarsPerRow_);
      X_localBlocks_.emplace_back(X_local_, blockBounds, Kokkos::ALL());
      Y_localBlocks_.emplace_back(Y_local_, blockBounds, Kokkos::ALL());
    }
  }

  const ArrayView<const local_ordinal_type> localRows = this->getBlockRows(blockIndex);

  if(this->scalarsPerRow_ == 1)
    mvgs.gatherViewToView (X_localBlocks_[blockIndex], X, localRows);
  else
    mvgs.gatherViewToViewBlock (X_localBlocks_[blockIndex], X, localRows, this->scalarsPerRow_);

  // We must gather the contents of the output multivector Y even on
  // input to solveBlock(), since the inverse operator might use it as
  // an initial guess for a linear solve.  We have no way of knowing
  // whether it does or does not.

  if(this->scalarsPerRow_ == 1)
    mvgs.gatherViewToView (Y_localBlocks_[blockIndex], Y, localRows);
  else
    mvgs.gatherViewToViewBlock (Y_localBlocks_[blockIndex], Y, localRows, this->scalarsPerRow_);

  // Apply the local operator:
  // Y_local := beta*Y_local + alpha*M^{-1}*X_local
  this->solveBlock (X_localBlocks_[blockIndex], Y_localBlocks_[blockIndex], blockIndex, mode,
                   as<local_scalar_type>(alpha), as<local_scalar_type>(beta));

  // Scatter the permuted subset output vector Y_local back into the
  // original output multivector Y.
  if(this->scalarsPerRow_ == 1)
    mvgs.scatterViewToView (Y, Y_localBlocks_[blockIndex], localRows);
  else
    mvgs.scatterViewToViewBlock (Y, Y_localBlocks_[blockIndex], localRows, this->scalarsPerRow_);
}

template<class MatrixType, class LocalScalarType>
void ContainerImpl<MatrixType, LocalScalarType>::
weightedApply(HostView& X,
              HostView& Y,
              HostView& D,
              int blockIndex,
              Teuchos::ETransp mode,
              scalar_type alpha,
              scalar_type beta) const
{
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::Range1D;
  using Teuchos::Ptr;
  using Teuchos::ptr;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using std::endl;
  typedef Teuchos::ScalarTraits<scalar_type> STS;

  // The local operator template parameter might have a different
  // Scalar type than MatrixType.  This means that we might have to
  // convert X and Y to the Tpetra::MultiVector specialization that
  // the local operator wants.  This class' X_ and Y_ internal fields
  // are of the right type for the local operator, so we can use those
  // as targets.

  const char prefix[] = "Ifpack2::Container::weightedApply: ";
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! this->IsComputed_, std::runtime_error, prefix << "You must have called the "
    "compute() method before you may call this method.  You may call "
    "weightedApply() as many times as you want after calling compute() once, "
    "but you must have called compute() at least once first.");
  
  //bmk 7-2019: BlockRelaxation already checked this, but if that changes...
  TEUCHOS_TEST_FOR_EXCEPTION(
    this->scalarsPerRow_ > 1, std::logic_error, prefix << "Use of block rows isn't allowed "
    "in overlapping Jacobi (the only method that uses weightedApply");

  const size_t numVecs = X.extent(1);

  TEUCHOS_TEST_FOR_EXCEPTION(
    X.extent(1) != Y.extent(1), std::runtime_error,
    prefix << "X and Y have different numbers of vectors (columns).  X has "
    << X.extent(1) << ", but Y has " << Y.extent(1) << ".");

  if(numVecs == 0) {
    return; // done! nothing to do
  }

  const size_t numRows = this->blockSizes_[blockIndex];

  // The local operator works on a permuted subset of the local parts
  // of X and Y.  The subset and permutation are defined by the index
  // array returned by getBlockRows().  If the permutation is trivial
  // and the subset is exactly equal to the local indices, then we
  // could use the local parts of X and Y exactly, without needing to
  // permute.  Otherwise, we have to use temporary storage to permute
  // X and Y.  For now, we always use temporary storage.
  //
  // Create temporary permuted versions of the input and output.
  // (Re)allocate X_ and/or Y_ only if necessary.  We'll use them to
  // store the permuted versions of X resp. Y.  Note that X_local has
  // the domain Map of the operator, which may be a permuted subset of
  // the local Map corresponding to X.getMap().  Similarly, Y_local
  // has the range Map of the operator, which may be a permuted subset
  // of the local Map corresponding to Y.getMap().  numRows_ here
  // gives the number of rows in the row Map of the local operator.
  //
  // FIXME (mfh 20 Aug 2013) There might be an implicit assumption
  // here that the row Map and the range Map of that operator are
  // the same.
  //
  // FIXME (mfh 20 Aug 2013) This "local permutation" functionality
  // really belongs in Tpetra.
  if(X_localBlocks_.size() == 0 || X.extent(1) != X_local_.extent(1))
  {
    //need to resize (or create for the first time) the three scratch arrays
    X_localBlocks_.clear();
    Y_localBlocks_.clear();
    X_localBlocks_.reserve(this->numBlocks_);
    Y_localBlocks_.reserve(this->numBlocks_);

    X_local_ = HostViewLocal("X_local", this->blockRows_.size() * this->scalarsPerRow_, numVecs);
    Y_local_ = HostViewLocal("Y_local", this->blockRows_.size() * this->scalarsPerRow_, numVecs);

    //create all X_local and Y_local managed Views at once, are
    //reused in subsequent apply() calls
    for(int i = 0; i < this->numBlocks_; i++)
    {
      auto blockBounds = std::make_pair(this->blockOffsets_[i] * this->scalarsPerRow_,
          (this->blockOffsets_[i] + this->blockSizes_[i]) * this->scalarsPerRow_);
      X_localBlocks_.emplace_back(X_local_, blockBounds, Kokkos::ALL());
      Y_localBlocks_.emplace_back(Y_local_, blockBounds, Kokkos::ALL());
    }
  }
  if((int) weightedApplyScratch_.extent(0) != 3 * this->maxBlockSize_ ||
      weightedApplyScratch_.extent(1) != numVecs)
  {
    weightedApplyScratch_ = HostViewLocal("weightedApply scratch", 3 * this->maxBlockSize_, numVecs);
  }

  ArrayView<const local_ordinal_type> localRows = this->getBlockRows(blockIndex);

  Details::MultiVectorLocalGatherScatter<mv_type, local_mv_type> mvgs;

  //note: BlockCrs w/ weighted Jacobi isn't allowed, so no need to use block gather/scatter
  mvgs.gatherViewToView (X_localBlocks_[blockIndex], X, localRows);
  // We must gather the output multivector Y even on input to
  // solveBlock(), since the local operator might use it as an initial
  // guess for a linear solve.  We have no way of knowing whether it
  // does or does not.

  mvgs.gatherViewToView (Y_localBlocks_[blockIndex], Y, localRows);

  // Apply the diagonal scaling D to the input X.  It's our choice
  // whether the result has the original input Map of X, or the
  // permuted subset Map of X_local.  If the latter, we also need to
  // gather D into the permuted subset Map.  We choose the latter, to
  // save memory and computation.  Thus, we do the following:
  //
  // 1. Gather D into a temporary vector D_local.
  // 2. Create a temporary X_scaled to hold diag(D_local) * X_local.
  // 3. Compute X_scaled := diag(D_loca) * X_local.
  auto maxBS = this->maxBlockSize_;
  auto bs = this->blockSizes_[blockIndex] * this->scalarsPerRow_;

  HostSubview D_local(weightedApplyScratch_, std::make_pair(0, bs), std::make_pair(0, 1));
  mvgs.gatherViewToView (D_local, D, localRows);
  HostSubview X_scaled(weightedApplyScratch_, std::make_pair(maxBS, maxBS + bs), Kokkos::ALL());
  for(size_t j = 0; j < numVecs; j++)
    for(size_t i = 0; i < numRows; i++)
      X_scaled(i, j) = X_localBlocks_[blockIndex](i, j) * D_local(i, 0);

  HostSubview Y_temp(weightedApplyScratch_, std::make_pair(maxBS * 2, maxBS * 2 + bs), Kokkos::ALL());
  // Apply the local operator: Y_temp := M^{-1} * X_scaled
  this->solveBlock (X_scaled, Y_temp, blockIndex, mode, STS::one(), STS::zero());
  // Y_local := beta * Y_local + alpha * diag(D_local) * Y_temp.
  //
  // Note that we still use the permuted subset scaling D_local here,
  // because Y_temp has the same permuted subset Map.  That's good, in
  // fact, because it's a subset: less data to read and multiply.
  local_impl_scalar_type a = alpha;
  local_impl_scalar_type b = beta;
  for(size_t j = 0; j < numVecs; j++)
    for(size_t i = 0; i < numRows; i++)
      Y_localBlocks_[blockIndex](i, j) = b * Y_localBlocks_[blockIndex](i, j) + a * Y_temp(i, j) * D_local(i, 0);

  // Copy the permuted subset output vector Y_local into the original
  // output multivector Y.
  mvgs.scatterViewToView (Y, Y_localBlocks_[blockIndex], localRows);
}

template<class MatrixType>
void Container<MatrixType>::
clearBlocks()
{
  numBlocks_ = 0;
  blockRows_.clear();
  blockSizes_.clear();
  blockOffsets_.clear();
  Diag_ = Teuchos::null;      //Diag_ will be recreated if needed
}

template<class MatrixType, class LocalScalarType>
void ContainerImpl<MatrixType, LocalScalarType>::
clearBlocks()
{
  X_localBlocks_.clear();
  Y_localBlocks_.clear();
  X_local_ = HostViewLocal();
  Y_local_ = HostViewLocal();
  Container<MatrixType>::clearBlocks();
}

} // namespace Ifpack2

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
