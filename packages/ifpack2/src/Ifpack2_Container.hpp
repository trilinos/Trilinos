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
#include <Tpetra_Map.hpp>
#include <Tpetra_Experimental_BlockCrsMatrix_decl.hpp>
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
class Container : public Teuchos::Describable {
  //! @name Internal typedefs (protected)
  //@{
protected:
  typedef typename MatrixType::scalar_type scalar_type;
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;
  typedef typename MatrixType::node_type node_type;
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> mv_type;
  typedef Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> vector_type;
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type> import_type;
  typedef Partitioner<Tpetra::RowGraph<local_ordinal_type, global_ordinal_type, node_type> > partitioner_type;
  typedef Tpetra::Experimental::BlockCrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> block_crs_matrix_type;
  typedef Tpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> row_matrix_type;

  static_assert(std::is_same<MatrixType, row_matrix_type>::value,
                "Ifpack2::Container: Please use MatrixType = Tpetra::RowMatrix.");

  //! Internal representation of Scalar in Kokkos::View
  typedef typename Kokkos::Details::ArithTraits<scalar_type>::val_type impl_scalar_type;
  //@}

public:
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
             const Teuchos::RCP<const import_type>& importer,
             int OverlapLevel,
             scalar_type DampingFactor) :
    inputMatrix_ (matrix),
    OverlapLevel_ (OverlapLevel),
    DampingFactor_ (DampingFactor),
    Importer_ (importer)
  {
    using Teuchos::Ptr;
    using Teuchos::RCP;
    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::Comm;
    // typedef Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type> import_type; // unused
    NumLocalRows_ = inputMatrix_->getNodeNumRows();
    NumGlobalRows_ = inputMatrix_->getGlobalNumRows();
    NumGlobalNonzeros_ = inputMatrix_->getGlobalNumEntries();
    IsParallel_ = inputMatrix_->getRowMap()->getComm()->getSize() != 1;
    auto global_bcrs =
      Teuchos::rcp_dynamic_cast<const block_crs_matrix_type>(inputMatrix_);
    hasBlockCrs_ = !global_bcrs.is_null();
    if(hasBlockCrs_)
      bcrsBlockSize_ = global_bcrs->getBlockSize();
    else
      bcrsBlockSize_ = 1;
    setBlockSizes(partitions);
  }

  /// \brief Constructor for single block.
  ///
  /// \brief matrix [in] The original input matrix.  This Container
  ///   will construct a local diagonal block from the rows given by
  ///   <tt>localRows</tt>.
  ///
  /// \param localRows [in] The set of (local) rows assigned to this
  ///   container.  <tt>localRows[i] == j</tt>, where i (from 0 to
  ///   <tt>getNumRows() - 1</tt>) indicates the Container's row, and
  ///   j indicates the local row in the calling process.  Subclasses
  ///   must always pass along these indices to the base class.
  Container (const Teuchos::RCP<const row_matrix_type>& matrix,
             const Teuchos::Array<local_ordinal_type>& localRows) :
    inputMatrix_ (matrix),
    numBlocks_ (1),
    partitions_ (localRows.size()),
    blockRows_ (1),
    partitionIndices_ (1),
    OverlapLevel_ (0),
    DampingFactor_ (STS::one()),
    Importer_ (Teuchos::null)
  {
    NumLocalRows_ = inputMatrix_->getNodeNumRows();
    NumGlobalRows_ = inputMatrix_->getGlobalNumRows();
    NumGlobalNonzeros_ = inputMatrix_->getGlobalNumEntries();
    IsParallel_ = inputMatrix_->getRowMap()->getComm()->getSize() > 1;
    blockRows_[0] = localRows.size();
    partitionIndices_[0] = 0;
    const block_crs_matrix_type* global_bcrs =
      Teuchos::rcp_dynamic_cast<const block_crs_matrix_type>(inputMatrix_).get();
    hasBlockCrs_ = global_bcrs;
    if(hasBlockCrs_)
      bcrsBlockSize_ = global_bcrs->getBlockSize();
    else
      bcrsBlockSize_ = 1;
    for(local_ordinal_type i = 0; i < localRows.size(); i++)
      partitions_[i] = localRows[i];
  }

  //! Destructor.
  virtual ~Container() {};

  /// \brief Local indices of the rows of the input matrix that belong to this block.
  ///
  /// The set of (local) rows assigned to this Container is defined by
  /// passing in a set of indices <tt>localRows[i] = j</tt> to the
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
  Teuchos::ArrayView<const local_ordinal_type> getLocalRows(int blockIndex) const
  {
    return Teuchos::ArrayView<const local_ordinal_type>
      (&partitions_[partitionIndices_[blockIndex]], blockRows_[blockIndex]);
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
    blockRows_.resize(numBlocks_);
    partitionIndices_.resize(numBlocks_);
    for(int i = 0; i < numBlocks_; i++)
    {
      local_ordinal_type rowsInBlock = partitions[i].size();
      blockRows_[i] = rowsInBlock;
      partitionIndices_[i] = totalBlockRows;
      totalBlockRows += rowsInBlock;
    }
    partitions_.resize(totalBlockRows);
    //set partitions_: each entry is the partition/block of the row
    local_ordinal_type iter = 0;
    for(int i = 0; i < numBlocks_; i++)
    {
      for(int j = 0; j < blockRows_[i]; j++)
      {
        partitions_[iter++] = partitions[i][j];
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

  /// \brief Extract the local diagonal block and prepare the solver.
  ///
  /// If any entries' values in the input matrix have changed, you
  /// must call this method before you may call apply() or
  /// weightedApply().
  virtual void compute () = 0;

  void DoJacobi(HostView& X, HostView& Y, int stride) const;
  void DoOverlappingJacobi(HostView& X, HostView& Y, HostView& W, int stride) const;
  void DoGaussSeidel(HostView& X, HostView& Y, HostView& Y2, int stride) const;
  void DoSGS(HostView& X, HostView& Y, HostView& Y2, int stride) const;

  //! Set parameters.
  virtual void setParameters (const Teuchos::ParameterList& List) = 0;

  //! Return \c true if the container has been successfully initialized.
  virtual bool isInitialized () const = 0;

  //! Return \c true if the container has been successfully computed.
  virtual bool isComputed () const = 0;

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
  apply (HostView& X,
         HostView& Y,
         int blockIndex,
         int stride,
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
  virtual void applyInverseJacobi (const mv_type& X, mv_type& Y,
                                   bool zeroStartingSolution = false,
                                   int numSweeps = 1) const
  { TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Not implemented."); }

  //! Wrapper for apply with MVs, used in unit tests (never called by BlockRelaxation)
  void applyMV (mv_type& X, mv_type& Y) const
  {
    HostView XView = X.template getLocalView<Kokkos::HostSpace>();
    HostView YView = Y.template getLocalView<Kokkos::HostSpace>();
    this->apply (XView, YView, 0, X.getStride());
  }

  /// \brief Compute <tt>Y := alpha * diag(D) * M^{-1} (diag(D) * X) + beta*Y</tt>.
  ///
  /// X is in the domain Map of the original matrix (the argument to
  /// compute()), and Y is in the range Map of the original matrix.
  /// This method only reads resp. modifies the permuted subset of
  /// entries of X resp. Y related to the diagonal block M.  That
  /// permuted subset is defined by the indices passed into the
  /// constructor.  The D scaling vector must have the same number of
  /// entries on each process as X and Y, but otherwise need not have
  /// the same Map.  (For example, D could be locally replicated, or
  /// could be a different object on each process with a local (\c
  /// MPI_COMM_SELF) communicator.)
  ///
  /// This method supports overlap techniques, such as those used in
  /// Schwarz methods.
  ///
  /// This method is marked \c const by analogy with apply(), which
  /// itself is marked \c const for compatibility with
  /// Tpetra::Operator's method of the same name.  This might require
  /// subclasses to mark some of their instance data as \c mutable.
  virtual void
  weightedApply (HostView& X,
                 HostView& Y,
                 HostView& W,
                 int blockIndex,
                 int stride,
                 Teuchos::ETransp mode = Teuchos::NO_TRANS,
                 scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
                 scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const = 0;

  //! Wrapper for weightedApply with MVs, used in unit tests (never called by BlockRelaxation)
  void weightedApplyMV (mv_type& X,
                        mv_type& Y,
                        vector_type& W)
  {
    HostView XView = X.template getLocalView<Kokkos::HostSpace>();
    HostView YView = Y.template getLocalView<Kokkos::HostSpace>();
    HostView WView = W.template getLocalView<Kokkos::HostSpace>();
    weightedApply (XView, YView, WView, 0, X.getStride());
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
  //! The input matrix to the constructor.
  Teuchos::RCP<const row_matrix_type> inputMatrix_;

  //! The number of blocks (partitions) in the container.
  int numBlocks_;
  //! Local indices of the rows of the input matrix that belong to this block.
  Teuchos::Array<local_ordinal_type> partitions_;      //size: total # of local rows (in all local blocks)
  //! Number of rows in each block.
  Teuchos::Array<local_ordinal_type> blockRows_;     //size: # of blocks
  //! Starting index in partitions_ of local row indices for each block.
  Teuchos::Array<local_ordinal_type> partitionIndices_;
  //! Diagonal elements.
  mutable Teuchos::RCP<vector_type> Diag_;
  //! Whether the problem is distributed across multiple MPI processes.
  bool IsParallel_;
  //! Number of rows of overlap for adjacent blocks.
  int OverlapLevel_;
  //! Damping factor, passed to apply() as alpha.
  scalar_type DampingFactor_;
  //! Importer for importing off-process elements of MultiVectors.
  Teuchos::RCP<const Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type>> Importer_;
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
};

//! Print information about the given Container to the output stream \c os.
template <class MatrixType>
inline std::ostream&
operator<< (std::ostream& os, const Ifpack2::Container<MatrixType>& obj)
{
  return obj.print (os);
}

template <class MatrixType>
void Container<MatrixType>::DoJacobi(HostView& X, HostView& Y, int stride) const
{
  const scalar_type one = STS::one();
  // Note: Flop counts copied naively from Ifpack.
  // use partitions_ and blockRows_
  size_t numVecs = X.extent(1);
  // Non-overlapping Jacobi
  for (local_ordinal_type i = 0; i < numBlocks_; i++)
  {
    // may happen that a partition is empty
    if(blockRows_[i] != 1 || hasBlockCrs_)
    {
      if(blockRows_[i] == 0 )
        continue;
      apply(X, Y, i, stride, Teuchos::NO_TRANS, DampingFactor_, one);
    }
    else    // singleton, can't access Containers_[i] as it was never filled and may be null.
    {
      local_ordinal_type LRID = partitions_[partitionIndices_[i]];
      getMatDiag();
      HostView diagView = Diag_->template getLocalView<Kokkos::HostSpace>();
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
void Container<MatrixType>::DoOverlappingJacobi(HostView& X, HostView& Y, HostView& W, int stride) const
{
  // Overlapping Jacobi
  for(local_ordinal_type i = 0; i < numBlocks_; i++)
  {
    // may happen that a partition is empty
    if(blockRows_[i] == 0)
      continue;
    if(blockRows_[i] != 1)
      weightedApply(X, Y, W, i, stride, Teuchos::NO_TRANS, DampingFactor_, STS::one());
  }
}

template<class MatrixType>
void Container<MatrixType>::DoGaussSeidel(HostView& X, HostView& Y, HostView& Y2, int stride) const
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::Ptr;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  // Note: Flop counts copied naively from Ifpack.
  const scalar_type one = STS::one();
  const size_t Length = inputMatrix_->getNodeMaxNumRowEntries();
  auto numVecs = X.extent(1);
  Array<scalar_type> Values;
  Array<local_ordinal_type> Indices;
  Indices.resize(Length);

  Values.resize(bcrsBlockSize_ * bcrsBlockSize_ * Length);    //note: if A was not a BlockCRS, bcrsBlockSize_ is 1
  // I think I decided I need two extra vectors:
  // One to store the sum of the corrections (initialized to zero)
  // One to store the temporary residual (doesn't matter if it is zeroed or not)
  // My apologies for making the names clear and meaningful. (X=RHS, Y=guess?! Nice.)
  HostView Resid("", X.extent(0), X.extent(1));
  for(local_ordinal_type i = 0; i < numBlocks_; i++)
  {
    if(blockRows_[i] > 1 || hasBlockCrs_)
    {
      if (blockRows_[i] == 0)
        continue;
      // update from previous block
      ArrayView<const local_ordinal_type> localRows = getLocalRows (i);
      const size_t localNumRows = blockRows_[i];
      for(size_t j = 0; j < localNumRows; j++)
      {
        const local_ordinal_type LID = localRows[j]; // Containers_[i]->ID (j);
        size_t NumEntries;
        inputMatrix_->getLocalRowCopy(LID, Indices(), Values(), NumEntries);
        for(size_t m = 0; m < numVecs; m++)
        {
          if(hasBlockCrs_)
          {
            for (int localR = 0; localR < bcrsBlockSize_; localR++)
              Resid(LID * bcrsBlockSize_ + localR, m) = X(LID * bcrsBlockSize_ + localR, m);
            for (size_t k = 0; k < NumEntries; ++k)
            {
              const local_ordinal_type col = Indices[k];
              for (int localR = 0; localR < bcrsBlockSize_; localR++)
              {
                for(int localC = 0; localC < bcrsBlockSize_; localC++)
                {
                  Resid(LID * bcrsBlockSize_ + localR, m) -=
                    Values[k * bcrsBlockSize_ * bcrsBlockSize_ + localR + localC * bcrsBlockSize_]
                    * Y2(col * bcrsBlockSize_ + localC, m);
                }
              }
            }
          }
          else
          {
            Resid(LID, m) = X(LID, m);
            for (size_t k = 0; k < NumEntries; ++k)
            {
              const local_ordinal_type col = Indices[k];
              Resid(LID, m) -= Values[k] * Y2(col, m);
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
      apply(Resid, Y2, i, stride, Teuchos::NO_TRANS, DampingFactor_, one);
    }
    else if(blockRows_[i] == 1)
    {
      // singleton, can't access Containers_[i] as it was never filled and may be null.
      // a singleton calculation is exact, all residuals should be zero.
      local_ordinal_type LRID = partitionIndices_[i];  // by definition, a singleton 1 row in block.
      getMatDiag();
      HostView diagView = Diag_->template getLocalView<Kokkos::HostSpace>();
      impl_scalar_type d = (impl_scalar_type) one / diagView(LRID, 0);
      for(size_t m = 0; m < numVecs; m++)
      {
        impl_scalar_type x = X(LRID, m);
        impl_scalar_type newy = x * d;
        Y2(LRID, m) = newy;
     }
    } // end else
  } // end for numBlocks_
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
void Container<MatrixType>::DoSGS(HostView& X, HostView& Y, HostView& Y2, int stride) const
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::Ptr;
  using Teuchos::ptr;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  const scalar_type one = STS::one();
  auto numVecs = X.extent(1);
  const size_t Length = inputMatrix_->getNodeMaxNumRowEntries();
  Array<scalar_type> Values;
  Array<local_ordinal_type> Indices(Length);
  Values.resize(bcrsBlockSize_ * bcrsBlockSize_ * Length);
  // I think I decided I need two extra vectors:
  // One to store the sum of the corrections (initialized to zero)
  // One to store the temporary residual (doesn't matter if it is zeroed or not)
  // My apologies for making the names clear and meaningful. (X=RHS, Y=guess?! Nice.)
  HostView Resid("", X.extent(0), X.extent(1));
  // Forward Sweep
  for(local_ordinal_type i = 0; i < numBlocks_; i++)
  {
    if(blockRows_[i] != 1 || hasBlockCrs_)
    {
      if(blockRows_[i] == 0)
        continue; // Skip empty partitions
      // update from previous block
      ArrayView<const local_ordinal_type> localRows = getLocalRows(i);
      for(local_ordinal_type j = 0; j < blockRows_[i]; j++)
      {
        const local_ordinal_type LID = localRows[j]; // Containers_[i]->ID (j);
        size_t NumEntries;
        inputMatrix_->getLocalRowCopy(LID, Indices(), Values(), NumEntries);
        //set tmpresid = initresid - A*correction
        for(size_t m = 0; m < numVecs; m++)
        {
          if(hasBlockCrs_)
          {
            for(int localR = 0; localR < bcrsBlockSize_; localR++)
              Resid(LID * bcrsBlockSize_ + localR, m) = X(LID * bcrsBlockSize_ + localR, m);
            for(size_t k = 0; k < NumEntries; ++k)
            {
              const local_ordinal_type col = Indices[k];
              for (int localR = 0; localR < bcrsBlockSize_; localR++)
              {
                for(int localC = 0; localC < bcrsBlockSize_; localC++)
                  Resid(LID * bcrsBlockSize_ + localR, m) -=
                    Values[k * (bcrsBlockSize_ * bcrsBlockSize_) + (localR + localC * bcrsBlockSize_)]
                    * Y2(col * bcrsBlockSize_ + localC, m);
              }
            }
          }
          else
          {
            Resid(LID, m) = X(LID, m);
            for(size_t k = 0; k < NumEntries; k++)
            {
              local_ordinal_type col = Indices[k];
              Resid(LID, m) -= Values[k] * Y2(col, m);
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
      apply(Resid, Y2, i, stride, Teuchos::NO_TRANS, DampingFactor_, one);
      // operations for all getrow's
    }
    else // singleton, can't access Containers_[i] as it was never filled and may be null.
    {
      local_ordinal_type LRID  = partitions_[partitionIndices_[i]];
      getMatDiag();
      HostView diagView = Diag_->template getLocalView<Kokkos::HostSpace>();
      impl_scalar_type d = (impl_scalar_type) one / diagView(LRID, 0);
      for(size_t m = 0; m < numVecs; m++)
      {
        impl_scalar_type x = X(LRID, m);
        impl_scalar_type newy = x * d;
        Y2(LRID, m) = newy;
      }
    } // end else
  } // end forward sweep over NumLocalBlocks
  // Reverse Sweep
  //
  // mfh 12 July 2013: The unusual iteration bounds, and the use of
  // i-1 rather than i in the loop body, ensure correctness even if
  // local_ordinal_type is unsigned.  "i = numBlocks_-1; i >= 0;
  // i--" will loop forever if local_ordinal_type is unsigned, because
  // unsigned integers are (trivially) always nonnegative.
  for(local_ordinal_type i = numBlocks_; i > 0; --i)
  {
    if(hasBlockCrs_ || blockRows_[i-1] != 1)
    {
      if(blockRows_[i - 1] == 0)
        continue;
      // update from previous block
      ArrayView<const local_ordinal_type> localRows = getLocalRows(i-1);
      for(local_ordinal_type j = 0; j < blockRows_[i-1]; j++)
      {
        const local_ordinal_type LID = localRows[j]; // Containers_[i-1]->ID (j);
        size_t NumEntries;
        inputMatrix_->getLocalRowCopy(LID, Indices(), Values(), NumEntries);
        //set tmpresid = initresid - A*correction
        for (size_t m = 0; m < numVecs; m++)
        {
          if(hasBlockCrs_)
          {
            for(int localR = 0; localR < bcrsBlockSize_; localR++)
              Resid(LID * bcrsBlockSize_ + localR, m) = X(LID * bcrsBlockSize_ + localR, m);
            for(size_t k = 0; k < NumEntries; ++k)
            {
              const local_ordinal_type col = Indices[k];
              for(int localR = 0; localR < bcrsBlockSize_; localR++)
              {
                for(int localC = 0; localC < bcrsBlockSize_; localC++)
                  Resid(LID*bcrsBlockSize_+localR, m) -=
                    Values[k * (bcrsBlockSize_ * bcrsBlockSize_) + (localR + localC * bcrsBlockSize_)]
                    * Y2(col * bcrsBlockSize_ + localC, m);
              }
            }
          }
          else
          {
            Resid(LID, m) = X(LID, m);
            for(size_t k = 0; k < NumEntries; ++k)
            {
              local_ordinal_type col = Indices[k];
              Resid(LID, m) -= Values[k] * Y2(col, m);
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
      apply(Resid, Y2, i - 1, stride, Teuchos::NO_TRANS, DampingFactor_, one);
      // operations for all getrow's
    } // end  Partitioner_->numRowsInPart (i) != 1 )
    // else do nothing, as by definition with a singleton, the residuals are zero.
  } //end reverse sweep
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
void Container<MatrixType>::clearBlocks()
{
  numBlocks_ = 0;
  partitions_.clear();
  blockRows_.clear();
  partitionIndices_.clear();
  Diag_ = Teuchos::null;      //Diag_ will be recreated if needed
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
      TypeNameTraits<MatrixType>::name () + " >";
  }

  static std::string concreteName (const ::Ifpack2::Container<MatrixType>&) {
    return name ();
  }
};

} // namespace Teuchos

#endif // IFPACK2_CONTAINER_HPP
