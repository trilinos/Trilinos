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

#ifndef IFPACK2_TRIDICONTAINER_DEF_HPP
#define IFPACK2_TRIDICONTAINER_DEF_HPP

#include "Ifpack2_TriDiContainer_decl.hpp"
#include "Teuchos_LAPACK.hpp"

#ifdef HAVE_MPI
#  include <mpi.h>
#  include "Teuchos_DefaultMpiComm.hpp"
#else
#  include "Teuchos_DefaultSerialComm.hpp"
#endif // HAVE_MPI


namespace Ifpack2 {

template<class MatrixType, class LocalScalarType>
TriDiContainer<MatrixType, LocalScalarType, true>::
TriDiContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
                const Teuchos::RCP<const import_type>& importer,
                int OverlapLevel,
                scalar_type DampingFactor) :
  Container<MatrixType> (matrix, partitions, importer, OverlapLevel,
                         DampingFactor),
  ipiv_ (this->partitions_.size(), 0),
  IsInitialized_ (false),
  IsComputed_ (false),
  scalars_ (nullptr),
  scalarOffsets_ (this->numBlocks_)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::toString;
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! matrix->hasColMap (), std::invalid_argument, "Ifpack2::TriDiContainer: "
    "The constructor's input matrix must have a column Map.");

  // Check whether the input set of local row indices is correct.
  const map_type& rowMap = * (matrix->getRowMap ());
  {
    for(int i = 0; i < this->numBlocks_; i++)
    {
      Teuchos::ArrayView<const local_ordinal_type> localRows = this->getLocalRows(i);
      for(local_ordinal_type j = 0; j < this->blockRows_[i]; j++)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(
          !rowMap.isNodeLocalElement(this->partitions_[this->partitionIndices_[i] + j]),
          std::invalid_argument, "Ifpack2::TriDiContainer: "
          "On process " << rowMap.getComm()->getRank() << " of "
          << rowMap.getComm()->getSize() << ", in the given set of local row "
          "indices localRows = " << Teuchos::toString(localRows) << ", the following "
          "entries is not valid local row index on the calling process: "
          << localRows[j] << ".");
      }
    }
  }

  // FIXME (mfh 25 Aug 2013) What if the matrix's row Map has a
  // different index base than zero?
  //compute scalar array offsets (probably different from partitionIndices_)
  local_ordinal_type scalarTotal = 0;
  for(local_ordinal_type i = 0; i < this->numBlocks_; i++)
  {
    scalarOffsets_[i] = scalarTotal;
    if(this->blockRows_[i] == 1)
      scalarTotal++;
    else
      scalarTotal += 4 * (this->blockRows_[i] - 1);
  }
  //Allocate scalar arrays
  scalars_ = new local_scalar_type[scalarTotal];
  diagBlocks_.reserve(this->numBlocks_);
  for(int i = 0; i < this->numBlocks_; i++)
    diagBlocks_.emplace_back(Teuchos::View, scalars_ + scalarOffsets_[i], this->blockRows_[i]);
}

template<class MatrixType, class LocalScalarType>
TriDiContainer<MatrixType, LocalScalarType, true>::
TriDiContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                const Teuchos::Array<local_ordinal_type>& localRows) :
  Container<MatrixType> (matrix, localRows),
  ipiv_ (this->partitions_.size(), 0),
  IsInitialized_ (false),
  IsComputed_ (false),
  scalars_ (nullptr)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::toString;
  TEUCHOS_TEST_FOR_EXCEPTION(
    !matrix->hasColMap(), std::invalid_argument, "Ifpack2::TriDiContainer: "
    "The constructor's input matrix must have a column Map.");

  // Check whether the input set of local row indices is correct.
  const map_type& rowMap = *(matrix->getRowMap());
  {
    for(local_ordinal_type j = 0; j < this->blockRows_[0]; j++)
    {
      //Check that all rows in only block are valid and locally owned
      TEUCHOS_TEST_FOR_EXCEPTION(
        !rowMap.isNodeLocalElement(this->partitions_[this->partitionIndices_[0] + j]),
        std::invalid_argument, "Ifpack2::TriDiContainer: "
        "On process " << rowMap.getComm ()->getRank () << " of "
        << rowMap.getComm ()->getSize () << ", in the given set of local row "
        "indices localRows = " << Teuchos::toString (localRows) << ", the following "
        "entries is not valid local row index on the calling process: "
        << localRows[j] << ".");
    }
  }
  // FIXME (mfh 25 Aug 2013) What if the matrix's row Map has a
  // different index base than zero?
  //for single block, let the SerialTriDiMat own the scalar memory, as there would be no speed gain
  diagBlocks_.emplace_back(this->blockRows_[0], this->blockRows_[0], true);
}

template<class MatrixType, class LocalScalarType>
TriDiContainer<MatrixType, LocalScalarType, true>::~TriDiContainer ()
{
  if(scalars_)
    delete[] scalars_;
}

template<class MatrixType, class LocalScalarType>
bool TriDiContainer<MatrixType, LocalScalarType, true>::isInitialized () const
{
  return IsInitialized_;
}

template<class MatrixType, class LocalScalarType>
bool TriDiContainer<MatrixType, LocalScalarType, true>::isComputed () const
{
  return IsComputed_;
}

template<class MatrixType, class LocalScalarType>
void TriDiContainer<MatrixType, LocalScalarType, true>::
setParameters (const Teuchos::ParameterList& /* List */)
{
  // the solver doesn't currently take any parameters
}

template<class MatrixType, class LocalScalarType>
void TriDiContainer<MatrixType, LocalScalarType, true>::initialize ()
{
  for(int i = 0; i < this->numBlocks_; i++)
    diagBlocks_[i].putScalar(Teuchos::ScalarTraits<local_scalar_type>::zero());
  std::fill(ipiv_.begin(), ipiv_.end(), 0);
  IsInitialized_ = true;
  // We assume that if you called this method, you intend to recompute
  // everything.
  IsComputed_ = false;
}

template<class MatrixType, class LocalScalarType>
void TriDiContainer<MatrixType, LocalScalarType, true>::compute ()
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    ipiv_.size () != this->partitions_.size(), std::logic_error,
    "Ifpack2::TriDiContainer::compute: ipiv_ array has the wrong size.  "
    "Please report this bug to the Ifpack2 developers.");

  IsComputed_ = false;
  if (! this->isInitialized ()) {
    this->initialize ();
  }

  // Extract the submatrix.
  extract ();
  factor (); // factor the submatrix

  IsComputed_ = true;
}

template<class MatrixType, class LocalScalarType>
void TriDiContainer<MatrixType, LocalScalarType, true>::clearBlocks ()
{
  std::vector<HostViewLocal> empty1;
  std::swap(empty1, X_local);
  std::vector<HostViewLocal> empty2;
  std::swap(empty2, Y_local);
  Container<MatrixType>::clearBlocks();
}

template<class MatrixType, class LocalScalarType>
void TriDiContainer<MatrixType, LocalScalarType, true>::factor ()
{
  for(int i = 0; i < this->numBlocks_; i++)
  {
    Teuchos::LAPACK<int, local_scalar_type> lapack;
    int INFO = 0;
    int* blockIpiv = (int*) ipiv_.getRawPtr() + this->partitionIndices_[i];
    lapack.GTTRF (diagBlocks_[i].numRowsCols (),
                  diagBlocks_[i].DL(),
                  diagBlocks_[i].D(),
                  diagBlocks_[i].DU(),
                  diagBlocks_[i].DU2(),
                  blockIpiv, &INFO);
    // INFO < 0 is a bug.
    TEUCHOS_TEST_FOR_EXCEPTION(
      INFO < 0, std::logic_error, "Ifpack2::TriDiContainer::factor: "
      "LAPACK's _GTTRF (LU factorization with partial pivoting) was called "
      "incorrectly.  INFO = " << INFO << " < 0.  "
      "Please report this bug to the Ifpack2 developers.");
    // INFO > 0 means the matrix is singular.  This is probably an issue
    // either with the choice of rows the rows we extracted, or with the
    // input matrix itself.
    TEUCHOS_TEST_FOR_EXCEPTION(
      INFO > 0, std::runtime_error, "Ifpack2::TriDiContainer::factor: "
      "LAPACK's _GTTRF (LU factorization with partial pivoting) reports that the "
      "computed U factor is exactly singular.  U(" << INFO << "," << INFO << ") "
      "(one-based index i) is exactly zero.  This probably means that the input "
      "matrix has a singular diagonal block.");
  }
}

template<class MatrixType, class LocalScalarType>
void TriDiContainer<MatrixType, LocalScalarType, true>::
applyImpl (HostViewLocal& X,
           HostViewLocal& Y,
           int blockIndex,
           int stride,
           Teuchos::ETransp mode,
           local_scalar_type alpha,
           local_scalar_type beta) const
{
  typedef Teuchos::ScalarTraits<local_scalar_type> STS;
  auto zero = STS::zero();
  size_t numVecs = X.extent(1);
  size_t numRows = X.extent(0);

  TEUCHOS_TEST_FOR_EXCEPTION(
    X.extent (0) != Y.extent (0),
    std::logic_error, "Ifpack2::TriDiContainer::applyImpl: X and Y have "
    "incompatible dimensions (" << X.extent (0) << " resp. "
    << Y.extent (0) << ").  Please report this bug to "
    "the Ifpack2 developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.extent (0) != static_cast<size_t> (diagBlocks_[blockIndex].numRowsCols()),
    std::logic_error, "Ifpack2::TriDiContainer::applyImpl: The input "
    "multivector X has incompatible dimensions from those of the "
    "inverse operator (" << X.extent (0) << " vs. "
    << (mode == Teuchos::NO_TRANS ? diagBlocks_[blockIndex].numRowsCols () : diagBlocks_[blockIndex].numRowsCols())
    << ").  Please report this bug to the Ifpack2 developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    Y.extent (0) != static_cast<size_t> (diagBlocks_[blockIndex].numRowsCols()),
    std::logic_error, "Ifpack2::TriDiContainer::applyImpl: The output "
    "multivector Y has incompatible dimensions from those of the "
    "inverse operator (" << Y.extent (0) << " vs. "
    << (mode == Teuchos::NO_TRANS ? diagBlocks_[blockIndex].numRowsCols() : diagBlocks_[blockIndex].numRowsCols ())
    << ").  Please report this bug to the Ifpack2 developers.");

  if(alpha == zero) { // don't need to solve the linear system
    if(beta == zero) {
      // Use BLAS AXPY semantics for beta == 0: overwrite, clobbering
      // any Inf or NaN values in Y (rather than multiplying them by
      // zero, resulting in NaN values).
      for(size_t j = 0; j < Y.extent(1); j++)
        for(size_t i = 0; i < Y.extent(0); i++)
          Y(i, j) = zero;
    }
    else { // beta != 0
      for(size_t j = 0; j < Y.extent(1); j++)
        for(size_t i = 0; i < Y.extent(0); i++)
          Y(i, j) *= beta;
    }
  }
  else { // alpha != 0; must solve the linear system
    Teuchos::LAPACK<int, local_scalar_type> lapack;
    // If beta is nonzero or Y is not constant stride, we have to use
    // a temporary output multivector.  It gets a copy of X, since
    // GETRS overwrites its (multi)vector input with its output.
    HostViewLocal Y_tmp("", numRows, numVecs);
    Kokkos::deep_copy(Y_tmp, X);
    scalar_type* Y_ptr = Y_tmp.data();
    int INFO = 0;
    const char trans =
      (mode == Teuchos::CONJ_TRANS ? 'C' : (mode == Teuchos::TRANS ? 'T' : 'N'));
    int* blockIpiv = (int*) ipiv_.getRawPtr() + this->partitionIndices_[blockIndex];
    lapack.GTTRS (trans,
                  diagBlocks_[blockIndex].numRowsCols(),
                  numVecs,
                  diagBlocks_[blockIndex].DL(),
                  diagBlocks_[blockIndex].D(),
                  diagBlocks_[blockIndex].DU(),
                  diagBlocks_[blockIndex].DU2(),
                  blockIpiv,
                  Y_ptr,
                  stride,
                  &INFO);
    TEUCHOS_TEST_FOR_EXCEPTION(
      INFO != 0, std::runtime_error, "Ifpack2::TriDiContainer::applyImpl: "
      "LAPACK's _GETRS (solve using LU factorization with partial pivoting) "
      "failed with INFO = " << INFO << " != 0.");

    if (beta != STS::zero ()) {
      for(size_t j = 0; j < Y.extent(1); j++)
      {
        for(size_t i = 0; i < Y.extent(0); i++)
        {
          Y(i, j) *= beta;
          Y(i, j) += alpha * Y_tmp(i, j);
        }
      }
    }
    else {
      for(size_t j = 0; j < Y.extent(1); j++)
      {
        for(size_t i = 0; i < Y.extent(0); i++)
          Y(i, j) = Y_tmp(i, j);
      }
    }
  }
}

template<class MatrixType, class LocalScalarType>
void TriDiContainer<MatrixType, LocalScalarType, true>::
apply (HostView& X,
       HostView& Y,
       int blockIndex,
       int stride,
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
    ! IsComputed_, std::runtime_error, "Ifpack2::TriDiContainer::apply: "
    "You must have called the compute() method before you may call apply().  "
    "You may call the apply() method as many times as you want after calling "
    "compute() once, but you must have called compute() at least once.");

  const size_t numVecs = X.extent(1);

  if(numVecs == 0) {
    return; // done! nothing to do
  }

  // The local operator works on a permuted subset of the local parts
  // of X and Y.  The subset and permutation are defined by the index
  // array returned by getLocalRows().  If the permutation is trivial
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

  if(X_local.size() == 0)
  {
    //create all X_local and Y_local managed Views at once, are
    //reused in subsequent apply() calls
    for(int i = 0; i < this->numBlocks_; i++)
    {
      X_local.emplace_back("", this->blockRows_[i], numVecs);
    }
    for(int i = 0; i < this->numBlocks_; i++)
    {
      Y_local.emplace_back("", this->blockRows_[i], numVecs);
    }
  }

  const ArrayView<const local_ordinal_type> localRows = this->getLocalRows(blockIndex);

  mvgs.gatherViewToView (X_local[blockIndex], X, localRows);

  // We must gather the contents of the output multivector Y even on
  // input to applyImpl(), since the inverse operator might use it as
  // an initial guess for a linear solve.  We have no way of knowing
  // whether it does or does not.

  mvgs.gatherViewToView (Y_local[blockIndex], Y, localRows);

  // Apply the local operator:
  // Y_local := beta*Y_local + alpha*M^{-1}*X_local
  this->applyImpl (X_local[blockIndex], Y_local[blockIndex], blockIndex, stride, mode,
                   as<local_scalar_type>(alpha), as<local_scalar_type>(beta));

  // Scatter the permuted subset output vector Y_local back into the
  // original output multivector Y.
  mvgs.scatterViewToView (Y, Y_local[blockIndex], localRows);
}

template<class MatrixType, class LocalScalarType>
void TriDiContainer<MatrixType, LocalScalarType, true>::
weightedApply (HostView& X,
               HostView& Y,
               HostView& D,
               int blockIndex,
               int stride,
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
  using std::cerr;
  using std::endl;
  typedef Teuchos::ScalarTraits<scalar_type> STS;

  // The local operator template parameter might have a different
  // Scalar type than MatrixType.  This means that we might have to
  // convert X and Y to the Tpetra::MultiVector specialization that
  // the local operator wants.  This class' X_ and Y_ internal fields
  // are of the right type for the local operator, so we can use those
  // as targets.

  Details::MultiVectorLocalGatherScatter<mv_type, local_mv_type> mvgs;

  size_t numRows = this->blockRows_[blockIndex];
  size_t numVecs = X.extent(1);

  if(numVecs == 0) {
    return; // done! nothing to do
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! IsComputed_, std::runtime_error, "Ifpack2::TriDiContainer::"
    "weightedApply: You must have called the compute() method before you may "
    "call apply().  You may call the apply() method as many times as you want "
    "after calling compute() once, but you must have called compute() at least "
    "once.");

  // The local operator works on a permuted subset of the local parts
  // of X and Y.  The subset and permutation are defined by the index
  // array returned by getLocalRows().  If the permutation is trivial
  // and the subset is exactly equal to the local indices, then we
  // could use the local parts of X and Y exactly, without needing to
  // permute.  Otherwise, we have to use temporary storage to permute
  // X and Y.  For now, we always use temporary storage.
  //
  // Ensure we have temporary permuted versions of the input and output.
  // Initialize X_ and/or Y_ if necessary.  We'll use them to
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
  if(X_local.size() == 0)
  {
    //create all X_local and Y_local managed Views at once, are
    //reused in subsequent apply() calls
    for(int i = 0; i < this->numBlocks_; i++)
    {
      X_local.emplace_back("", this->blockRows_[i], numVecs);
    }
    for(int i = 0; i < this->numBlocks_; i++)
    {
      Y_local.emplace_back("", this->blockRows_[i], numVecs);
    }
  }

  ArrayView<const local_ordinal_type> localRows = this->getLocalRows(blockIndex);

  mvgs.gatherViewToView (X_local[blockIndex], X, localRows);

  // We must gather the output multivector Y even on input to
  // applyImpl(), since the local operator might use it as an initial
  // guess for a linear solve.  We have no way of knowing whether it
  // does or does not.

  mvgs.gatherViewToView (Y_local[blockIndex], Y, localRows);

  // Apply the diagonal scaling D to the input X.  It's our choice
  // whether the result has the original input Map of X, or the
  // permuted subset Map of X_local.  If the latter, we also need to
  // gather D into the permuted subset Map.  We choose the latter, to
  // save memory and computation.  Thus, we do the following:
  //
  // 1. Gather D into a temporary vector D_local.
  // 2. Create a temporary X_scaled to hold diag(D_local) * X_local.
  // 3. Compute X_scaled := diag(D_loca) * X_local.

  HostViewLocal D_local("", numVecs, numRows);

  mvgs.gatherViewToView (D_local, D, localRows);

  HostViewLocal X_scaled("", numVecs, numRows);

  for(size_t i = 0; i < X_scaled.extent(0); i++) {
    for(size_t j = 0; j < X_scaled.extent(1); j++) {
      X_scaled(i, j) = X_local[blockIndex](i, j) * D_local(0, j);
    }
  }

  // Y_temp will hold the result of M^{-1}*X_scaled.  If beta == 0, we
  // can write the result of Inverse_->apply() directly to Y_local, so
  // Y_temp may alias Y_local.  Otherwise, if beta != 0, we need
  // temporary storage for M^{-1}*X_scaled, so Y_temp must be
  // different than Y_local.
  HostViewLocal Y_temp("", Y.extent(0), Y.extent(1));

  // Apply the local operator: Y_tmp := M^{-1} * X_scaled
  applyImpl(X_scaled, Y_temp, blockIndex, stride, mode, STS::one(), STS::zero());
  // Y_local := beta * Y_local + alpha * diag(D_local) * Y_temp.
  //
  // Note that we still use the permuted subset scaling D_local here,
  // because Y_temp has the same permuted subset Map.  That's good, in
  // fact, because it's a subset: less data to read and multiply.
  for(size_t i = 0; i < Y.extent(0); i++) {
    for(size_t j = 0; j < Y.extent(1); j++) {
      Y_local[blockIndex](i, j) *= beta;
      Y_local[blockIndex](i, j) += alpha * D_local(i, 0) * Y_temp(i, j);
    }
  }

  // Copy the permuted subset output vector Y_local into the original
  // output multivector Y.
  mvgs.scatterViewToView (Y, Y_local[blockIndex], localRows);
}

template<class MatrixType, class LocalScalarType>
std::ostream& TriDiContainer<MatrixType, LocalScalarType, true>::print(std::ostream& os) const
{
  Teuchos::FancyOStream fos(Teuchos::rcp(&os,false));
  fos.setOutputToRootOnly(0);
  describe(fos);
  return(os);
}

template<class MatrixType, class LocalScalarType>
std::string TriDiContainer<MatrixType, LocalScalarType, true>::description() const
{
  std::ostringstream oss;
  oss << Teuchos::Describable::description();
  if (isInitialized()) {
    if (isComputed()) {
      oss << "{status = initialized, computed";
    }
    else {
      oss << "{status = initialized, not computed";
    }
  }
  else {
    oss << "{status = not initialized, not computed";
  }

  oss << "}";
  return oss.str();
}

template<class MatrixType, class LocalScalarType>
void
TriDiContainer<MatrixType, LocalScalarType, true>::
describe (Teuchos::FancyOStream& os,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using std::endl;
  if(verbLevel==Teuchos::VERB_NONE) return;
  os << "================================================================================" << endl;
  os << "Ifpack2::TriDiContainer" << endl;
  os << "Number of blocks        = " << this->numBlocks_ << endl;
  os << "isInitialized()         = " << IsInitialized_ << endl;
  os << "isComputed()            = " << IsComputed_ << endl;
  os << "================================================================================" << endl;
  os << endl;
}

template<class MatrixType, class LocalScalarType>
void
TriDiContainer<MatrixType, LocalScalarType, true>::
extract ()
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::toString;
  auto& A = *this->inputMatrix_;
  const size_t inputMatrixNumRows = A.getNodeNumRows();
  // We only use the rank of the calling process and the number of MPI
  // processes for generating error messages.  Extraction itself is
  // entirely local to each participating MPI process.
  const int myRank = A.getRowMap()->getComm()->getRank();
  const int numProcs = A.getRowMap()->getComm()->getSize();

  // Sanity check that the local row indices to extract fall within
  // the valid range of local row indices for the input matrix.
  for(int i = 0; i < this->numBlocks_; i++)
  {
    const local_ordinal_type numRows_ = this->blockRows_[i];
    Teuchos::ArrayView<const local_ordinal_type> localRows = this->getLocalRows(i);
    for(local_ordinal_type j = 0; j < numRows_; j++)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(
        localRows[j] < 0 ||
        static_cast<size_t> (localRows[j]) >= inputMatrixNumRows,
        std::runtime_error, "Ifpack2::TriDiContainer::extract: On process " <<
        myRank << " of " << numProcs << ", localRows[j=" << j << "] = " <<
        localRows[j] << ", which is out of the valid range of local row indices "
        "indices [0, " << (inputMatrixNumRows - 1) << "] for the input matrix.");
    }

    // Convert the local row indices we want into local column indices.
    // For every local row ii_local = localRows[i] we take, we also want
    // to take the corresponding column.  To find the corresponding
    // column, we use the row Map to convert the local row index
    // ii_local into a global index ii_global, and then use the column
    // Map to convert ii_global into a local column index jj_local.  If
    // the input matrix doesn't have a column Map, we need to be using
    // global indices anyway...

    // We use the domain Map to exclude off-process global entries.
    const map_type& globalRowMap = *(A.getRowMap());
    const map_type& globalColMap = *(A.getColMap());
    const map_type& globalDomMap = *(A.getDomainMap());

    bool rowIndsValid = true;
    bool colIndsValid = true;
    Array<local_ordinal_type> localCols (numRows_);
    // For error messages, collect the sets of invalid row indices and
    // invalid column indices.  They are otherwise not useful.
    Array<local_ordinal_type> invalidLocalRowInds;
    Array<global_ordinal_type> invalidGlobalColInds;
    for (local_ordinal_type j = 0; j < numRows_; j++)
    {
      // ii_local is the (local) row index we want to look up.
      const local_ordinal_type ii_local = localRows[j];
      // Find the global index jj_global corresponding to ii_local.
      // Global indices are the same (rather, are required to be the
      // same) in all three Maps, which is why we use jj (suggesting a
      // column index, which is how we will use it below).
      const global_ordinal_type jj_global = globalRowMap.getGlobalElement(ii_local);
      if(jj_global == Teuchos::OrdinalTraits<global_ordinal_type>::invalid())
      {
        // If ii_local is not a local index in the row Map on the
        // calling process, that means localRows is incorrect.  We've
        // already checked for this in the constructor, but we might as
        // well check again here, since it's cheap to do so (just an
        // integer comparison, since we need jj_global anyway).
        rowIndsValid = false;
        invalidLocalRowInds.push_back(ii_local);
        break;
      }
      // Exclude "off-process" entries: that is, those in the column Map
      // on this process that are not in the domain Map on this process.
      if(globalDomMap.isNodeGlobalElement (jj_global))
      {
        // jj_global is not an off-process entry.  Look up its local
        // index in the column Map; we want to extract this column index
        // from the input matrix.  If jj_global is _not_ in the column
        // Map on the calling process, that could mean that the column
        // in question is empty on this process.  That would be bad for
        // solving linear systems with the extract submatrix.  We could
        // solve the resulting singular linear systems in a minimum-norm
        // least-squares sense, but for now we simply raise an exception.
        const local_ordinal_type jj_local = globalColMap.getLocalElement(jj_global);
        if(jj_local == Teuchos::OrdinalTraits<local_ordinal_type>::invalid())
        {
          colIndsValid = false;
          invalidGlobalColInds.push_back(jj_global);
          break;
        }
        localCols[j] = jj_local;
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
      !rowIndsValid, std::logic_error, "Ifpack2::TriDiContainer::extract: "
      "On process " << myRank << ", at least one row index in the set of local "
      "row indices given to the constructor is not a valid local row index in "
      "the input matrix's row Map on this process.  This should be impossible "
      "because the constructor checks for this case.  Here is the complete set "
      "of invalid local row indices: " << toString (invalidLocalRowInds) << ".  "
      "Please report this bug to the Ifpack2 developers.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      !colIndsValid, std::runtime_error, "Ifpack2::TriDiContainer::extract: "
      "On process " << myRank << ", "
      "At least one row index in the set of row indices given to the constructor "
      "does not have a corresponding column index in the input matrix's column "
      "Map.  This probably means that the column(s) in question is/are empty on "
      "this process, which would make the submatrix to extract structurally "
      "singular.  Here is the compete set of invalid global column indices: "
      << toString (invalidGlobalColInds) << ".");

    diagBlocks_[i].putScalar(Teuchos::ScalarTraits<local_scalar_type>::zero());
    const size_t maxNumEntriesInRow = A.getNodeMaxNumRowEntries();
    Array<scalar_type> val(maxNumEntriesInRow);
    Array<local_ordinal_type> ind(maxNumEntriesInRow);

    const local_ordinal_type INVALID = Teuchos::OrdinalTraits<local_ordinal_type>::invalid();
    for(local_ordinal_type j = 0; j < numRows_; j++)
    {
      const local_ordinal_type localRow = localRows[j];
      size_t numEntries;
      A.getLocalRowCopy(localRow, ind(), val(), numEntries);

      for(size_t k = 0; k < numEntries; k++)
      {
        const local_ordinal_type localCol = ind[k];

        // Skip off-process elements
        //
        // FIXME (mfh 24 Aug 2013) This assumes the following:
        //
        // 1. The column and row Maps begin with the same set of
        //    on-process entries, in the same order.  That is,
        //    on-process row and column indices are the same.
        // 2. All off-process indices in the column Map of the input
        //    matrix occur after that initial set.
        if(localCol >= 0 && static_cast<size_t> (localCol) < inputMatrixNumRows)
        {
          // for local column IDs, look for each ID in the list
          // of columns hosted by this object
          local_ordinal_type jj = INVALID;
          for (local_ordinal_type kk = 0; kk < numRows_; kk++)
          {
            if(localRows[kk] == localCol)
              jj = kk;
          }
          if (jj != INVALID)
            diagBlocks_[i](j, jj) += val[k];
        }
      }
    }
  }
}

template<class MatrixType, class LocalScalarType>
std::string TriDiContainer<MatrixType, LocalScalarType, true>::getName()
{
  return "TriDi";
}

template<class MatrixType, class LocalScalarType>
TriDiContainer<MatrixType, LocalScalarType, false>::
TriDiContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
                const Teuchos::RCP<const import_type>& importer,
                int OverlapLevel,
                scalar_type DampingFactor) :
  Container<MatrixType> (matrix, partitions, importer, OverlapLevel,
                         DampingFactor)
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Ifpack2::TriDiContainer: Not implemented for "
     "LocalScalarType = " << Teuchos::TypeNameTraits<LocalScalarType>::name ()
     << ".");
}

template<class MatrixType, class LocalScalarType>
TriDiContainer<MatrixType, LocalScalarType, false>::
TriDiContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                const Teuchos::Array<local_ordinal_type>& localRows) :
  Container<MatrixType> (matrix, localRows)
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Ifpack2::TriDiContainer: Not implemented for "
     "LocalScalarType = " << Teuchos::TypeNameTraits<LocalScalarType>::name ()
     << ".");
}

template<class MatrixType, class LocalScalarType>
TriDiContainer<MatrixType, LocalScalarType, false>::~TriDiContainer () {}

template<class MatrixType, class LocalScalarType>
bool TriDiContainer<MatrixType, LocalScalarType, false>::isInitialized () const
{
  return false;
}

template<class MatrixType, class LocalScalarType>
bool TriDiContainer<MatrixType, LocalScalarType, false>::isComputed () const
{
  return false;
}

template<class MatrixType, class LocalScalarType>
void TriDiContainer<MatrixType, LocalScalarType, false>::
setParameters (const Teuchos::ParameterList& /* List */) {}

template<class MatrixType, class LocalScalarType>
void TriDiContainer<MatrixType, LocalScalarType, false>::initialize () {}

template<class MatrixType, class LocalScalarType>
void TriDiContainer<MatrixType, LocalScalarType, false>::compute () {}

template<class MatrixType, class LocalScalarType>
void TriDiContainer<MatrixType, LocalScalarType, false>::clearBlocks () {}

template<class MatrixType, class LocalScalarType>
void TriDiContainer<MatrixType, LocalScalarType, false>::factor () {}

template<class MatrixType, class LocalScalarType>
void TriDiContainer<MatrixType, LocalScalarType, false>::
applyImpl (HostViewLocal& X,
           HostViewLocal& Y,
           int blockIndex,
           int stride,
           Teuchos::ETransp mode,
           local_scalar_type alpha,
           local_scalar_type beta) const {}

template<class MatrixType, class LocalScalarType>
void TriDiContainer<MatrixType, LocalScalarType, false>::
apply (HostView& X,
       HostView& Y,
       int blockIndex,
       int stride,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const {}

template<class MatrixType, class LocalScalarType>
void TriDiContainer<MatrixType, LocalScalarType, false>::
weightedApply (HostView& X,
               HostView& Y,
               HostView& D,
               int blockIndex,
               int stride,
               Teuchos::ETransp mode,
               scalar_type alpha,
               scalar_type beta) const {}

template<class MatrixType, class LocalScalarType>
std::ostream& TriDiContainer<MatrixType, LocalScalarType, false>::print(std::ostream& os) const
{
  return os;
}

template<class MatrixType, class LocalScalarType>
std::string TriDiContainer<MatrixType, LocalScalarType, false>::description() const
{
  return "";
}

template<class MatrixType, class LocalScalarType>
void
TriDiContainer<MatrixType, LocalScalarType, false>::
describe (Teuchos::FancyOStream& os,
          const Teuchos::EVerbosityLevel verbLevel) const {}

template<class MatrixType, class LocalScalarType>
void
TriDiContainer<MatrixType, LocalScalarType, false>::
extract () {}

template<class MatrixType, class LocalScalarType>
std::string TriDiContainer<MatrixType, LocalScalarType, false>::getName()
{
  return "";
}

#define IFPACK2_TRIDICONTAINER_INSTANT(S,LO,GO,N) \
  template class Ifpack2::TriDiContainer< Tpetra::RowMatrix<S, LO, GO, N>, S >;

} // namespace Ifpack2

#endif
