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

#ifndef IFPACK2_BANDEDCONTAINER_DEF_HPP
#define IFPACK2_BANDEDCONTAINER_DEF_HPP

#include "Teuchos_LAPACK.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include <iostream>
#include <sstream>

#ifdef HAVE_MPI
#  include <mpi.h>
#  include "Teuchos_DefaultMpiComm.hpp"
#else
#  include "Teuchos_DefaultSerialComm.hpp"
#endif // HAVE_MPI

namespace Ifpack2 {

template<class MatrixType, class LocalScalarType>
BandedContainer<MatrixType, LocalScalarType, true>::
BandedContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                 const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
                 const Teuchos::RCP<const import_type>& importer,
                 int OverlapLevel,
                 scalar_type DampingFactor) :
  Container<MatrixType>(matrix, partitions, importer, OverlapLevel, DampingFactor),
  ipiv_(this->partitions_.size()),
  kl_(this->numBlocks_, -1),
  ku_(this->numBlocks_, -1),
  scalars_(nullptr),
  scalarOffsets_(this->numBlocks_)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! matrix->hasColMap (), std::invalid_argument, "Ifpack2::BandedContainer: "
    "The constructor's input matrix must have a column Map.");

  // Check whether the input set of local row indices is correct.
  const map_type& rowMap = * (matrix->getRowMap ());
  for(int i = 0; i < this->numBlocks_; i++)
  {
    Teuchos::ArrayView<const local_ordinal_type> localRows = this->getLocalRows(i);
    for(local_ordinal_type j = 0; j < this->blockRows_[i]; j++)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(
        !rowMap.isNodeLocalElement(localRows[j]),
        std::invalid_argument, "Ifpack2::BandedContainer: "
        "On process " << rowMap.getComm ()->getRank () << " of "
        << rowMap.getComm ()->getSize () << ", in the given set of local row "
        "indices localRows = " << Teuchos::toString (localRows) << ", the following "
        "entry is not valid local row indices on the calling process: "
        << localRows[j] << ".");
    }
  }
  IsInitialized_ = false;
  IsComputed_ = false;
}

template<class MatrixType, class LocalScalarType>
BandedContainer<MatrixType, LocalScalarType, true>::
BandedContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                 const Teuchos::Array<local_ordinal_type>& localRows) :
  Container<MatrixType>(matrix, localRows),
  ipiv_(this->blockRows_[0]),
  kl_(1, -1),
  ku_(1, -1),
  scalars_(nullptr),
  scalarOffsets_(1, 0)
{
  TEUCHOS_TEST_FOR_EXCEPTION(!matrix->hasColMap(), std::invalid_argument, "Ifpack2::BandedContainer: "
    "The constructor's input matrix must have a column Map.");

  // Check whether the input set of local row indices is correct.
  const map_type& rowMap = *(matrix->getRowMap());
  for(local_ordinal_type j = 0; j < this->blockRows_[0]; j++)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      !rowMap.isNodeLocalElement(localRows[j]),
      std::invalid_argument, "Ifpack2::BandedContainer: "
      "On process " << rowMap.getComm()->getRank() << " of "
      << rowMap.getComm()->getSize() << ", in the given set of local row "
      "indices localRows = " << Teuchos::toString(localRows) << ", the following "
      "entry is not valid local row indices on the calling process: "
      << localRows[j] << ".");
  }
  IsInitialized_ = false;
  IsComputed_ = false;
}

template<class MatrixType, class LocalScalarType>
BandedContainer<MatrixType, LocalScalarType, true>::
~BandedContainer ()
{
  if(scalars_)
    delete[] scalars_;
}

template<class MatrixType, class LocalScalarType>
void BandedContainer<MatrixType, LocalScalarType, true>::
setParameters (const Teuchos::ParameterList& List)
{
  typedef typename Teuchos::ArrayView<const local_ordinal_type>::size_type size_type;
  if(List.isParameter("relaxation: banded container superdiagonals"))
    kl_[0] = List.get<int>("relaxation: banded container superdiagonals");
  if(List.isParameter("relaxation: banded container subdiagonals"))
    ku_[0] = List.get<int>("relaxation: banded container subdiagonals");

  for(local_ordinal_type b = 1; b < this->numBlocks_; b++)
  {
    kl_[b] = kl_[0];
    ku_[b] = ku_[0];
  }

  // The user provided insufficient information. If this is the case we check for the optimal values.
  // User information may be overwritten only if necessary.
  for(local_ordinal_type b = 0; b < this->numBlocks_; b++)
  {
    if (ku_[b] == -1 || kl_[b] == -1)
    {
      const Teuchos::ArrayView<const local_ordinal_type> localRows = this->getLocalRows(b);
      const size_type numRows = localRows.size();

      // loop over local rows in current block
      for(size_type i = 0; i < numRows; ++i)
      {
        Teuchos::ArrayView<const local_ordinal_type> indices;
        Teuchos::ArrayView<const scalar_type> values;
        this->inputMatrix_->getLocalRowView(localRows[i], indices, values);

        size_type min_col_it = numRows > 0 ? numRows - 1 : 0; // just a guess
        size_type max_col_it = 0;

        size_type cntCols = 0;

        // loop over all column entries
        for(size_type c = 0; c < indices.size(); c++)
        {
          const local_ordinal_type lColIdx = indices[c]; // current column idx
          // check whether lColIdx is contained in localRows[]
          for(size_type j = 0; j < numRows; j++)
          {
            if (localRows[j] == lColIdx)
            {
              if(localRows[min_col_it] > lColIdx)
                min_col_it = j;
              if(localRows[max_col_it] < lColIdx)
                max_col_it = j;
              cntCols++;
            }
          }
          if(cntCols == numRows)
            break; // skip remaining entries in column
        }
        ku_[b] = std::max(ku_[b], Teuchos::as<local_ordinal_type>(max_col_it - i));
        kl_[b] = std::max(kl_[b], Teuchos::as<local_ordinal_type>(i - min_col_it));
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION
      (kl_[b] == -1 || ku_[b] == -1, std::invalid_argument,
       "Ifpack2::BandedContainer::setParameters: the user must provide the number"
       " of sub- and superdiagonals in the 'kl' and 'ku' parameters.");
  }
}

template<class MatrixType, class LocalScalarType>
void
BandedContainer<MatrixType, LocalScalarType, true>::
initialize ()
{
  using Teuchos::null;
  using Teuchos::rcp;
  for(local_ordinal_type b = 0; b < this->numBlocks_; b++)
  {
    TEUCHOS_TEST_FOR_EXCEPTION
      (kl_[b] == -1 || ku_[b] == -1, std::invalid_argument,
       "Ifpack2::BandedContainer::initialize: the user must provide the number of"
       " sub- and superdiagonals in the 'kl' and 'ku' parameters. Make sure that "
       "you call BandedContainer<T>::setParameters!");
  }
  global_ordinal_type totalScalars = 0;
  for(local_ordinal_type b = 0; b < this->numBlocks_; b++)
  {
    local_ordinal_type stride = 2 * kl_[b] + ku_[b] + 1;
    scalarOffsets_[b] = totalScalars;
    totalScalars += stride * this->blockRows_[b];
  }
  scalars_ = new local_scalar_type[totalScalars];
  for(int b = 0; b < this->numBlocks_; b++)
  {
    local_ordinal_type nrows = this->blockRows_[b];
    diagBlocks_.emplace_back(Teuchos::View, scalars_ + scalarOffsets_[b], 2 * kl_[b] + ku_[b] + 1, nrows, nrows, kl_[b], kl_[b] + ku_[b]);
    diagBlocks_[b].putScalar(Teuchos::ScalarTraits<local_scalar_type>::zero());
  }
  // We assume that if you called this method, you intend to recompute
  // everything.
  IsInitialized_ = false;
  IsComputed_ = false;
  std::fill (ipiv_.begin (), ipiv_.end (), 0);
  IsInitialized_ = true;
}

template<class MatrixType, class LocalScalarType>
void
BandedContainer<MatrixType, LocalScalarType, true>::
compute ()
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    ipiv_.size () != this->partitions_.size(), std::logic_error,
    "Ifpack2::BandedContainer::compute: ipiv_ array has the wrong size.  "
    "Please report this bug to the Ifpack2 developers.");

  IsComputed_ = false;
  if (! this->isInitialized ()) {
    this->initialize ();
  }

  // Extract the submatrices from input matrix.
  extract ();
  factor (); // factor the submatrix

  IsComputed_ = true;
}

template<class MatrixType, class LocalScalarType>
void
BandedContainer<MatrixType, LocalScalarType, true>::
clearBlocks ()
{
  std::vector<HostViewLocal> empty1;
  std::swap(empty1, X_local);
  std::vector<HostViewLocal> empty2;
  std::swap(empty2, Y_local);
  Container<MatrixType>::clearBlocks ();
}

template<class MatrixType, class LocalScalarType>
void
BandedContainer<MatrixType, LocalScalarType, true>::
factor ()
{
  Teuchos::LAPACK<int, local_scalar_type> lapack;
  int INFO = 0;

  for(int i = 0; i < this->numBlocks_; i++)
  {
    // Plausibility checks for matrix
    TEUCHOS_TEST_FOR_EXCEPTION(diagBlocks_[i].values()==0, std::invalid_argument,
                       "BandedContainer<T>::factor: Diagonal block is an empty SerialBandDenseMatrix<T>!");
    TEUCHOS_TEST_FOR_EXCEPTION(diagBlocks_[i].upperBandwidth() < diagBlocks_[i].lowerBandwidth(), std::invalid_argument,
                       "BandedContainer<T>::factor: Diagonal block needs kl additional superdiagonals for factorization! However, the number of superdiagonals is smaller than the number of subdiagonals!");
    int* blockIpiv = &ipiv_[this->partitionIndices_[i]];
    lapack.GBTRF (diagBlocks_[i].numRows(),
        diagBlocks_[i].numCols(),
        diagBlocks_[i].lowerBandwidth(),
        diagBlocks_[i].upperBandwidth() - diagBlocks_[i].lowerBandwidth(), /* enter the real number of superdiagonals (see Teuchos_SerialBandDenseSolver)*/
        diagBlocks_[i].values(),
        diagBlocks_[i].stride(),
        blockIpiv,
        &INFO);

    // INFO < 0 is a bug.
    TEUCHOS_TEST_FOR_EXCEPTION(
      INFO < 0, std::logic_error, "Ifpack2::BandedContainer::factor: "
      "LAPACK's _GBTRF (LU factorization with partial pivoting) was called "
      "incorrectly.  INFO = " << INFO << " < 0.  "
      "Please report this bug to the Ifpack2 developers.");
    // INFO > 0 means the matrix is singular.  This is probably an issue
    // either with the choice of rows the rows we extracted, or with the
    // input matrix itself.
    TEUCHOS_TEST_FOR_EXCEPTION(
      INFO > 0, std::runtime_error, "Ifpack2::BandedContainer::factor: "
      "LAPACK's _GBTRF (LU factorization with partial pivoting) reports that the "
      "computed U factor is exactly singular.  U(" << INFO << "," << INFO << ") "
      "(one-based index i) is exactly zero.  This probably means that the input "
      "matrix has a singular diagonal block.");
  }
}

template<class MatrixType, class LocalScalarType>
void
BandedContainer<MatrixType, LocalScalarType, true>::
applyImpl (HostViewLocal& X,
           HostViewLocal& Y,
           int blockIndex,
           int stride,
           Teuchos::ETransp mode,
           const local_scalar_type alpha,
           const local_scalar_type beta) const
{
  using Teuchos::ArrayRCP;
  using Teuchos::Ptr;
  using Teuchos::ptr;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;

  TEUCHOS_TEST_FOR_EXCEPTION(
    X.extent (0) != Y.extent (0),
    std::logic_error, "Ifpack2::BandedContainer::applyImpl: X and Y have "
    "incompatible dimensions (" << X.extent (0) << " resp. "
    << Y.extent (0) << ").  Please report this bug to "
    "the Ifpack2 developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.extent (0) != static_cast<size_t> (mode == Teuchos::NO_TRANS ? diagBlocks_[blockIndex].numCols() : diagBlocks_[blockIndex].numRows()),
    std::logic_error, "Ifpack2::BandedContainer::applyImpl: The input "
    "multivector X has incompatible dimensions from those of the "
    "inverse operator (" << X.extent (0) << " vs. "
    << (mode == Teuchos::NO_TRANS ? diagBlocks_[blockIndex].numCols() : diagBlocks_[blockIndex].numRows())
    << ").  Please report this bug to the Ifpack2 developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    Y.extent (0) != static_cast<size_t> (mode == Teuchos::NO_TRANS ? diagBlocks_[blockIndex].numRows() : diagBlocks_[blockIndex].numCols()),
    std::logic_error, "Ifpack2::BandedContainer::applyImpl: The output "
    "multivector Y has incompatible dimensions from those of the "
    "inverse operator (" << Y.extent (0) << " vs. "
    << (mode == Teuchos::NO_TRANS ? diagBlocks_[blockIndex].numRows() : diagBlocks_[blockIndex].numCols())
    << ").  Please report this bug to the Ifpack2 developers.");

  size_t numVecs = (int) X.extent (1);

  auto zero = Teuchos::ScalarTraits<scalar_type>::zero ();
  if (alpha == zero) { // don't need to solve the linear system
    if (beta == zero) {
      // Use BLAS AXPY semantics for beta == 0: overwrite, clobbering
      // any Inf or NaN values in Y (rather than multiplying them by
      // zero, resulting in NaN values).
      for(size_t j = 0; j < Y.extent(0); j++)
        for(size_t i = 0; i < Y.extent(1); i++)
          Y(i, j) = zero;
    }
    else { // beta != 0
      for(size_t j = 0; j < Y.extent(0); j++)
        for(size_t i = 0; i < Y.extent(1); i++)
          Y(i, j) = beta * (local_impl_scalar_type) Y(i, j);
    }
  }
  else { // alpha != 0; must solve the linear system
    Teuchos::LAPACK<int, local_scalar_type> lapack;
    // If beta is nonzero or Y is not constant stride, we have to use
    // a temporary output multivector.  It gets a copy of X, since
    // GBTRS overwrites its (multi)vector input with its output.
    Ptr<HostViewLocal> Y_tmp;
    bool deleteYT = false;
    if(beta == zero) {
      Y = X;
      Y_tmp = ptr(&Y);
    }
    else {
      Y_tmp = ptr (new HostViewLocal ("", X.extent (0), X.extent (1))); // constructor copies X
      deleteYT = true;
      Kokkos::deep_copy(*Y_tmp, X);
    }

    local_scalar_type* const Y_ptr = (local_scalar_type*) Y_tmp->data();

    int INFO = 0;
    const char trans =(mode == Teuchos::CONJ_TRANS ? 'C' : (mode == Teuchos::TRANS ? 'T' : 'N'));

    const int* blockIpiv = &ipiv_[this->partitionIndices_[blockIndex]];
    lapack.GBTRS(trans,
        diagBlocks_[blockIndex].numCols(),
        diagBlocks_[blockIndex].lowerBandwidth(),
        /* enter the real number of superdiagonals (see Teuchos_SerialBandDenseSolver)*/
        diagBlocks_[blockIndex].upperBandwidth() - diagBlocks_[blockIndex].lowerBandwidth(),
        numVecs,
        diagBlocks_[blockIndex].values(),
        diagBlocks_[blockIndex].stride(),
        blockIpiv,
        Y_ptr, stride, &INFO);

    TEUCHOS_TEST_FOR_EXCEPTION(
      INFO != 0, std::runtime_error, "Ifpack2::BandedContainer::applyImpl: "
      "LAPACK's _GBTRS (solve using LU factorization with partial pivoting) "
      "failed with INFO = " << INFO << " != 0.");

    if (beta != zero) {
      for(size_t j = 0; j < Y.extent(1); j++)
        for(size_t i = 0; i < Y.extent(0); i++)
          Y(i, j) = beta * Y(i, j) + alpha * (*Y_tmp)(i, j);
    }
    if(deleteYT)
      delete Y_tmp.get();
  }
}

template<class MatrixType, class LocalScalarType>
void
BandedContainer<MatrixType, LocalScalarType, true>::
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

  // Tpetra::MultiVector specialization corresponding to MatrixType.
  Details::MultiVectorLocalGatherScatter<mv_type, local_mv_type> mvgs;
  const size_t numVecs = X.extent(1);

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! IsComputed_, std::runtime_error, "Ifpack2::BandedContainer::apply: "
    "You must have called the compute() method before you may call apply().  "
    "You may call the apply() method as many times as you want after calling "
    "compute() once, but you must have called compute() at least once.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.extent(1) != Y.extent(1), std::runtime_error,
    "Ifpack2::BandedContainer::apply: X and Y have different numbers of "
    "vectors.  X has " << X.extent(1)
    << ", but Y has " << Y.extent(1) << ".");

  if (numVecs == 0) {
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

  mvgs.gatherViewToView(X_local[blockIndex], X, localRows);

  // We must gather the contents of the output multivector Y even on
  // input to applyImpl(), since the inverse operator might use it as
  // an initial guess for a linear solve.  We have no way of knowing
  // whether it does or does not.

  mvgs.gatherViewToView (Y_local[blockIndex], Y, localRows);

  // Apply the local operator:
  // Y_local := beta*Y_local + alpha*M^{-1}*X_local
  this->applyImpl (X_local[blockIndex], Y_local[blockIndex], blockIndex, stride, mode, as<local_scalar_type>(alpha),
                   as<local_scalar_type>(beta));

  // Scatter the permuted subset output vector Y_local back into the
  // original output multivector Y.
  mvgs.scatterViewToView(Y, Y_local[blockIndex], localRows);
}

template<class MatrixType, class LocalScalarType>
void
BandedContainer<MatrixType, LocalScalarType, true>::
weightedApply (HostView& /* X */,
               HostView& /* Y */,
               HostView& /* D */,
               int /* blockIndex */,
               int /* stride */,
               Teuchos::ETransp /* mode */,
               scalar_type /* alpha */,
               scalar_type /* beta */) const
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

  TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::runtime_error, "Ifpack2::BandedContainer::"
      "weightedApply: This code is not tested and not used. Expect bugs.");

  // mfh 09 Aug 2017: Comment out unreachable code to prevent compiler warning.
#if 0

  // The local operator template parameter might have a different
  // Scalar type than MatrixType.  This means that we might have to
  // convert X and Y to the Tpetra::MultiVector specialization that
  // the local operator wants.  This class' X_ and Y_ internal fields
  // are of the right type for the local operator, so we can use those
  // as targets.

  auto zero = Teuchos::ScalarTraits<scalar_type>::zero ();
  auto one = Teuchos::ScalarTraits<scalar_type>::one ();
  // typedef Tpetra::Vector<local_scalar_type, local_ordinal_type, global_ordinal_type, node_type> LV; // unused

  Details::MultiVectorLocalGatherScatter<mv_type, local_mv_type> mvgs;
  const size_t numVecs = X.extent(1);

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! IsComputed_, std::runtime_error, "Ifpack2::BandedContainer::"
    "weightedApply: You must have called the compute() method before you may "
    "call apply().  You may call the apply() method as many times as you want "
    "after calling compute() once, but you must have called compute() at least "
    "once.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    numVecs != Y.extent(1), std::runtime_error,
    "Ifpack2::BandedContainer::weightedApply: X and Y have different numbers "
    "of vectors.  X has " << X.extent(1) << ", but Y has "
    << Y.extent(1) << ".");

  if (numVecs == 0) {
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
  // gives the number of rows in the row Map of the local operator.
  //
  // FIXME (mfh 20 Aug 2013) There might be an implicit assumption
  // here that the row Map and the range Map of that operator are
  // the same.
  //
  // FIXME (mfh 20 Aug 2013) This "local permutation" functionality
  // really belongs in Tpetra.

  const size_t numRows = this->blockRows_[blockIndex];

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

  HostViewLocal D_local("", numRows, 1);
  HostViewLocal X_scaled("", numRows, numVecs);

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

  mvgs.gatherViewToView (D_local, D, localRows);

  for(size_t j = 0; j < numVecs; j++)
    for(size_t i = 0; i < numRows; i++)
      X_scaled(i, j) = X_local[blockIndex](i, j) * D_local(i, 0);

  // Y_temp will hold the result of M^{-1}*X_scaled.  If beta == 0, we
  // can write the result of Inverse_->apply() directly to Y_local, so
  // Y_temp may alias Y_local.  Otherwise, if beta != 0, we need
  // temporary storage for M^{-1}*X_scaled, so Y_temp must be
  // different than Y_local.
  Ptr<HostViewLocal> Y_temp;
  bool deleteYT = false;
  if(beta == zero)
    Y_temp = ptr(&Y_local[blockIndex]);
  else
  {
    Y_temp = ptr(new HostViewLocal("", numRows, numVecs));
    deleteYT = true;
  }

  // Apply the local operator: Y_temp := M^{-1} * X_scaled
  applyImpl(X_scaled, *Y_temp, blockIndex, stride, mode, one, one);
  // Y_local := beta * Y_local + alpha * diag(D_local) * Y_temp.
  //
  // Note that we still use the permuted subset scaling D_local here,
  // because Y_temp has the same permuted subset Map.  That's good, in
  // fact, because it's a subset: less data to read and multiply.

  for(size_t j = 0; j < numVecs; j++)
    for(size_t i = 0; i < numRows; i++)
      Y_local[blockIndex](i, j) = Y_local[blockIndex](i, j) * (local_impl_scalar_type) beta + (local_impl_scalar_type) alpha * (*Y_temp)(i, j) * D_local(i, 0);

  if(deleteYT)
    delete Y_temp.get();

  // Copy the permuted subset output vector Y_local into the original
  // output multivector Y.
  mvgs.scatterViewToView (Y, Y_local[blockIndex], localRows);

#endif // 0
}

template<class MatrixType, class LocalScalarType>
std::ostream&
BandedContainer<MatrixType, LocalScalarType, true>::
print (std::ostream& os) const
{
  Teuchos::FancyOStream fos (Teuchos::rcpFromRef (os));
  fos.setOutputToRootOnly (0);
  describe (fos);
  return os;
}

template<class MatrixType, class LocalScalarType>
std::string
BandedContainer<MatrixType, LocalScalarType, true>::
description () const
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
BandedContainer<MatrixType, LocalScalarType, true>::
describe (Teuchos::FancyOStream& os,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  if(verbLevel==Teuchos::VERB_NONE) return;
  os << "================================================================================" << std::endl;
  os << "Ifpack2::BandedContainer" << std::endl;
  for(int i = 0; i < this->numBlocks_; i++)
  {
    os << "Block " << i << ": Number of rows           = " << this->blockRows_[i] << std::endl;
    os << "Block " << i << ": Number of subdiagonals   = " << diagBlocks_[i].lowerBandwidth() << std::endl;
    os << "Block " << i << ": Number of superdiagonals = " << diagBlocks_[i].upperBandwidth() << std::endl;
  }
  os << "isInitialized()          = " << IsInitialized_ << std::endl;
  os << "isComputed()             = " << IsComputed_ << std::endl;
  os << "================================================================================" << std::endl;
  os << std::endl;
}

template<class MatrixType, class LocalScalarType>
void
BandedContainer<MatrixType, LocalScalarType, true>::
extract ()
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::toString;
  auto& A = *this->inputMatrix_;
  const size_t inputMatrixNumRows = A.getNodeNumRows ();
  // We only use the rank of the calling process and the number of MPI
  // processes for generating error messages.  Extraction itself is
  // entirely local to each participating MPI process.
  const int myRank = A.getRowMap()->getComm()->getRank();
  const int numProcs = A.getRowMap()->getComm()->getSize();

  for(int blockIndex = 0; blockIndex < this->numBlocks_; blockIndex++)
  {
    const local_ordinal_type numRows_ = this->blockRows_[blockIndex];
    // Sanity check that the local row indices to extract fall within
    // the valid range of local row indices for the input matrix.
    ArrayView<const local_ordinal_type> localRows = this->getLocalRows(blockIndex);
    for(local_ordinal_type j = 0; j < numRows_; j++)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(
        localRows[j] < 0 ||
        static_cast<size_t> (localRows[j]) >= inputMatrixNumRows,
        std::runtime_error, "Ifpack2::BandedContainer::extract: On process " <<
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
    const map_type& globalRowMap = *(A.getRowMap ());
    const map_type& globalColMap = *(A.getColMap ());
    const map_type& globalDomMap = *(A.getDomainMap ());

    bool rowIndsValid = true;
    bool colIndsValid = true;
    Array<local_ordinal_type> localCols (numRows_);
    // For error messages, collect the sets of invalid row indices and
    // invalid column indices.  They are otherwise not useful.
    Array<local_ordinal_type> invalidLocalRowInds;
    Array<global_ordinal_type> invalidGlobalColInds;
    for(local_ordinal_type i = 0; i < numRows_; i++)
    {
      // ii_local is the (local) row index we want to look up.
      const local_ordinal_type ii_local = localRows[i];
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
      if(globalDomMap.isNodeGlobalElement(jj_global))
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
        localCols[i] = jj_local;
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! rowIndsValid, std::logic_error, "Ifpack2::BandedContainer::extract: "
      "On process " << myRank << ", at least one row index in the set of local "
      "row indices given to the constructor is not a valid local row index in "
      "the input matrix's row Map on this process.  This should be impossible "
      "because the constructor checks for this case.  Here is the complete set "
      "of invalid local row indices: " << toString(invalidLocalRowInds) << ".  "
      "Please report this bug to the Ifpack2 developers.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! colIndsValid, std::runtime_error, "Ifpack2::BandedContainer::extract: "
      "On process " << myRank << ", "
      "At least one row index in the set of row indices given to the constructor "
      "does not have a corresponding column index in the input matrix's column "
      "Map.  This probably means that the column(s) in question is/are empty on "
      "this process, which would make the submatrix to extract structurally "
      "singular.  Here is the compete set of invalid global column indices: "
      << toString(invalidGlobalColInds) << ".");

    const size_t maxNumEntriesInRow = A.getNodeMaxNumRowEntries();
    Array<scalar_type> val(maxNumEntriesInRow);
    Array<local_ordinal_type> ind(maxNumEntriesInRow);

    const local_ordinal_type INVALID = Teuchos::OrdinalTraits<local_ordinal_type>::invalid();
    for (local_ordinal_type i = 0; i < numRows_; i++)
    {
      const local_ordinal_type localRow = this->partitions_[this->partitionIndices_[blockIndex] + i];
      size_t numEntries;
      A.getLocalRowCopy(localRow, ind(), val(), numEntries);
      for (size_t k = 0; k < numEntries; ++k)
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
        if(localCol >= 0 && static_cast<size_t>(localCol) < inputMatrixNumRows)
        {
          // for local column IDs, look for each ID in the list
          // of columns hosted by this object
          local_ordinal_type jj = INVALID;
          for (size_t kk = 0; kk < (size_t) numRows_; kk++)
          {
            if(localRows[kk] == localCol)
              jj = kk;
          }
          if (jj != INVALID)
            diagBlocks_[blockIndex](i, jj) += val[k]; // ???
        }
      }
    }
  }
}

template<class MatrixType, class LocalScalarType>
std::string BandedContainer<MatrixType, LocalScalarType, true>::getName()
{
  return "Banded";
}

template<class MatrixType, class LocalScalarType>
BandedContainer<MatrixType, LocalScalarType, false>::
BandedContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                 const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
                 const Teuchos::RCP<const import_type>& importer,
                 int OverlapLevel,
                 scalar_type DampingFactor) :
  Container<MatrixType>(matrix, partitions, importer, OverlapLevel, DampingFactor)
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Ifpack2::BandedContainer: Not implemented for "
     "LocalScalarType = " << Teuchos::TypeNameTraits<LocalScalarType>::name ()
     << ".");
}

template<class MatrixType, class LocalScalarType>
BandedContainer<MatrixType, LocalScalarType, false>::
BandedContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                 const Teuchos::Array<local_ordinal_type>& localRows) :
  Container<MatrixType>(matrix, localRows)
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::logic_error, "Ifpack2::BandedContainer: Not implemented for "
     "LocalScalarType = " << Teuchos::TypeNameTraits<LocalScalarType>::name ()
     << ".");
}

template<class MatrixType, class LocalScalarType>
BandedContainer<MatrixType, LocalScalarType, false>::
~BandedContainer () {}

template<class MatrixType, class LocalScalarType>
void BandedContainer<MatrixType, LocalScalarType, false>::
setParameters (const Teuchos::ParameterList& List) {}

template<class MatrixType, class LocalScalarType>
void
BandedContainer<MatrixType, LocalScalarType, false>::
initialize () {}

template<class MatrixType, class LocalScalarType>
void
BandedContainer<MatrixType, LocalScalarType, false>::
compute () {}

template<class MatrixType, class LocalScalarType>
void
BandedContainer<MatrixType, LocalScalarType, false>::
clearBlocks () {}

template<class MatrixType, class LocalScalarType>
void
BandedContainer<MatrixType, LocalScalarType, false>::
factor () {}

template<class MatrixType, class LocalScalarType>
void
BandedContainer<MatrixType, LocalScalarType, false>::
applyImpl (HostViewLocal& X,
           HostViewLocal& Y,
           int blockIndex,
           int stride,
           Teuchos::ETransp mode,
           const local_scalar_type alpha,
           const local_scalar_type beta) const {}

template<class MatrixType, class LocalScalarType>
void
BandedContainer<MatrixType, LocalScalarType, false>::
apply (HostView& X,
       HostView& Y,
       int blockIndex,
       int stride,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const {}

template<class MatrixType, class LocalScalarType>
void
BandedContainer<MatrixType, LocalScalarType, false>::
weightedApply (HostView& X,
               HostView& Y,
               HostView& D,
               int blockIndex,
               int stride,
               Teuchos::ETransp mode,
               scalar_type alpha,
               scalar_type beta) const {}

template<class MatrixType, class LocalScalarType>
std::ostream&
BandedContainer<MatrixType, LocalScalarType, false>::
print (std::ostream& os) const
{
  return os;
}

template<class MatrixType, class LocalScalarType>
std::string
BandedContainer<MatrixType, LocalScalarType, false>::
description () const
{
  return "";
}

template<class MatrixType, class LocalScalarType>
void
BandedContainer<MatrixType, LocalScalarType, false>::
describe (Teuchos::FancyOStream& os,
          const Teuchos::EVerbosityLevel verbLevel) const {}

template<class MatrixType, class LocalScalarType>
void
BandedContainer<MatrixType, LocalScalarType, false>::
extract () {}

template<class MatrixType, class LocalScalarType>
std::string BandedContainer<MatrixType, LocalScalarType, false>::getName()
{
  return "";
}

} // namespace Ifpack2

#define IFPACK2_BANDEDCONTAINER_INSTANT(S,LO,GO,N) \
  template class Ifpack2::BandedContainer< Tpetra::RowMatrix<S, LO, GO, N>, S >;

#endif // IFPACK2_BANDEDCONTAINER_HPP
