/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
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

#ifndef IFPACK2_DENSECONTAINER_DEF_HPP
#define IFPACK2_DENSECONTAINER_DEF_HPP

#include "Ifpack2_DenseContainer_decl.hpp"
#include "Teuchos_LAPACK.hpp"

#ifdef HAVE_MPI
#  include <mpi.h>
#  include "Teuchos_DefaultMpiComm.hpp"
#else
#  include "Teuchos_DefaultSerialComm.hpp"
#endif // HAVE_MPI


namespace Ifpack2 {

//==============================================================================
template<class MatrixType, class LocalScalarType>
DenseContainer<MatrixType, LocalScalarType>::
DenseContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                const Teuchos::ArrayView<const local_ordinal_type>& localRows) :
  Container<MatrixType> (matrix, localRows),
  numRows_ (localRows.size ()),
  diagBlock_ (numRows_, numRows_),
  ipiv_ (numRows_, 0)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::toString;
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;
  typedef typename ArrayView<const local_ordinal_type>::size_type size_type;
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! matrix->hasColMap (), std::invalid_argument, "Ifpack2::DenseContainer: "
    "The constructor's input matrix must have a column Map.");

  // Check whether the input set of local row indices is correct.
  const map_type& rowMap = * (matrix->getRowMap ());
  const size_type numRows = localRows.size ();
  bool rowIndicesValid = true;
  Array<local_ordinal_type> invalidLocalRowIndices;
  for (size_type i = 0; i < numRows; ++i) {
    if (! rowMap.isNodeLocalElement (localRows[i])) {
      rowIndicesValid = false;
      invalidLocalRowIndices.push_back (localRows[i]);
      break;
    }
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! rowIndicesValid, std::invalid_argument, "Ifpack2::DenseContainer: "
    "On process " << rowMap.getComm ()->getRank () << " of "
    << rowMap.getComm ()->getSize () << ", in the given set of local row "
    "indices localRows = " << toString (localRows) << ", the following "
    "entries are not valid local row indices on the calling process: "
    << toString (invalidLocalRowIndices) << ".");

#ifdef HAVE_MPI
  RCP<const Teuchos::Comm<int> > localComm =
    rcp (new Teuchos::MpiComm<int> (MPI_COMM_SELF));
#else
  RCP<const Teuchos::Comm<int> > localComm =
    rcp (new Teuchos::SerialComm<int> ());
#endif // HAVE_MPI

  // FIXME (mfh 25 Aug 2013) What if the matrix's row Map has a
  // different index base than zero?
  const global_ordinal_type indexBase = 0;
  localMap_ = rcp (new map_type (numRows_, indexBase, localComm));
}

//==============================================================================
template<class MatrixType, class LocalScalarType>
DenseContainer<MatrixType, LocalScalarType>::~DenseContainer()
{}

//==============================================================================
template<class MatrixType, class LocalScalarType>
size_t DenseContainer<MatrixType, LocalScalarType>::getNumRows () const
{
  return numRows_;
}

//==============================================================================
template<class MatrixType, class LocalScalarType>
bool DenseContainer<MatrixType, LocalScalarType>::isInitialized () const
{
  return IsInitialized_;
}

//==============================================================================
template<class MatrixType, class LocalScalarType>
bool DenseContainer<MatrixType, LocalScalarType>::isComputed() const
{
  return IsComputed_;
}

//==============================================================================
template<class MatrixType, class LocalScalarType>
void DenseContainer<MatrixType, LocalScalarType>::
setParameters (const Teuchos::ParameterList& List)
{
  (void) List; // the solver doesn't currently take any parameters
}

//==============================================================================
template<class MatrixType, class LocalScalarType>
void DenseContainer<MatrixType, LocalScalarType>::initialize ()
{
  using Teuchos::null;
  using Teuchos::rcp;

  // We assume that if you called this method, you intend to recompute
  // everything.
  IsInitialized_ = false;
  IsComputed_ = false;

  // Fill the diagonal block and LU permutation array with zeros.
  diagBlock_.putScalar (Teuchos::ScalarTraits<local_scalar_type>::zero ());
  std::fill (ipiv_.begin (), ipiv_.end (), 0);

  IsInitialized_ = true;
}

//==============================================================================
template<class MatrixType, class LocalScalarType>
void DenseContainer<MatrixType, LocalScalarType>::compute ()
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    static_cast<size_t> (ipiv_.size ()) != numRows_, std::logic_error,
    "Ifpack2::DenseContainer::compute: ipiv_ array has the wrong size.  "
    "Please report this bug to the Ifpack2 developers.");

  IsComputed_ = false;
  if (! this->isInitialized ()) {
    this->initialize ();
  }

  // Extract the submatrix.
  extract (this->getMatrix ()); // extract the submatrix
  factor (); // factor the submatrix

  IsComputed_ = true;
}

template<class MatrixType, class LocalScalarType>
void DenseContainer<MatrixType, LocalScalarType>::factor ()
{
  Teuchos::LAPACK<int, local_scalar_type> lapack;
  int INFO = 0;
  lapack.GETRF (diagBlock_.numRows (), diagBlock_.numCols (),
                diagBlock_.values (), diagBlock_.stride (),
                ipiv_.getRawPtr (), &INFO);
  // INFO < 0 is a bug.
  TEUCHOS_TEST_FOR_EXCEPTION(
    INFO < 0, std::logic_error, "Ifpack2::DenseContainer::factor: "
    "LAPACK's _GETRF (LU factorization with partial pivoting) was called "
    "incorrectly.  INFO = " << INFO << " < 0.  "
    "Please report this bug to the Ifpack2 developers.");
  // INFO > 0 means the matrix is singular.  This is probably an issue
  // either with the choice of rows the rows we extracted, or with the
  // input matrix itself.
  TEUCHOS_TEST_FOR_EXCEPTION(
    INFO > 0, std::runtime_error, "Ifpack2::DenseContainer::factor: "
    "LAPACK's _GETRF (LU factorization with partial pivoting) reports that the "
    "computed U factor is exactly singular.  U(" << INFO << "," << INFO << ") "
    "(one-based index i) is exactly zero.  This probably means that the input "
    "matrix has a singular diagonal block.");
}

//==============================================================================
template<class MatrixType, class LocalScalarType>
void DenseContainer<MatrixType, LocalScalarType>::
applyImpl (const local_mv_type& X,
           local_mv_type& Y,
           Teuchos::ETransp mode,
           LocalScalarType alpha,
           LocalScalarType beta) const
{
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;

  TEUCHOS_TEST_FOR_EXCEPTION(
    X.getLocalLength () != Y.getLocalLength (),
    std::logic_error, "Ifpack2::DenseContainer::applyImpl: X and Y have "
    "incompatible dimensions (" << X.getLocalLength () << " resp. "
    << Y.getLocalLength () << ").  Please report this bug to "
    "the Ifpack2 developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    localMap_->getNodeNumElements () != X.getLocalLength (),
    std::logic_error, "Ifpack2::DenseContainer::applyImpl: The inverse "
    "operator and X have incompatible dimensions (" <<
    localMap_->getNodeNumElements () << " resp. "
    << X.getLocalLength () << ").  Please report this bug to "
    "the Ifpack2 developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    localMap_->getNodeNumElements () != Y.getLocalLength (),
    std::logic_error, "Ifpack2::DenseContainer::applyImpl: The inverse "
    "operator and Y have incompatible dimensions (" <<
    localMap_->getNodeNumElements () << " resp. "
    << Y.getLocalLength () << ").  Please report this bug to "
    "the Ifpack2 developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.getLocalLength () != static_cast<size_t> (diagBlock_.numRows ()),
    std::logic_error, "Ifpack2::DenseContainer::applyImpl: The input "
    "multivector X has incompatible dimensions from those of the "
    "inverse operator (" << X.getLocalLength () << " vs. "
    << (mode == Teuchos::NO_TRANS ? diagBlock_.numCols () : diagBlock_.numRows ())
    << ").  Please report this bug to the Ifpack2 developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.getLocalLength () != static_cast<size_t> (diagBlock_.numRows ()),
    std::logic_error, "Ifpack2::DenseContainer::applyImpl: The output "
    "multivector Y has incompatible dimensions from those of the "
    "inverse operator (" << Y.getLocalLength () << " vs. "
    << (mode == Teuchos::NO_TRANS ? diagBlock_.numRows () : diagBlock_.numCols ())
    << ").  Please report this bug to the Ifpack2 developers.");

  typedef Teuchos::ScalarTraits<local_scalar_type> STS;
  const int numVecs = static_cast<int> (X.getNumVectors ());
  if (alpha == STS::zero ()) { // don't need to solve the linear system
    if (beta == STS::zero ()) {
      // Use BLAS AXPY semantics for beta == 0: overwrite, clobbering
      // any Inf or NaN values in Y (rather than multiplying them by
      // zero, resulting in NaN values).
      Y.putScalar (STS::zero ());
    }
    else { // beta != 0
      Y.scale (beta);
    }
  }
  else { // alpha != 0; must solve the linear system
    Teuchos::LAPACK<int, local_scalar_type> lapack;
    // If beta is nonzero or Y is not constant stride, we have to use
    // a temporary output multivector.  It gets a (deep) copy of X,
    // since GETRS overwrites its (multi)vector input with its output.
    RCP<local_mv_type> Y_tmp;
    if (beta == STS::zero () ){
      Tpetra::deep_copy (Y, X);
      Y_tmp = rcpFromRef (Y);
    }
    else {
      Y_tmp = rcp (new local_mv_type (X, Teuchos::Copy));
    }
    const int Y_stride = static_cast<int> (Y_tmp->getStride ());
    ArrayRCP<local_scalar_type> Y_view = Y_tmp->get1dViewNonConst ();
    local_scalar_type* const Y_ptr = Y_view.getRawPtr ();

    int INFO = 0;
    const char trans =
      (mode == Teuchos::CONJ_TRANS ? 'C' : (mode == Teuchos::TRANS ? 'T' : 'N'));
    lapack.GETRS (trans, diagBlock_.numRows (), numVecs,
                  diagBlock_.values (), diagBlock_.stride (),
                  ipiv_.getRawPtr (), Y_ptr, Y_stride, &INFO);
    TEUCHOS_TEST_FOR_EXCEPTION(
      INFO != 0, std::runtime_error, "Ifpack2::DenseContainer::applyImpl: "
      "LAPACK's _GETRS (solve using LU factorization with partial pivoting) "
      "failed with INFO = " << INFO << " != 0.");

    if (beta != STS::zero ()) {
      Y.update (alpha, *Y_tmp, beta);
    }
    else if (! Y.isConstantStride ()) {
      Tpetra::deep_copy (Y, *Y_tmp);
    }
  }
}

//==============================================================================
template<class MatrixType, class LocalScalarType>
void DenseContainer<MatrixType, LocalScalarType>::
apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
       Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
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

  // Tpetra::MultiVector specialization corresponding to the global operator.
  typedef Tpetra::MultiVector<scalar_type,
    local_ordinal_type, global_ordinal_type, node_type> global_mv_type;

  const char prefix[] = "Ifpack2::DenseContainer::weightedApply: ";
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! IsComputed_, std::runtime_error, prefix << "You must have called the "
    "compute() method before you may call this method.  You may call "
    "apply() as many times as you want after calling compute() once, "
    "but you must have called compute() at least once first.");
  const size_t numVecs = X.getNumVectors ();
  TEUCHOS_TEST_FOR_EXCEPTION(
    numVecs != Y.getNumVectors (), std::runtime_error,
    prefix << "X and Y have different numbers of vectors (columns).  X has "
    << X.getNumVectors () << ", but Y has " << X.getNumVectors () << ".");

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
  // gives the number of rows in the row Map of the local Inverse_
  // operator.
  //
  // FIXME (mfh 20 Aug 2013) There might be an implicit assumption
  // here that the row Map and the range Map of that operator are
  // the same.
  //
  // FIXME (mfh 20 Aug 2013) This "local permutation" functionality
  // really belongs in Tpetra.

  if (X_.is_null ()) {
    X_ = rcp (new local_mv_type (localMap_, numVecs));
  }
  RCP<local_mv_type> X_local = X_;
  TEUCHOS_TEST_FOR_EXCEPTION(
    X_local->getLocalLength () != numRows_, std::logic_error,
    "Ifpack2::DenseContainer::apply: "
    "X_local has length " << X_local->getLocalLength () << ", which does "
    "not match numRows_ = " << numRows_ << ".  Please report this bug to "
    "the Ifpack2 developers.");
  ArrayView<const local_ordinal_type> localRows = this->getLocalRows ();

  Details::MultiVectorLocalGatherScatter<global_mv_type, local_mv_type> mvgs;
  mvgs.gather (*X_local, X, localRows);

  // We must gather the contents of the output multivector Y even on
  // input to applyImpl(), since the inverse operator might use it as
  // an initial guess for a linear solve.  We have no way of knowing
  // whether it does or does not.

  if (Y_.is_null ()) {
    Y_ = rcp (new local_mv_type (localMap_, numVecs));
  }
  RCP<local_mv_type> Y_local = Y_;
  TEUCHOS_TEST_FOR_EXCEPTION(
    Y_local->getLocalLength () != numRows_, std::logic_error,
    "Ifpack2::DenseContainer::apply: "
    "Y_local has length " << X_local->getLocalLength () << ", which does "
    "not match numRows_ = " << numRows_ << ".  Please report this bug to "
    "the Ifpack2 developers.");
  mvgs.gather (*Y_local, Y, localRows);

  // Apply the local operator:
  // Y_local := beta*Y_local + alpha*M^{-1}*X_local
  this->applyImpl (*X_local, *Y_local, mode, as<local_scalar_type> (alpha),
                   as<local_scalar_type> (beta));

  // Scatter the permuted subset output vector Y_local back into the
  // original output multivector Y.
  mvgs.scatter (Y, *Y_local, localRows);
}

//==============================================================================
template<class MatrixType, class LocalScalarType>
void DenseContainer<MatrixType,LocalScalarType>::
weightedApply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
               Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
               const Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& D,
               Teuchos::ETransp mode,
               scalar_type alpha,
               scalar_type beta) const
{
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::Range1D;
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

  // Tpetra::MultiVector specialization corresponding to the global operator.
  typedef Tpetra::MultiVector<scalar_type,
    local_ordinal_type, global_ordinal_type, node_type> global_mv_type;
  // Tpetra::Vector specialization corresponding to the local
  // operator.  We will use this for the subset permutation of the
  // diagonal scaling D.
  typedef Tpetra::Vector<local_scalar_type, local_ordinal_type,
                         global_ordinal_type, node_type> local_vec_type;

  const char prefix[] = "Ifpack2::DenseContainer::weightedApply: ";
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! IsComputed_, std::runtime_error, prefix << "You must have called the "
    "compute() method before you may call this method.  You may call "
    "weightedApply() as many times as you want after calling compute() once, "
    "but you must have called compute() at least once first.");
  const size_t numVecs = X.getNumVectors ();
  TEUCHOS_TEST_FOR_EXCEPTION(
    numVecs != Y.getNumVectors (), std::runtime_error,
    prefix << "X and Y have different numbers of vectors (columns).  X has "
    << X.getNumVectors () << ", but Y has " << X.getNumVectors () << ".");

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

  if (X_.is_null ()) {
    X_ = rcp (new local_mv_type (localMap_, numVecs));
  }
  RCP<local_mv_type> X_local = X_;
  TEUCHOS_TEST_FOR_EXCEPTION(
    X_local->getLocalLength () != numRows_, std::logic_error,
    "Ifpack2::DenseContainer::weightedApply: "
    "X_local has length " << X_local->getLocalLength () << ", which does "
    "not match numRows_ = " << numRows_ << ".  Please report this bug to "
    "the Ifpack2 developers.");
  ArrayView<const local_ordinal_type> localRows = this->getLocalRows ();

  Details::MultiVectorLocalGatherScatter<global_mv_type, local_mv_type> mvgs;
  mvgs.gather (*X_local, X, localRows);

  // We must gather the output multivector Y even on input to
  // applyImpl(), since the local operator might use it as an initial
  // guess for a linear solve.  We have no way of knowing whether it
  // does or does not.

  if (Y_.is_null ()) {
    Y_ = rcp (new local_mv_type (localMap_, numVecs));
  }
  RCP<local_mv_type> Y_local = Y_;
  TEUCHOS_TEST_FOR_EXCEPTION(
    Y_local->getLocalLength () != numRows_, std::logic_error,
    "Ifpack2::DenseContainer::weightedApply: "
    "Y_local has length " << X_local->getLocalLength () << ", which does "
    "not match numRows_ = " << numRows_ << ".  Please report this bug to "
    "the Ifpack2 developers.");
  mvgs.gather (*Y_local, Y, localRows);

  // Apply the diagonal scaling D to the input X.  It's our choice
  // whether the result has the original input Map of X, or the
  // permuted subset Map of X_local.  If the latter, we also need to
  // gather D into the permuted subset Map.  We choose the latter, to
  // save memory and computation.  Thus, we do the following:
  //
  // 1. Gather D into a temporary vector D_local.
  // 2. Create a temporary X_scaled to hold diag(D_local) * X_local.
  // 3. Compute X_scaled := diag(D_loca) * X_local.

  local_vec_type D_local (localMap_);
  TEUCHOS_TEST_FOR_EXCEPTION(
    D_local.getLocalLength () != numRows_, std::logic_error,
    "Ifpack2::DenseContainer::weightedApply: "
    "D_local has length " << D_local.getLocalLength () << ", which does "
    "not match numRows_ = " << numRows_ << ".  Please report this bug to "
    "the Ifpack2 developers.");
  mvgs.gather (D_local, D, localRows);
  local_mv_type X_scaled (localMap_, numVecs);
  X_scaled.elementWiseMultiply (STS::one (), D_local, *X_local, STS::zero ());

  // Y_temp will hold the result of M^{-1}*X_scaled.  If beta == 0, we
  // can write the result of Inverse_->apply() directly to Y_local, so
  // Y_temp may alias Y_local.  Otherwise, if beta != 0, we need
  // temporary storage for M^{-1}*X_scaled, so Y_temp must be
  // different than Y_local.
  RCP<local_mv_type> Y_temp;
  if (beta == STS::zero ()) {
    Y_temp = Y_local;
  } else {
    Y_temp = rcp (new local_mv_type (localMap_, numVecs));
  }

  // Apply the local operator: Y_temp := M^{-1} * X_scaled
  this->applyImpl (X_scaled, *Y_temp, mode, STS::one (), STS::zero ());
  // Y_local := beta * Y_local + alpha * diag(D_local) * Y_temp.
  //
  // Note that we still use the permuted subset scaling D_local here,
  // because Y_temp has the same permuted subset Map.  That's good, in
  // fact, because it's a subset: less data to read and multiply.
  Y_local->elementWiseMultiply (alpha, D_local, *Y_temp, beta);

  // Copy the permuted subset output vector Y_local into the original
  // output multivector Y.
  mvgs.scatter (Y, *Y_local, localRows);
}

//==============================================================================
template<class MatrixType, class LocalScalarType>
std::ostream& DenseContainer<MatrixType,LocalScalarType>::
print (std::ostream& os) const
{
  Teuchos::FancyOStream fos (Teuchos::rcpFromRef (os));
  fos.setOutputToRootOnly (0);
  this->describe (fos);
  return os;
}

//==============================================================================
template<class MatrixType, class LocalScalarType>
std::string DenseContainer<MatrixType,LocalScalarType>::description() const
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

//==============================================================================
template<class MatrixType, class LocalScalarType>
void DenseContainer<MatrixType,LocalScalarType>::describe(Teuchos::FancyOStream &os, const Teuchos::EVerbosityLevel verbLevel) const
{
  using std::endl;
  if(verbLevel==Teuchos::VERB_NONE) return;
  os << "================================================================================" << endl;
  os << "Ifpack2::DenseContainer" << endl;
  os << "Number of rows          = " << numRows_ << endl;
  os << "isInitialized()         = " << IsInitialized_ << endl;
  os << "isComputed()            = " << IsComputed_ << endl;
  os << "================================================================================" << endl;
  os << endl;
}

//==============================================================================
template<class MatrixType, class LocalScalarType>
void DenseContainer<MatrixType,LocalScalarType>::
extract (const Teuchos::RCP<const row_matrix_type>& globalMatrix)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::toString;
  typedef local_ordinal_type LO;
  typedef global_ordinal_type GO;
  typedef Tpetra::Map<LO, GO, node_type> map_type;
  const size_t inputMatrixNumRows = globalMatrix->getNodeNumRows ();
  // We only use the rank of the calling process and the number of MPI
  // processes for generating error messages.  Extraction itself is
  // entirely local to each participating MPI process.
  const int myRank = globalMatrix->getRowMap ()->getComm ()->getRank ();
  const int numProcs = globalMatrix->getRowMap ()->getComm ()->getSize ();

  // Sanity check that the local row indices to extract fall within
  // the valid range of local row indices for the input matrix.
  ArrayView<const LO> localRows = this->getLocalRows ();
  for (size_t j = 0; j < numRows_; ++j) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      localRows[j] < 0 ||
      static_cast<size_t> (localRows[j]) >= inputMatrixNumRows,
      std::runtime_error, "Ifpack2::DenseContainer::extract: On process " <<
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
  const map_type& globalRowMap = * (globalMatrix->getRowMap ());
  const map_type& globalColMap = * (globalMatrix->getColMap ());
  const map_type& globalDomMap = * (globalMatrix->getDomainMap ());

  bool rowIndsValid = true;
  bool colIndsValid = true;
  Array<LO> localCols (numRows_);
  // For error messages, collect the sets of invalid row indices and
  // invalid column indices.  They are otherwise not useful.
  Array<LO> invalidLocalRowInds;
  Array<GO> invalidGlobalColInds;
  for (size_t i = 0; i < numRows_; ++i) {
    // ii_local is the (local) row index we want to look up.
    const LO ii_local = localRows[i];
    // Find the global index jj_global corresponding to ii_local.
    // Global indices are the same (rather, are required to be the
    // same) in all three Maps, which is why we use jj (suggesting a
    // column index, which is how we will use it below).
    const GO jj_global = globalRowMap.getGlobalElement (ii_local);
    if (jj_global == Teuchos::OrdinalTraits<GO>::invalid ()) {
      // If ii_local is not a local index in the row Map on the
      // calling process, that means localRows is incorrect.  We've
      // already checked for this in the constructor, but we might as
      // well check again here, since it's cheap to do so (just an
      // integer comparison, since we need jj_global anyway).
      rowIndsValid = false;
      invalidLocalRowInds.push_back (ii_local);
      break;
    }
    // Exclude "off-process" entries: that is, those in the column Map
    // on this process that are not in the domain Map on this process.
    if (globalDomMap.isNodeGlobalElement (jj_global)) {
      // jj_global is not an off-process entry.  Look up its local
      // index in the column Map; we want to extract this column index
      // from the input matrix.  If jj_global is _not_ in the column
      // Map on the calling process, that could mean that the column
      // in question is empty on this process.  That would be bad for
      // solving linear systems with the extract submatrix.  We could
      // solve the resulting singular linear systems in a minimum-norm
      // least-squares sense, but for now we simply raise an exception.
      const LO jj_local = globalColMap.getLocalElement (jj_global);
      if (jj_local == Teuchos::OrdinalTraits<local_ordinal_type>::invalid ()) {
        colIndsValid = false;
        invalidGlobalColInds.push_back (jj_global);
        break;
      }
      localCols[i] = jj_local;
    }
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! rowIndsValid, std::logic_error, "Ifpack2::DenseContainer::extract: "
    "On process " << myRank << ", at least one row index in the set of local "
    "row indices given to the constructor is not a valid local row index in "
    "the input matrix's row Map on this process.  This should be impossible "
    "because the constructor checks for this case.  Here is the complete set "
    "of invalid local row indices: " << toString (invalidLocalRowInds) << ".  "
    "Please report this bug to the Ifpack2 developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! colIndsValid, std::runtime_error, "Ifpack2::DenseContainer::extract: "
    "On process " << myRank << ", "
    "At least one row index in the set of row indices given to the constructor "
    "does not have a corresponding column index in the input matrix's column "
    "Map.  This probably means that the column(s) in question is/are empty on "
    "this process, which would make the submatrix to extract structurally "
    "singular.  Here is the compete set of invalid global column indices: "
    << toString (invalidGlobalColInds) << ".");

  diagBlock_.putScalar (Teuchos::ScalarTraits<local_scalar_type>::zero ());

  const size_t maxNumEntriesInRow = globalMatrix->getNodeMaxNumRowEntries ();
  Array<scalar_type> val (maxNumEntriesInRow);
  Array<local_ordinal_type> ind (maxNumEntriesInRow);

  const local_ordinal_type INVALID =
    Teuchos::OrdinalTraits<local_ordinal_type>::invalid ();
  for (size_t i = 0; i < numRows_; ++i) {
    const local_ordinal_type localRow = localRows[i];
    size_t numEntries;
    globalMatrix->getLocalRowCopy (localRow, ind (), val (), numEntries);

    for (size_t k = 0; k < numEntries; ++k) {
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
      if (localCol >= 0 && static_cast<size_t> (localCol) < inputMatrixNumRows) {
        // for local column IDs, look for each ID in the list
        // of columns hosted by this object
        local_ordinal_type jj = INVALID;
        for (size_t kk = 0; kk < numRows_; ++kk) {
          if (localRows[kk] == localCol) {
            jj = kk;
          }
        }
        if (jj != INVALID) {
          diagBlock_ (i, jj) += val[k]; // ???
        }
      }
    }
  }
}

} // namespace Ifpack2

// FIXME (mfh 16 Sep 2014) We should really only use RowMatrix here!
// There's no need to instantiate for CrsMatrix too.  All Ifpack2
// preconditioners can and should do dynamic casts if they need a type
// more specific than RowMatrix.

#define IFPACK2_DENSECONTAINER_INSTANT(S,LO,GO,N) \
  template class Ifpack2::DenseContainer< Tpetra::RowMatrix<S, LO, GO, N>, S >; \
  template class Ifpack2::DenseContainer< Tpetra::CrsMatrix<S, LO, GO, N>, S >;

#endif // IFPACK2_SPARSECONTAINER_HPP
