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

#ifndef IFPACK2_EXPERIMENTAL_CRSRBILUK_DEF_HPP
#define IFPACK2_EXPERIMENTAL_CRSRBILUK_DEF_HPP

#include <Tpetra_Experimental_BlockMultiVector.hpp>
#include <Ifpack2_Experimental_RBILUK.hpp>
#include <Ifpack2_RILUK.hpp>

namespace Ifpack2 {

namespace Experimental {

template<class MatrixType>
RBILUK<MatrixType>::RBILUK (const Teuchos::RCP<const block_crs_matrix_type>& Matrix_in)
  : RILUK<row_matrix_type>(Teuchos::rcp_dynamic_cast<const row_matrix_type>(Matrix_in) ),
    A_block_(Matrix_in)
{}


template<class MatrixType>
RBILUK<MatrixType>::~RBILUK() {}


template<class MatrixType>
void
RBILUK<MatrixType>::setMatrix (const Teuchos::RCP<const block_crs_matrix_type>& A)
{
  // It's legal for A to be null; in that case, you may not call
  // initialize() until calling setMatrix() with a nonnull input.
  // Regardless, setting the matrix invalidates any previous
  // factorization.
  if (A.getRawPtr () != A_block_.getRawPtr ())
  {
    this->isAllocated_ = false;
    this->isInitialized_ = false;
    this->isComputed_ = false;
    A_local_block_crs_ = Teuchos::null;
    this->Graph_ = Teuchos::null;
    L_block_ = Teuchos::null;
    U_block_ = Teuchos::null;
    D_block_ = Teuchos::null;
    A_block_ = A;
  }
}



template<class MatrixType>
const typename RBILUK<MatrixType>::block_crs_matrix_type&
RBILUK<MatrixType>::getLBlock () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    L_block_.is_null (), std::runtime_error, "Ifpack2::RILUK::getL: The L factor "
    "is null.  Please call initialize() (and preferably also compute()) "
    "before calling this method.  If the input matrix has not yet been set, "
    "you must first call setMatrix() with a nonnull input matrix before you "
    "may call initialize() or compute().");
  return *L_block_;
}


template<class MatrixType>
const typename RBILUK<MatrixType>::block_crs_matrix_type&
RBILUK<MatrixType>::getDBlock () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    D_block_.is_null (), std::runtime_error, "Ifpack2::RILUK::getD: The D factor "
    "(of diagonal entries) is null.  Please call initialize() (and "
    "preferably also compute()) before calling this method.  If the input "
    "matrix has not yet been set, you must first call setMatrix() with a "
    "nonnull input matrix before you may call initialize() or compute().");
  return *D_block_;
}


template<class MatrixType>
const typename RBILUK<MatrixType>::block_crs_matrix_type&
RBILUK<MatrixType>::getUBlock () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    U_block_.is_null (), std::runtime_error, "Ifpack2::RILUK::getU: The U factor "
    "is null.  Please call initialize() (and preferably also compute()) "
    "before calling this method.  If the input matrix has not yet been set, "
    "you must first call setMatrix() with a nonnull input matrix before you "
    "may call initialize() or compute().");
  return *U_block_;
}

template<class MatrixType>
void RBILUK<MatrixType>::allocate_L_and_U_blocks ()
{
  using Teuchos::null;
  using Teuchos::rcp;

  if (! this->isAllocated_) {
    // Deallocate any existing storage.  This avoids storing 2x
    // memory, since RCP op= does not deallocate until after the
    // assignment.
    this->L_ = null;
    this->U_ = null;
    this->D_ = null;
    L_block_ = null;
    U_block_ = null;
    D_block_ = null;

    // Allocate Matrix using ILUK graphs
    const local_ordinal_type blockSize = A_block_->getBlockSize ();
    L_block_ = rcp(new block_crs_matrix_type(*this->Graph_->getL_Graph (), blockSize) );
    U_block_ = rcp(new block_crs_matrix_type(*this->Graph_->getU_Graph (), blockSize) );
    if (!A_block_->isComputedDiagonalGraph() ) Teuchos::rcp_const_cast<block_crs_matrix_type>(A_block_)->computeDiagonalGraph();
    D_block_ = rcp(new block_crs_matrix_type(*(A_block_->getDiagonalGraph ()), blockSize) );
    L_block_->setAllToScalar (STM::zero ()); // Zero out L and U matrices
    U_block_->setAllToScalar (STM::zero ());
    D_block_->setAllToScalar (STM::zero ());

  }
  this->isAllocated_ = true;
}

template<class MatrixType>
Teuchos::RCP<const typename RBILUK<MatrixType>::block_crs_matrix_type>
RBILUK<MatrixType>::getBlockMatrix () const {
  return A_block_;
}

template<class MatrixType>
void RBILUK<MatrixType>::initialize ()
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcp_implicit_cast;

  TEUCHOS_TEST_FOR_EXCEPTION(
    A_block_.is_null (), std::runtime_error, "Ifpack2::Experimental::RBILUK::initialize: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");

  Teuchos::Time timer ("RBILUK::initialize");
  { // Start timing
    Teuchos::TimeMonitor timeMon (timer);

    // Calling initialize() means that the user asserts that the graph
    // of the sparse matrix may have changed.  We must not just reuse
    // the previous graph in that case.
    //
    // Regarding setting isAllocated_ to false: Eventually, we may want
    // some kind of clever memory reuse strategy, but it's always
    // correct just to blow everything away and start over.
    this->isInitialized_ = false;
    this->isAllocated_ = false;
    this->isComputed_ = false;
    this->Graph_ = Teuchos::null;

    typedef Tpetra::CrsGraph<local_ordinal_type,
                             global_ordinal_type,
                             node_type> crs_graph_type;

    Teuchos::RCP<const crs_graph_type> matrixCrsGraph = Teuchos::rcpFromRef(A_block_->getCrsGraph() );
    this->Graph_ = rcp (new Ifpack2::IlukGraph<crs_graph_type> (matrixCrsGraph,
        this->LevelOfFill_, 0));

    this->Graph_->initialize ();
    allocate_L_and_U_blocks ();
    initAllValues (*A_block_);
  } // Stop timing

  this->isInitialized_ = true;
  this->numInitialize_ += 1;
  this->initializeTime_ += timer.totalElapsedTime ();
}


template<class MatrixType>
void
RBILUK<MatrixType>::
initAllValues (const block_crs_matrix_type& A)
{
  using Teuchos::ArrayRCP;
  using Teuchos::Comm;
  using Teuchos::ptr;
  using Teuchos::RCP;
  using Teuchos::REDUCE_SUM;
  using Teuchos::reduceAll;
  typedef Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> map_type;

  local_ordinal_type NumIn = 0, NumL = 0, NumU = 0;
  bool DiagFound = false;
  size_t NumNonzeroDiags = 0;
  size_t MaxNumEntries = A.getNodeMaxNumRowEntries();
  local_ordinal_type blockSize = A.getBlockSize();
  local_ordinal_type blockMatSize = blockSize*blockSize;

  // First check that the local row map ordering is the same as the local portion of the column map.
  // The extraction of the strictly lower/upper parts of A, as well as the factorization,
  // implicitly assume that this is the case.
  Teuchos::ArrayView<const global_ordinal_type> rowGIDs = A.getRowMap()->getNodeElementList();
  Teuchos::ArrayView<const global_ordinal_type> colGIDs = A.getColMap()->getNodeElementList();
  bool gidsAreConsistentlyOrdered=true;
  global_ordinal_type indexOfInconsistentGID=0;
  for (global_ordinal_type i=0; i<rowGIDs.size(); ++i) {
    if (rowGIDs[i] != colGIDs[i]) {
      gidsAreConsistentlyOrdered=false;
      indexOfInconsistentGID=i;
      break;
    }
  }
  TEUCHOS_TEST_FOR_EXCEPTION(gidsAreConsistentlyOrdered==false, std::runtime_error,
                             "The ordering of the local GIDs in the row and column maps is not the same"
                             << std::endl << "at index " << indexOfInconsistentGID
                             << ".  Consistency is required, as all calculations are done with"
                             << std::endl << "local indexing.");

  // Allocate temporary space for extracting the strictly
  // lower and upper parts of the matrix A.
  Teuchos::Array<local_ordinal_type> LI(MaxNumEntries);
  Teuchos::Array<local_ordinal_type> UI(MaxNumEntries);
  Teuchos::Array<scalar_type> LV(MaxNumEntries*blockMatSize);
  Teuchos::Array<scalar_type> UV(MaxNumEntries*blockMatSize);

  Teuchos::Array<scalar_type> diagValues(blockMatSize);

  L_block_->setAllToScalar (STM::zero ()); // Zero out L and U matrices
  U_block_->setAllToScalar (STM::zero ());
  D_block_->setAllToScalar (STM::zero ()); // Set diagonal values to zero

  RCP<const map_type> rowMap = L_block_->getRowMap ();

  // First we copy the user's matrix into L and U, regardless of fill level.
  // It is important to note that L and U are populated using local indices.
  // This means that if the row map GIDs are not monotonically increasing
  // (i.e., permuted or gappy), then the strictly lower (upper) part of the
  // matrix is not the one that you would get if you based L (U) on GIDs.
  // This is ok, as the *order* of the GIDs in the rowmap is a better
  // expression of the user's intent than the GIDs themselves.

  Teuchos::ArrayView<const global_ordinal_type> nodeGIDs = rowMap->getNodeElementList();
  for (size_t myRow=0; myRow<A.getNodeNumRows(); ++myRow) {
    local_ordinal_type local_row = myRow;

    //TODO JJH 4April2014 An optimization is to use getLocalRowView.  Not all matrices support this,
    //                    we'd need to check via the Tpetra::RowMatrix method supportsRowViews().
    const local_ordinal_type * InI = 0;
    scalar_type * InV = 0;
    A.getLocalRowView(local_row, InI, InV, NumIn);

    // Split into L and U (we don't assume that indices are ordered).

    NumL = 0;
    NumU = 0;
    DiagFound = false;

    for (local_ordinal_type j = 0; j < NumIn; ++j) {
      const local_ordinal_type k = InI[j];
      const local_ordinal_type blockOffset = blockMatSize*j;

      if (k == local_row) {
        DiagFound = true;
        // Store perturbed diagonal in Tpetra::Vector D_
        for (local_ordinal_type jj = 0; jj < blockMatSize; ++jj)
          diagValues[jj] = this->Rthresh_ * InV[blockOffset+jj] + IFPACK2_SGN(InV[blockOffset+jj]) * this->Athresh_;
        D_block_->replaceLocalValues(local_row, &InI[j], diagValues.getRawPtr(), 1);
      }
      else if (k < 0) { // Out of range
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::runtime_error, "Ifpack2::RILUK::initAllValues: current "
          "GID k = " << k << " < 0.  I'm not sure why this is an error; it is "
          "probably an artifact of the undocumented assumptions of the "
          "original implementation (likely copied and pasted from Ifpack).  "
          "Nevertheless, the code I found here insisted on this being an error "
          "state, so I will throw an exception here.");
      }
      else if (k < local_row) {
        LI[NumL] = k;
        const local_ordinal_type LBlockOffset = NumL*blockMatSize;
        for (local_ordinal_type jj = 0; jj < blockMatSize; ++jj)
          LV[LBlockOffset+jj] = InV[blockOffset+jj];
        NumL++;
      }
      else if (Teuchos::as<size_t>(k) <= rowMap->getNodeNumElements()) {
        UI[NumU] = k;
        const local_ordinal_type UBlockOffset = NumU*blockMatSize;
        for (local_ordinal_type jj = 0; jj < blockMatSize; ++jj)
          UV[UBlockOffset+jj] = InV[blockOffset+jj];
        NumU++;
      }
    }

    // Check in things for this row of L and U

    if (DiagFound) {
      ++NumNonzeroDiags;
    } else
    {
      for (local_ordinal_type jj = 0; jj < blockSize; ++jj)
        diagValues[jj*(blockSize+1)] = this->Athresh_;
      D_block_->replaceLocalValues(local_row, &local_row, diagValues.getRawPtr(), 1);
    }

    if (NumL) {
      L_block_->replaceLocalValues(local_row, &LI[0], &LV[0], NumL);
    }

    if (NumU) {
      U_block_->replaceLocalValues(local_row, &UI[0], &UV[0], NumU);
    }
  }

  this->isInitialized_ = true;
}

template<class MatrixType>
void RBILUK<MatrixType>::compute ()
{
  // initialize() checks this too, but it's easier for users if the
  // error shows them the name of the method that they actually
  // called, rather than the name of some internally called method.
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_block_.is_null (), std::runtime_error, "Ifpack2::Experimental::RBILUK::compute: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");

  if (! this->isInitialized ()) {
    initialize (); // Don't count this in the compute() time
  }

  typedef typename GetLapackType<impl_scalar_type>::lapack_scalar_type LST;
  typedef typename GetLapackType<impl_scalar_type>::lapack_type lapack_type;

  lapack_type lapack;

  Teuchos::Time timer ("RBILUK::compute");
  { // Start timing
    this->isComputed_ = false;

    // MinMachNum should be officially defined, for now pick something a little
    // bigger than IEEE underflow value

//    const scalar_type MinDiagonalValue = STS::rmin ();
//    const scalar_type MaxDiagonalValue = STS::one () / MinDiagonalValue;
    initAllValues (*A_block_);

    size_t NumIn;
    local_ordinal_type NumL, NumU;

    // Get Maximum Row length
    const size_t MaxNumEntries =
      L_block_->getNodeMaxNumRowEntries () + U_block_->getNodeMaxNumRowEntries () + 1;

    const local_ordinal_type blockSize = A_block_->getBlockSize();
    const local_ordinal_type blockMatSize = blockSize*blockSize;

    Teuchos::Array<int> ipiv(blockSize);
    Teuchos::Array<LST> work(1);

    size_t num_cols = U_block_->getColMap()->getNodeNumElements();
    Teuchos::Array<int> colflag(num_cols);

    Teuchos::Array<impl_scalar_type> diagMod(blockMatSize);
    little_block_type diagModBlock(&diagMod[0], blockSize, 1, blockSize);
    Teuchos::Array<impl_scalar_type> matTmpArray(blockMatSize);
    little_block_type matTmp(&matTmpArray[0], blockSize, 1, blockSize);
    Teuchos::Array<impl_scalar_type> multiplierArray(blockMatSize);
    little_block_type multiplier(&multiplierArray[0], blockSize, 1, blockSize);

//    Teuchos::ArrayRCP<scalar_type> DV = D_->get1dViewNonConst(); // Get view of diagonal

    // Now start the factorization.

    // Need some integer workspace and pointers
    local_ordinal_type NumUU;
    Teuchos::ArrayView<const local_ordinal_type> UUI;
    Teuchos::ArrayView<const scalar_type> UUV;
    for (size_t j = 0; j < num_cols; ++j) {
      colflag[j] = -1;
    }
    Teuchos::Array<local_ordinal_type> InI(MaxNumEntries);
    Teuchos::Array<scalar_type> InV(MaxNumEntries*blockMatSize);

    for (size_t i = 0; i < L_block_->getNodeNumRows (); ++i) {
      local_ordinal_type local_row = i;

      // Fill InV, InI with current row of L, D and U combined

      NumIn = MaxNumEntries;
      const local_ordinal_type * colValsL;
      scalar_type * valsL;

      L_block_->getLocalRowView(local_row, colValsL, valsL, NumL);
      for (local_ordinal_type j = 0; j < NumL; ++j)
      {
        const local_ordinal_type matOffset = blockMatSize*j;
        little_block_type lmat(&valsL[matOffset],blockSize,1,blockSize);
        little_block_type lmatV(&InV[matOffset],blockSize,1,blockSize);
        lmatV.assign(lmat);
        InI[j] = colValsL[j];
      }

      little_block_type dmat = D_block_->getLocalBlock(local_row, local_row);
      little_block_type dmatV(&InV[NumL*blockMatSize], blockSize, 1, blockSize);
      dmatV.assign(dmat);

      const local_ordinal_type * colValsU;
      scalar_type * valsU;
      U_block_->getLocalRowView(local_row, colValsU, valsU, NumU);
      for (local_ordinal_type j = 0; j < NumU; ++j)
      {
        InI[NumL+1+j] = colValsU[j];
        const local_ordinal_type matOffset = blockMatSize*(NumL+1+j);
        little_block_type umat(&valsU[blockMatSize*j], blockSize, 1, blockSize);
        little_block_type umatV(&InV[matOffset], blockSize, 1, blockSize);
        umatV.assign(umat);
      }
      NumIn = NumL+NumU+1;

      // Set column flags
      for (size_t j = 0; j < NumIn; ++j) {
        colflag[InI[j]] = j;
      }


      scalar_type diagmod = STM::zero (); // Off-diagonal accumulator
      diagModBlock.fill(diagmod);

      for (local_ordinal_type jj = 0; jj < NumL; ++jj) {
        local_ordinal_type j = InI[jj];
        little_block_type currentVal(&InV[jj*blockMatSize], blockSize, 1, blockSize); // current_mults++;
        multiplier.assign(currentVal);

        const little_block_type dmatInverse = D_block_->getLocalBlock(j,j);
        square_matrix_matrix_multiply(currentVal.getRawPtr(), dmatInverse.getRawPtr(), matTmp.getRawPtr(), blockSize);
        currentVal.assign(matTmp);

        const local_ordinal_type * UUI;
        scalar_type * UUV;
        U_block_->getLocalRowView(j, UUI, UUV, NumUU);

        if (this->RelaxValue_ == STM::zero ()) {
          for (local_ordinal_type k = 0; k < NumUU; ++k) {
            const int kk = colflag[UUI[k]];
            if (kk > -1) {
              little_block_type kkval(&InV[kk*blockMatSize], blockSize, 1, blockSize);
              little_block_type uumat(&UUV[k*blockMatSize], blockSize, 1, blockSize);
              square_matrix_matrix_multiply(multiplier.getRawPtr(), uumat.getRawPtr(), kkval.getRawPtr(), blockSize, -STM::one(), STM::one());
            }
          }
        }
        else {
          for (local_ordinal_type k = 0; k < NumUU; ++k) {
            const int kk = colflag[UUI[k]];
            little_block_type uumat(&UUV[k*blockMatSize], blockSize, 1, blockSize);
            if (kk > -1) {
              little_block_type kkval(&InV[kk*blockMatSize], blockSize, 1, blockSize);
              square_matrix_matrix_multiply(multiplier.getRawPtr(), uumat.getRawPtr(), kkval.getRawPtr(), blockSize, -STM::one(), STM::one());
            }
            else {
              square_matrix_matrix_multiply(multiplier.getRawPtr(), uumat.getRawPtr(), diagModBlock.getRawPtr(), blockSize, -STM::one(), STM::one());
            }
          }
        }
      }
      if (NumL) {
        // Replace current row of L
        L_block_->replaceLocalValues(local_row, InI.getRawPtr(), InV.getRawPtr(), NumL);
      }

      dmat.assign(dmatV);

      if (this->RelaxValue_ != STM::zero ()) {
        dmat.update(this->RelaxValue_, diagModBlock);
      }

//      if (STS::magnitude (DV[i]) > STS::magnitude (MaxDiagonalValue)) {
//        if (STS::real (DV[i]) < STM::zero ()) {
//          DV[i] = -MinDiagonalValue;
//        }
//        else {
//          DV[i] = MinDiagonalValue;
//        }
//      }
//      else
      {

        LST* const d_raw = reinterpret_cast<LST*> (dmat.getRawPtr());
        int lapackInfo;
        for (int i = 0; i < blockSize; ++i)
          ipiv[i] = 0;

        lapack.GETRF(blockSize, blockSize, d_raw, blockSize, ipiv.getRawPtr(), &lapackInfo);
        TEUCHOS_TEST_FOR_EXCEPTION(
          lapackInfo != 0, std::runtime_error, "Ifpack2::Experimental::RBILUK::compute: "
          "lapackInfo = " << lapackInfo << " which indicates an error in the factorization GETRF.");

        int lwork = -1;
        lapack.GETRI(blockSize, d_raw, blockSize, ipiv.getRawPtr(), work.getRawPtr(), lwork, &lapackInfo);
        TEUCHOS_TEST_FOR_EXCEPTION(
          lapackInfo != 0, std::runtime_error, "Ifpack2::Experimental::RBILUK::compute: "
          "lapackInfo = " << lapackInfo << " which indicates an error in the matrix inverse GETRI.");

        typedef typename Kokkos::Details::ArithTraits<impl_scalar_type>::mag_type ImplMagnitudeType;
        ImplMagnitudeType worksize = Kokkos::Details::ArithTraits<impl_scalar_type>::magnitude(work[0]);
        lwork = static_cast<int>(worksize);
        work.resize(lwork);
        lapack.GETRI(blockSize, d_raw, blockSize, ipiv.getRawPtr(), work.getRawPtr(), lwork, &lapackInfo);
        TEUCHOS_TEST_FOR_EXCEPTION(
          lapackInfo != 0, std::runtime_error, "Ifpack2::Experimental::RBILUK::compute: "
          "lapackInfo = " << lapackInfo << " which indicates an error in the matrix inverse GETRI.");
      }

      for (local_ordinal_type j = 0; j < NumU; ++j) {
        little_block_type currentVal(&InV[(NumL+1+j)*blockMatSize], blockSize, 1, blockSize); // current_mults++;
        // scale U by the diagonal inverse
        square_matrix_matrix_multiply(dmat.getRawPtr(), currentVal.getRawPtr(), matTmp.getRawPtr(), blockSize);
        currentVal.assign(matTmp);
      }

      if (NumU) {
        // Replace current row of L and U
        U_block_->replaceLocalValues (local_row, &InI[NumL+1], &InV[blockMatSize*(NumL+1)], NumU);
      }

      // Reset column flags
      for (size_t j = 0; j < NumIn; ++j) {
        colflag[InI[j]] = -1;
      }
    }

  } // Stop timing

  this->isComputed_ = true;
  this->numCompute_ += 1;
  this->computeTime_ += timer.totalElapsedTime ();
}


template<class MatrixType>
void
RBILUK<MatrixType>::
apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
       Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  using Teuchos::RCP;
  using Teuchos::rcpFromRef;
  typedef Tpetra::Experimental::BlockMultiVector<scalar_type,
    local_ordinal_type, global_ordinal_type, node_type> BMV;
  typedef Tpetra::MultiVector<scalar_type,
    local_ordinal_type, global_ordinal_type, node_type> MV;

  TEUCHOS_TEST_FOR_EXCEPTION(
    A_block_.is_null (), std::runtime_error, "Ifpack2::Experimental::RBILUK::apply: The matrix is "
    "null.  Please call setMatrix() with a nonnull input, then initialize() "
    "and compute(), before calling this method.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! this->isComputed (), std::runtime_error,
    "Ifpack2::Experimental::RBILUK::apply: If you have not yet called compute(), "
    "you must call compute() before calling this method.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.getNumVectors () != Y.getNumVectors (), std::invalid_argument,
    "Ifpack2::Experimental::RBILUK::apply: X and Y do not have the same number of columns.  "
    "X.getNumVectors() = " << X.getNumVectors ()
    << " != Y.getNumVectors() = " << Y.getNumVectors () << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(
    STS::isComplex && mode == Teuchos::CONJ_TRANS, std::logic_error,
    "Ifpack2::Experimental::RBILUK::apply: mode = Teuchos::CONJ_TRANS is not implemented for "
    "complex Scalar type.  Please talk to the Ifpack2 developers to get this "
    "fixed.  There is a FIXME in this file about this very issue.");

  const local_ordinal_type blockSize = A_block_->getBlockSize();
  const local_ordinal_type blockMatSize = blockSize*blockSize;

  BMV yBlock (Y, * (A_block_->getGraph ()->getDomainMap ()), blockSize);
  const BMV xBlock (X, * (A_block_->getColMap ()), blockSize);

  Teuchos::Array<scalar_type> lclarray(blockSize);
  little_vec_type lclvec(&lclarray[0], blockSize, 1);
  const scalar_type one = STM::one ();
  const scalar_type zero = STM::zero ();

  Teuchos::Time timer ("RBILUK::apply");
  { // Start timing
    Teuchos::TimeMonitor timeMon (timer);
    if (alpha == one && beta == zero) {
      if (mode == Teuchos::NO_TRANS) { // Solve L (D (U Y)) = X for Y.
        // Start by solving L C = X for C.  C must have the same Map
        // as D.  We have to use a temp multivector, since
        // localSolve() does not allow its input and output to alias
        // one another.
        //
        // FIXME (mfh 24 Jan 2014) Cache this temp multivector.
        const local_ordinal_type numVectors = xBlock.getNumVectors();
        BMV cBlock (* (A_block_->getGraph ()->getDomainMap ()), blockSize, numVectors);
        BMV rBlock (* (A_block_->getGraph ()->getDomainMap ()), blockSize, numVectors);
        for (local_ordinal_type imv = 0; imv < numVectors; ++imv)
        {
          for (size_t i = 0; i < D_block_->getNodeNumRows(); ++i)
          {
            local_ordinal_type local_row = i;
            little_vec_type xval = xBlock.getLocalBlock(local_row,imv);
            little_vec_type cval = cBlock.getLocalBlock(local_row,imv);
            cval.assign(xval);

            local_ordinal_type NumL;
            const local_ordinal_type * colValsL;
            scalar_type * valsL;

            L_block_->getLocalRowView(local_row, colValsL, valsL, NumL);

            for (local_ordinal_type j = 0; j < NumL; ++j)
            {
              local_ordinal_type col = colValsL[j];
              little_vec_type prevVal = cBlock.getLocalBlock(col, imv);

              const local_ordinal_type matOffset = blockMatSize*j;
              little_block_type lij(&valsL[matOffset],blockSize,1,blockSize);

              cval.matvecUpdate(-one, lij, prevVal);
            }
          }
        }

        // Solve D R = C. Note that D has been replaced by D^{-1} at this point.
        D_block_->applyBlock(cBlock, rBlock);

        // Solve U Y = R.
        for (local_ordinal_type imv = 0; imv < numVectors; ++imv)
        {
          const local_ordinal_type numRows = D_block_->getNodeNumRows();
          for (local_ordinal_type i = 0; i < numRows; ++i)
          {
            local_ordinal_type local_row = (numRows-1)-i;
            little_vec_type rval = rBlock.getLocalBlock(local_row,imv);
            little_vec_type yval = yBlock.getLocalBlock(local_row,imv);
            yval.assign(rval);

            local_ordinal_type NumU;
            const local_ordinal_type * colValsU;
            scalar_type * valsU;

            U_block_->getLocalRowView(local_row, colValsU, valsU, NumU);

            for (local_ordinal_type j = 0; j < NumU; ++j)
            {
              local_ordinal_type col = colValsU[NumU-1-j];
              little_vec_type prevVal = yBlock.getLocalBlock(col, imv);

              const local_ordinal_type matOffset = blockMatSize*(NumU-1-j);
              little_block_type uij(&valsU[matOffset], blockSize, 1, blockSize);

              yval.matvecUpdate(-one, uij, prevVal);
            }
          }
        }
      }
      else { // Solve U^P (D^P (L^P Y)) = X for Y (where P is * or T).
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::runtime_error,
          "Ifpack2::Experimental::RBILUK::apply: transpose apply is not implemented for the block algorithm. ");
      }
    }
    else { // alpha != 1 or beta != 0
      if (alpha == zero) {
        if (beta == zero) {
          Y.putScalar (zero);
        } else {
          Y.scale (beta);
        }
      } else { // alpha != zero
        MV Y_tmp (Y.getMap (), Y.getNumVectors ());
        apply (X, Y_tmp, mode);
        Y.update (alpha, Y_tmp, beta);
      }
    }
  } // Stop timing

  this->numApply_ += 1;
  this->applyTime_ = timer.totalElapsedTime ();
}


template<class MatrixType>
void RBILUK<MatrixType>::
multiply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
          Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
          const Teuchos::ETransp mode) const
{
  const scalar_type zero = STM::zero ();
  const scalar_type one = STM::one ();

  if (mode != Teuchos::NO_TRANS) {
    U_block_->apply (X, Y, mode); //
    Y.update (one, X, one); // Y = Y + X (account for implicit unit diagonal)

    MV Y_tmp (Y, Teuchos::Copy); // Need a temp copy of Y1
    D_block_->apply (Y, Y_tmp, mode, one, one);

    L_block_->apply (Y_tmp, Y, mode, one, one);
    Y.update (one, Y_tmp, one); // (account for implicit unit diagonal)
  }
  else {
    L_block_->apply (X, Y, mode); //
    Y.update (one, X, one); // Y = Y + X (account for implicit unit diagonal)

    MV Y_tmp (Y, Teuchos::Copy); // Need a temp copy of Y1
    D_block_->apply (Y, Y_tmp, mode, one, one);

    U_block_->apply (Y_tmp, Y, mode, one, one);
    Y.update (one, Y_tmp, one); // (account for implicit unit diagonal)
  }
}

template<class MatrixType>
std::string RBILUK<MatrixType>::description () const
{
  std::ostringstream os;

  // Output is a valid YAML dictionary in flow style.  If you don't
  // like everything on a single line, you should call describe()
  // instead.
  os << "\"Ifpack2::Experimental::RBILUK\": {";
  os << "Initialized: " << (this->isInitialized () ? "true" : "false") << ", "
     << "Computed: " << (this->isComputed () ? "true" : "false") << ", ";

  os << "Level-of-fill: " << this->getLevelOfFill() << ", ";

  if (A_block_.is_null ()) {
    os << "Matrix: null";
  }
  else {
    os << "Global matrix dimensions: ["
       << A_block_->getGlobalNumRows () << ", " << A_block_->getGlobalNumCols () << "]"
       << ", Global nnz: " << A_block_->getGlobalNumEntries();
  }

  os << "}";
  return os.str ();
}

} // namespace Experimental

} // namespace Ifpack2

#define IFPACK2_EXPERIMENTAL_RBILUK_INSTANT(S,LO,GO,N)                            \
  template class Ifpack2::Experimental::RBILUK< Tpetra::BlockCrsMatrix<S, LO, GO, N> >;

#endif
