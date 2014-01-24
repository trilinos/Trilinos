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

#ifndef IFPACK2_CRSRILUK_DEF_HPP
#define IFPACK2_CRSRILUK_DEF_HPP

namespace Ifpack2 {

template<class MatrixType>
RILUK<MatrixType>::RILUK (const Teuchos::RCP<const row_matrix_type>& Matrix_in)
  : A_ (Matrix_in),
    isOverlapped_ (false),
    LevelOfFill_ (0),
    LevelOfOverlap_ (0),
    isAllocated_ (false),
    isInitialized_ (false),
    isComputed_ (false),
    numInitialize_ (0),
    numCompute_ (0),
    numApply_ (0),
    RelaxValue_ (Teuchos::ScalarTraits<magnitude_type>::zero ()),
    Athresh_ (Teuchos::ScalarTraits<magnitude_type>::zero ()),
    Rthresh_ (Teuchos::ScalarTraits<magnitude_type>::one ()),
    Condest_ (-Teuchos::ScalarTraits<magnitude_type>::one ()),
    OverlapMode_ (Tpetra::REPLACE)
{}


template<class MatrixType>
RILUK<MatrixType>::RILUK (const Teuchos::RCP<const crs_matrix_type>& Matrix_in)
  : A_ (Matrix_in),
    isOverlapped_ (false),
    LevelOfFill_ (0),
    LevelOfOverlap_ (0),
    isAllocated_ (false),
    isInitialized_ (false),
    isComputed_ (false),
    numInitialize_ (0),
    numCompute_ (0),
    numApply_ (0),
    RelaxValue_ (Teuchos::ScalarTraits<magnitude_type>::zero ()),
    Athresh_ (Teuchos::ScalarTraits<magnitude_type>::zero ()),
    Rthresh_ (Teuchos::ScalarTraits<magnitude_type>::one ()),
    Condest_ (-Teuchos::ScalarTraits<magnitude_type>::one ()),
    OverlapMode_ (Tpetra::REPLACE)
{}


template<class MatrixType>
RILUK<MatrixType>::~RILUK() {}


template<class MatrixType>
void RILUK<MatrixType>::allocate_L_and_U()
{
  using Teuchos::null;
  using Teuchos::rcp;

  if (! isAllocated_) {
    // Deallocate any existing storage.  This avoids storing 2x
    // memory, since RCP op= does not deallocate until after the
    // assignment.
    L_ = null;
    U_ = null;
    D_ = null;

    // Allocate Matrix using ILUK graphs
    L_ = rcp (new crs_matrix_type (Graph_->getL_Graph ()));
    U_ = rcp (new crs_matrix_type (Graph_->getU_Graph ()));
    L_->setAllToScalar (STS::zero ()); // Zero out L and U matrices
    U_->setAllToScalar (STS::zero ());

    // FIXME (mfh 24 Jan 2014) This assumes domain == range Map for L and U.
    L_->fillComplete ();
    U_->fillComplete ();

    D_ = rcp (new vec_type (Graph_->getL_Graph ()->getRowMap ()));
  }
  isAllocated_ = true;
}


template<class MatrixType>
void
RILUK<MatrixType>::
setParameters (const Teuchos::ParameterList& params)
{
  using Teuchos::as;
  using Teuchos::Exceptions::InvalidParameterName;
  using Teuchos::Exceptions::InvalidParameterType;

  // Default values of the various parameters.
  int fillLevel = 0;
  int overlapLevel = 0;
  magnitude_type absThresh = STM::zero ();
  magnitude_type relThresh = STM::one ();
  magnitude_type relaxValue = STM::zero ();

  //
  // "fact: iluk level-of-fill" parsing is more complicated, because
  // we want to allow as many types as make sense.  int is the native
  // type, but we also want to accept magnitude_type (for
  // compatibility with ILUT) and double (for backwards compatibilty
  // with ILUT).
  //

  bool gotFillLevel = false;
  try {
    fillLevel = params.get<int> ("fact: iluk level-of-fill");
    gotFillLevel = true;
  }
  catch (InvalidParameterType&) {
    // Throwing again in the catch block would just unwind the stack.
    // Instead, we do nothing here, and check the Boolean outside to
    // see if we got the value.
  }
  catch (InvalidParameterName&) {
    gotFillLevel = true; // Accept the default value.
  }

  if (! gotFillLevel) {
    try {
      // Try magnitude_type, for compatibility with ILUT.
      // The cast from magnitude_type to int must succeed.
      fillLevel = as<int> (params.get<magnitude_type> ("fact: iluk level-of-fill"));
      gotFillLevel = true;
    }
    catch (InvalidParameterType&) {
      // Try double next.
    }
    // Don't catch InvalidParameterName here; we've already done that above.
  }

  if (! gotFillLevel) {
    try {
      // Try double, for compatibility with ILUT.
      // The cast from double to int must succeed.
      fillLevel = as<int> (params.get<double> ("fact: iluk level-of-fill"));
      gotFillLevel = true;
    }
    catch (InvalidParameterType& e) {
      // We're out of options.  The user gave us the parameter, but it
      // doesn't have the right type.  The best thing for us to do in
      // that case is to throw, telling the user to use the right
      // type.
      throw e;
    }
    // Don't catch InvalidParameterName here; we've already done that above.
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! gotFillLevel,
    std::logic_error,
    "Ifpack2::RILUK::setParameters: We should never get here!  "
    "The method should either have read the \"fact: iluk level-of-fill\"  "
    "parameter by this point, or have thrown an exception.  "
    "Please let the Ifpack2 developers know about this bug.");

  //
  // overlapLevel was always int.  ILUT doesn't have this parameter.
  // However, some tests (as of 28 Nov 2012:
  // Ifpack2_RILUK_small_belos_MPI_1) depend on being able to set this
  // as a double instead of an int.  Thus, we go through the same
  // procedure as above with fill level.
  //

  bool gotOverlapLevel = false;
  try {
    overlapLevel = params.get<int> ("fact: iluk level-of-overlap");
    gotOverlapLevel = true;
  }
  catch (InvalidParameterType&) {
    // Throwing again in the catch block would just unwind the stack.
    // Instead, we do nothing here, and check the Boolean outside to
    // see if we got the value.
  }
  catch (InvalidParameterName&) {
    gotOverlapLevel = true; // Accept the default value.
  }

  if (! gotOverlapLevel) {
    try {
      // Try magnitude_type, for compatibility with ILUT.
      // The cast from magnitude_type to int must succeed.
      overlapLevel = as<int> (params.get<magnitude_type> ("fact: iluk level-of-overlap"));
      gotOverlapLevel = true;
    }
    catch (InvalidParameterType&) {
      // Try double next.
    }
    // Don't catch InvalidParameterName here; we've already done that above.
  }

  if (! gotOverlapLevel) {
    try {
      // Try double, for compatibility with ILUT.
      // The cast from double to int must succeed.
      overlapLevel = as<int> (params.get<double> ("fact: iluk level-of-overlap"));
      gotOverlapLevel = true;
    }
    catch (InvalidParameterType& e) {
      // We're out of options.  The user gave us the parameter, but it
      // doesn't have the right type.  The best thing for us to do in
      // that case is to throw, telling the user to use the right
      // type.
      throw e;
    }
    // Don't catch InvalidParameterName here; we've already done that above.
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! gotOverlapLevel,
    std::logic_error,
    "Ifpack2::RILUK::setParameters: We should never get here!  "
    "The method should either have read the \"fact: iluk level-of-overlap\"  "
    "parameter by this point, or have thrown an exception.  "
    "Please let the Ifpack2 developers know about this bug.");

  //
  // For the other parameters, we prefer magnitude_type, but allow
  // double for backwards compatibility.
  //

  try {
    absThresh = params.get<magnitude_type> ("fact: absolute threshold");
  }
  catch (InvalidParameterType&) {
    // Try double, for backwards compatibility.
    // The cast from double to magnitude_type must succeed.
    absThresh = as<magnitude_type> (params.get<double> ("fact: absolute threshold"));
  }
  catch (InvalidParameterName&) {
    // Accept the default value.
  }

  try {
    relThresh = params.get<magnitude_type> ("fact: relative threshold");
  }
  catch (InvalidParameterType&) {
    // Try double, for backwards compatibility.
    // The cast from double to magnitude_type must succeed.
    relThresh = as<magnitude_type> (params.get<double> ("fact: relative threshold"));
  }
  catch (InvalidParameterName&) {
    // Accept the default value.
  }

  try {
    relaxValue = params.get<magnitude_type> ("fact: relax value");
  }
  catch (InvalidParameterType&) {
    // Try double, for backwards compatibility.
    // The cast from double to magnitude_type must succeed.
    relaxValue = as<magnitude_type> (params.get<double> ("fact: relax value"));
  }
  catch (InvalidParameterName&) {
    // Accept the default value.
  }

  // "Commit" the values only after validating all of them.  This
  // ensures that there are no side effects if this routine throws an
  // exception.

  LevelOfFill_ = fillLevel;
  LevelOfOverlap_ = overlapLevel;

  // mfh 28 Nov 2012: The previous code would not assign Athresh_,
  // Rthresh_, or RelaxValue_, if the read-in value was -1.  I don't
  // know if keeping this behavior is correct, but I'll keep it just
  // so as not to change previous behavior.

  if (absThresh != -STM::one ()) {
    Athresh_ = absThresh;
  }
  if (relThresh != -STM::one ()) {
    Rthresh_ = relThresh;
  }
  if (relaxValue != -STM::one ()) {
    RelaxValue_ = relaxValue;
  }
}


template<class MatrixType>
Teuchos::RCP<const typename RILUK<MatrixType>::row_matrix_type>
RILUK<MatrixType>::getMatrix () const {
  return Teuchos::rcp_implicit_cast<const row_matrix_type> (A_);
}


template<class MatrixType>
Teuchos::RCP<const typename RILUK<MatrixType>::crs_matrix_type>
RILUK<MatrixType>::getCrsMatrix () const {
  return Teuchos::rcp_dynamic_cast<const crs_matrix_type> (A_, true);
}



template<class MatrixType>
void RILUK<MatrixType>::initialize ()
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  typedef Tpetra::CrsGraph<local_ordinal_type,
                           global_ordinal_type,
                           node_type> crs_graph_type;

  // Calling initialize() means that the user asserts that the graph
  // of the sparse matrix may have changed.  We must not just reuse
  // the previous graph in that case.
  //
  // Regarding setting isAllocated_ to false: Eventually, we may want
  // some kind of clever memory reuse strategy, but it's always
  // correct just to blow everything away and start over.
  isInitialized_ = false;
  isAllocated_ = false;
  isComputed_ = false;
  Graph_ = Teuchos::null;

  // mfh 13 Dec 2013: If it's a Tpetra::CrsMatrix, just give Graph_
  // its Tpetra::CrsGraph.  Otherwise, we might need to rewrite
  // IlukGraph to handle a Tpetra::RowGraph.  Just throw an exception
  // for now.
  {
    RCP<const crs_matrix_type> A_crs =
      rcp_dynamic_cast<const crs_matrix_type> (A_);
    if (A_crs.is_null ()) {
      // FIXME (mfh 24 Jan 2014) It would be more efficient to rewrite
      // RILUK so that it works with any RowMatrix input, not just
      // CrsMatrix.  However, to make it work for now, we can copy the
      // input matrix and hope for the best.

      // FIXME (mfh 24 Jan 2014) It would be smarter to count up the
      // number of elements in each row of A_, so that we can create
      // A_crs_nc using static profile.  The code below is correct but
      // potentially slow.
      RCP<crs_matrix_type> A_crs_nc =
        rcp (new crs_matrix_type (A_->getRowMap (), A_->getColMap (), 0));

      // FIXME (mfh 24 Jan 2014) This Import approach will only work
      // if A_ has a one-to-one row Map.  This is generally the case
      // with matrices given to Ifpack2.
      typedef Tpetra::Import<local_ordinal_type, global_ordinal_type,
        node_type> import_type;

      // Source and destination Maps are the same in this case.
      // That way, the Import just implements a copy.
      import_type import (A_->getRowMap (), A_->getRowMap ());
      A_crs_nc->doImport (*A_, import, Tpetra::REPLACE);
      A_crs_nc->fillComplete (A_->getDomainMap (), A_->getRangeMap ());
      A_crs = rcp_const_cast<const crs_matrix_type> (A_crs_nc);
    }
    A_crs_ = A_crs;
    Graph_ = rcp (new Ifpack2::IlukGraph<crs_graph_type> (A_crs->getCrsGraph (),
                                                          LevelOfFill_,
                                                          LevelOfOverlap_));
  }

  Graph_->initialize ();
  allocate_L_and_U ();

  // FIXME (mfh 24 Jan 2014) RILUK should not be responsible for
  // overlap; AdditiveSchwarz should be.
  if (isOverlapped_) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "Ifpack2::RILUK::initialize: This class does not implement overlapping.  "
      "You should use AdditiveSchwarz as the outer preconditioner if you want "
      "overlapping.");
//    OverlapA = Teuchos::rcp( new MatrixType(Graph_->OverlapGraph()) );
//    EPETRA_CHK_ERR(OverlapA->Import(A, *Graph_->OverlapImporter(), Insert));
//    EPETRA_CHK_ERR(OverlapA->FillComplete());
  }
  else {
    initAllValues (*A_crs_);
  }

  isInitialized_ = true;
  ++numInitialize_;
}


template<class MatrixType>
void
RILUK<MatrixType>::
initAllValues (const row_matrix_type& OverlapA)
{
  using Teuchos::ArrayRCP;
  using Teuchos::Comm;
  using Teuchos::ptr;
  using Teuchos::RCP;
  using Teuchos::REDUCE_SUM;
  using Teuchos::reduceAll;
  typedef Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> map_type;

  size_t NumIn = 0, NumL = 0, NumU = 0;
  bool DiagFound = false;
  size_t NumNonzeroDiags = 0;
  size_t MaxNumEntries = OverlapA.getGlobalMaxNumRowEntries();

  Teuchos::Array<global_ordinal_type> InI(MaxNumEntries); // Allocate temp space
  Teuchos::Array<global_ordinal_type> LI(MaxNumEntries);
  Teuchos::Array<global_ordinal_type> UI(MaxNumEntries);
  Teuchos::Array<scalar_type> InV(MaxNumEntries);
  Teuchos::Array<scalar_type> LV(MaxNumEntries);
  Teuchos::Array<scalar_type> UV(MaxNumEntries);

  // Check if values should be inserted or replaced
  const bool ReplaceValues = L_->isStaticGraph () || L_->isLocallyIndexed ();

  L_->resumeFill ();
  U_->resumeFill ();
  if (ReplaceValues) {
    L_->setAllToScalar (STS::zero ()); // Zero out L and U matrices
    U_->setAllToScalar (STS::zero ());
  }

  D_->putScalar (STS::zero ()); // Set diagonal values to zero
  ArrayRCP<scalar_type> DV = D_->get1dViewNonConst (); // Get view of diagonal

  RCP<const map_type> rowMap = L_->getRowMap ();

  // First we copy the user's matrix into L and U, regardless of fill level

  // FIXME (mfh 24 Jan 2014) This assumes that the row Map's global
  // indices are contiguous on the calling process!
  for (global_ordinal_type i = rowMap->getMinGlobalIndex ();
       i <= rowMap->getMaxGlobalIndex (); ++i) {
    global_ordinal_type global_row = i;
    local_ordinal_type local_row = rowMap->getLocalElement (global_row);

    OverlapA.getGlobalRowCopy (global_row, InI(), InV(), NumIn); // Get Values and Indices

    // Split into L and U (we don't assume that indices are ordered).

    NumL = 0;
    NumU = 0;
    DiagFound = false;

    for (size_t j = 0; j < NumIn; ++j) {
      const global_ordinal_type k = InI[j];

      if (k==i) {
        DiagFound = true;
        DV[local_row] += Rthresh_ * InV[j] + IFPACK2_SGN(InV[j]) * Athresh_; // Store perturbed diagonal in Tpetra::Vector D_
      }
      else if (k < 0) { // Out of range
        throw std::runtime_error("out of range (k<0) in Ifpack2::RILUK::initAllValues");
      }
      else if (k < i) {
        LI[NumL] = k;
        LV[NumL] = InV[j];
        NumL++;
      }
      else if (k <= rowMap->getMaxGlobalIndex()) {
        UI[NumU] = k;
        UV[NumU] = InV[j];
        NumU++;
      }
//      else {
//        throw std::runtime_error("out of range in Ifpack2::RILUK::initAllValues");
//      }
    }

    // Check in things for this row of L and U

    if (DiagFound) {
      ++NumNonzeroDiags;
    } else {
      DV[local_row] = Athresh_;
    }

    if (NumL) {
      if (ReplaceValues) {
        L_->replaceGlobalValues(global_row, LI(0, NumL), LV(0,NumL));
      } else {
        L_->insertGlobalValues(global_row, LI(0,NumL), LV(0,NumL));
      }
    }

    if (NumU) {
      if (ReplaceValues) {
        U_->replaceGlobalValues(global_row, UI(0,NumU), UV(0,NumU));
      } else {
        U_->insertGlobalValues(global_row, UI(0,NumU), UV(0,NumU));
      }
    }
  }

  // The domain of L and the range of U are exactly their own row maps
  // (there is no communication).  The domain of U and the range of L
  // must be the same as those of the original matrix, However if the
  // original matrix is a VbrMatrix, these two latter maps are
  // translation from a block map to a point map.
  L_->fillComplete (L_->getColMap (), A_->getRangeMap ());
  U_->fillComplete (A_->getDomainMap (), U_->getRowMap ());

  // At this point L and U have the values of A in the structure of L
  // and U, and diagonal vector D

  isInitialized_ = true;
}


template<class MatrixType>
void RILUK<MatrixType>::compute ()
{
  if (! isInitialized ()) {
    initialize ();
  }
  isComputed_ = false;

  L_->resumeFill ();
  U_->resumeFill ();

  // MinMachNum should be officially defined, for now pick something a little
  // bigger than IEEE underflow value

  const scalar_type MinDiagonalValue = STS::rmin ();
  const scalar_type MaxDiagonalValue = STS::one () / MinDiagonalValue;

  size_t NumIn, NumL, NumU;

  // Get Maximum Row length
  const size_t MaxNumEntries =
    L_->getNodeMaxNumRowEntries () + U_->getNodeMaxNumRowEntries () + 1;

  Teuchos::Array<local_ordinal_type> InI(MaxNumEntries); // Allocate temp space
  Teuchos::Array<scalar_type> InV(MaxNumEntries);
  size_t num_cols = U_->getColMap()->getNodeNumElements();
  Teuchos::Array<int> colflag(num_cols);

  Teuchos::ArrayRCP<scalar_type> DV = D_->get1dViewNonConst(); // Get view of diagonal

  size_t current_madds = 0; // We will count multiply-add as they happen

  // Now start the factorization.

  // Need some integer workspace and pointers
  size_t NumUU;
  Teuchos::ArrayView<const local_ordinal_type> UUI;
  Teuchos::ArrayView<const scalar_type> UUV;
  for (size_t j = 0; j < num_cols; ++j) {
    colflag[j] = -1;
  }

  for (size_t i = 0; i < L_->getNodeNumRows (); ++i) {
    local_ordinal_type local_row = i;

 // Fill InV, InI with current row of L, D and U combined

    NumIn = MaxNumEntries;
    L_->getLocalRowCopy (local_row, InI (), InV (), NumL);

    InV[NumL] = DV[i]; // Put in diagonal
    InI[NumL] = local_row;

    U_->getLocalRowCopy (local_row, InI (NumL+1, MaxNumEntries-NumL-1),
                         InV (NumL+1, MaxNumEntries-NumL-1), NumU);
    NumIn = NumL+NumU+1;

    // Set column flags
    for (size_t j = 0; j < NumIn; ++j) {
      colflag[InI[j]] = j;
    }

    scalar_type diagmod = STS::zero (); // Off-diagonal accumulator

    for (size_t jj = 0; jj < NumL; ++jj) {
      local_ordinal_type j = InI[jj];
      scalar_type multiplier = InV[jj]; // current_mults++;

      InV[jj] *= DV[j];

      U_->getLocalRowView(j, UUI, UUV); // View of row above
      NumUU = UUI.size();

      if (RelaxValue_ == STM::zero ()) {
        for (size_t k = 0; k < NumUU; ++k) {
          const int kk = colflag[UUI[k]];
          // FIXME (mfh 23 Dec 2013) Wait a second, we just set
          // colflag above using size_t (which is generally unsigned),
          // but now we're querying it using int (which is signed).
          if (kk > -1) {
            InV[kk] -= multiplier * UUV[k];
            current_madds++;
          }
        }
      }
      else {
        for (size_t k = 0; k < NumUU; ++k) {
          // FIXME (mfh 23 Dec 2013) Wait a second, we just set
          // colflag above using size_t (which is generally unsigned),
          // but now we're querying it using int (which is signed).
          const int kk = colflag[UUI[k]];
          if (kk > -1) {
            InV[kk] -= multiplier*UUV[k];
          }
          else {
            diagmod -= multiplier*UUV[k];
          }
          current_madds++;
        }
      }
    }
    if (NumL) {
      // Replace current row of L
      L_->replaceLocalValues (local_row, InI (0, NumL), InV (0, NumL));
    }

    DV[i] = InV[NumL]; // Extract Diagonal value

    if (RelaxValue_ != STM::zero ()) {
      DV[i] += RelaxValue_*diagmod; // Add off diagonal modifications
      // current_madds++;
    }

    if (STS::magnitude (DV[i]) > STS::magnitude (MaxDiagonalValue)) {
      if (STS::real (DV[i]) < STM::zero ()) {
        DV[i] = -MinDiagonalValue;
      }
      else {
        DV[i] = MinDiagonalValue;
      }
    }
    else {
      DV[i] = STS::one () / DV[i]; // Invert diagonal value
    }

    for (size_t j = 0; j < NumU; ++j) {
      InV[NumL+1+j] *= DV[i]; // Scale U by inverse of diagonal
    }

    if (NumU) {
      // Replace current row of L and U
      U_->replaceLocalValues (local_row, InI (NumL+1, NumU), InV (NumL+1, NumU));
    }

    // Reset column flags
    for (size_t j = 0; j < NumIn; ++j) {
      colflag[InI[j]] = -1;
    }
  }

  // FIXME (mfh 23 Dec 2013) Do we know that the column Map of L_ is
  // always one-to-one?
  L_->fillComplete (L_->getColMap (), A_->getRangeMap ());
  U_->fillComplete (A_->getDomainMap (), U_->getRowMap ());

  // Validate that the L and U factors are actually lower and upper triangular

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! L_->isLowerTriangular (), std::runtime_error,
    "Ifpack2::RILUK::compute: L isn't lower triangular.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! U_->isUpperTriangular (), std::runtime_error,
    "Ifpack2::RILUK::compute: U isn't lower triangular.");

  // Add up flops

  double current_flops = 2.0 * current_madds;
  double total_flops = 0.0;

  // Get total madds across all PEs
  Teuchos::reduceAll (* (L_->getRowMap ()->getComm ()), Teuchos::REDUCE_SUM,
                      1, &current_flops, &total_flops);

  // Now count the rest
  total_flops += static_cast<double> (L_->getGlobalNumEntries ()); // multiplier above
  total_flops += static_cast<double> (D_->getGlobalLength ()); // reciprocal of diagonal
  if (RelaxValue_ != STM::zero ()) {
    total_flops += 2 * static_cast<double> (D_->getGlobalLength ()); // relax update of diag
  }

  //UpdateFlops(total_flops); // Update flop count

  isComputed_ = true;
  ++numCompute_;
}


template<class MatrixType>
void
RILUK<MatrixType>::
apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
       Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! isComputed (), std::runtime_error,
    "Ifpack2::RILUK::apply: If you have not yet called compute(), "
    "you must call compute() before calling this method.");

  if (alpha == STS::one () && beta == STS::zero ()) {
    //
    // This function finds Y such that
    // LDU Y = X, or
    // U(trans) D L(trans) Y = X
    // for multiple RHS
    //
    // First generate X and Y as needed for this function
    Teuchos::RCP<const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > X1;
    Teuchos::RCP<Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > Y1;
    generateXY (mode, X, Y, X1, Y1);

    const scalar_type one = Teuchos::ScalarTraits<scalar_type>::one();
    const scalar_type zero = Teuchos::ScalarTraits<scalar_type>::zero();
    if (mode == Teuchos::NO_TRANS) {
      L_->localSolve (*X1, *Y1,mode);
      Y1->elementWiseMultiply (one, *D_, *Y1, zero); // y = D*y (D_ has inverse of diagonal)
      U_->localSolve (*Y1, *Y1, mode); // Solve Uy = y
      if (isOverlapped_) {
        // Export computed Y values if needed
        Y.doExport (*Y1, *L_->getGraph ()->getExporter (), OverlapMode_);
      }
    }
    else {
      U_->localSolve (*X1, *Y1, mode); // Solve Uy = y
      Y1->elementWiseMultiply (one, *D_, *Y1, zero); // y = D*y (D_ has inverse of diagonal)
      L_->localSolve (*Y1, *Y1,mode);
      if (isOverlapped_) {
        // Export computed Y values if needed
        Y.doExport (*Y1, *U_->getGraph ()->getImporter (), OverlapMode_);
      }
    }
    ++numApply_; // inside the 'if', so it's not called twice
  }
  else { // alpha != 1 or beta != 0
    MV Y_tmp (Y.getMap (), Y.getNumVectors ());
    apply (X, Y_tmp, mode);
    Y.update (alpha, Y_tmp, beta);
  }
}

//=============================================================================
template<class MatrixType>
int RILUK<MatrixType>::Multiply(const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
            Teuchos::ETransp mode) const
{
//
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
//

  // First generate X and Y as needed for this function
  Teuchos::RCP<const MV> X1;
  Teuchos::RCP<MV> Y1;
  generateXY (mode, X, Y, X1, Y1);

//  Epetra_Flops * counter = this->GetFlopCounter();
//  if (counter!=0) {
//    L_->SetFlopCounter(*counter);
//    Y1->SetFlopCounter(*counter);
//    U_->SetFlopCounter(*counter);
//  }

  const scalar_type zero = Teuchos::ScalarTraits<scalar_type>::zero ();
  const scalar_type one = Teuchos::ScalarTraits<scalar_type>::one ();

  if (!mode == Teuchos::NO_TRANS) {
    U_->apply (*X1, *Y1,mode); //
    Y1->update (one, *X1, one); // Y1 = Y1 + X1 (account for implicit unit diagonal)
    Y1->elementWiseMultiply (one, *D_, *Y1, zero); // y = D*y (D_ has inverse of diagonal)
    MV Y1temp (*Y1); // Need a temp copy of Y1
    L_->apply (Y1temp, *Y1,mode);
    Y1->update (one, Y1temp, one); // (account for implicit unit diagonal)
    if (isOverlapped_) {
      Y.doExport (*Y1, *L_->getGraph ()->getExporter (), OverlapMode_);
    } // Export computed Y values if needed
  }
  else {
    L_->apply (*X1, *Y1,mode);
    Y1->update (one, *X1, one); // Y1 = Y1 + X1 (account for implicit unit diagonal)
    Y1->elementWiseMultiply (one, *D_, *Y1, zero); // y = D*y (D_ has inverse of diagonal)
    MV Y1temp (*Y1); // Need a temp copy of Y1
    U_->apply (Y1temp, *Y1,mode);
    Y1->update (one, Y1temp, one); // (account for implicit unit diagonal)
    if (isOverlapped_) {
      Y.doExport(*Y1, *L_->getGraph ()->getExporter (), OverlapMode_);
    }
  }
  return 0;
}


template<class MatrixType>
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType
RILUK<MatrixType>::computeCondEst (Teuchos::ETransp mode) const
{
  if (Condest_ != -Teuchos::ScalarTraits<magnitude_type>::one ()) {
    return Condest_;
  }
  // Create a vector with all values equal to one
  vec_type ones (U_->getDomainMap ());
  vec_type onesResult (L_->getRangeMap ());
  ones.putScalar (Teuchos::ScalarTraits<scalar_type>::one ());

  apply (ones, onesResult, mode); // Compute the effect of the solve on the vector of ones
  onesResult.abs (onesResult); // Make all values non-negative
  Teuchos::Array<magnitude_type> norms (1);
  onesResult.normInf (norms ());
  Condest_ = norms[0];
  return Condest_;
}


template<class MatrixType>
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType
RILUK<MatrixType>::
computeCondEst (CondestType CT,
                local_ordinal_type MaxIters,
                magnitude_type Tol,
                const Teuchos::Ptr<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >& Matrix)
{
  // Forestall "unused variable" compiler warnings.
  (void) CT;
  (void) MaxIters;
  (void) Tol;
  (void) Matrix;

  if (Condest_ != -Teuchos::ScalarTraits<magnitude_type>::one() ) {
    return Condest_;
  }
  // Create a vector with all values equal to one
  vec_type ones (U_->getDomainMap ());
  vec_type onesResult (L_->getRangeMap ());
  ones.putScalar (Teuchos::ScalarTraits<scalar_type>::one ());

  // Compute the effect of the solve on the vector of ones
  apply (ones, onesResult, Teuchos::NO_TRANS);
  onesResult.abs (onesResult); // Make all values non-negative
  Teuchos::Array<magnitude_type> norms (1);
  onesResult.normInf (norms ());
  Condest_ = norms[0];
  return Condest_;
}


template<class MatrixType>
void
RILUK<MatrixType>::
generateXY (Teuchos::ETransp mode,
            const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Xin,
            const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Yin,
            Teuchos::RCP<const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >& Xout,
            Teuchos::RCP<Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >& Yout) const
{
  // mfh 05 Dec 2013: I just cleaned up the implementation of this
  // method to make it more legible; I didn't otherwise try to improve
  // it.  I don't know if it works or not.

  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  // Generate an X and Y suitable for solve() and multiply()

  TEUCHOS_TEST_FOR_EXCEPTION(
    Xin.getNumVectors () != Yin.getNumVectors (), std::runtime_error,
    "Ifpack2::RILUK::generateXY: X and Y do not have the same number of "
    "columns.  X has " << Xin.getNumVectors () << " columns, but Y has "
    << Yin.getNumVectors () << " columns.");

  Xout = rcpFromRef (Xin);
  Yout = rcpFromRef (const_cast<MV&> (Yin));
  if (! isOverlapped_) {
    return; // Nothing more to do
  }

  if (isOverlapped_) {
    // Make sure the number of vectors in the multivector is the same as before.
    if (! OverlapX_.is_null ()) {
      if (OverlapX_->getNumVectors () != Xin.getNumVectors ()) {
        OverlapX_ = Teuchos::null;
        OverlapY_ = Teuchos::null;
      }
    }
    if (OverlapX_.is_null ()) { // Need to allocate space for overlap X and Y
      OverlapX_ = rcp (new MV (U_->getColMap (), Xout->getNumVectors ()));
      OverlapY_ = rcp (new MV (L_->getRowMap (), Yout->getNumVectors ()));
    }

    // Import X values for solve
    if (mode == Teuchos::NO_TRANS) {
      OverlapX_->doImport (*Xout, *U_->getGraph ()->getImporter (), Tpetra::INSERT);
    } else {
      OverlapX_->doImport (*Xout, *L_->getGraph ()->getExporter (), Tpetra::INSERT);
    }
    Xout = OverlapX_;
    Yout = OverlapY_; // Set pointers for Xout and Yout to point to overlap space
  }
}

}//namespace Ifpack2

#endif

