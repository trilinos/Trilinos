//@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2010) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef IFPACK2_CRSRILUK_DEF_HPP
#define IFPACK2_CRSRILUK_DEF_HPP

namespace Ifpack2 {

//==============================================================================
template<class MatrixType>
RILUK<MatrixType>::RILUK(const Teuchos::RCP<const MatrixType>& Matrix_in)
  : isOverlapped_(false),
    Graph_(),
    A_(Matrix_in),
    UseTranspose_(false),
    LevelOfFill_(0),
    LevelOfOverlap_(0),
    NumMyDiagonals_(0),
    isAllocated_(false),
    isInitialized_(false),
    numInitialize_(0),
    numCompute_(0),
    numApply_(0),
    Factored_(false),
    RelaxValue_(0.0),
    Athresh_(0.0),
    Rthresh_(1.0),
    Condest_(-1.0),
    OverlapMode_(Tpetra::REPLACE)
{
}

//==============================================================================
//template<class MatrixType>
//RILUK<MatrixType>::RILUK(const RILUK<MatrixType>& src)
//  : isOverlapped_(src.isOverlapped_),
//    Graph_(src.Graph_),
//    UseTranspose_(src.UseTranspose_),
//    LevelOfFill_(src.LevelOfFill_),
//    LevelOfOverlap_(src.LevelOfOverlap_),
//    NumMyDiagonals_(src.NumMyDiagonals_),
//    isAllocated_(src.isAllocated_),
//    isInitialized_(src.isInitialized_),
//    Factored_(src.Factored_),
//    RelaxValue_(src.RelaxValue_),
//    Athresh_(src.Athresh_),
//    Rthresh_(src.Rthresh_),
//    Condest_(src.Condest_),
//    OverlapMode_(src.OverlapMode_)
//{
//  L_ = Teuchos::rcp( new MatrixType(src.getL()) );
//  U_ = Teuchos::rcp( new MatrixType(src.getU()) );
//  D_ = Teuchos::rcp( new Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>(src.getD()) );
//}

//==============================================================================
template<class MatrixType>
RILUK<MatrixType>::~RILUK() {
}

//==============================================================================
template<class MatrixType>
void RILUK<MatrixType>::allocate_L_and_U() {

  // Allocate Matrix using ILUK graphs
  L_ = Teuchos::rcp( new MatrixType(Graph_->getL_Graph()) );
  U_ = Teuchos::rcp( new MatrixType(Graph_->getU_Graph()) );
  L_->setAllToScalar(0.0); // Zero out L and U matrices
  U_->setAllToScalar(0.0);
  L_->fillComplete();
  U_->fillComplete();
  bool isLower = L_->isLowerTriangular();
  bool isUpper = U_->isUpperTriangular();
  if (!isLower || !isUpper) {
    std::cout << "error in triangular detection" << std::endl;
  }

  D_ = Teuchos::rcp( new Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>(Graph_->getL_Graph()->getRowMap()) );
  setAllocated(true);
}

//==========================================================================
template<class MatrixType>
void
RILUK<MatrixType>::
setParameters (const Teuchos::ParameterList& params)
{
  using Teuchos::as;
  using Teuchos::Exceptions::InvalidParameterName;
  using Teuchos::Exceptions::InvalidParameterType;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;

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

//==========================================================================
template<class MatrixType>
void RILUK<MatrixType>::initialize() {

  if (Graph_ != Teuchos::null) return;

  Graph_ = Teuchos::rcp(new Ifpack2::IlukGraph<local_ordinal_type,global_ordinal_type,node_type>(A_->getCrsGraph(), LevelOfFill_, LevelOfOverlap_));

  Graph_->constructFilledGraph();

  setInitialized(true);

  if (!isAllocated()) allocate_L_and_U();

  if (isOverlapped_) {
    throw std::runtime_error("Error, Ifpack2::RILUK::initialize: overlapping not yet supported.");
//    OverlapA = Teuchos::rcp( new MatrixType(Graph_->OverlapGraph()) );
//    EPETRA_CHK_ERR(OverlapA->Import(A, *Graph_->OverlapImporter(), Insert));
//    EPETRA_CHK_ERR(OverlapA->FillComplete());
  }
  else {
    initAllValues(*A_);
  }

  ++numInitialize_;
}

//==========================================================================
template<class MatrixType>
void RILUK<MatrixType>::initAllValues(const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> & OverlapA) {

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

  bool ReplaceValues = (L_->isStaticGraph() || L_->isLocallyIndexed()); // Check if values should be inserted or replaced

  L_->resumeFill();
  U_->resumeFill();
  if (ReplaceValues) {
    L_->setAllToScalar(0.0); // Zero out L and U matrices
    U_->setAllToScalar(0.0);
  }

  D_->putScalar(0.0); // Set diagonal values to zero
  Teuchos::ArrayRCP<scalar_type> DV = D_->get1dViewNonConst(); // Get view of diagonal

  const Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> >& rowMap =
    L_->getRowMap();

  // First we copy the user's matrix into L and U, regardless of fill level

  for (global_ordinal_type i=rowMap->getMinGlobalIndex(); i<=rowMap->getMaxGlobalIndex(); ++i) {
    global_ordinal_type global_row = i;
    local_ordinal_type local_row = rowMap->getLocalElement(global_row);

    OverlapA.getGlobalRowCopy(global_row, InI(), InV(), NumIn); // Get Values and Indices

    // Split into L and U (we don't assume that indices are ordered).

    NumL = 0;
    NumU = 0;
    DiagFound = false;

    for (size_t j=0; j< NumIn; j++) {
      global_ordinal_type k = InI[j];

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

    if (DiagFound) ++NumNonzeroDiags;
    else DV[local_row] = Athresh_;

    if (NumL) {
      if (ReplaceValues) {
        L_->replaceGlobalValues(global_row, LI(0, NumL), LV(0,NumL));
      }
      else {
        L_->insertGlobalValues(global_row, LI(0,NumL), LV(0,NumL));
      }
    }

    if (NumU) {
      if (ReplaceValues) {
        U_->replaceGlobalValues(global_row, UI(0,NumU), UV(0,NumU));
      }
      else {
        U_->insertGlobalValues(global_row, UI(0,NumU), UV(0,NumU));
      }
    }

  }

  // The domain of L and the range of U are exactly their own row maps (there is no communication).
  // The domain of U and the range of L must be the same as those of the original matrix,
  // However if the original matrix is a VbrMatrix, these two latter maps are translation from
  // a block map to a point map.
  L_->fillComplete(L_->getColMap(), A_->getRangeMap());
  U_->fillComplete(A_->getDomainMap(), U_->getRowMap());

  // At this point L and U have the values of A in the structure of L and U, and diagonal vector D

  setInitialized(true);
  setFactored(false);

  size_t TotalNonzeroDiags = 0;
  Teuchos::reduceAll(*L_->getRowMap()->getComm(),Teuchos::REDUCE_SUM,
                     1,&NumNonzeroDiags,&TotalNonzeroDiags);
  NumMyDiagonals_ = NumNonzeroDiags;
  if (NumNonzeroDiags != U_->getNodeNumRows()) {
    throw std::runtime_error("Error in Ifpack2::RILUK::initAllValues, wrong number of diagonals.");
  }
}

//==========================================================================
template<class MatrixType>
void RILUK<MatrixType>::compute() {

  L_->resumeFill();
  U_->resumeFill();

  TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized(), std::runtime_error,
      "Ifpack2::RILUK::compute() ERROR: isInitialized() must be true.");
  TEUCHOS_TEST_FOR_EXCEPTION(isComputed() == true, std::runtime_error,
      "Ifpack2::RILUK::compute() ERROR: Can't have already computed factors.");

  // MinMachNum should be officially defined, for now pick something a little
  // bigger than IEEE underflow value

  scalar_type MinDiagonalValue = Teuchos::ScalarTraits<scalar_type>::rmin();
  scalar_type MaxDiagonalValue = Teuchos::ScalarTraits<scalar_type>::one()/MinDiagonalValue;

  size_t NumIn, NumL, NumU;

  // Get Maximun Row length
  size_t MaxNumEntries = L_->getNodeMaxNumRowEntries() + U_->getNodeMaxNumRowEntries() + 1;

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
  for (size_t j=0; j<num_cols; j++) colflag[j] = - 1;

  for(size_t i=0; i<L_->getNodeNumRows(); i++) {
    local_ordinal_type local_row = i;

 // Fill InV, InI with current row of L, D and U combined

    NumIn = MaxNumEntries;
    L_->getLocalRowCopy(local_row, InI(), InV(), NumL);

    InV[NumL] = DV[i]; // Put in diagonal
    InI[NumL] = local_row;

    U_->getLocalRowCopy(local_row, InI(NumL+1,MaxNumEntries-NumL-1), InV(NumL+1,MaxNumEntries-NumL-1), NumU);
    NumIn = NumL+NumU+1;

    // Set column flags
    for (size_t j=0; j<NumIn; j++) colflag[InI[j]] = j;

    scalar_type diagmod = 0.0; // Off-diagonal accumulator

    for (size_t jj=0; jj<NumL; jj++) {
      local_ordinal_type j = InI[jj];
      scalar_type multiplier = InV[jj]; // current_mults++;

      InV[jj] *= DV[j];

      U_->getLocalRowView(j, UUI, UUV); // View of row above
      NumUU = UUI.size();

      if (RelaxValue_==0.0) {
        for (size_t k=0; k<NumUU; k++) {
          int kk = colflag[UUI[k]];
          if (kk>-1) {
            InV[kk] -= multiplier*UUV[k];
            current_madds++;
          }
        }
      }
      else {
        for (size_t k=0; k<NumUU; k++) {
          int kk = colflag[UUI[k]];
          if (kk>-1) InV[kk] -= multiplier*UUV[k];
          else diagmod -= multiplier*UUV[k];
          current_madds++;
        }
      }
    }
    if (NumL) {
      L_->replaceLocalValues(local_row, InI(0,NumL), InV(0,NumL));  // Replace current row of L
    }

    DV[i] = InV[NumL]; // Extract Diagonal value

    if (RelaxValue_!=0.0) {
      DV[i] += RelaxValue_*diagmod; // Add off diagonal modifications
      // current_madds++;
    }

    if (Teuchos::ScalarTraits<scalar_type>::magnitude(DV[i]) > Teuchos::ScalarTraits<scalar_type>::magnitude(MaxDiagonalValue)) {
      if (Teuchos::ScalarTraits<scalar_type>::real(DV[i]) < 0) DV[i] = - MinDiagonalValue;
      else DV[i] = MinDiagonalValue;
    }
    else
      DV[i] = Teuchos::ScalarTraits<scalar_type>::one()/DV[i]; // Invert diagonal value

    for (size_t j=0; j<NumU; j++) InV[NumL+1+j] *= DV[i]; // Scale U by inverse of diagonal

    if (NumU) {
      U_->replaceLocalValues(local_row, InI(NumL+1,NumU), InV(NumL+1,NumU));  // Replace current row of L and U
    }

    // Reset column flags
    for (size_t j=0; j<NumIn; j++) colflag[InI[j]] = -1;
  }

  L_->fillComplete(L_->getColMap(), A_->getRangeMap());
  U_->fillComplete(A_->getDomainMap(), U_->getRowMap());

  // Validate that the L and U factors are actually lower and upper triangular

  if( !L_->isLowerTriangular() )
    throw std::runtime_error("Ifpack2::RILUK::compute() ERROR, L isn't lower triangular.");
  if( !U_->isUpperTriangular() )
    throw std::runtime_error("Ifpack2::RILUK::compute() ERROR, U isn't lower triangular.");

  // Add up flops

  double current_flops = 2 * current_madds;
  double total_flops = 0;

  // Get total madds across all PEs
  Teuchos::reduceAll(*L_->getRowMap()->getComm(),Teuchos::REDUCE_SUM,
                     1,&current_flops,&total_flops);

  // Now count the rest
  total_flops += (double) L_->getGlobalNumEntries(); // Accounts for multiplier above
  total_flops += (double) D_->getGlobalLength(); // Accounts for reciprocal of diagonal
  if (RelaxValue_!=0.0) total_flops += 2 * (double)D_->getGlobalLength(); // Accounts for relax update of diag

  //UpdateFlops(total_flops); // Update flop count

  setFactored(true);
  ++numCompute_;
}

//=============================================================================
template<class MatrixType>
void RILUK<MatrixType>::apply(
       const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
             Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
             Teuchos::ETransp mode, scalar_type alpha, scalar_type beta) const
{
  typedef Teuchos::ScalarTraits<scalar_type> STS;

  TEUCHOS_TEST_FOR_EXCEPTION(!isComputed(), std::runtime_error,
    "Ifpack2::RILUK::apply() ERROR, compute() hasn't been called yet.");

  TEUCHOS_TEST_FOR_EXCEPTION(
    alpha != STS::one (), 
    std::logic_error,
    "Ifpack2::RILUK::apply() does not currently allow alpha != 1.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    beta != STS::zero (), 
    std::logic_error,
    "Ifpack2::RILUK::apply() does not currently allow zero != 0.");

//
// This function finds Y such that
// LDU Y = X, or
// U(trans) D L(trans) Y = X
// for multiple RHS
//

  // First generate X and Y as needed for this function
  Teuchos::RCP<const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > X1;
  Teuchos::RCP<Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > Y1;
  generateXY(mode, X, Y, X1, Y1);

  scalar_type one = Teuchos::ScalarTraits<scalar_type>::one();
  scalar_type zero = Teuchos::ScalarTraits<scalar_type>::zero();

  if (mode == Teuchos::NO_TRANS) {

    L_->localSolve(*X1, *Y1,mode);
    Y1->elementWiseMultiply(one, *D_, *Y1, zero); // y = D*y (D_ has inverse of diagonal)
    U_->localSolve(*Y1, *Y1,mode); // Solve Uy = y
    if (isOverlapped_) {
      // Export computed Y values if needed
      Y.doExport(*Y1,*L_->getGraph()->getExporter(), OverlapMode_);
    }
  }
  else {
    U_->localSolve(*X1, *Y1,mode); // Solve Uy = y
    Y1->elementWiseMultiply(one, *D_, *Y1, zero); // y = D*y (D_ has inverse of diagonal)
    L_->localSolve(*Y1, *Y1,mode);
    if (isOverlapped_) {Y.doExport(*Y1,*U_->getGraph()->getImporter(), OverlapMode_);} // Export computed Y values if needed
  }

  ++numApply_;
}

//=============================================================================
template<class MatrixType>
int RILUK<MatrixType>::Multiply(const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                              Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
            Teuchos::ETransp mode) const {
//
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
//

  // First generate X and Y as needed for this function
  Teuchos::RCP<const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > X1;
  Teuchos::RCP<Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > Y1;
  generateXY(mode, X, Y, X1, Y1);

//  Epetra_Flops * counter = this->GetFlopCounter();
//  if (counter!=0) {
//    L_->SetFlopCounter(*counter);
//    Y1->SetFlopCounter(*counter);
//    U_->SetFlopCounter(*counter);
//  }

  if (!mode == Teuchos::NO_TRANS) {
    U_->apply(*X1, *Y1,mode); //
    Y1->update(1.0, *X1, 1.0); // Y1 = Y1 + X1 (account for implicit unit diagonal)
    Y1->elementWiseMultiply(1.0, *D_, *Y1, 0.0); // y = D*y (D_ has inverse of diagonal)
    Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> Y1temp(*Y1); // Need a temp copy of Y1
    L_->apply(Y1temp, *Y1,mode);
    Y1->update(1.0, Y1temp, 1.0); // (account for implicit unit diagonal)
    if (isOverlapped_) {Y.doExport(*Y1,*L_->getGraph()->getExporter(), OverlapMode_);} // Export computed Y values if needed
  }
  else {

    L_->apply(*X1, *Y1,mode);
    Y1->update(1.0, *X1, 1.0); // Y1 = Y1 + X1 (account for implicit unit diagonal)
    Y1->elementWiseMultiply(1, *D_, *Y1, 0); // y = D*y (D_ has inverse of diagonal)
    Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> Y1temp(*Y1); // Need a temp copy of Y1
    U_->apply(Y1temp, *Y1,mode);
    Y1->update(1.0, Y1temp, 1.0); // (account for implicit unit diagonal)
    if (isOverlapped_) {Y.doExport(*Y1,*L_->getGraph()->getExporter(), OverlapMode_);}
  }
  return(0);
}

//=============================================================================
template<class MatrixType>
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType
RILUK<MatrixType>::computeCondEst(Teuchos::ETransp mode) const {

  if (Condest_>=0.0) {
    return Condest_;
  }
  // Create a vector with all values equal to one
  Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> Ones(U_->getDomainMap());
  Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> OnesResult(L_->getRangeMap());
  Ones.putScalar(1.0);

  apply(Ones, OnesResult,mode); // Compute the effect of the solve on the vector of ones
  OnesResult.abs(OnesResult); // Make all values non-negative
  Teuchos::Array<magnitude_type> norms(1);
  OnesResult.normInf(norms());
  Condest_ = norms[0];
  return Condest_;
}


//=========================================================================
template<class MatrixType>
void RILUK<MatrixType>::generateXY(Teuchos::ETransp mode,
    const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Xin,
    const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Yin,
    Teuchos::RCP<const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >& Xout,
    Teuchos::RCP<Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >& Yout) const {

  // Generate an X and Y suitable for performing Solve() and Multiply() methods

  TEUCHOS_TEST_FOR_EXCEPTION(Xin.getNumVectors()!=Yin.getNumVectors(), std::runtime_error,
       "Ifpack2::RILUK::GenerateXY ERROR: X and Y not the same size");

  //cout << "Xin = " << Xin << endl;
  Xout = Teuchos::rcp( (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> *) &Xin, false );
  Yout = Teuchos::rcp( (Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> *) &Yin, false );
  if (!isOverlapped_) return; // Nothing more to do

  if (isOverlapped_) {
    // Make sure the number of vectors in the multivector is the same as before.
    if (OverlapX_!=Teuchos::null) {
      if (OverlapX_->getNumVectors()!=Xin.getNumVectors()) {
        OverlapX_ = Teuchos::null;
        OverlapY_ = Teuchos::null;
      }
    }
    if (OverlapX_==Teuchos::null) { // Need to allocate space for overlap X and Y
      OverlapX_ = Teuchos::rcp( new Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>(U_->getColMap(), Xout->getNumVectors()) );
      OverlapY_ = Teuchos::rcp( new Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>(L_->getRowMap(), Yout->getNumVectors()) );
    }
    if (mode == Teuchos::NO_TRANS) {
      OverlapX_->doImport(*Xout,*U_->getGraph()->getImporter(), Tpetra::INSERT); // Import X values for solve
    }
    else {
      OverlapX_->doImport(*Xout,*L_->getGraph()->getExporter(), Tpetra::INSERT); // Import X values for solve
    }
    Xout = OverlapX_;
    Yout = OverlapY_; // Set pointers for Xout and Yout to point to overlap space
    //cout << "OverlapX_ = " << *OverlapX_ << endl;
  }
}

}//namespace Ifpack2

#endif

