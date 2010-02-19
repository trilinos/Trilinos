//@HEADER
// ***********************************************************************
// 
//       Tifpack: Templated Object-Oriented Algebraic Preconditioner Package
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

#ifndef TIFPACK_CRSRILUK_DEF_HPP
#define TIFPACK_CRSRILUK_DEF_HPP

namespace Tifpack {

//==============================================================================
template<class MatrixType>
CrsRiluk<MatrixType>::CrsRiluk(const Teuchos::RCP<const MatrixType>& Matrix_in) 
  : UserMatrixIsCrs_(false),
    Graph_(),
    A_(Matrix_in),
    UseTranspose_(false),
    LevelOfFill_(0),
    LevelOfOverlap_(0),
    NumMyDiagonals_(0),
    Allocated_(false),
    isInitialized_(false),
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
//CrsRiluk<MatrixType>::CrsRiluk(const CrsRiluk<MatrixType>& src) 
//  : UserMatrixIsCrs_(src.UserMatrixIsCrs_),
//    IsOverlapped_(src.IsOverlapped_),
//    Graph_(src.Graph_),
//    UseTranspose_(src.UseTranspose_),
//    LevelOfFill_(src.LevelOfFill_),
//    LevelOfOverlap_(src.LevelOfOverlap_),
//    NumMyDiagonals_(src.NumMyDiagonals_),
//    Allocated_(src.Allocated_),
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
//  D_ = Teuchos::rcp( new Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(src.getD()) );
//}

//==============================================================================
template<class MatrixType>
CrsRiluk<MatrixType>::~CrsRiluk() {
}

//==============================================================================
template<class MatrixType>
void CrsRiluk<MatrixType>::AllocateCrs() {

  // Allocate Matrix using ILUK graphs
  L_ = Teuchos::rcp( new MatrixType(Graph_->getL_Graph()) );
  U_ = Teuchos::rcp( new MatrixType(Graph_->getU_Graph()) );
  D_ = Teuchos::rcp( new Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(Graph_->getL_Graph()->getRowMap()) );
  SetAllocated(true);
}

//==========================================================================
template<class MatrixType>
void CrsRiluk<MatrixType>::setParameters(const Teuchos::ParameterList& parameterlist) {
  Tifpack::GetParameter(parameterlist, "fact: iluk level-of-fill", LevelOfFill_);
  Tifpack::GetParameter(parameterlist, "fact: iluk level-of-overlap", LevelOfOverlap_);
  double tmp = -1;
  Tifpack::GetParameter(parameterlist, "fact: absolute threshold", tmp);
  if (tmp != -1) Athresh_ = tmp;
  tmp = -1;
  Tifpack::GetParameter(parameterlist, "fact: relative threshold", tmp);
  if (tmp != -1) Rthresh_ = tmp;
  tmp = -1;
  Tifpack::GetParameter(parameterlist, "fact: relax value", tmp);
  if (tmp != -1) RelaxValue_ = tmp;
}

//==========================================================================
template<class MatrixType>
void CrsRiluk<MatrixType>::initialize() {

  if (Graph_ != Teuchos::null) return;

  Graph_ = Teuchos::rcp(new Tifpack::IlukGraph<LocalOrdinal,GlobalOrdinal,Node>(A_->getCrsGraph(), LevelOfFill_, LevelOfOverlap_));

  Graph_->constructFilledGraph();

  setInitialized(true);

  if (!Allocated()) AllocateCrs();

//  Teuchos::RCP<MatrixType> OverlapA = Teuchos::rcp( (Epetra_CrsMatrix *) &A, false );

  if (IsOverlapped_) {
  
//    OverlapA = Teuchos::rcp( new MatrixType(Graph_->OverlapGraph()) );
//    EPETRA_CHK_ERR(OverlapA->Import(A, *Graph_->OverlapImporter(), Insert));
//    EPETRA_CHK_ERR(OverlapA->FillComplete());
  }
  
  // Get Maximum Row length
  //int MaxNumEntries = OverlapA->MaxNumEntries();

  // Do the rest using generic Epetra_RowMatrix interface

//  EPETRA_CHK_ERR(InitAllValues(*OverlapA, MaxNumEntries));

}

//==========================================================================
template<class MatrixType>
bool CrsRiluk<MatrixType>::isInitialized() const {
  return isInitialized_;
}

//==========================================================================
template<class MatrixType>
void CrsRiluk<MatrixType>::InitAllValues(const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> & OverlapA, int MaxNumEntries) {

  size_t NumIn, NumL, NumU;
  bool DiagFound;
  size_t NumNonzeroDiags = 0;


  Teuchos::Array<LocalOrdinal> InI(MaxNumEntries); // Allocate temp space
  Teuchos::Array<LocalOrdinal> LI(MaxNumEntries);
  Teuchos::Array<LocalOrdinal> UI(MaxNumEntries);
  Teuchos::Array<Scalar> InV(MaxNumEntries);
  Teuchos::Array<Scalar> LV(MaxNumEntries);
  Teuchos::Array<Scalar> UV(MaxNumEntries);

  bool ReplaceValues = false ; //(L_->isStaticGraph() || L_->isLocallyIndexed()); // Check if values should be inserted or replaced

  if (ReplaceValues) {
    L_->setAllToScalar(0.0); // Zero out L and U matrices
    U_->setAllToScalar(0.0);
  }

  D_->putScalar(0.0); // Set diagonal values to zero
  Teuchos::ArrayRCP<Scalar> DV = D_->get1dViewNonConst(); // Get view of diagonal
    

  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap =
    L_->getRowMap();

  // First we copy the user's matrix into L and U, regardless of fill level

  for (LocalOrdinal i=rowMap->getMinLocalIndex(); i<=rowMap->getMaxLocalIndex(); ++i) {
    LocalOrdinal local_row = i;

    OverlapA.getLocalRowCopy(local_row, InI(), InV(), NumIn); // Get Values and Indices
    
    // Split into L and U (we don't assume that indices are ordered).
    
    NumL = 0; 
    NumU = 0; 
    DiagFound = false;
    
    for (size_t j=0; j< NumIn; j++) {
      LocalOrdinal k = InI[j];

      if (k==i) {
        DiagFound = true;
        DV[i] += Rthresh_ * InV[j] + TIFPACK_SGN(InV[j]) * Athresh_; // Store perturbed diagonal in Tpetra::Vector D_
      }

      else if (k < 0) { // Out of range
        throw std::runtime_error("out of range in Tifpack::CrsRiluk::InitAllValues");
      }

      else if (k < i) {
        LI[NumL] = k;
        LV[NumL] = InV[j];
        NumL++;
      }
      else if (k<rowMap->getMaxLocalIndex()) {
        UI[NumU] = k;
        UV[NumU] = InV[j];
        NumU++;
      }
    }
    
    // Check in things for this row of L and U

    if (DiagFound) NumNonzeroDiags++;
    else DV[i] = Athresh_;

    if (NumL) {
//      if (ReplaceValues) {
//        //replaceGlobalValues works with Global row and Local column-indices:
//        GlobalOrdinal global_row = rowMap->getGlobalElement(local_row);
//        L_->replaceGlobalValues(global_row, LI(0, NumL), LV(0,NumL));
//      }
//      else {
        L_->insertLocalValues(local_row, LI(0,NumL), LV(0,NumL));
//      }
    }

    if (NumU) {
//      if (ReplaceValues) {
//        //replaceGlobalValues works with Global row and Local column-indices:
//        GlobalOrdinal global_row = rowMap->getGlobalElement(local_row);
//        U_->replaceGlobalValues(global_row, UI(0,NumU), UV(0,NumU));
//      }
//      else {
        U_->insertLocalValues(local_row, UI(0,NumU), UV(0,NumU));
//      }
    }

  }

  if (!ReplaceValues) {
    // The domain of L and the range of U are exactly their own row maps (there is no communication).
    // The domain of U and the range of L must be the same as those of the original matrix,
    // However if the original matrix is a VbrMatrix, these two latter maps are translation from
    // a block map to a point map.
    L_->fillComplete(L_->getColMap(), A_->getRangeMap());
    U_->fillComplete(A_->getDomainMap(), U_->getRowMap());
  }

  // At this point L and U have the values of A in the structure of L and U, and diagonal vector D

  setInitialized(true);
  SetFactored(false);

  size_t TotalNonzeroDiags = 0;
  Teuchos::reduceAll(*L_->getRowMap()->getComm(),Teuchos::REDUCE_SUM,
                     1,&NumNonzeroDiags,&TotalNonzeroDiags);
  NumMyDiagonals_ = NumNonzeroDiags;
  if (NumNonzeroDiags != U_->getNodeNumRows()) {
    throw std::runtime_error("Error in Tifpack::CrsRiluk::InitAllValues, wrong number of diagonals.");
  }
}

//==========================================================================
template<class MatrixType>
void CrsRiluk<MatrixType>::compute() {

  TEST_FOR_EXCEPTION(!isInitialized(), std::runtime_error,
      "Tifpack::CrsRiluk::compute() ERROR: isInitialized() must be true.");
  TEST_FOR_EXCEPTION(isComputed() == true, std::runtime_error,
      "Tifpack::CrsRiluk::compute() ERROR: Can't have already computed factors.");

  setInitialized(false);

  // MinMachNum should be officially defined, for now pick something a little 
  // bigger than IEEE underflow value

  Scalar MinDiagonalValue = Teuchos::ScalarTraits<Scalar>::rmin();
  Scalar MaxDiagonalValue = Teuchos::ScalarTraits<Scalar>::one()/MinDiagonalValue;

  size_t NumIn, NumL, NumU;

  // Get Maximun Row length
  size_t MaxNumEntries = L_->getNodeMaxNumRowEntries() + U_->getNodeMaxNumRowEntries() + 1;

  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap =
    L_->getRowMap();

  Teuchos::Array<GlobalOrdinal> InI(MaxNumEntries); // Allocate temp space
  Teuchos::Array<Scalar> InV(MaxNumEntries);
  Teuchos::Array<int> colflag(MaxNumEntries);

  Teuchos::ArrayRCP<Scalar> DV = D_->get1dViewNonConst(); // Get view of diagonal

  int current_madds = 0; // We will count multiply-add as they happen

  // Now start the factorization.

  // Need some integer workspace and pointers
  size_t NumUU; 
  Teuchos::ArrayRCP<const LocalOrdinal> UUI;
  Teuchos::ArrayRCP<const Scalar> UUV;
  for (size_t j=0; j<MaxNumEntries; j++) colflag[j] = - 1;

  for(size_t i=0; i<L_->getNodeNumRows(); i++) {
    LocalOrdinal local_row = i;

 // Fill InV, InI with current row of L, D and U combined

    NumIn = MaxNumEntries;
    L_->getGlobalRowCopy(local_row, InI(), InV(), NumL);

    InV[NumL] = DV[i]; // Put in diagonal
    InI[NumL] = local_row;
    
    U_->getGlobalRowCopy(local_row, InI(NumL+1,MaxNumEntries-NumL-1), InV(NumL+1,MaxNumEntries-NumL-1), NumU);
    NumIn = NumL+NumU+1;

    // Set column flags
    for (size_t j=0; j<NumIn; j++) colflag[InI[j]] = j;

    Scalar diagmod = 0.0; // Off-diagonal accumulator

    for (size_t jj=0; jj<NumL; jj++) {
      LocalOrdinal j = InI[jj];
      Scalar multiplier = InV[jj]; // current_mults++;

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
      GlobalOrdinal global_row = rowMap->getGlobalElement(i);
      L_->replaceGlobalValues(global_row, InI(0,NumL), InV(0,NumL));  // Replace current row of L
    }

    DV[i] = InV[NumL]; // Extract Diagonal value

    if (RelaxValue_!=0.0) {
      DV[i] += RelaxValue_*diagmod; // Add off diagonal modifications
      // current_madds++;
    }

    if (Teuchos::ScalarTraits<Scalar>::magnitude(DV[i]) > MaxDiagonalValue) {
      if (DV[i] < 0) DV[i] = - MinDiagonalValue;
      else DV[i] = MinDiagonalValue;
    }
    else
      DV[i] = 1.0/DV[i]; // Invert diagonal value

    for (size_t j=0; j<NumU; j++) InV[NumL+1+j] *= DV[i]; // Scale U by inverse of diagonal

    if (NumU) {
      GlobalOrdinal global_row = rowMap->getGlobalElement(i);
      U_->replaceGlobalValues(global_row, InI(NumL+1,NumU), InV(NumL+1,NumU));  // Replace current row of L and U
    }

    // Reset column flags
    for (size_t j=0; j<NumIn; j++) colflag[InI[j]] = -1;
  }

  // Validate that the L and U factors are actually lower and upper triangular

  if( !L_->isLowerTriangular() ) 
    throw std::runtime_error("Tifpack::CrsRiluk::compute() ERROR, L isn't lower triangular.");
  if( !U_->isUpperTriangular() ) 
    throw std::runtime_error("Tifpack::CrsRiluk::compute() ERROR, U isn't lower triangular.");
  
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

  SetFactored(true);
}

//=============================================================================
template<class MatrixType>
void CrsRiluk<MatrixType>::apply(
       const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
             Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
             Teuchos::ETransp mode, Scalar alpha, Scalar beta) const {
//
// This function finds Y such that
// LDU Y = X, or
// U(trans) D L(trans) Y = X
// for multiple RHS
//

  // First generate X and Y as needed for this function
  Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > X1;
  Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Y1;
  generateXY(mode, X, Y, X1, Y1);

//  Epetra_Flops * counter = this->GetFlopCounter();
//  if (counter!=0) {
//    L_->SetFlopCounter(*counter);
//    Y1->SetFlopCounter(*counter);
//    U_->SetFlopCounter(*counter);
//  }

  Scalar one = Teuchos::ScalarTraits<Scalar>::one();
  Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();

  if (!mode == Teuchos::NO_TRANS) {

    L_->solve(*X1, *Y1,mode);
    Y1->elementWiseMultiply(one, *D_, *Y1, zero); // y = D*y (D_ has inverse of diagonal)
    U_->solve(*Y1, *Y1,mode); // Solve Uy = y
    if (IsOverlapped_) {Y.doExport(*Y1,*L_->getGraph()->getExporter(), OverlapMode_);} // Export computed Y values if needed
  }
  else {
    U_->solve(*X1, *Y1,mode); // Solve Uy = y
    Y1->elementWiseMultiply(one, *D_, *Y1, zero); // y = D*y (D_ has inverse of diagonal)
    L_->solve(*Y1, *Y1,mode);
    if (IsOverlapped_) {Y.doExport(*Y1,*U_->getGraph()->getImporter(), OverlapMode_);} // Export computed Y values if needed
  } 
}

//=============================================================================
template<class MatrixType>
int CrsRiluk<MatrixType>::Multiply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
			      Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
            Teuchos::ETransp mode) const {
//
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
//
    
  // First generate X and Y as needed for this function
  Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > X1;
  Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Y1;
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
    Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Y1temp(*Y1); // Need a temp copy of Y1
    L_->apply(Y1temp, *Y1,mode);
    Y1->update(1.0, Y1temp, 1.0); // (account for implicit unit diagonal)
    if (IsOverlapped_) {Y.doExport(*Y1,*L_->getGraph()->getExporter(), OverlapMode_);} // Export computed Y values if needed
  }
  else {

    L_->apply(*X1, *Y1,mode);
    Y1->update(1.0, *X1, 1.0); // Y1 = Y1 + X1 (account for implicit unit diagonal)
    Y1->elementWiseMultiply(1, *D_, *Y1, 0); // y = D*y (D_ has inverse of diagonal)
    Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Y1temp(*Y1); // Need a temp copy of Y1
    U_->apply(Y1temp, *Y1,mode);
    Y1->update(1.0, Y1temp, 1.0); // (account for implicit unit diagonal)
    if (IsOverlapped_) {Y.doExport(*Y1,*L_->getGraph()->getExporter(), OverlapMode_);}
  } 
  return(0);
}

//=============================================================================
template<class MatrixType>
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType
CrsRiluk<MatrixType>::computeCondEst(Teuchos::ETransp mode) const {

  if (Condest_>=0.0) {
    return Condest_;
  }
  // Create a vector with all values equal to one
  Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> Ones(U_->getDomainMap());
  Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> OnesResult(L_->getRangeMap());
  Ones.putScalar(1.0);

  apply(Ones, OnesResult,mode); // Compute the effect of the solve on the vector of ones
  OnesResult.abs(OnesResult); // Make all values non-negative
  Condest_ = OnesResult.normInf(); // Get the maximum value across all processors
  return Condest_;
}


//=========================================================================
template<class MatrixType>
void CrsRiluk<MatrixType>::generateXY(Teuchos::ETransp mode, 
    const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Xin,
    const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Yin,
    Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Xout, 
    Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Yout) const {

  // Generate an X and Y suitable for performing Solve() and Multiply() methods

  TEST_FOR_EXCEPTION(Xin.getNumVectors()!=Yin.getNumVectors(), std::runtime_error,
       "Tifpack::CrsRiluk::GenerateXY ERROR: X and Y not the same size");

  //cout << "Xin = " << Xin << endl;
  Xout = Teuchos::rcp( (Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> *) &Xin, false );
  Yout = Teuchos::rcp( (Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> *) &Yin, false );
  if (!IsOverlapped_ && UserMatrixIsCrs_) return; // Nothing more to do

  if (IsOverlapped_) {
    // Make sure the number of vectors in the multivector is the same as before.
    if (OverlapX_!=Teuchos::null) {
      if (OverlapX_->getNumVectors()!=Xin.getNumVectors()) {
        OverlapX_ = Teuchos::null;
        OverlapY_ = Teuchos::null;
      }
    }
    if (OverlapX_==Teuchos::null) { // Need to allocate space for overlap X and Y
      OverlapX_ = Teuchos::rcp( new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(U_->getColMap(), Xout->getNumVectors()) );
      OverlapY_ = Teuchos::rcp( new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(L_->getRowMap(), Yout->getNumVectors()) );
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

}//namespace Tifpack

#endif

