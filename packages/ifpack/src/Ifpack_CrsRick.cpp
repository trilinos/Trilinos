/*@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
*/

#include "Ifpack_CrsRick.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

#ifdef HAVE_IFPACK_TEUCHOS
#include <Teuchos_ParameterList.hpp>
#include <ifp_parameters.h>
#endif

//==============================================================================
Ifpack_CrsRick::Ifpack_CrsRick(const Epetra_CrsMatrix &A, const Ifpack_IlukGraph & Graph) 
  : A_(A),
    Graph_(Graph),
    UseTranspose_(false),
    Allocated_(false),
    ValuesInitialized_(false),
    Factored_(false),
    RelaxValue_(0.0),
    Condest_(-1.0),
    Athresh_(0.0),
    Rthresh_(1.0),
    OverlapX_(0),
    OverlapY_(0),
    OverlapMode_(Zero)
{
  int ierr = Allocate();
}

//==============================================================================
Ifpack_CrsRick::Ifpack_CrsRick(const Ifpack_CrsRick & FactoredMatrix) 
  : A_(FactoredMatrix.A_),
    Graph_(FactoredMatrix.Graph_),
    UseTranspose_(FactoredMatrix.UseTranspose_),
    Allocated_(FactoredMatrix.Allocated_),
    ValuesInitialized_(FactoredMatrix.ValuesInitialized_),
    Factored_(FactoredMatrix.Factored_),
    RelaxValue_(FactoredMatrix.RelaxValue_),
    Condest_(FactoredMatrix.Condest_),
    Athresh_(FactoredMatrix.Athresh_),
    Rthresh_(FactoredMatrix.Rthresh_),
    OverlapX_(0),
    OverlapY_(0),
    OverlapMode_(FactoredMatrix.OverlapMode_)
{
  U_ = new Epetra_CrsMatrix(FactoredMatrix.U());
  D_ = new Epetra_Vector(Graph_.L_Graph().RowMap());
  
}

//==============================================================================
int Ifpack_CrsRick::Allocate() {

  // Allocate Epetra_CrsMatrix using ILUK graphs
  U_ = new Epetra_CrsMatrix(Copy, Graph_.U_Graph());
  D_ = new Epetra_Vector(Graph_.L_Graph().RowMap());
  
  
    SetAllocated(true);
    return(0);
}
//==============================================================================
Ifpack_CrsRick::~Ifpack_CrsRick(){


  delete U_;
  delete D_; // Diagonal is stored separately.  We store the inverse.

  if (OverlapX_!=0) delete OverlapX_;
  if (OverlapY_!=0) delete OverlapY_;

  ValuesInitialized_ = false;
  Factored_ = false;
  Allocated_ = false;
}

#ifdef HAVE_IFPACK_TEUCHOS
//==========================================================================
int Ifpack_CrsRick::SetParameters(const Teuchos::ParameterList& parameterlist,
				  bool cerr_warning_if_unused)
{
  Ifpack::param_struct params;
  params.double_params[Ifpack::relax_value] = RelaxValue_;
  params.double_params[Ifpack::absolute_threshold] = Athresh_;
  params.double_params[Ifpack::relative_threshold] = Rthresh_;
  params.overlap_mode = OverlapMode_;

  Ifpack::set_parameters(parameterlist, params, cerr_warning_if_unused);

  RelaxValue_ = params.double_params[Ifpack::relax_value];
  Athresh_ = params.double_params[Ifpack::absolute_threshold];
  Rthresh_ = params.double_params[Ifpack::relative_threshold];
  OverlapMode_ = params.overlap_mode;

  return(0);
}
#endif

//==========================================================================
int Ifpack_CrsRick::InitValues() {

  // if (!Allocated()) return(-1); // This test is not needed at this time.  All constructors allocate.

  int ierr = 0;
  int i, j;
  int * InI=0, * LI=0, * UI = 0;
  double * InV=0, * LV=0, * UV = 0;
  int NumIn, NumL, NumU;
  bool DiagFound;
  int NumNonzeroDiags = 0;

  Epetra_CrsMatrix * OverlapA = (Epetra_CrsMatrix *) &A_;

  if (Graph_.LevelOverlap()>0 && Graph_.L_Graph().DomainMap().DistributedGlobal()) {
  
  OverlapA = new Epetra_CrsMatrix(Copy, *Graph_.OverlapGraph());
  OverlapA->Import(A_, *Graph_.OverlapImporter(), Insert);
  }

  // Get Maximun Row length
  int MaxNumEntries = OverlapA->MaxNumEntries();

  InI = new int[MaxNumEntries]; // Allocate temp space
  UI = new int[MaxNumEntries];
  InV = new double[MaxNumEntries];
  UV = new double[MaxNumEntries];

  double *DV;
  ierr = D_->ExtractView(&DV); // Get view of diagonal
    

  // First we copy the user's matrix into diagonal vector and U, regardless of fill level

  for (i=0; i< NumMyRows(); i++) {

    OverlapA->ExtractMyRowCopy(i, MaxNumEntries, NumIn, InV, InI); // Get Values and Indices
    
    // Split into L and U (we don't assume that indices are ordered).
    
    NumL = 0; 
    NumU = 0; 
    DiagFound = false;
    
    for (j=0; j< NumIn; j++) {
      int k = InI[j];

      if (k==i) {
	DiagFound = true;
	DV[i] += Rthresh_ * InV[j] + EPETRA_SGN(InV[j]) * Athresh_; // Store perturbed diagonal in Epetra_Vector D_
      }

      else if (k < 0) return(-1); // Out of range
      else if (k<NumMyRows()) {
	UI[NumU] = k;
	UV[NumU] = InV[j];
	NumU++;
      }
    }
    
    // Check in things for this row of L and U

    if (DiagFound) NumNonzeroDiags++;
    if (NumU) U_->ReplaceMyValues(i, NumU, UV, UI);
    
  }

  delete [] UI;
  delete [] UV;
  delete [] InI;
  delete [] InV;

  if (Graph_.LevelOverlap()>0 && Graph_.L_Graph().DomainMap().DistributedGlobal()) delete OverlapA;


  U_->TransformToLocal();

  // At this point L and U have the values of A in the structure of L and U, and diagonal vector D

  SetValuesInitialized(true);
  SetFactored(false);

  int TotalNonzeroDiags = 0;
  Graph_.U_Graph().RowMap().Comm().SumAll(&NumNonzeroDiags, &TotalNonzeroDiags, 1);
  if (Graph_.LevelOverlap()==0 &&
      ((TotalNonzeroDiags!=NumGlobalRows()) || 
       (TotalNonzeroDiags!=NumGlobalDiagonals()))) ierr = 1;
  if (NumNonzeroDiags != NumMyDiagonals()) ierr = 1; // Diagonals are not right, warn user

  return(ierr);
}

//==========================================================================
int Ifpack_CrsRick::Factor() {

  // if (!Allocated()) return(-1); // This test is not needed at this time.  All constructors allocate.
  if (!ValuesInitialized()) return(-2); // Must have values initialized.
  if (Factored()) return(-3); // Can't have already computed factors.

  SetValuesInitialized(false);

  // MinMachNum should be officially defined, for now pick something a little 
  // bigger than IEEE underflow value

  double MinDiagonalValue = Epetra_MinDouble;
  double MaxDiagonalValue = 1.0/MinDiagonalValue;

  int ierr = 0;
  int i, j, k;
  int * UI = 0;
  double * UV = 0;
  int NumIn, NumU;

  // Get Maximun Row length
  int MaxNumEntries = U_->MaxNumEntries() + 1;

  int * InI = new int[MaxNumEntries]; // Allocate temp space
  double * InV = new double[MaxNumEntries];
  int * colflag = new int[NumMyCols()];

  double *DV;
  ierr = D_->ExtractView(&DV); // Get view of diagonal

  int current_madds = 0; // We will count multiply-add as they happen

  // Now start the factorization.

  // Need some integer workspace and pointers
  int NumUU; 
  int * UUI;
  double * UUV;
  for (j=0; j<NumMyCols(); j++) colflag[j] = - 1;

  for(i=0; i<NumMyRows(); i++) {

 // Fill InV, InI with current row of L, D and U combined

    NumIn = MaxNumEntries;
    assert(L_->ExtractMyRowCopy(i, NumIn, NumL, InV, InI)==0);
    LV = InV;
    LI = InI;

    InV[NumL] = DV[i]; // Put in diagonal
    InI[NumL] = i;
    
    assert(U_->ExtractMyRowCopy(i, NumIn-NumL-1, NumU, InV+NumL+1, InI+NumL+1)==0);
    NumIn = NumL+NumU+1;
    UV = InV+NumL+1;
    UI = InI+NumL+1;

    // Set column flags
    for (j=0; j<NumIn; j++) colflag[InI[j]] = j;

    double diagmod = 0.0; // Off-diagonal accumulator

    for (int jj=0; jj<NumL; jj++) {
      j = InI[jj];
      double multiplier = InV[jj]; // current_mults++;

      InV[jj] *= DV[j];
      
      assert(U_->ExtractMyRowView(j, NumUU, UUV, UUI)==0); // View of row above

      if (RelaxValue_==0.0) {
	for (k=0; k<NumUU; k++) {
	  int kk = colflag[UUI[k]];
	  if (kk>-1) {
	    InV[kk] -= multiplier*UUV[k];
	    current_madds++;
	  }
	}
      }
      else {
	for (k=0; k<NumUU; k++) {
	  int kk = colflag[UUI[k]];
	  if (kk>-1) InV[kk] -= multiplier*UUV[k];
	  else diagmod -= multiplier*UUV[k];
	  current_madds++;
	}
      }
     }
    if (NumL) assert(L_->ReplaceMyValues(i, NumL, LV, LI)==0);  // Replace current row of L

    DV[i] = InV[NumL]; // Extract Diagonal value

    if (RelaxValue_!=0.0) {
      DV[i] += RelaxValue_*diagmod; // Add off diagonal modifications
      // current_madds++;
    }

    if (fabs(DV[i]) > MaxDiagonalValue) {
      if (DV[i] < 0) DV[i] = - MinDiagonalValue;
      else DV[i] = MinDiagonalValue;
    }
    else
      DV[i] = 1.0/DV[i]; // Invert diagonal value

    for (j=0; j<NumU; j++) UV[j] *= DV[i]; // Scale U by inverse of diagonal

    if (NumU) assert(U_->ReplaceMyValues(i, NumU, UV, UI)==0);  // Replace current row of L and U


    // Reset column flags
    for (j=0; j<NumIn; j++) colflag[InI[j]] = -1;
  }

  
  // Add up flops
 
  double current_flops = 2 * current_madds;
  double total_flops = 0;
    
  Graph_.L_Graph().RowMap().Comm().SumAll(&current_flops, &total_flops, 1); // Get total madds across all PEs

  // Now count the rest
  total_flops += (double) L_->NumGlobalNonzeros(); // Accounts for multiplier above
  total_flops += (double) D_->GlobalLength(); // Accounts for reciprocal of diagonal
  if (RelaxValue_!=0.0) total_flops += 2 * (double)D_->GlobalLength(); // Accounts for relax update of diag

  UpdateFlops(total_flops); // Update flop count

  delete [] InI;
  delete [] InV;
  delete [] colflag;
  
  SetFactored(true);

  return(ierr);

}

//=============================================================================
int Ifpack_CrsRick::Solve(bool Trans, const Epetra_Vector& X, 
				Epetra_Vector& Y) const {
//
// This function finds Y such that LDU Y = X or U(trans) D L(trans) Y = X for a single RHS
//

  bool Upper = true;
  bool Lower = false;
  bool UnitDiagonal = true;

  Epetra_Vector * X1 = (Epetra_Vector *) &X;
  Epetra_Vector * Y1 = (Epetra_Vector *) &Y;

  if (Graph_.LevelOverlap()>0 && Graph_.L_Graph().DomainMap().DistributedGlobal()) {
    if (OverlapX_==0) { // Need to allocate space for overlap X and Y
      OverlapX_ = new Epetra_Vector(Graph_.OverlapGraph()->RowMap());
      OverlapY_ = new Epetra_Vector(Graph_.OverlapGraph()->RowMap());
    }
    OverlapX_->Import(X,*Graph_.OverlapImporter(), Insert); // Import X values for solve
    X1 = (Epetra_Vector *) OverlapX_;
    Y1 = (Epetra_Vector *) OverlapY_; // Set pointers for X1 and Y1 to point to overlap space
  }

  Epetra_Flops * counter = this->GetFlopCounter();
  if (counter!=0) {
    L_->SetFlopCounter(*counter);
    Y1->SetFlopCounter(*counter);
    U_->SetFlopCounter(*counter);
  }

  if (!Trans) {

    L_->Solve(Lower, Trans, UnitDiagonal, *X1, *Y1);
    Y1->Multiply(1.0, *D_, *Y1, 0.0); // y = D*y (D_ has inverse of diagonal)
    U_->Solve(Upper, Trans, UnitDiagonal, *Y1, *Y1); // Solve Uy = y
  }
  else
    {
      U_->Solve(Upper, Trans, UnitDiagonal, *X1, *Y1); // Solve Uy = y
      Y1->Multiply(1.0, *D_, *Y1, 0.0); // y = D*y (D_ has inverse of diagonal)
      L_->Solve(Lower, Trans, UnitDiagonal, *Y1, *Y1);
      
    } 

  // Export computed Y values as directed
  if (Graph_.LevelOverlap()>0 && Graph_.L_Graph().DomainMap().DistributedGlobal())
    Y.Export(*OverlapY_,*Graph_.OverlapImporter(), OverlapMode_);
  return(0);
}


//=============================================================================
int Ifpack_CrsRick::Solve(bool Trans, const Epetra_MultiVector& X, 
				Epetra_MultiVector& Y) const {
//
// This function finds Y such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
//

  if (X.NumVectors()!=Y.NumVectors()) EPETRA_CHK_ERR(-1); // Return error: X and Y not the same size

  bool Upper = true;
  bool Lower = false;
  bool UnitDiagonal = true;

  Epetra_MultiVector * X1 = (Epetra_MultiVector *) &X;
  Epetra_MultiVector * Y1 = (Epetra_MultiVector *) &Y;

  if (Graph_.LevelOverlap()>0 && Graph_.L_Graph().DomainMap().DistributedGlobal()) {
    // Make sure the number of vectors in the multivector is the same as before.
    if (OverlapX_!=0) {
      if (OverlapX_->NumVectors()!=X.NumVectors()) {
	delete OverlapX_; OverlapX_ = 0;
	delete OverlapY_; OverlapY_ = 0;
      }
    }
    if (OverlapX_==0) { // Need to allocate space for overlap X and Y
      OverlapX_ = new Epetra_MultiVector(Graph_.OverlapGraph()->RowMap(), X.NumVectors());
      OverlapY_ = new Epetra_MultiVector(Graph_.OverlapGraph()->RowMap(), Y.NumVectors());
    }
    OverlapX_->Import(X,*Graph_.OverlapImporter(), Insert); // Import X values for solve
    X1 = OverlapX_;
    Y1 = OverlapY_; // Set pointers for X1 and Y1 to point to overlap space
  }

  Epetra_Flops * counter = this->GetFlopCounter();
  if (counter!=0) {
    L_->SetFlopCounter(*counter);
    Y1->SetFlopCounter(*counter);
    U_->SetFlopCounter(*counter);
  }

  if (!Trans) {

    L_->Solve(Lower, Trans, UnitDiagonal, *X1, *Y1);
    Y1->Multiply(1.0, *D_, *Y1, 0.0); // y = D*y (D_ has inverse of diagonal)
    U_->Solve(Upper, Trans, UnitDiagonal, *Y1, *Y1); // Solve Uy = y
  }
  else
    {
      U_->Solve(Upper, Trans, UnitDiagonal, *X1, *Y1); // Solve Uy = y
      Y1->Multiply(1.0, *D_, *Y1, 0.0); // y = D*y (D_ has inverse of diagonal)
      L_->Solve(Lower, Trans, UnitDiagonal, *Y1, *Y1);
      
    } 

  // Export computed Y values as directed
  if (Graph_.LevelOverlap()>0 && Graph_.L_Graph().DomainMap().DistributedGlobal())
    Y.Export(*OverlapY_,*Graph_.OverlapImporter(), OverlapMode_);
  return(0);
}
//=============================================================================
int Ifpack_CrsRick::Multiply(bool Trans, const Epetra_MultiVector& X, 
				Epetra_MultiVector& Y) const {
//
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
//

  if (X.NumVectors()!=Y.NumVectors()) EPETRA_CHK_ERR(-1); // Return error: X and Y not the same size

  bool Upper = true;
  bool Lower = false;
  bool UnitDiagonal = true;

  Epetra_MultiVector * X1 = (Epetra_MultiVector *) &X;
  Epetra_MultiVector * Y1 = (Epetra_MultiVector *) &Y;

  if (Graph_.LevelOverlap()>0 && Graph_.L_Graph().DomainMap().DistributedGlobal()) {
    // Make sure the number of vectors in the multivector is the same as before.
    if (OverlapX_!=0) {
      if (OverlapX_->NumVectors()!=X.NumVectors()) {
	delete OverlapX_; OverlapX_ = 0;
	delete OverlapY_; OverlapY_ = 0;
      }
    }
    if (OverlapX_==0) { // Need to allocate space for overlap X and Y
      OverlapX_ = new Epetra_MultiVector(Graph_.OverlapGraph()->RowMap(), X.NumVectors());
      OverlapY_ = new Epetra_MultiVector(Graph_.OverlapGraph()->RowMap(), Y.NumVectors());
    }
    OverlapX_->Import(X,*Graph_.OverlapImporter(), Insert); // Import X values for solve
    X1 = OverlapX_;
    Y1 = OverlapY_; // Set pointers for X1 and Y1 to point to overlap space
  }

  Epetra_Flops * counter = this->GetFlopCounter();
  if (counter!=0) {
    L_->SetFlopCounter(*counter);
    Y1->SetFlopCounter(*counter);
    U_->SetFlopCounter(*counter);
  }

  if (Trans) {

    L_->Multiply(Trans, *X1, *Y1);
    Y1->Update(1.0, *X1, 1.0); // Y1 = Y1 + X1 (account for implicit unit diagonal)
    Y1->ReciprocalMultiply(1.0, *D_, *Y1, 0.0); // y = D*y (D_ has inverse of diagonal)
    Epetra_MultiVector Y1temp(*Y1); // Need a temp copy of Y1
    U_->Multiply(Trans, Y1temp, *Y1);
    Y1->Update(1.0, Y1temp, 1.0); // (account for implicit unit diagonal)
  }
  else {
    U_->Multiply(Trans, *X1, *Y1); // 
    Y1->Update(1.0, *X1, 1.0); // Y1 = Y1 + X1 (account for implicit unit diagonal)
    Y1->ReciprocalMultiply(1.0, *D_, *Y1, 0.0); // y = D*y (D_ has inverse of diagonal)
    Epetra_MultiVector Y1temp(*Y1); // Need a temp copy of Y1
    L_->Multiply(Trans, Y1temp, *Y1);
    Y1->Update(1.0, Y1temp, 1.0); // (account for implicit unit diagonal)
    } 

  // Export computed Y values as directed
  if (Graph_.LevelOverlap()>0 && Graph_.L_Graph().DomainMap().DistributedGlobal())
    Y.Export(*OverlapY_,*Graph_.OverlapImporter(), OverlapMode_);
  return(0);
}
//=============================================================================
int Ifpack_CrsRick::Condest(bool Trans, double & ConditionNumberEstimate) const {

  if (Condest_>=0.0) {
    ConditionNumberEstimate = Condest_;
    return(0);
  }
  // Create a vector with all values equal to one
  Epetra_Vector Ones(A_.RowMap());
  Epetra_Vector OnesResult(Ones);
  Ones.PutScalar(1.0);

  EPETRA_CHK_ERR(Solve(Trans, Ones, OnesResult)); // Compute the effect of the solve on the vector of ones
  EPETRA_CHK_ERR(OnesResult.Abs(OnesResult)); // Make all values non-negative
  EPETRA_CHK_ERR(OnesResult.MaxValue(&ConditionNumberEstimate)); // Get the maximum value across all processors
  Condest_ = ConditionNumberEstimate; // Save value for possible later calls
  return(0);
}
//=============================================================================
// Non-member functions

ostream& operator << (ostream& os, const Ifpack_CrsRick& A)
{
/*  Epetra_fmtflags olda = os.setf(ios::right,ios::adjustfield);
  Epetra_fmtflags oldf = os.setf(ios::scientific,ios::floatfield);
  int oldp = os.precision(12); */
  int LevelFill = A.Graph().LevelFill();
  int LevelOverlap = A.Graph().LevelOverlap();
  Epetra_CrsMatrix & L = (Epetra_CrsMatrix &) A.L();
  Epetra_CrsMatrix & U = (Epetra_CrsMatrix &) A.U();
  Epetra_Vector & D = (Epetra_Vector &) A.D();

  os.width(14);
  os << endl;
  os <<  "     Level of Fill = "; os << LevelFill;
  os << endl;
  os.width(14);
  os <<  "     Level of Overlap = "; os << LevelOverlap;
  os << endl;

  os.width(14);
  os <<  "     Lower Triangle = ";
  os << endl;
  os << L; // Let Epetra_CrsMatrix handle the rest.
  os << endl;

  os.width(14);
  os <<  "     Inverse of Diagonal = ";
  os << endl;
  os << D; // Let Epetra_Vector handle the rest.
  os << endl;

  os.width(14);
  os <<  "     Upper Triangle = ";
  os << endl;
  os << U; // Let Epetra_CrsMatrix handle the rest.
  os << endl;
 
  // Reset os flags

/*  os.setf(olda,ios::adjustfield);
  os.setf(oldf,ios::floatfield);
  os.precision(oldp); */

  return os;
}
