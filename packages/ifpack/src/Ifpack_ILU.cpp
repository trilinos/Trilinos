//@HEADER
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

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_CondestType.h"
#include "Ifpack_ILU.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

using Teuchos::RefCountPtr;
using Teuchos::rcp;

#ifdef IFPACK_TEUCHOS_TIME_MONITOR
#  include "Teuchos_TimeMonitor.hpp"
#endif

//==============================================================================
Ifpack_ILU::Ifpack_ILU(Epetra_RowMatrix* Matrix_in) :
  A_(rcp(Matrix_in,false)),
  Comm_(Matrix_in->Comm()),
  UseTranspose_(false),
  NumMyDiagonals_(0),
  RelaxValue_(0.0),
  Athresh_(0.0),
  Rthresh_(1.0),
  Condest_(-1.0),
  LevelOfFill_(0),
  IsInitialized_(false),
  IsComputed_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyInverseTime_(0.0),
  ComputeFlops_(0.0),
  ApplyInverseFlops_(0.0),
  Time_(Comm())
{
  Teuchos::ParameterList List;
  SetParameters(List);
}

//==============================================================================
void Ifpack_ILU::Destroy()
{
  // reset pointers to already allocated stuff
  U_DomainMap_ = 0;
  L_RangeMap_ = 0;
}

//==========================================================================
int Ifpack_ILU::SetParameters(Teuchos::ParameterList& List)
{
  RelaxValue_ = List.get("fact: relax value", RelaxValue_);
  Athresh_ = List.get("fact: absolute threshold", Athresh_);
  Rthresh_ = List.get("fact: relative threshold", Rthresh_);
  LevelOfFill_ = List.get("fact: level-of-fill", LevelOfFill_);

  // set label
  sprintf(Label_, "IFPACK ILU (fill=%d, relax=%f, athr=%f, rthr=%f)",
	  LevelOfFill(), RelaxValue(), AbsoluteThreshold(), 
	  RelativeThreshold());
  return(0);
}

//==========================================================================
int Ifpack_ILU::ComputeSetup() 
{

#ifdef IFPACK_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Ifpack_ILU::ComputeSetup");
#endif

  L_ = rcp(new Epetra_CrsMatrix(Copy, Graph().L_Graph()));
  U_ = rcp(new Epetra_CrsMatrix(Copy, Graph().U_Graph()));
  D_ = rcp(new Epetra_Vector(Graph().L_Graph().RowMap()));
  if ((L_.get() == 0) || (U_.get() == 0) || (D_.get() == 0))
    IFPACK_CHK_ERR(-5);

  // Get Maximun Row length
  int MaxNumEntries = Matrix().MaxNumEntries();

  // Set L range map and U domain map
  U_DomainMap_ = &(Matrix().OperatorDomainMap());
  L_RangeMap_ = &(Matrix().OperatorRangeMap());
 
  // this is the old InitAllValues()
  int ierr = 0;
  int i, j;
  int NumIn, NumL, NumU;
  bool DiagFound;
  int NumNonzeroDiags = 0;

  vector<int> InI(MaxNumEntries); // Allocate temp space
  vector<int> LI(MaxNumEntries);
  vector<int> UI(MaxNumEntries);
  vector<double> InV(MaxNumEntries);
  vector<double> LV(MaxNumEntries);
  vector<double> UV(MaxNumEntries);

  bool ReplaceValues = (L_->StaticGraph() || L_->IndicesAreLocal()); // Check if values should be inserted or replaced

  if (ReplaceValues) {
    L_->PutScalar(0.0); // Zero out L and U matrices
    U_->PutScalar(0.0);
  }

  D_->PutScalar(0.0); // Set diagonal values to zero
  double *DV;
  IFPACK_CHK_ERR(D_->ExtractView(&DV)); // Get view of diagonal

  // First we copy the user's matrix into L and U, regardless of fill level

  for (i=0; i< NumMyRows(); i++) {

    IFPACK_CHK_ERR(Matrix().ExtractMyRowCopy(i, MaxNumEntries, NumIn, &InV[0], &InI[0])); // Get Values and Indices
    
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

      else if (k < 0) {IFPACK_CHK_ERR(-4);} // Out of range

      else if (k < i) {
	LI[NumL] = k;
	LV[NumL] = InV[j];
	NumL++;
      }
      else if (k<NumMyRows()) {
	UI[NumU] = k;
	UV[NumU] = InV[j];
	NumU++;
      }
    }
    
    // Check in things for this row of L and U

    if (DiagFound) NumNonzeroDiags++;
    else DV[i] = Athresh_;

    if (NumL) {
      if (ReplaceValues) {
	(L_->ReplaceMyValues(i, NumL, &LV[0], &LI[0]));
      }
      else {
	IFPACK_CHK_ERR(L_->InsertMyValues(i, NumL, &LV[0], &LI[0]));
      }
    }

    if (NumU) {
      if (ReplaceValues) {
	(U_->ReplaceMyValues(i, NumU, &UV[0], &UI[0]));
      }
      else {
	IFPACK_CHK_ERR(U_->InsertMyValues(i, NumU, &UV[0], &UI[0]));
      }
    }
    
  }

  if (!ReplaceValues) {
    // The domain of L and the range of U are exactly their own row maps (there is no communication).
    // The domain of U and the range of L must be the same as those of the original matrix,
    // However if the original matrix is a VbrMatrix, these two latter maps are translation from
    // a block map to a point map.
    IFPACK_CHK_ERR(L_->FillComplete((L_->RowMatrixColMap()), *L_RangeMap_));
    IFPACK_CHK_ERR(U_->FillComplete(*U_DomainMap_, U_->RowMatrixRowMap()));
  }

  // At this point L and U have the values of A in the structure of L and U, and diagonal vector D
  int TotalNonzeroDiags = 0;
  IFPACK_CHK_ERR(Graph().L_Graph().RowMap().Comm().SumAll(&NumNonzeroDiags, &TotalNonzeroDiags, 1));
  NumMyDiagonals_ = NumNonzeroDiags;
  if (NumNonzeroDiags != NumMyRows()) ierr = 1; // Diagonals are not right, warn user

  IFPACK_CHK_ERR(ierr);
  return(ierr);
}

//==========================================================================
int Ifpack_ILU::Initialize() 
{

#ifdef IFPACK_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Ifpack_ILU::Initialize");
#endif

  Time_.ResetStartTime();
  IsInitialized_ = false;

  // reset this object
  Destroy();

  Epetra_CrsMatrix* CrsMatrix;
  CrsMatrix = dynamic_cast<Epetra_CrsMatrix*>(&*A_);
  if (CrsMatrix == 0) {
    // this means that we have to create
    // the graph from a given Epetra_RowMatrix. Note
    // that at this point we are ignoring any possible
    // graph coming from VBR matrices.
    int size = A_->MaxNumEntries();
    CrsGraph_ = rcp(new Epetra_CrsGraph(Copy,A_->RowMatrixRowMap(), size));
    if (CrsGraph_.get() == 0)
      IFPACK_CHK_ERR(-5); // memory allocation error

    vector<int> Indices(size);
    vector<double> Values(size);

    // extract each row at-a-time, and insert it into
    // the graph, ignore all off-process entries
    for (int i = 0 ; i < A_->NumMyRows() ; ++i) {
      int NumEntries;
      int GlobalRow = A_->RowMatrixRowMap().GID(i);
      IFPACK_CHK_ERR(A_->ExtractMyRowCopy(i, size, NumEntries, 
					  &Values[0], &Indices[0]));
      // convert to global indices
      for (int j = 0 ; j < NumEntries ; ++j) {
	Indices[j] = A_->RowMatrixColMap().GID(Indices[j]); 
      }
      IFPACK_CHK_ERR(CrsGraph_->InsertGlobalIndices(GlobalRow,NumEntries,
						   &Indices[0]));
    }
    
    IFPACK_CHK_ERR(CrsGraph_->FillComplete(A_->RowMatrixRowMap(),
					  A_->RowMatrixRowMap()));

    // always overlap zero, wider overlap will be handled
    // by the AdditiveSchwarz preconditioner.
    Graph_ = rcp(new Ifpack_IlukGraph(*CrsGraph_, LevelOfFill_, 0));

  }
  else {
    // see comment above for the overlap.
    Graph_ = rcp(new Ifpack_IlukGraph(CrsMatrix->Graph(), LevelOfFill_, 0));
  }

  if (Graph_.get() == 0)
    IFPACK_CHK_ERR(-5); // memory allocation error
  IFPACK_CHK_ERR(Graph_->ConstructFilledGraph());

  IsInitialized_ = true;
  NumInitialize_++;
  InitializeTime_ += Time_.ElapsedTime();

  return(0);
}

//==========================================================================
int Ifpack_ILU::Compute() 
{

#ifdef IFPACK_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Ifpack_ILU::Compute");
#endif

  if (!IsInitialized()) 
    IFPACK_CHK_ERR(Initialize());

  Time_.ResetStartTime();
  IsComputed_ = false;

  // convert Matrix() into L and U factors.
  IFPACK_CHK_ERR(ComputeSetup());

  // MinMachNum should be officially defined, for now pick something a little 
  // bigger than IEEE underflow value

  double MinDiagonalValue = Epetra_MinDouble;
  double MaxDiagonalValue = 1.0/MinDiagonalValue;

  int ierr = 0;
  int i, j, k;
  int *LI, *UI;
  double *LV, *UV;
  int NumIn, NumL, NumU;

  // Get Maximun Row length
  int MaxNumEntries = L_->MaxNumEntries() + U_->MaxNumEntries() + 1;

  vector<int> InI(MaxNumEntries+1);    // Allocate temp space, pad by one to 
  vector<double> InV(MaxNumEntries+1); // to avoid debugger complaints for pathological cases
  vector<int> colflag(NumMyCols());

  double *DV;
  ierr = D_->ExtractView(&DV); // Get view of diagonal

#ifdef IFPACK_FLOPCOUNTERS
  int current_madds = 0; // We will count multiply-add as they happen
#endif

  // =========================== //
  // Now start the factorization //
  // =========================== //

  // Need some integer workspace and pointers
  int NumUU; 
  int * UUI;
  double * UUV;
  for (j = 0; j < NumMyCols(); ++j) colflag[j] = - 1;

  for (i = 0; i < NumMyRows(); ++i) {

    // Fill InV, InI with current row of L, D and U combined

    NumIn = MaxNumEntries;
    IFPACK_CHK_ERR(L_->ExtractMyRowCopy(i, NumIn, NumL, &InV[0], &InI[0]));
    LV = &InV[0];
    LI = &InI[0];

    InV[NumL] = DV[i]; // Put in diagonal
    InI[NumL] = i;
    
    IFPACK_CHK_ERR(U_->ExtractMyRowCopy(i, NumIn-NumL-1, NumU, &InV[NumL+1], &InI[NumL+1]));
    NumIn = NumL+NumU+1;
    UV = &InV[NumL+1];
    UI = &InI[NumL+1];

    // Set column flags
    for (j=0; j<NumIn; j++) colflag[InI[j]] = j;

    double diagmod = 0.0; // Off-diagonal accumulator

    for (int jj=0; jj<NumL; jj++) {
      j = InI[jj];
      double multiplier = InV[jj]; // current_mults++;

      InV[jj] *= DV[j];
      
      IFPACK_CHK_ERR(U_->ExtractMyRowView(j, NumUU, UUV, UUI)); // View of row above

      if (RelaxValue_==0.0) {
	for (k=0; k<NumUU; k++) {
	  int kk = colflag[UUI[k]];
	  if (kk>-1) {
	    InV[kk] -= multiplier*UUV[k];
#ifdef IFPACK_FLOPCOUNTERS
	    current_madds++;
#endif
	  }
	}
      }
      else {
	for (k=0; k<NumUU; k++) {
	  int kk = colflag[UUI[k]];
	  if (kk>-1) InV[kk] -= multiplier*UUV[k];
	  else diagmod -= multiplier*UUV[k];
#ifdef IFPACK_FLOPCOUNTERS
	  current_madds++;
#endif
	}
      }
     }
    if (NumL) {
      IFPACK_CHK_ERR(L_->ReplaceMyValues(i, NumL, LV, LI));  // Replace current row of L
    }

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

    if (NumU) {
      IFPACK_CHK_ERR(U_->ReplaceMyValues(i, NumU, UV, UI));  // Replace current row of L and U
    }

    // Reset column flags
    for (j=0; j<NumIn; j++) colflag[InI[j]] = -1;
  }

  // Validate that the L and U factors are actually lower and upper triangular

  if (!L_->LowerTriangular()) 
    IFPACK_CHK_ERR(-4);
  if (!U_->UpperTriangular()) 
    IFPACK_CHK_ERR(-4);
  
#ifdef IFPACK_FLOPCOUNTERS
  // Add up flops
 
  double current_flops = 2 * current_madds;
  double total_flops = 0;
    
  IFPACK_CHK_ERR(Graph().L_Graph().RowMap().Comm().SumAll(&current_flops, &total_flops, 1)); // Get total madds across all PEs

  // Now count the rest
  total_flops += (double) L_->NumGlobalNonzeros(); // Accounts for multiplier above
  total_flops += (double) D_->GlobalLength(); // Accounts for reciprocal of diagonal
  if (RelaxValue_!=0.0) total_flops += 2 * (double)D_->GlobalLength(); // Accounts for relax update of diag

  // add to this object's counter
  ComputeFlops_ += total_flops;
#endif
  
  IsComputed_ = true;
  NumCompute_++;
  ComputeTime_ += Time_.ElapsedTime();

  return(ierr);

}

//=============================================================================
// This function finds Y such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int Ifpack_ILU::Solve(bool Trans, const Epetra_MultiVector& X, 
                      Epetra_MultiVector& Y) const 
{

#ifdef IFPACK_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Ifpack_ILU::ApplyInverse - Solve");
#endif

  // in this function the overlap is always zero
  bool Upper = true;
  bool Lower = false;
  bool UnitDiagonal = true;

  if (!Trans) {

    IFPACK_CHK_ERR(L_->Solve(Lower, Trans, UnitDiagonal, X, Y));
    // y = D*y (D_ has inverse of diagonal)
    IFPACK_CHK_ERR(Y.Multiply(1.0, *D_, Y, 0.0)); 
    // Solve Uy = y
    IFPACK_CHK_ERR(U_->Solve(Upper, Trans, UnitDiagonal, Y, Y)); 
  }
  else {
    // Solve Uy = y
    IFPACK_CHK_ERR(U_->Solve(Upper, Trans, UnitDiagonal, X, Y)); 
    // y = D*y (D_ has inverse of diagonal)
    IFPACK_CHK_ERR(Y.Multiply(1.0, *D_, Y, 0.0)); 
    IFPACK_CHK_ERR(L_->Solve(Lower, Trans, UnitDiagonal, Y, Y));
  } 


  return(0);
}

//=============================================================================
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int Ifpack_ILU::Multiply(bool Trans, const Epetra_MultiVector& X, 
				Epetra_MultiVector& Y) const 
{

#ifdef IFPACK_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Ifpack_ILU::Multiply");
#endif

  if (!IsComputed())
    IFPACK_CHK_ERR(-3);

  if (!Trans) {
    IFPACK_CHK_ERR(U_->Multiply(Trans, X, Y)); 
    // Y1 = Y1 + X1 (account for implicit unit diagonal)
    IFPACK_CHK_ERR(Y.Update(1.0, X, 1.0)); 
    // y = D*y (D_ has inverse of diagonal)
    IFPACK_CHK_ERR(Y.ReciprocalMultiply(1.0, *D_, Y, 0.0)); 
    Epetra_MultiVector Y1temp(Y); // Need a temp copy of Y1
    IFPACK_CHK_ERR(L_->Multiply(Trans, Y1temp, Y));
    // (account for implicit unit diagonal)
    IFPACK_CHK_ERR(Y.Update(1.0, Y1temp, 1.0)); 
  }
  else {

    IFPACK_CHK_ERR(L_->Multiply(Trans, X, Y));
    // Y1 = Y1 + X1 (account for implicit unit diagonal)
    IFPACK_CHK_ERR(Y.Update(1.0, X, 1.0)); 
    // y = D*y (D_ has inverse of diagonal)
    IFPACK_CHK_ERR(Y.ReciprocalMultiply(1.0, *D_, Y, 0.0)); 
    Epetra_MultiVector Y1temp(Y); // Need a temp copy of Y1
    IFPACK_CHK_ERR(U_->Multiply(Trans, Y1temp, Y));
    // (account for implicit unit diagonal)
    IFPACK_CHK_ERR(Y.Update(1.0, Y1temp, 1.0)); 
  } 

  return(0);
}

//=============================================================================
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int Ifpack_ILU::ApplyInverse(const Epetra_MultiVector& X, 
                             Epetra_MultiVector& Y) const
{

#ifdef IFPACK_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Ifpack_ILU::ApplyInverse");
#endif

  if (!IsComputed())
    IFPACK_CHK_ERR(-3);

  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-2);

  Time_.ResetStartTime();

  // AztecOO gives X and Y pointing to the same memory location,
  // need to create an auxiliary vector, Xcopy
  Teuchos::RefCountPtr< const Epetra_MultiVector > Xcopy;
  if (X.Pointers()[0] == Y.Pointers()[0])
    Xcopy = Teuchos::rcp( new Epetra_MultiVector(X) );
  else
    Xcopy = Teuchos::rcp( &X, false );

  IFPACK_CHK_ERR(Solve(Ifpack_ILU::UseTranspose(), *Xcopy, Y));

  // approx is the number of nonzeros in L and U
#ifdef IFPACK_FLOPCOUNTERS
  ApplyInverseFlops_ += X.NumVectors() * 4 * 
    (L_->NumGlobalNonzeros() + U_->NumGlobalNonzeros());
#endif

  ++NumApplyInverse_;
  ApplyInverseTime_ += Time_.ElapsedTime();

  return(0);

}

//=============================================================================
double Ifpack_ILU::Condest(const Ifpack_CondestType CT, 
                           const int MaxIters, const double Tol,
                              Epetra_RowMatrix* Matrix_in)
{

#ifdef IFPACK_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Ifpack_ILU::Condest");
#endif

  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  Condest_ = Ifpack_Condest(*this, CT, MaxIters, Tol, Matrix_in);

  return(Condest_);
}

//=============================================================================
std::ostream&
Ifpack_ILU::Print(std::ostream& os) const
{
  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack_ILU: " << Label() << endl << endl;
    os << "Level-of-fill      = " << LevelOfFill() << endl;
    os << "Absolute threshold = " << AbsoluteThreshold() << endl;
    os << "Relative threshold = " << RelativeThreshold() << endl;
    os << "Relax value        = " << RelaxValue() << endl;
    os << "Condition number estimate = " << Condest() << endl;
    os << "Global number of rows            = " << A_->NumGlobalRows() << endl;
    if (IsComputed_) {
      os << "Number of rows of L, D, U       = " << L_->NumGlobalRows() << endl;
      os << "Number of nonzeros of L + U     = " << NumGlobalNonzeros() << endl;
      os << "nonzeros / rows                 = " 
        << 1.0 * NumGlobalNonzeros() / U_->NumGlobalRows() << endl;
    }
    os << endl;
    os << "Phase           # calls   Total Time (s)       Total MFlops     MFlops/s" << endl;
    os << "-----           -------   --------------       ------------     --------" << endl;
    os << "Initialize()    "   << std::setw(5) << NumInitialize() 
       << "  " << std::setw(15) << InitializeTime() 
       << "               0.0            0.0" << endl;
    os << "Compute()       "   << std::setw(5) << NumCompute() 
       << "  " << std::setw(15) << ComputeTime()
       << "  " << std::setw(15) << 1.0e-6 * ComputeFlops();
    if (ComputeTime() != 0.0) 
      os << "  " << std::setw(15) << 1.0e-6 * ComputeFlops() / ComputeTime() << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "ApplyInverse()  "   << std::setw(5) << NumApplyInverse() 
       << "  " << std::setw(15) << ApplyInverseTime()
       << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops();
    if (ApplyInverseTime() != 0.0) 
      os << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops() / ApplyInverseTime() << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "================================================================================" << endl;
    os << endl;
  }

  return(os);
}
