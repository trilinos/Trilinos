#include "Ifpack_CrsRiluk.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

//==============================================================================
Ifpack_CrsRiluk::Ifpack_CrsRiluk(const Ifpack_IlukGraph & Graph) 
  : UserMatrixIsVbr_(false),
    UserMatrixIsCrs_(false),
    Graph_(Graph),
    DomainMap_(Graph.DomainMap()),
    RangeMap_(Graph.RangeMap()),
    IlukRowMap_(0),
    Comm_(Graph.Comm()),
    UseTranspose_(false),
    NumMyDiagonals_(0),
    Allocated_(false),
    ValuesInitialized_(false),
    Factored_(false),
    RelaxValue_(0.0),
    Athresh_(0.0),
    Rthresh_(1.0),
    Condest_(-1.0),
    OverlapX_(0),
    OverlapY_(0),
    OverlapMode_(Zero)
{
  // Test for non-trivial overlap here so we can use it later.
  IsOverlapped_ = (Graph_.LevelOverlap()>0 && DomainMap_.DistributedGlobal());
}

//==============================================================================
Ifpack_CrsRiluk::Ifpack_CrsRiluk(const Ifpack_CrsRiluk & FactoredMatrix) 
  : UserMatrixIsVbr_(FactoredMatrix.UserMatrixIsVbr_),
    UserMatrixIsCrs_(FactoredMatrix.UserMatrixIsCrs_),
    IsOverlapped_(FactoredMatrix.IsOverlapped_),
    Graph_(FactoredMatrix.Graph_),
    DomainMap_(FactoredMatrix.DomainMap_),
    RangeMap_(FactoredMatrix.RangeMap_),
    IlukRowMap_(FactoredMatrix.IlukRowMap_),
    Comm_(FactoredMatrix.Comm_),
    UseTranspose_(FactoredMatrix.UseTranspose_),
    NumMyDiagonals_(FactoredMatrix.NumMyDiagonals_),
    Allocated_(FactoredMatrix.Allocated_),
    ValuesInitialized_(FactoredMatrix.ValuesInitialized_),
    Factored_(FactoredMatrix.Factored_),
    RelaxValue_(FactoredMatrix.RelaxValue_),
    Athresh_(FactoredMatrix.Athresh_),
    Rthresh_(FactoredMatrix.Rthresh_),
    Condest_(FactoredMatrix.Condest_),
    OverlapX_(0),
    OverlapY_(0),
    OverlapMode_(FactoredMatrix.OverlapMode_)
{
  L_ = new Epetra_CrsMatrix(FactoredMatrix.L());
  U_ = new Epetra_CrsMatrix(FactoredMatrix.U());
  D_ = new Epetra_Vector(FactoredMatrix.D());
  
}

//==============================================================================
int Ifpack_CrsRiluk::AllocateCrs() {

  // Allocate Epetra_CrsMatrix using ILUK graphs
  L_ = new Epetra_CrsMatrix(Copy, Graph_.L_Graph());
  U_ = new Epetra_CrsMatrix(Copy, Graph_.U_Graph());
  D_ = new Epetra_Vector(Graph_.L_Graph().RowMap());
  L_Graph_ = 0;
  U_Graph_ = 0;
  SetAllocated(true);
  return(0);
}
//==============================================================================
int Ifpack_CrsRiluk::AllocateVbr() {

  // First we need to create an Epetra_Map that has the same number of points  as the
  // Epetra_BlockMap associated with the Overlap Graph.  Since all ILU computations are local
  // we need not worry about off-processor communication and can choose a simple linear 
  // set of global IDs.

  int NumMyPoints = Graph_.L_Graph().RowMap().NumMyPoints(); // Number of points (equations) on this processor
  
  IlukRowMap_ = new Epetra_Map(-1, NumMyPoints, 0, Comm());

  if (Graph().LevelFill()) { // If there is fill, then pre-build the L and U structures from the Block version of L and U.
    L_Graph_ = new Epetra_CrsGraph(Copy, *IlukRowMap_, *IlukRowMap_, 0);
    U_Graph_ = new Epetra_CrsGraph(Copy, *IlukRowMap_, *IlukRowMap_, 0);
    EPETRA_CHK_ERR(BlockGraph2PointGraph(Graph_.L_Graph(), *L_Graph_, false));
    EPETRA_CHK_ERR(BlockGraph2PointGraph(Graph_.U_Graph(), *U_Graph_, true));
    L_Graph_->TransformToLocal();
    U_Graph_->TransformToLocal();

    L_ = new Epetra_CrsMatrix(Copy, *L_Graph_);
    U_ = new Epetra_CrsMatrix(Copy, *U_Graph_);
    D_ = new Epetra_Vector(*IlukRowMap_);
  }
  else {
    // Allocate Epetra_CrsMatrix using ILUK graphs
    L_ = new Epetra_CrsMatrix(Copy, *IlukRowMap_, *IlukRowMap_, 0);
    U_ = new Epetra_CrsMatrix(Copy, *IlukRowMap_, *IlukRowMap_, 0);
    D_ = new Epetra_Vector(*IlukRowMap_);
    L_Graph_ = 0;
    U_Graph_ = 0;
  }
  return(0);
}
//==============================================================================
int Ifpack_CrsRiluk::BlockGraph2PointGraph(const Epetra_CrsGraph & BG, Epetra_CrsGraph & PG, bool Upper) {

  if (!BG.IndicesAreLocal()) {EPETRA_CHK_ERR(-1);} // Must have done TransformToLocal on BG

  int * ColFirstPointInElementList = BG.RowMap().FirstPointInElementList();
  int * ColElementSizeList = BG.RowMap().ElementSizeList();
  if (BG.Importer()!=0) {
    ColFirstPointInElementList = BG.ImportMap().FirstPointInElementList();
    ColElementSizeList = BG.ImportMap().ElementSizeList();
  }

  int Length = (BG.MaxNumIndices()+1) * BG.ImportMap().MaxMyElementSize();
  int * tmpIndices = new int[Length];


  int BlockRow, BlockOffset, NumEntries;
  int NumBlockEntries;
  int * BlockIndices;

  int NumMyRows = PG.NumMyRows();

  for (int i=0; i<NumMyRows; i++) {
    EPETRA_CHK_ERR(BG.RowMap().FindLocalElementID(i, BlockRow, BlockOffset));
    EPETRA_CHK_ERR(BG.ExtractMyRowView(BlockRow, NumBlockEntries, BlockIndices));

    int * ptr = tmpIndices; // Set pointer to beginning of buffer

    int RowDim = BG.RowMap().ElementSize(BlockRow);
    NumEntries = 0;

    // This next line make sure that the off-diagonal entries in the block diagonal of the 
    // original block entry matrix are included in the nonzero pattern of the point graph
    if (Upper) {
      int jstart = i+1;
      int jstop = EPETRA_MIN(NumMyRows,i+RowDim-BlockOffset);
      for (int j= jstart; j< jstop; j++) {*ptr++ = j; NumEntries++;}
    }

    for (int j=0; j<NumBlockEntries; j++) {
      int ColDim = ColElementSizeList[BlockIndices[j]];
      NumEntries += ColDim;
      assert(NumEntries<=Length); // Sanity test
      int Index = ColFirstPointInElementList[BlockIndices[j]];
      for (int k=0; k < ColDim; k++) *ptr++ = Index++;
    }

    // This next line make sure that the off-diagonal entries in the block diagonal of the 
    // original block entry matrix are included in the nonzero pattern of the point graph
    if (!Upper) {
      int jstart = EPETRA_MAX(0,i-RowDim+1);
      int jstop = i;
      for (int j = jstart; j < jstop; j++) {*ptr++ = j; NumEntries++;}
    }

    EPETRA_CHK_ERR(PG.InsertMyIndices(i, NumEntries, tmpIndices));
  }
  delete [] tmpIndices;

  SetAllocated(true);

  return(0);
}
//==============================================================================
Ifpack_CrsRiluk::~Ifpack_CrsRiluk(){


  delete L_;
  delete U_;
  delete D_; // Diagonal is stored separately.  We store the inverse.

  if (OverlapX_!=0) delete OverlapX_;
  if (OverlapY_!=0) delete OverlapY_;
  if (L_Graph_!=0) delete L_Graph_;
  if (U_Graph_!=0) delete U_Graph_;
  if (IlukRowMap_!=0) delete IlukRowMap_;

  OverlapX_ = 0;
  OverlapY_ = 0;
  IlukRowMap_ = 0;
  
  ValuesInitialized_ = false;
  Factored_ = false;
  Allocated_ = false;
}

//==========================================================================
int Ifpack_CrsRiluk::InitValues(const Epetra_CrsMatrix & A) {

  UserMatrixIsCrs_ = true;

  if (!Allocated()) AllocateCrs();

  Epetra_CrsMatrix * OverlapA = (Epetra_CrsMatrix *) &A;

  if (IsOverlapped_) {
  
    OverlapA = new Epetra_CrsMatrix(Copy, *Graph_.OverlapGraph());
    EPETRA_CHK_ERR(OverlapA->Import(A, *Graph_.OverlapImporter(), Insert));
    EPETRA_CHK_ERR(OverlapA->TransformToLocal());
  }
  
  // Get Maximun Row length
  int MaxNumEntries = OverlapA->MaxNumEntries();

  // Do the rest using generic Epetra_RowMatrix interface

  EPETRA_CHK_ERR(InitAllValues(*OverlapA, MaxNumEntries));

  if (IsOverlapped_) delete OverlapA;

  return(0);
}

//==========================================================================
int Ifpack_CrsRiluk::InitValues(const Epetra_VbrMatrix & A) {

  UserMatrixIsVbr_ = true;

  if (!Allocated()) AllocateVbr();

  //cout << "Original Graph " << endl <<  A.Graph() << endl << flush;
  //A.Comm().Barrier(); 
  //if (A.Comm().MyPID()==0) cout << "*****************************************************" <<endl;
  //cout << "Original Matrix " << endl << A << endl << flush;
  //A.Comm().Barrier(); 
  //if (A.Comm().MyPID()==0) cout << "*****************************************************" <<endl;
  //cout << "Overlap Graph " << endl << *Graph_.OverlapGraph() << endl << flush;
  //A.Comm().Barrier(); 
  //if (A.Comm().MyPID()==0) cout << "*****************************************************" <<endl;

  Epetra_VbrMatrix * OverlapA = (Epetra_VbrMatrix *) &A;

  if (IsOverlapped_) {
  
    OverlapA = new Epetra_VbrMatrix(Copy, *Graph_.OverlapGraph());
    EPETRA_CHK_ERR(OverlapA->Import(A, *Graph_.OverlapImporter(), Insert));
    EPETRA_CHK_ERR(OverlapA->TransformToLocal());
  }
  
  //cout << "Overlap Matrix " << endl << *OverlapA << endl << flush;

  // Get Maximun Row length
  int MaxNumEntries = OverlapA->MaxNumNonzeros();

  // Do the rest using generic Epetra_RowMatrix interface

  EPETRA_CHK_ERR(InitAllValues(*OverlapA, MaxNumEntries));

  if (IsOverlapped_) delete OverlapA;

  return(0);
}
//==========================================================================

int Ifpack_CrsRiluk::InitAllValues(const Epetra_RowMatrix & OverlapA, int MaxNumEntries) {

  int ierr = 0;
  int i, j;
  int * InI=0, * LI=0, * UI = 0;
  double * InV=0, * LV=0, * UV = 0;
  int NumIn, NumL, NumU;
  bool DiagFound;
  int NumNonzeroDiags = 0;


  InI = new int[MaxNumEntries]; // Allocate temp space
  LI = new int[MaxNumEntries];
  UI = new int[MaxNumEntries];
  InV = new double[MaxNumEntries];
  LV = new double[MaxNumEntries];
  UV = new double[MaxNumEntries];

  bool ReplaceValues = (L_->StaticGraph() || L_->IndicesAreLocal()); // Check if values should be inserted or replaced

  if (ReplaceValues) {
    L_->PutScalar(0.0); // Zero out L and U matrices
    U_->PutScalar(0.0);
  }


  double *DV;
  EPETRA_CHK_ERR(D_->ExtractView(&DV)); // Get view of diagonal
    

  // First we copy the user's matrix into L and U, regardless of fill level

  for (i=0; i< NumMyRows(); i++) {

    EPETRA_CHK_ERR(OverlapA.ExtractMyRowCopy(i, MaxNumEntries, NumIn, InV, InI)); // Get Values and Indices
    
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

      else if (k < 0) {EPETRA_CHK_ERR(-1);} // Out of range

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
	EPETRA_CHK_ERR(L_->ReplaceMyValues(i, NumL, LV, LI));
      }
      else {
	EPETRA_CHK_ERR(L_->InsertMyValues(i, NumL, LV, LI));
      }
    }

    if (NumU) {
      if (ReplaceValues) {
	EPETRA_CHK_ERR(U_->ReplaceMyValues(i, NumU, UV, UI));
      }
      else {
	EPETRA_CHK_ERR(U_->InsertMyValues(i, NumU, UV, UI));
      }
    }
    
  }

  delete [] LI;
  delete [] UI;
  delete [] LV;
  delete [] UV;
  delete [] InI;
  delete [] InV;



  if (!ReplaceValues) {
    EPETRA_CHK_ERR(L_->TransformToLocal());
    EPETRA_CHK_ERR(U_->TransformToLocal());
  }

  // At this point L and U have the values of A in the structure of L and U, and diagonal vector D

  SetValuesInitialized(true);
  SetFactored(false);

  int TotalNonzeroDiags = 0;
  EPETRA_CHK_ERR(Graph_.L_Graph().RowMap().Comm().SumAll(&NumNonzeroDiags, &TotalNonzeroDiags, 1));
  NumMyDiagonals_ = NumNonzeroDiags;
  if (NumNonzeroDiags != NumMyRows()) ierr = 1; // Diagonals are not right, warn user

  return(ierr);
}

//==========================================================================
int Ifpack_CrsRiluk::Factor() {

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
  int * LI=0, * UI = 0;
  double * LV=0, * UV = 0;
  int NumIn, NumL, NumU;

  // Get Maximun Row length
  int MaxNumEntries = L_->MaxNumEntries() + U_->MaxNumEntries() + 1;

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
    EPETRA_CHK_ERR(L_->ExtractMyRowCopy(i, NumIn, NumL, InV, InI));
    LV = InV;
    LI = InI;

    InV[NumL] = DV[i]; // Put in diagonal
    InI[NumL] = i;
    
    EPETRA_CHK_ERR(U_->ExtractMyRowCopy(i, NumIn-NumL-1, NumU, InV+NumL+1, InI+NumL+1));
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
      
      EPETRA_CHK_ERR(U_->ExtractMyRowView(j, NumUU, UUV, UUI)); // View of row above

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
    if (NumL) {
      EPETRA_CHK_ERR(L_->ReplaceMyValues(i, NumL, LV, LI));  // Replace current row of L
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
      EPETRA_CHK_ERR(U_->ReplaceMyValues(i, NumU, UV, UI));  // Replace current row of L and U
    }

    // Reset column flags
    for (j=0; j<NumIn; j++) colflag[InI[j]] = -1;
  }

  
  // Add up flops
 
  double current_flops = 2 * current_madds;
  double total_flops = 0;
    
  EPETRA_CHK_ERR(Graph_.L_Graph().RowMap().Comm().SumAll(&current_flops, &total_flops, 1)); // Get total madds across all PEs

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
int Ifpack_CrsRiluk::Solve(bool Trans, const Epetra_Vector& X, 
				Epetra_Vector& Y) const {
//
// This function finds Y such that LDU Y = X or U(trans) D L(trans) Y = X for a single RHS
//

  bool Upper = true;
  bool Lower = false;
  bool UnitDiagonal = true;

  Epetra_Vector * X1 = (Epetra_Vector *) &X;
  Epetra_Vector * Y1 = (Epetra_Vector *) &Y;

  if (IsOverlapped_) {
    if (OverlapX_==0) { // Need to allocate space for overlap X and Y
      OverlapX_ = new Epetra_Vector(Graph_.OverlapGraph()->RowMap());
      OverlapY_ = new Epetra_Vector(Graph_.OverlapGraph()->RowMap());
    }
    EPETRA_CHK_ERR(OverlapX_->Import(X,*Graph_.OverlapImporter(), Insert)); // Import X values for solve
    X1 = (Epetra_Vector *) OverlapX_;
    Y1 = (Epetra_Vector *) OverlapY_; // Set pointers for X1 and Y1 to point to overlap space
  }

  // This code block checks to see if the user matrix, and therefore the user vectors, overlap graph, etc. 
  // are using a block map.  If so, we want make a view of these vectors using a map that is compatible with
  // L and U (which have simple maps).

  Epetra_Vector * X2 = 0;
  Epetra_Vector * Y2 = 0;

  if (UserMatrixIsVbr_) {
    double * xtmp, * ytmp;
    X1->ExtractView(&xtmp);
    Y1->ExtractView(&ytmp);
    X2 = new Epetra_Vector(View, L_->RowMap(), xtmp);
    Y2 = new Epetra_Vector(View, L_->RowMap(), ytmp);
  }
  else {
    X2 = X1;
    Y2 = Y1;
  }

  Epetra_Flops * counter = this->GetFlopCounter();
  if (counter!=0) {
    L_->SetFlopCounter(*counter);
    Y2->SetFlopCounter(*counter);
    U_->SetFlopCounter(*counter);
  }

  if (!Trans) {

    EPETRA_CHK_ERR(L_->Solve(Lower, Trans, UnitDiagonal, *X2, *Y2));
    EPETRA_CHK_ERR(Y2->Multiply(1.0, *D_, *Y2, 0.0)); // y = D*y (D_ has inverse of diagonal)
    EPETRA_CHK_ERR(U_->Solve(Upper, Trans, UnitDiagonal, *Y2, *Y2)); // Solve Uy = y
  }
  else
    {
      EPETRA_CHK_ERR(U_->Solve(Upper, Trans, UnitDiagonal, *X2, *Y2)); // Solve Uy = y
      EPETRA_CHK_ERR(Y2->Multiply(1.0, *D_, *Y2, 0.0)); // y = D*y (D_ has inverse of diagonal)
      EPETRA_CHK_ERR(L_->Solve(Lower, Trans, UnitDiagonal, *Y2, *Y2));
      
    } 

  if (UserMatrixIsVbr_) {
    delete X2;
    delete Y2;
  }

  // Export computed Y values as directed
  if (IsOverlapped_) {
    EPETRA_CHK_ERR(Y.Export(*OverlapY_,*Graph_.OverlapImporter(), OverlapMode_));
  }
  return(0);
}


//=============================================================================
int Ifpack_CrsRiluk::Solve(bool Trans, const Epetra_MultiVector& X, 
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

  if (IsOverlapped_) {
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
    EPETRA_CHK_ERR(OverlapX_->Import(X,*Graph_.OverlapImporter(), Insert)); // Import X values for solve
    X1 = OverlapX_;
    Y1 = OverlapY_; // Set pointers for X1 and Y1 to point to overlap space
  }

  // This code block checks to see if the user matrix, and therefore the user vectors, overlap graph, etc. 
  // are using a block map.  If so, we want make a view of these vectors using a map that is compatible with
  // L and U (which have simple maps).

  Epetra_MultiVector * X2 = 0;
  Epetra_MultiVector * Y2 = 0;

  if (UserMatrixIsVbr_) {
    double ** xtmp, ** ytmp;
    X1->ExtractView(&xtmp);
    Y1->ExtractView(&ytmp);
    X2 = new Epetra_MultiVector(View, L_->RowMap(), xtmp, X.NumVectors());
    Y2 = new Epetra_MultiVector(View, L_->RowMap(), ytmp, X.NumVectors());
  }
  else {
    X2 = X1;
    Y2 = Y1;
  }

  Epetra_Flops * counter = this->GetFlopCounter();
  if (counter!=0) {
    L_->SetFlopCounter(*counter);
    Y2->SetFlopCounter(*counter);
    U_->SetFlopCounter(*counter);
  }

  if (!Trans) {

    EPETRA_CHK_ERR(L_->Solve(Lower, Trans, UnitDiagonal, *X2, *Y2));
    EPETRA_CHK_ERR(Y2->Multiply(1.0, *D_, *Y2, 0.0)); // y = D*y (D_ has inverse of diagonal)
    EPETRA_CHK_ERR(U_->Solve(Upper, Trans, UnitDiagonal, *Y2, *Y2)); // Solve Uy = y
  }
  else
    {
      EPETRA_CHK_ERR(U_->Solve(Upper, Trans, UnitDiagonal, *X2, *Y2)); // Solve Uy = y
      EPETRA_CHK_ERR(Y2->Multiply(1.0, *D_, *Y2, 0.0)); // y = D*y (D_ has inverse of diagonal)
      EPETRA_CHK_ERR(L_->Solve(Lower, Trans, UnitDiagonal, *Y2, *Y2));
      
    } 

  if (UserMatrixIsVbr_) {
    delete X2;
    delete Y2;
  }

  // Export computed Y values as directed
  if (IsOverlapped_) {
    EPETRA_CHK_ERR(Y.Export(*OverlapY_,*Graph_.OverlapImporter(), OverlapMode_));
  }
  return(0);
}
//=============================================================================
int Ifpack_CrsRiluk::Multiply(bool Trans, const Epetra_MultiVector& X, 
				Epetra_MultiVector& Y) const {
//
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
//

  if (X.NumVectors()!=Y.NumVectors()) EPETRA_CHK_ERR(-1); // Return error: X and Y not the same size


  Epetra_MultiVector * X1 = (Epetra_MultiVector *) &X;
  Epetra_MultiVector * Y1 = (Epetra_MultiVector *) &Y;

  if (IsOverlapped_) {
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
    EPETRA_CHK_ERR(OverlapX_->Import(X,*Graph_.OverlapImporter(), Insert)); // Import X values for solve
    X1 = OverlapX_;
    Y1 = OverlapY_; // Set pointers for X1 and Y1 to point to overlap space
  }

  Epetra_MultiVector * X2 = 0;
  Epetra_MultiVector * Y2 = 0;

  if (UserMatrixIsVbr_) {
    double ** xtmp, ** ytmp;
    X1->ExtractView(&xtmp);
    Y1->ExtractView(&ytmp);
    X2 = new Epetra_MultiVector(View, L_->RowMap(), xtmp, X.NumVectors());
    Y2 = new Epetra_MultiVector(View, L_->RowMap(), ytmp, X.NumVectors());
  }
  else {
    X2 = X1;
    Y2 = Y1;
  }

  Epetra_Flops * counter = this->GetFlopCounter();
  if (counter!=0) {
    L_->SetFlopCounter(*counter);
    Y2->SetFlopCounter(*counter);
    U_->SetFlopCounter(*counter);
  }

  if (Trans) {

    EPETRA_CHK_ERR(L_->Multiply(Trans, *X2, *Y2));
    EPETRA_CHK_ERR(Y2->Update(1.0, *X2, 1.0)); // Y2 = Y2 + X2 (account for implicit unit diagonal)
    EPETRA_CHK_ERR(Y2->ReciprocalMultiply(1.0, *D_, *Y2, 0.0)); // y = D*y (D_ has inverse of diagonal)
    Epetra_MultiVector Y2temp(*Y2); // Need a temp copy of Y2
    EPETRA_CHK_ERR(U_->Multiply(Trans, Y2temp, *Y2));
    EPETRA_CHK_ERR(Y2->Update(1.0, Y2temp, 1.0)); // (account for implicit unit diagonal)
  }
  else {
    EPETRA_CHK_ERR(U_->Multiply(Trans, *X2, *Y2)); // 
    EPETRA_CHK_ERR(Y2->Update(1.0, *X2, 1.0)); // Y2 = Y2 + X2 (account for implicit unit diagonal)
    EPETRA_CHK_ERR(Y2->ReciprocalMultiply(1.0, *D_, *Y2, 0.0)); // y = D*y (D_ has inverse of diagonal)
    Epetra_MultiVector Y2temp(*Y2); // Need a temp copy of Y2
    EPETRA_CHK_ERR(L_->Multiply(Trans, Y2temp, *Y2));
    EPETRA_CHK_ERR(Y2->Update(1.0, Y2temp, 1.0)); // (account for implicit unit diagonal)
    } 

  if (UserMatrixIsVbr_) {
    delete X2;
    delete Y2;
  }

  // Export computed Y values as directed
  if (IsOverlapped_) {
    EPETRA_CHK_ERR(Y.Export(*OverlapY_,*Graph_.OverlapImporter(), OverlapMode_));
  }
  return(0);
}
//=============================================================================
int Ifpack_CrsRiluk::Condest(bool Trans, double & ConditionNumberEstimate) const {

  if (Condest_>=0.0) {
    ConditionNumberEstimate = Condest_;
    return(0);
  }
  // Create a vector with all values equal to one
  Epetra_Vector Ones(Graph_.DomainMap());
  Epetra_Vector OnesResult(Graph_.RangeMap());
  Ones.PutScalar(1.0);

  EPETRA_CHK_ERR(Solve(Trans, Ones, OnesResult)); // Compute the effect of the solve on the vector of ones
  EPETRA_CHK_ERR(OnesResult.Abs(OnesResult)); // Make all values non-negative
  EPETRA_CHK_ERR(OnesResult.MaxValue(&ConditionNumberEstimate)); // Get the maximum value across all processors
  Condest_ = ConditionNumberEstimate; // Save value for possible later calls
  return(0);
}
//=============================================================================
// Non-member functions

ostream& operator << (ostream& os, const Ifpack_CrsRiluk& A)
{
  Epetra_fmtflags olda = os.setf(ios::right,ios::adjustfield);
  Epetra_fmtflags oldf = os.setf(ios::scientific,ios::floatfield);
  int oldp = os.precision(12);
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

  os.setf(olda,ios::adjustfield);
  os.setf(oldf,ios::floatfield);
  os.precision(oldp);

  return os;
}
