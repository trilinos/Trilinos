#include "Ifpack_OverlapSolveObject.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Flops.h"

//==============================================================================
Ifpack_OverlapSolveObject::Ifpack_OverlapSolveObject(char * Label, const Epetra_Comm & Comm) 
  : Label_(Label),
    L_(0),
    UseLTrans_(false),
    D_(0),
    U_(0),
    UseUTrans_(false),
    UseTranspose_(false),
    Comm_(Comm),
    Condest_(-1.0),
    OverlapMode_(Zero)
{
}

//==============================================================================
Ifpack_OverlapSolveObject::Ifpack_OverlapSolveObject(const Ifpack_OverlapSolveObject & Source) 
  : Label_(Source.Label_),
    L_(Source.L_),
    UseLTrans_(Source.UseLTrans_),
    D_(Source.D_),
    U_(Source.U_),
    UseUTrans_(Source.UseUTrans_),
    UseTranspose_(Source.UseTranspose_),
    Comm_(Source.Comm_),
    Condest_(Source.Condest_),
    OverlapMode_(Source.OverlapMode_)
{
}
//==============================================================================
Ifpack_OverlapSolveObject::~Ifpack_OverlapSolveObject(){
}
//==============================================================================
//=============================================================================
int Ifpack_OverlapSolveObject::Solve(bool Trans, const Epetra_MultiVector& X, 
				Epetra_MultiVector& Y) const {
//
// This function finds Y such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
//

  // First generate X and Y as needed for this function
  Epetra_MultiVector * X1 = 0;
  Epetra_MultiVector * Y1 = 0;
  EPETRA_CHK_ERR(SetupXY(Trans, X, Y, X1, Y1));

  bool Upper = true;
  bool Lower = false;
  bool UnitDiagonal = true;

  if (!Trans) {

    EPETRA_CHK_ERR(L_->Solve(Lower, Trans, UnitDiagonal, *X1, *Y1));
    EPETRA_CHK_ERR(Y1->Multiply(1.0, *D_, *Y1, 0.0)); // y = D*y (D_ has inverse of diagonal)
    EPETRA_CHK_ERR(U_->Solve(Upper, Trans, UnitDiagonal, *Y1, *Y1)); // Solve Uy = y
    if (L_->Exporter()!=0) {EPETRA_CHK_ERR(Y.Export(*Y1,*L_->Exporter(), OverlapMode_));} // Export computed Y values if needed
  }
  else {
    EPETRA_CHK_ERR(U_->Solve(Upper, Trans, UnitDiagonal, *X1, *Y1)); // Solve Uy = y
    EPETRA_CHK_ERR(Y1->Multiply(1.0, *D_, *Y1, 0.0)); // y = D*y (D_ has inverse of diagonal)
    EPETRA_CHK_ERR(L_->Solve(Lower, Trans, UnitDiagonal, *Y1, *Y1));
    if (U_->Importer()!=0) {EPETRA_CHK_ERR(Y.Export(*Y1,*U_->Importer(), OverlapMode_));} // Export computed Y values if needed
  } 

  return(0);
}
//=============================================================================
int Ifpack_OverlapSolveObject::Multiply(bool Trans, const Epetra_MultiVector& X, 
				Epetra_MultiVector& Y) const {
//
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
//
    
  // First generate X and Y as needed for this function
  Epetra_MultiVector * X1 = 0;
  Epetra_MultiVector * Y1 = 0;
  EPETRA_CHK_ERR(SetupXY(Trans, X, Y, X1, Y1));

  if (!Trans) {
    EPETRA_CHK_ERR(U_->Multiply(Trans, *X1, *Y1)); // 
    EPETRA_CHK_ERR(Y1->Update(1.0, *X1, 1.0)); // Y1 = Y1 + X1 (account for implicit unit diagonal)
    EPETRA_CHK_ERR(Y1->ReciprocalMultiply(1.0, *D_, *Y1, 0.0)); // y = D*y (D_ has inverse of diagonal)
    Epetra_MultiVector Y1temp(*Y1); // Need a temp copy of Y1
    EPETRA_CHK_ERR(L_->Multiply(Trans, Y1temp, *Y1));
    EPETRA_CHK_ERR(Y1->Update(1.0, Y1temp, 1.0)); // (account for implicit unit diagonal)
    if (L_->Exporter()!=0) {EPETRA_CHK_ERR(Y.Export(*Y1,*L_->Exporter(), OverlapMode_));} // Export computed Y values if needed
  }
  else {

    EPETRA_CHK_ERR(L_->Multiply(Trans, *X1, *Y1));
    EPETRA_CHK_ERR(Y1->Update(1.0, *X1, 1.0)); // Y1 = Y1 + X1 (account for implicit unit diagonal)
    EPETRA_CHK_ERR(Y1->ReciprocalMultiply(1.0, *D_, *Y1, 0.0)); // y = D*y (D_ has inverse of diagonal)
    Epetra_MultiVector Y1temp(*Y1); // Need a temp copy of Y1
    EPETRA_CHK_ERR(U_->Multiply(Trans, Y1temp, *Y1));
    EPETRA_CHK_ERR(Y1->Update(1.0, Y1temp, 1.0)); // (account for implicit unit diagonal)
    if (L_->Exporter()!=0) {EPETRA_CHK_ERR(Y.Export(*Y1,*L_->Exporter(), OverlapMode_));}
  } 
  return(0);
} 
//=========================================================================
int Ifpack_OverlapSolveObject::SetupXY(bool Trans, 
				       const Epetra_MultiVector& Xin, const Epetra_MultiVector& Yin,
				       Epetra_MultiVector * & Xout, Epetra_MultiVector * & Yout) const {

  // Generate an X and Y suitable for performing Solve() and Multiply() methods
  
  if (Xin.NumVectors()!=Yin.NumVectors()) EPETRA_CHK_ERR(-1); // Return error: X and Y not the same size
  
  //cout << "Xin = " << Xin << endl;
  Xout = (Epetra_MultiVector *) &Xin;
  Yout = (Epetra_MultiVector *) &Yin;
  return(0);
}
//=============================================================================
int Ifpack_OverlapSolveObject::Condest(bool Trans, double & ConditionNumberEstimate) const {

  if (Condest_>=0.0) {
    ConditionNumberEstimate = Condest_;
    return(0);
  }
  // Create a vector with all values equal to one
  Epetra_Vector Ones(U_->DomainMap());
  Epetra_Vector OnesResult(L_->RangeMap());
  Ones.PutScalar(1.0);

  EPETRA_CHK_ERR(Solve(Trans, Ones, OnesResult)); // Compute the effect of the solve on the vector of ones
  EPETRA_CHK_ERR(OnesResult.Abs(OnesResult)); // Make all values non-negative
  EPETRA_CHK_ERR(OnesResult.MaxValue(&ConditionNumberEstimate)); // Get the maximum value across all processors
  Condest_ = ConditionNumberEstimate; // Save value for possible later calls
  return(0);
}
