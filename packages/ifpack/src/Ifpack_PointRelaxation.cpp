#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_PointRelaxation.h"
#include "Ifpack_Utils.h"
#include "Ifpack_Condest.h"

static const int IFPACK_JACOBI = 0;
static const int IFPACK_GS = 1;
static const int IFPACK_SGS = 2;
static const int IFPACK_SOR = 3;
static const int IFPACK_SSOR = 4;

//==============================================================================
Ifpack_PointRelaxation::
Ifpack_PointRelaxation(const Epetra_RowMatrix* Matrix) :
  Diagonal_(0),
  ZeroStartingSolution_(true),
  IsInitialized_(false),
  IsComputed_(false),
  NumSweeps_(1),
  DampingFactor_(1.0),
  NumMyRows_(0),
  Matrix_(Matrix),
  UseTranspose_(false),
  PrintFrequency_(0),
  Condest_(-1.0),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyInverseTime_(0.0),
  Time_(0),
  ComputeFlops_(0),
  ApplyInverseFlops_(0),
  PrecType_(IFPACK_JACOBI),
  MinDiagonalValue_(0.0)
{
}

//==============================================================================
Ifpack_PointRelaxation::
Ifpack_PointRelaxation(const Ifpack_PointRelaxation& rhs) :
  Diagonal_(0),
  ZeroStartingSolution_(rhs.ZeroStartingSolution()),
  IsInitialized_(rhs.IsInitialized()),
  IsComputed_(rhs.IsComputed()),
  NumSweeps_(rhs.NumSweeps()),
  DampingFactor_(rhs.DampingFactor()),
  NumMyRows_(Matrix_->NumMyRows()),
  Matrix_(&rhs.Matrix()),
  UseTranspose_(rhs.UseTranspose()),
  PrintFrequency_(rhs.PrintFrequency()),
  Condest_(rhs.Condest()),
  NumInitialize_(rhs.NumInitialize()),
  NumCompute_(rhs.NumCompute()),
  NumApplyInverse_(rhs.NumApplyInverse()),
  InitializeTime_(rhs.InitializeTime()),
  ComputeTime_(rhs.ComputeTime()),
  ApplyInverseTime_(rhs.ApplyInverseTime()),
  Time_(0),
  ComputeFlops_(rhs.ComputeFlops()),
  ApplyInverseFlops_(rhs.ApplyInverseFlops()),
  PrecType_(rhs.PrecType()),
  MinDiagonalValue_(rhs.MinDiagonalValue())
{
  Diagonal_ = new Epetra_Vector((*rhs.Diagonal()));
}

//==============================================================================
Ifpack_PointRelaxation::~Ifpack_PointRelaxation()
{
  if (Diagonal_)
    delete Diagonal_;
  if (Time_)
    delete Time_;
}

//==============================================================================
Ifpack_PointRelaxation& 
Ifpack_PointRelaxation::operator=(const Ifpack_PointRelaxation& rhs)
{
  if (this == &rhs)
    return(*this);

  if (Diagonal_)
    delete Diagonal_;

  Matrix_ = &rhs.Matrix();
  UseTranspose_ = rhs.UseTranspose();
  NumSweeps_ = rhs.NumSweeps();
  DampingFactor_ = rhs.DampingFactor();
  PrintFrequency_ = rhs.PrintFrequency();
  IsInitialized_ = rhs.IsInitialized();
  IsComputed_ = rhs.IsComputed();
  Condest_ = rhs.Condest();
  ZeroStartingSolution_ = rhs.ZeroStartingSolution();
  PrecType_ = rhs.PrecType();
  NumMyRows_ = Matrix_->NumMyRows();
  MinDiagonalValue_ = rhs.MinDiagonalValue();
  Diagonal_ = new Epetra_Vector((*rhs.Diagonal()));

  NumInitialize_ = rhs.NumInitialize();
  NumCompute_ = rhs.NumCompute();
  NumApplyInverse_ = rhs.NumApplyInverse();
  InitializeTime_ = rhs.InitializeTime();
  ComputeTime_ = rhs.ComputeTime();
  ApplyInverseTime_ = rhs.ApplyInverseTime();
  ComputeFlops_ = rhs.ComputeFlops();
  ApplyInverseFlops_ = rhs.ApplyInverseFlops();

  return(*this);
}

//==============================================================================
#ifdef HAVE_IFPACK_TEUCHOS
int Ifpack_PointRelaxation::SetParameters(Teuchos::ParameterList& List)
{

  string PT;
  if (PrecType_ == IFPACK_JACOBI)
    PT = "Jacobi";
  else if (PrecType_ == IFPACK_GS)
    PT = "Gauss-Seidel";
  else if (PrecType_ == IFPACK_SGS)
    PT = "symmetric Gauss-Seidel";
  else if (PrecType_ == IFPACK_SOR)
    PT = "SOR";
  else if (PrecType_ == IFPACK_SSOR)
    PT = "SSOR";

  PT = List.get("point: type", PT);

  if (PT == "Jacobi")
    PrecType_ = IFPACK_JACOBI;
  else if (PT == "Gauss-Seidel")
    PrecType_ = IFPACK_GS;
  else if (PT == "symmetric Gauss-Seidel")
    PrecType_ = IFPACK_SGS;
  else if (PT == "SOR")
    PrecType_ = IFPACK_SOR;
  else if (PT == "SSOR")
    PrecType_ = IFPACK_SSOR;
  
  NumSweeps_ = List.get("point: sweeps",NumSweeps_);
  DampingFactor_ = List.get("point: damping factor", DampingFactor_);
  PrintFrequency_ = List.get("point: print frequency", PrintFrequency_);
  ZeroStartingSolution_ = List.get("point: zero starting solution", 
				   ZeroStartingSolution_);
  MinDiagonalValue_ = List.get("point: min diagonal value", MinDiagonalValue_);

  SetLabel();

  return(0);
}
#endif

//==============================================================================
const Epetra_Comm& Ifpack_PointRelaxation::Comm() const
{
  return(Matrix_->Comm());
}

//==============================================================================
const Epetra_Map& Ifpack_PointRelaxation::OperatorDomainMap() const
{
  return(Matrix_->OperatorDomainMap());
}

//==============================================================================
const Epetra_Map& Ifpack_PointRelaxation::OperatorRangeMap() const
{
  return(Matrix_->OperatorRangeMap());
}

//==============================================================================
int Ifpack_PointRelaxation::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

  if (IsComputed() == false)
    IFPACK_CHK_ERR(-4);

  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-3);

  IFPACK_CHK_ERR(Matrix_->Multiply(UseTranspose(),X,Y));
  return(0);
}

//==============================================================================
int Ifpack_PointRelaxation::Initialize()
{
  IsInitialized_ = false;

  if (Matrix_ == 0)
    IFPACK_CHK_ERR(-1);

  if (Time_ == 0)
    Time_ = new Epetra_Time(Comm());

  if (Comm().NumProc() != 1) {
    cerr << "This class must be used with communicators containing" << endl;
    cerr << "only one process. More general uses are allowed through" << endl;
    cerr << "class Ifpack_AdditiveSchwarz" << endl;
    exit(EXIT_FAILURE);
  }

  if (Matrix().NumGlobalRows() != Matrix().NumGlobalCols())
    IFPACK_CHK_ERR(-3); // only square matrices

  NumMyRows_ = Matrix_->NumMyRows();

  ++NumInitialize_;
  InitializeTime_ += Time_->ElapsedTime();
  IsInitialized_ = true;
  return(0);
}

//==============================================================================
int Ifpack_PointRelaxation::Compute()
{
  if (!IsInitialized())
    IFPACK_CHK_ERR(Initialize());

  Time_->ResetStartTime();

  // reset values
  IsComputed_ = false;
  Condest_ = -1.0;

  if (NumSweeps() <= 0)
    IFPACK_CHK_ERR(-3); // at least one application
  
  if (Diagonal_)
    delete Diagonal_;
  Diagonal_ = new Epetra_Vector(Matrix().RowMatrixRowMap());

  if (Diagonal_ == 0)
    IFPACK_CHK_ERR(-11);

  IFPACK_CHK_ERR(Matrix().ExtractDiagonalCopy(*Diagonal_));

  // check diagonal elements, replace zeros with 1.0
  for (int i = 0 ; i < NumMyRows_ ; ++i) {
    double diag = (*Diagonal_)[i];
    if (IFPACK_ABS(diag) < MinDiagonalValue_)
      (*Diagonal_)[i] = MinDiagonalValue_;
  }

  // some methods require the inverse of the diagonal, compute it
  // now
  if ((PrecType_ == IFPACK_JACOBI) || (PrecType_ == IFPACK_GS) 
      || (PrecType_ == IFPACK_SOR)) {
    Diagonal_->Reciprocal(*Diagonal_);
  }

  ++NumCompute_;
  ComputeTime_ += Time_->ElapsedTime();
  IsComputed_ = true;

  return(0);
}

//==============================================================================
ostream& Ifpack_PointRelaxation::Print(ostream & os) const
{

  double MinVal, MeanVal, MaxVal;

  if (IsComputed_) {
    Diagonal_->MinValue(&MinVal);
    Diagonal_->MeanValue(&MeanVal);
    Diagonal_->MaxValue(&MaxVal);
  }

  os << endl;
  os << "*** " << Label() << endl << endl;
  os << "Number of rows    = " << NumMyRows_ << endl;
  os << "Number of sweeps  = " << NumSweeps_ << endl;
  os << "Damping Factor    = " << DampingFactor_ << endl;
  os << "Print frequency   = " << PrintFrequency_ << endl;
  os << "IsInitialized()   = " << IsInitialized_ << endl;
  os << "IsComputed()      = " << IsComputed_ << endl;
  if (IsComputed_) {
    os << "Minimum value on stored diagonal = " << MinVal << endl;
    os << "Maximum value on stored diagonal = " << MaxVal << endl;
    os << "Average value on stored diagonal = " << MeanVal << endl;
  }
  os << endl;

  return(os);
}

//==============================================================================
double Ifpack_PointRelaxation::
Condest(const Ifpack_CondestType CT, 
        const int MaxIters, const double Tol,
	Epetra_RowMatrix* Matrix)
{
  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  if (Condest_ == -1.0)
    Condest_ = Ifpack_Condest(*this, CT, MaxIters, Tol, Matrix);

  return(Condest_);
}

//==============================================================================
void Ifpack_PointRelaxation::SetLabel()
{
  string PT;
  if (PrecType_ == IFPACK_JACOBI)
    PT = "Jacobi";
  else if (PrecType_ == IFPACK_GS)
    PT = "Gauss-Seidel";
  else if (PrecType_ == IFPACK_SGS)
    PT = "sym Gauss-Seidel";
  else if (PrecType_ == IFPACK_SOR)
    PT = "SOR";
  else if (PrecType_ == IFPACK_SSOR)
    PT = "SSOR";

  sprintf(Label_, "IFPACK %s (sweeps=%d, damping=%f)",
          PT.c_str(), NumSweeps(), DampingFactor());
}

//==============================================================================
int Ifpack_PointRelaxation::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (IsComputed() == false)
    IFPACK_CHK_ERR(-4);

  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-3);

  Time_->ResetStartTime();

  switch (PrecType_) {
  case IFPACK_JACOBI:
    IFPACK_CHK_ERR(ApplyInverseJacobi(X,Y));
    break;
  case IFPACK_GS:
    IFPACK_CHK_ERR(ApplyInverseGS(X,Y));
    break;
  case IFPACK_SGS:
    IFPACK_CHK_ERR(ApplyInverseSGS(X,Y));
    break;
  case IFPACK_SOR:
    IFPACK_CHK_ERR(ApplyInverseSOR(X,Y));
    break;
  case IFPACK_SSOR:
    IFPACK_CHK_ERR(ApplyInverseSSOR(X,Y));
    break;
  }

  ++NumApplyInverse_;
  ApplyInverseTime_ += Time_->ElapsedTime();
  return(0);
}

//==============================================================================
int Ifpack_PointRelaxation::
ApplyInverseJacobi(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // ------------ //
  // single sweep //
  // ------------ //

  if (NumSweeps_ == 1 && ZeroStartingSolution_
      && (PrintFrequency() == 0)) {
    IFPACK_CHK_ERR(Y.Multiply(DampingFactor_,X,*Diagonal_,0.0));
    return(0);
  }

  // --------------------- //
  // general case (solver) //
  // --------------------- //
  // NOTE: I suppose that X and Y do NOT point to the same
  // memory. This means that this method will not work
  // directly with AztecOO. Users are not supposed to use
  // this class directly, but only through Ifpack_AdditiveSchwarz
  // (which handles the case of X and Y pointing to the same
  // memory space).

  Epetra_MultiVector AX(X);

  if (ZeroStartingSolution_)
    Y.PutScalar(0.0);

  if (PrintFrequency())
    Ifpack_PrintResidual(Label(),Matrix(),Y,X);

  for (int j = 0; j < NumSweeps_ ; j++) {

    // compute the residual
    if (j) {
      IFPACK_CHK_ERR(Apply(Y,AX));
      AX.Update(1.0,X,-1.0);
    }
    else {
      AX = X;
    }

    // apply the inverse of the diagonal
    AX.Multiply(1.0, AX, *Diagonal_, 0.0);

    // update the solution at step `j'
    Y.Update(DampingFactor_, AX, 1.0);

    if (PrintFrequency() && (j != 0) && (j % PrintFrequency() == 0))
      Ifpack_PrintResidual(j,Matrix(),Y,X);

  }

  if (PrintFrequency())
    Ifpack_PrintResidual(NumSweeps_,*Matrix_,Y,X);

  return(0);

}

//==============================================================================
int Ifpack_PointRelaxation::
ApplyInverseGS(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int NumVectors = X.NumVectors();

  int Length = Matrix().MaxNumEntries();
  vector<int> Indices;
  vector<double> Values;
  Indices.resize(Length);
  Values.resize(Length);

  // ============ //
  // single sweep //
  // ------------ //

  if ((NumSweeps() == 1)&& ZeroStartingSolution_
      && (PrintFrequency() != 0)) {

    for (int i = 0 ; i < NumMyRows_ ; ++i) {

      int NumEntries;
      int col;
      IFPACK_CHK_ERR(Matrix().ExtractMyRowCopy(i, Length,NumEntries,
						&Values[0], &Indices[0]));

      for (int m = 0 ; m < NumVectors ; ++m) {

	double dtemp = 0.0;
	for (int k = 0 ; k < NumEntries ; ++k) {

	  col = Indices[k];
	  if (col >= i || (col >= NumMyRows_)) continue;

	  dtemp += Values[k] * Y[m][col];
	}

	Y[m][i] = DampingFactor() * ((*Diagonal_)[i])
	  * (Y[m][i] - dtemp);
      }
    }

    return(0);
  }

  // --------------------- //
  // general case (solver) //
  // --------------------- //

  // need an additional vector for AztecOO preconditioning
  // (as X and Y both point to the same memory space)
  Epetra_MultiVector Xtmp(X);

  if (ZeroStartingSolution_)
    Y.PutScalar(0.0);
  
  if (PrintFrequency())
    Ifpack_PrintResidual(Label(),Matrix(),Y,Xtmp);

  for (int j = 0; j < NumSweeps() ; j++) {

    for (int i = 0 ; i < NumMyRows_ ; ++i) {

      int NumEntries;
      int col;
      IFPACK_CHK_ERR(Matrix().ExtractMyRowCopy(i, Length,NumEntries,
						&Values[0], &Indices[0]));

      for (int m = 0 ; m < NumVectors ; ++m) {

	double dtemp = 0.0;
	for (int k = 0 ; k < NumEntries ; ++k) {

	  col = Indices[k];

	  if ((col == i) || (col >= NumMyRows_)) continue;

	  dtemp += Values[k] * Y[m][col];
	}

	Y[m][i] = DampingFactor() * ((*Diagonal_)[i])
	  * (Xtmp[m][i] - dtemp);
      }
    }

    if (PrintFrequency() && (j != 0) && (j % PrintFrequency() == 0))
      Ifpack_PrintResidual(j,Matrix(),Y,Xtmp);
  }

  if (PrintFrequency())
    Ifpack_PrintResidual(NumSweeps(),Matrix(),Y,Xtmp);

  return(0);
}

//==============================================================================
int Ifpack_PointRelaxation::
ApplyInverseSGS(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

  // ------------ //
  // single sweep //
  // ------------ //
 
  if ((NumSweeps() == 1) && ZeroStartingSolution_
      && (PrintFrequency() != 0)) {
    IFPACK_CHK_ERR(ApplyInverseSGS2(Y));
    return(0);
  }
  
  // --------------------- //
  // general case (solver) //
  // --------------------- //
  
  Epetra_MultiVector Xtmp(X);
  Epetra_MultiVector AX(X);

  if (ZeroStartingSolution_)
    Y.PutScalar(0.0);

  if (PrintFrequency())
    Ifpack_PrintResidual(Label(),Matrix(),Y,Xtmp);
  
  for (int j = 0; j < NumSweeps() ; j++) {

    // compute the residual
    IFPACK_CHK_ERR(Apply(Y,AX));
    AX.Update(1.0,Xtmp,-1.0);

    // apply the lower triangular part of A
    IFPACK_CHK_ERR(ApplyInverseSGS2(AX));

    // update the residual
    Y.Update(DampingFactor(), AX, 1.0);

    if (PrintFrequency() && (j != 0) && (j % PrintFrequency() == 0))
      Ifpack_PrintResidual(j,Matrix(),Y,Xtmp);
    
  }

  if (PrintFrequency())
    Ifpack_PrintResidual(NumSweeps(),Matrix(),Y,Xtmp);

  return(0);

}

//==============================================================================
int Ifpack_PointRelaxation::
ApplyInverseSGS2(Epetra_MultiVector& X) const
{
  int NumVectors = X.NumVectors();
  int Length = Matrix().MaxNumEntries();
  vector<int> Indices;
  vector<double> Values;
  Indices.resize(Length);
  Values.resize(Length);

  // Apply (D - \omega E)^{-1}

  for (int i = 0 ; i < NumMyRows_ ; ++i) {

    int NumEntries;
    int col;
    double diag = 1.0;
    double dtemp = 0.0;

    IFPACK_CHK_ERR(Matrix().ExtractMyRowCopy(i, Length,NumEntries,
					     &Values[0], &Indices[0]));

    for (int m = 0 ; m < NumVectors ; ++m) {

      for (int k = 0 ; k < NumEntries ; ++k) {

	col = Indices[k];
	if (col > i) 
	  continue;

	if (col == i) 
	  diag = Values[k];
	else
	  dtemp += Values[k] * X[m][col];
      }

      X[m][i] = (X[m][i] - dtemp) / diag;

    }
  }

  // Apply D

  X.Multiply(1.0,X,*Diagonal_,0.0);

  // Apply (D - \omega F)^{-1}

  for (int i = NumMyRows_  - 1 ; i > -1 ; --i) {

    int NumEntries;
    int col;
    double dtemp = 0.0;

    IFPACK_CHK_ERR(Matrix().ExtractMyRowCopy(i, Length,NumEntries,
					     &Values[0], &Indices[0]));

    for (int m = 0 ; m < NumVectors ; ++m) {

      for (int k = 0 ; k < NumEntries ; ++k) {

	col = Indices[k];
	
	if (col <= i) 
	  continue;

	dtemp += Values[k] * X[m][col];
      }

      X[m][i] = (X[m][i] - dtemp) / ((*Diagonal_)[i]);

    }
  }
  return(0);

}

//==============================================================================
int Ifpack_PointRelaxation::
ApplyInverseSOR(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int NumVectors = X.NumVectors();

  int Length = Matrix().MaxNumEntries();
  vector<int> Indices;
  vector<double> Values;
  Indices.resize(Length);
  Values.resize(Length);

  // one iteration of SOR is as damped Gauss-Seidel.
  // Here I consider the general case only.

  Epetra_MultiVector Xtmp(X);

  if (ZeroStartingSolution_)  
    Y.PutScalar(0.0);

  if (PrintFrequency())
    Ifpack_PrintResidual(Label(),Matrix(),Y,Xtmp);

  for (int j = 0; j < NumSweeps() ; j++) {

    for (int i = 0 ; i < NumMyRows_ ; ++i) {

      int NumEntries;
      int col;
      IFPACK_CHK_ERR(Matrix().ExtractMyRowCopy(i, Length,NumEntries,
						&Values[0], &Indices[0]));

      for (int m = 0 ; m < NumVectors ; ++m) {

	double dtemp = 0.0;
	for (int k = 0 ; k < NumEntries ; ++k) {

	  col = Indices[k];

	  if ((col == i) || (col >= NumMyRows_)) continue;

	  dtemp += Values[k] * Y[m][col];
	}

	Y[m][i] = (1 - DampingFactor()) * Y[m][i] + 
	  DampingFactor() * ((*Diagonal_)[i]) * (Xtmp[m][i] - dtemp);
      }
    }

    if (PrintFrequency() && (j != 0) && (j % PrintFrequency() == 0))
      Ifpack_PrintResidual(j,Matrix(),Y,Xtmp);

  }

  if (PrintFrequency())
    Ifpack_PrintResidual(NumSweeps(),Matrix(),Y,Xtmp);

  return(0);
}

//==============================================================================
int Ifpack_PointRelaxation::
ApplyInverseSSOR(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

  // ---------------- //
  // single sweep can //
  // ---------------- //
 
  if ((NumSweeps() == 1) && ZeroStartingSolution_
      && (PrintFrequency() != 0)) {
    IFPACK_CHK_ERR(ApplyInverseSSOR2(Y));
    return(0);
  }
  
  // --------------------- //
  // general case (solver) //
  // --------------------- //
  
  Epetra_MultiVector Xtmp(X);
  Epetra_MultiVector AX(X);

  if (ZeroStartingSolution_)
    Y.PutScalar(0.0);

  if (PrintFrequency())
    Ifpack_PrintResidual(Label(),Matrix(),Y,Xtmp);
  
  for (int j = 0; j < NumSweeps() ; j++) {

    // compute the residual
    IFPACK_CHK_ERR(Apply(Y,AX));
    AX.Update(1.0,Xtmp,-1.0);

    // apply the lower triangular part of A
    IFPACK_CHK_ERR(ApplyInverseSSOR2(AX));

    // update the residual
    Y.Update(DampingFactor(), AX, 1.0);

    if (PrintFrequency() && (j != 0) && (j % PrintFrequency() == 0))
      Ifpack_PrintResidual(j,Matrix(),Y,Xtmp);
    
  }

  if (PrintFrequency())
    Ifpack_PrintResidual(NumSweeps(),Matrix(),Y,Xtmp);

  return(0);

}

//==============================================================================
int Ifpack_PointRelaxation::
ApplyInverseSSOR2(Epetra_MultiVector& X) const
{

  int NumVectors = X.NumVectors();
  int Length = Matrix().MaxNumEntries();
  vector<int> Indices;
  vector<double> Values;
  Indices.resize(Length);
  Values.resize(Length);

  // Apply (D - \omega E)^{-1}

  for (int i = 0 ; i < NumMyRows_ ; ++i) {

    int NumEntries;
    int col;
    double dtemp = 0.0;

    IFPACK_CHK_ERR(Matrix().ExtractMyRowCopy(i, Length,NumEntries,
					     &Values[0], &Indices[0]));

    for (int m = 0 ; m < NumVectors ; ++m) {

      for (int k = 0 ; k < NumEntries ; ++k) {

	col = Indices[k];
	if (col >= i) 
	  continue;

        dtemp += Values[k] * X[m][col];
      }

      X[m][i] = (X[m][i] - DampingFactor() * dtemp) / (*Diagonal_)[i];

    }
  }

  // Apply D

  X.Multiply(1.0,X,*Diagonal_,0.0);

  // Apply (D - \omega F)^{-1}

  for (int i = NumMyRows_  - 1 ; i > -1 ; --i) {

    int NumEntries;
    int col;
    double dtemp = 0.0;

    IFPACK_CHK_ERR(Matrix().ExtractMyRowCopy(i, Length,NumEntries,
					     &Values[0], &Indices[0]));

    for (int m = 0 ; m < NumVectors ; ++m) {

      for (int k = 0 ; k < NumEntries ; ++k) {

	col = Indices[k];
	
	if (col <= i) 
	  continue;

	dtemp += Values[k] * X[m][col];
      }

      X[m][i] = (X[m][i] - DampingFactor() * dtemp) / ((*Diagonal_)[i]);

    }
  }
  return(0);
}
