#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_PointPreconditioner.h"

//==============================================================================
Ifpack_PointPreconditioner::
Ifpack_PointPreconditioner(Epetra_RowMatrix* Matrix) :
  Matrix_(Matrix),
  UseTranspose_(false),
  NumSweeps_(1),
  DampingFactor_(1.0),
  PrintFrequency_(0),
  Diagonal_(0),
  IsComputed_(false),
  Condest_(-1.0),
  CondestMaxIters_(1550),
  CondestTol_(1e-9),
  ZeroStartingSolution_(true)
{
}

//==============================================================================
Ifpack_PointPreconditioner::~Ifpack_PointPreconditioner()
{
  if (Diagonal_)
    delete Diagonal_;
}

//==============================================================================
#ifdef HAVE_IFPACK_TEUCHOS
int Ifpack_PointPreconditioner::SetParameters(Teuchos::ParameterList& List)
{

  SetNumSweeps(List.get("point: sweeps",NumSweeps()));
  SetDampingFactor(List.get("point: damping factor", DampingFactor()));
  SetPrintFrequency(List.get("point: print frequency", PrintFrequency()));
  ZeroStartingSolution_ = List.get("point: zero starting solution", 
				   ZeroStartingSolution_);

  CondestMaxIters_ = List.get("condest: max iters", CondestMaxIters_);
  CondestTol_ = List.get("condest: tolerance", CondestTol_);

  SetLabel();

  return(0);
}
#endif

//==============================================================================
const Epetra_Comm& Ifpack_PointPreconditioner::Comm() const
{
  return(Matrix_->Comm());
}

//==============================================================================
const Epetra_Map& Ifpack_PointPreconditioner::OperatorDomainMap() const
{
  return(Matrix_->OperatorDomainMap());
}

//==============================================================================
const Epetra_Map& Ifpack_PointPreconditioner::OperatorRangeMap() const
{
  return(Matrix_->OperatorRangeMap());
}

//==============================================================================
int Ifpack_PointPreconditioner::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

  if (IsComputed() == false)
    IFPACK_CHK_ERR(-4);

  // FIXME: other checks with OperatorDomainMap() and similar
  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-3);

  IFPACK_CHK_ERR(Matrix_->Multiply(UseTranspose(),X,Y));
  return(0);
}

//==============================================================================
int Ifpack_PointPreconditioner::Compute()
{
  if (!IsInitialized())
    IFPACK_CHK_ERR(Initialize());

  if (Matrix_ == 0)
    IFPACK_CHK_ERR(-1);

  if (Matrix().NumGlobalRows() != Matrix().NumGlobalCols())
    IFPACK_CHK_ERR(-3); // only square matrices

  if (NumSweeps() <= 0)
    IFPACK_CHK_ERR(-3); // at least one application
  
  Diagonal_ = new Epetra_Vector(Matrix().RowMatrixRowMap());

  if (Diagonal_ == 0)
    IFPACK_CHK_ERR(-11);

  IFPACK_CHK_ERR(Matrix().ExtractDiagonalCopy(*Diagonal_));

  IsComputed_ = true;

  return(0);
}

//==============================================================================
ostream& Ifpack_PointPreconditioner::Print(ostream & os) const
{

  double MinVal, MeanVal, MaxVal;

  if (IsComputed()) {
    Diagonal_->MinValue(&MinVal);
    Diagonal_->MeanValue(&MeanVal);
    Diagonal_->MaxValue(&MaxVal);
  }

  if (Comm().MyPID())
    return(os);

  os << endl;
  os << "*** " << Label() << endl << endl;
  os << "*** " << "Number of global rows = " 
       << Matrix().NumGlobalRows() << endl;
  os << "*** " << "Number of global cols = " 
       << Matrix().NumGlobalCols() << endl;
  os << "*** " << "Print frequency = " << PrintFrequency() << endl;
  os << "*** " << "IsComputed()    = " << IsComputed() << endl;
  os << "*** " << "Use zero starting solution  = " 
       << ZeroStartingSolution_ << endl;
  if (IsComputed()) {
    os << "*** " << "Minimum value on diagonal = " << MinVal << endl;
    os << "*** " << "Maximum value on diagonal = " << MaxVal << endl;
    os << "*** " << "Average value on diagonal = " << MeanVal << endl;
    os << "*** Note: Jacobi and Gauss-Seidel reported values refer" << endl;
    os << "***       to the inverse of the diagonal" << endl;
  }
  os << endl;

  return(os);
}

#include "Ifpack_Condest.h"
//==============================================================================
double Ifpack_PointPreconditioner::
Condest(const Ifpack_CondestType CT, 
	Epetra_RowMatrix* Matrix)
{
  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  if (Condest_ == -1.0)
    Condest_ = Ifpack_Condest(*this, CT, CondestMaxIters_, CondestTol_,
			      Matrix);

  return(Condest_);
}
