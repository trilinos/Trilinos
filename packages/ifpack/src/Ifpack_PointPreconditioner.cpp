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
Ifpack_PointPreconditioner(const Epetra_RowMatrix* Matrix) :
  Matrix_(Matrix),
  UseTranspose_(false),
  NumSweeps_(1),
  DampingFactor_(1.0),
  PrintLevel_(0),
  Diagonal_(0),
  IsComputed_(false)
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
  SetNumSweeps(List.get("sweeps",NumSweeps()));
  SetDampingFactor(List.get("damping factor", DampingFactor()));
  PrintLevel_ = List.get("print level", PrintLevel());

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
}

//==============================================================================
int Ifpack_PointPreconditioner::Compute()
{

  if (Matrix()->NumGlobalRows() != Matrix()->NumGlobalCols())
    IFPACK_CHK_ERR(-3); // only square matrices

  if (NumSweeps() <= 0)
    IFPACK_CHK_ERR(-3); // at least one application
  
  Diagonal_ = new Epetra_Vector(Matrix()->RowMatrixRowMap());

  if (Diagonal_ == 0)
    IFPACK_CHK_ERR(-11);

  IFPACK_CHK_ERR(Matrix()->ExtractDiagonalCopy(*Diagonal_));

  IsComputed_ = true;

  return(0);
}

