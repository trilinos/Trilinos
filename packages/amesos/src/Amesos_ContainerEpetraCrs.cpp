#include "Amesos_ConfigDefs.h"
#include "Amesos_Container.h"
#include "Amesos_ContainerEpetraCrs.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

//==============================================================================
Amesos_ContainerEpetraCrs::Amesos_ContainerEpetraCrs() :
  NumRows_(-1),
  NumVectors_(-1),
  Map_(0),
  Matrix_(0),
  Inverse_(0),
  LHS_(0),
  RHS_(0),
  Solver_(0),
  IsProblemShaped_(false),
  IsProblemComputed_(false),
  SerialComm_(0)
{

#ifdef HAVE_MPI
  SerialComm_ = new Epetra_MpiComm(MPI_COMM_SELF);
#else
  SerialComm_ = new Epetra_SerialComm;
#endif

}

//==============================================================================
Amesos_ContainerEpetraCrs::~Amesos_ContainerEpetraCrs()
{
  Destroy();
}

//==============================================================================
int Amesos_ContainerEpetraCrs::NumRows() const
{
  if (IsProblemShaped() == false)
    return(0);
  else
    return(NumRows_);
}

//==============================================================================
int Amesos_ContainerEpetraCrs::Shape(const int NumRows, const int NumVectors) 
{
  
  if (IsProblemShaped() == true)
    AMESOS_CHK_ERR(-1); // already computed

  IsProblemShaped_ = true;

  NumRows_ = NumRows;
  NumVectors_ = NumVectors;

  Map_ = new Epetra_Map(NumRows_,0,*SerialComm_);

  LHS_ = new Epetra_MultiVector(*Map_,NumVectors);
  RHS_ = new Epetra_MultiVector(*Map_,NumVectors);

  // FIXME: try View??
  Matrix_ = new Epetra_CrsMatrix(Copy,*Map_,0);
  Solver_ = new Epetra_LinearProblem;
  
  Solver_->SetOperator(Matrix_);
  Solver_->SetLHS(LHS_);
  Solver_->SetRHS(RHS_);

  return(0);
  
}

//==============================================================================
int Amesos_ContainerEpetraCrs::Reshape(const int NumRows, const int NumVector) 
{

  if (IsProblemShaped() == true)
    Destroy();

  Shape(NumRows);
  
}

//==============================================================================
double& Amesos_ContainerEpetraCrs::LHS(const int i, const int Vector)
{
  return(((*LHS_)(Vector))->Values()[i]);
}
  
//==============================================================================
double& Amesos_ContainerEpetraCrs::RHS(const int i, const int Vector)
{
  return(((*RHS_)(Vector))->Values()[i]);
}

//==============================================================================
int Amesos_ContainerEpetraCrs::
SetMatrixElement(const int row, const int col, double value)
{
  if (IsProblemShaped() == false)
    AMESOS_CHK_ERR(-5); // problem not shaped yet

  if ((row < 0) || (row >= NumRows())) {
    AMESOS_CHK_ERR(-2); // not in range
  }

  if ((col < 0) || (col >= NumRows())) {
    AMESOS_CHK_ERR(-2); // not in range
  }

  int ierr = Matrix_->InsertGlobalValues((int)row,1,(double*)&value,(int*)&col);
  if (ierr < 0) {
    ierr = Matrix_->SumIntoGlobalValues((int)row,1,(double*)&value,(int*)&col);
    if (ierr < 0)
      AMESOS_CHK_ERR(-1);
  }

  return(0);

}

//==============================================================================
int Amesos_ContainerEpetraCrs::Compute()
{
  if (IsProblemComputed() == true) {
    AMESOS_CHK_ERR(-6); // already computed
  }

  IsProblemComputed_ = true;

  AMESOS_CHK_ERR(Matrix_->FillComplete());

  return(0);
}

//==============================================================================
int Amesos_ContainerEpetraCrs::Apply()
{
  if (IsProblemShaped() == false) {
    AMESOS_CHK_ERR(-5); // not yet shaped
  }

  if (IsProblemComputed() == false) {
    AMESOS_CHK_ERR(-6); // not yet computed
  }
  
  AMESOS_CHK_ERR(Matrix_->Apply(*RHS_, *LHS_));

  return(0);
}

//==============================================================================
int Amesos_ContainerEpetraCrs::ApplyInverse()
{
  if (IsProblemShaped() == false) {
    AMESOS_CHK_ERR(-5); // not yet shaped
  }

  if (IsProblemComputed() == false) {
    AMESOS_CHK_ERR(-6); // not yet computed
  }
  
  if (Inverse_ == 0) {
    AMESOS_CHK_ERR(-7);
  }

  AMESOS_CHK_ERR(Inverse_->ApplyInverse(*RHS_, *LHS_));

  return(0);
}
 

//==============================================================================
int Amesos_ContainerEpetraCrs::Destroy()
{
  if (Map_)
    delete Map_;

  if (Matrix_)
    delete Matrix_;

  if (LHS_)
    delete LHS_;

  if (RHS_)
    delete RHS_;

  if (Solver_)
    delete Solver_;

  Map_ = 0;
  Matrix_ = 0;
  Inverse_ = 0;
  LHS_ = 0;
  RHS_ = 0;
  Solver_ = 0;

  IsProblemShaped_ = false;
  IsProblemComputed_ = false;
  return(0);
}

//==============================================================================
int Amesos_ContainerEpetraCrs::GetMatrixPointer(void** Matrix)
{
  if (IsProblemShaped() == false) {
    AMESOS_CHK_ERR(-1);
  }

  if (IsProblemComputed() == false) {
    AMESOS_CHK_ERR(-1);
  }

  *Matrix = (void*) Matrix_;

  return(0);
}

//==============================================================================
int Amesos_ContainerEpetraCrs::SetInversePointer(void* Inverse)
{
  Inverse_ = (Epetra_Operator*)Inverse;

  return(0);
}


