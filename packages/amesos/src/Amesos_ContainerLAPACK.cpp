#include "Amesos_ConfigDefs.h"
#include "Amesos_LocalLinearProblem.h"
#include "Amesos_ContainerLAPACK.h"
#include "Epetra_IntSerialDenseVector.h"

//==============================================================================
Amesos_ContainerLAPACK::Amesos_ContainerLAPACK() :
  NumRows_(-1),
  IsProblemShaped_(false),
  IsProblemComputed_(false)
{}


//==============================================================================
Amesos_ContainerLAPACK::~Amesos_ContainerLAPACK()
{
  Destroy();
}

//==============================================================================
int Amesos_ContainerLAPACK::NumRows() const
{
  if (IsProblemShaped() == false)
    return(0);
  else
    return(NumRows_);
}

//==============================================================================
int Amesos_ContainerLAPACK::NumVectors() const
{
  if (IsProblemShaped() == false)
    return(0);
  else
    return(NumVectors_);
}

//==============================================================================
int Amesos_ContainerLAPACK::Shape(const int NumRows, const int NumVectors) 
{
  
  if (IsProblemShaped() == true)
    AMESOS_CHK_ERR(-1); // already computed

  IsProblemShaped_ = true;

  NumRows_ = NumRows;
  NumVectors_ = NumVectors;

  AMESOS_CHK_ERR(LHS_.Reshape(NumRows_,NumVectors));
  AMESOS_CHK_ERR(RHS_.Reshape(NumRows_,NumVectors));
  AMESOS_CHK_ERR(GID_.Reshape(NumRows_,NumVectors));

  Matrix_.Reshape(NumRows_,NumRows_);
  // zero out matrix elements
  for (int i = 0 ; i < NumRows_ ; ++i)
    for (int j = 0 ; j < NumRows_ ; ++j)
      Matrix_(i,j) = 0.0;

  // zero out vector elements
  for (int i = 0 ; i < NumRows_ ; ++i)
    for (int j = 0 ; j < NumVectors_ ; ++j) {
      LHS_(i,j) = 0.0;
      RHS_(i,j) = 0.0;
    }

  AMESOS_CHK_ERR(Solver_.SetMatrix(Matrix_));
  AMESOS_CHK_ERR(Solver_.SetVectors(LHS_,RHS_));

  return(0);
  
}

//==============================================================================
int Amesos_ContainerLAPACK::Reshape(const int NumRows, const int NumVectors) 
{

  if (IsProblemShaped() == true)
    Destroy();

  Shape(NumRows,NumVectors);
  
}

//==============================================================================
double& Amesos_ContainerLAPACK::LHS(const int i, const int Vector)
{
  return(LHS_.A()[Vector * NumRows_ + i]);
}
  
//==============================================================================
double& Amesos_ContainerLAPACK::RHS(const int i, const int Vector)
{
  return(RHS_.A()[Vector * NumRows_ + i]);
}

//==============================================================================
int Amesos_ContainerLAPACK::
SetMatrixElement(const int row, const int col, const double value)
{
  if (IsProblemShaped() == false)
    AMESOS_CHK_ERR(-5); // problem not shaped yet

  if ((row < 0) || (row >= NumRows())) {
    AMESOS_CHK_ERR(-2); // not in range
  }

  if ((col < 0) || (col >= NumRows())) {
    AMESOS_CHK_ERR(-2); // not in range
  }

  Matrix_(row, col) = value;

  return(0);

}

//==============================================================================
int Amesos_ContainerLAPACK::ComputeInverse(char* Type,
					   Amesos_InverseFactory& Factory,
					   Teuchos::ParameterList& List)
{

  AMESOS_CHK_ERR(Solver_.Factor());

}

//==============================================================================
int Amesos_ContainerLAPACK::Compute()
{
  if (IsProblemComputed() == true) {
    AMESOS_CHK_ERR(-6); // already computed
  }

  IsProblemComputed_ = true;

  return(0);
}

//==============================================================================
int Amesos_ContainerLAPACK::Apply()
{
  AMESOS_CHK_ERR(-99);
}

//==============================================================================
int Amesos_ContainerLAPACK::ApplyInverse()
{
  if (IsProblemShaped() == false) {
    AMESOS_CHK_ERR(-5); // not yet shaped
  }

  if (IsProblemComputed() == false) {
    AMESOS_CHK_ERR(-6); // not yet computed
  }
  
  AMESOS_CHK_ERR(Solver_.Solve());

  return(0);
}
 

//==============================================================================
int Amesos_ContainerLAPACK::Destroy() {
  IsProblemShaped_ = false;
  IsProblemComputed_ = false;
  return(0);
}

//==============================================================================
int& Amesos_ContainerLAPACK::GID(const int i)
{
  return(GID_[i]);
}

