#ifndef IFPACK_DENSECONTAINER_H
#define IFPACK_DENSECONTAINER_H

#include "Ifpack_Container.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
class Epetra_RowMatrix;

class Ifpack_DenseContainer : public Ifpack_Container {

public:

  Ifpack_DenseContainer() :
    NumRows_(-1),
    NumVectors_(1),
    IsComputed_(false),
    IsShaped_(false)
  {}

  virtual ~Ifpack_DenseContainer()
  {
    Destroy();
  }

  virtual int Shape(const int NumRows, const int NumVectors = 1);

  virtual int Reshape(const int NumRows, const int NumVectors = 1);

  virtual int NumRows() const;

  virtual int NumVectors() const
  {
    return(NumVectors_);
  }

  virtual int SetNumVectors(const int NumVectors)
  {
    NumVectors_ = NumVectors;
    return(0);
  }

  virtual double& LHS(const int i, const int Vector = 0);
  
  virtual double& RHS(const int i, const int Vector = 0);

  virtual int& ID(const int i);

  virtual int Extract(const Epetra_RowMatrix* Matrix);

  virtual int SetMatrixElement(const int row, const int col,
			       const double value);

  virtual int SetParameters(Teuchos::ParameterList& List)
  {
    return(0);
  }

  virtual int Compute()
  {
    if (IsComputed() == true) {
      IFPACK_CHK_ERR(-6); // already computed
    }

    // factorize the matrix using LAPACK
    IFPACK_CHK_ERR(Solver_.Factor());
    IsComputed_ = true;
  }


  virtual bool IsShaped() const
  {
    return(IsShaped_);
  }

  virtual bool IsComputed() const
  {
    return(IsComputed_);
  }

  virtual int Apply()
  {
    IFPACK_CHK_ERR(-1);
  }

  virtual int ApplyInverse();

  virtual int Destroy()
  {
    IsShaped_ = false;
    IsComputed_ = false;
  }    

private:
  
  int NumRows_; 
  int NumVectors_; 
  Epetra_SerialDenseMatrix Matrix_;
  Epetra_SerialDenseMatrix LHS_;
  Epetra_SerialDenseMatrix RHS_;
  Epetra_SerialDenseSolver Solver_;
  Epetra_IntSerialDenseVector ID_;
  bool IsShaped_;
  bool IsComputed_;

};

//==============================================================================
int Ifpack_DenseContainer::NumRows() const
{
  if (IsShaped() == false)
    return(0);
  else
    return(NumRows_);
}

//==============================================================================
int Ifpack_DenseContainer::Shape(const int NumRows, const int NumVectors) 
{
  
  if (IsShaped() == true)
    IFPACK_CHK_ERR(-1); // already computed

  IsShaped_ = true;

  NumRows_ = NumRows;
  NumVectors_ = NumVectors;

  IFPACK_CHK_ERR(LHS_.Reshape(NumRows_,NumVectors));
  IFPACK_CHK_ERR(RHS_.Reshape(NumRows_,NumVectors));
  IFPACK_CHK_ERR(ID_.Reshape(NumRows_,NumVectors));

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

  IFPACK_CHK_ERR(Solver_.SetMatrix(Matrix_));
  IFPACK_CHK_ERR(Solver_.SetVectors(LHS_,RHS_));

  return(0);
  
}

//==============================================================================
int Ifpack_DenseContainer::Reshape(const int NumRows, const int NumVectors) 
{

  if (IsShaped() == true)
    Destroy();

  Shape(NumRows,NumVectors);
  
}

//==============================================================================
double& Ifpack_DenseContainer::LHS(const int i, const int Vector)
{
  return(LHS_.A()[Vector * NumRows_ + i]);
}
  
//==============================================================================
double& Ifpack_DenseContainer::RHS(const int i, const int Vector)
{
  return(RHS_.A()[Vector * NumRows_ + i]);
}

//==============================================================================
int Ifpack_DenseContainer::
SetMatrixElement(const int row, const int col, const double value)
{
  if (IsShaped() == false)
    IFPACK_CHK_ERR(-5); // problem not shaped yet

  if ((row < 0) || (row >= NumRows())) {
    IFPACK_CHK_ERR(-2); // not in range
  }

  if ((col < 0) || (col >= NumRows())) {
    IFPACK_CHK_ERR(-2); // not in range
  }

  Matrix_(row, col) = value;

  return(0);

}

//==============================================================================
int Ifpack_DenseContainer::ApplyInverse()
{
  if (IsShaped() == false) {
    IFPACK_CHK_ERR(-5); // not yet shaped
  }

  if (IsComputed() == false) {
    IFPACK_CHK_ERR(-6); // not yet computed
  }
  
  IFPACK_CHK_ERR(Solver_.Solve());

  return(0);
}

//==============================================================================
int& Ifpack_DenseContainer::ID(const int i)
{
  return(ID_[i]);
}

//==============================================================================
// FIXME: optimize performances of this guy...
int Ifpack_DenseContainer::Extract(const Epetra_RowMatrix* Matrix)
{

  int Length = Matrix->MaxNumEntries();
  vector<double> Values;
  Values.resize(Length);
  vector<int> Indices;
  Indices.resize(Length);

  for (int j = 0 ; j < NumRows_ ; ++j) {

    int LRID = ID(j);

    int NumEntries;

    int ierr = 
      Matrix->ExtractMyRowCopy(LRID, Length, NumEntries, 
			       &Values[0], &Indices[0]);
    IFPACK_CHK_ERR(ierr);

    for (int k = 0 ; k < NumEntries ; ++k) {

      int LCID = Indices[k];

      // skip off-processor elements
      if (LCID >= Matrix->NumMyRows()) 
	continue;

      // for local column IDs, look for each ID in the list
      // of columns hosted by this object
      // FIXME: use STL
      int jj = -1;
      for (int kk = 0 ; kk < NumRows_ ; ++kk)
	if (ID(kk) == LCID)
	  jj = kk;

      if (jj != -1)
	SetMatrixElement(j,jj,Values[k]);

    }
  }

  return(0);
}


#endif
