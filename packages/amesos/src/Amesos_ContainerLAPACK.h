#ifndef AMESOS_CONTAINERLAPACK
#define AMESOS_CONTAINERLAPACK

#include "Amesos_ConfigDefs.h"
#include "Amesos_Container.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_IntSerialDenseVector.h"

class Amesos_ContainerLAPACK : public Amesos_Container {

public:

  Amesos_ContainerLAPACK();

  virtual ~Amesos_ContainerLAPACK();

  virtual int Shape(const int NumRows, const int NumVectors = 1);

  virtual int Shape(const Epetra_BlockMap& Map, const int NumVectors = 1)
  {
    AMESOS_CHK_ERR(-99);
  }

  virtual int Reshape(const int NumRows, const int NumVectors = 1);

  virtual int Reshape(const Epetra_BlockMap& Map, const int NumVectors = 1)
  {
    AMESOS_CHK_ERR(-99);
  }

  virtual int NumRows() const;

  virtual int NumVectors() const;

  virtual double& LHS(const int i, const int Vector = 0);
  
  virtual double& RHS(const int i, const int Vector = 0);

  virtual int& GID(const int i);

  virtual int SetMatrixElement(const int row, const int col,
			       const double value);

  virtual int Compute();

  virtual int ComputeInverse(char* Type,
			     Amesos_InverseFactory& Factory,
			     Teuchos::ParameterList& List);
  
  virtual bool IsProblemShaped() const
  {
    return(IsProblemShaped_);
  }

  virtual bool IsProblemComputed() const
  {
    return(IsProblemComputed_);
  }

  virtual int Apply();

  virtual int ApplyInverse();

  virtual int Destroy();

private:
  
  int NumRows_; 
  int NumVectors_; 
  Epetra_SerialDenseMatrix Matrix_;
  Epetra_SerialDenseMatrix LHS_;
  Epetra_SerialDenseMatrix RHS_;
  Epetra_SerialDenseSolver Solver_;
  Epetra_IntSerialDenseVector GID_;
  bool IsProblemShaped_;
  bool IsProblemComputed_;

};

#endif // AMESOS_CONTAINERLAPACK
