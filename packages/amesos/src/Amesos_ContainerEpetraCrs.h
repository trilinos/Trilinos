#ifndef AMESOS_CONTAINEREPETRACRS_H
#define AMESOS_CONTAINEREPETRACRS_H

#include "Amesos_Container.h"
class Epetra_Comm;
class Epetra_Map;
class Epetra_RowMatrix;
class Epetra_Operator;
class Epetra_CrsMatrix;
class Epetra_LinearProblem;
class Epetra_MultiVector;
class Epetra_BlockMap;

class Amesos_ContainerEpetraCrs : public Amesos_Container {

public:

  Amesos_ContainerEpetraCrs();

  virtual ~Amesos_ContainerEpetraCrs();

  virtual int Shape(const int NumRows, const int NumVectors = 1);

  virtual int Shape(const Epetra_BlockMap& Map, const int NumVectors = 1)
  {
    AMESOS_CHK_ERR(-98);
  }

  virtual int Reshape(const int NumRows, const int NumVectors = 1);

  virtual int Reshape(const Epetra_BlockMap& Map, const int NumVectors = 1)
  {
    AMESOS_CHK_ERR(-98);
  }

  virtual int NumRows() const;

  virtual int NumVectors() const
  {
    return(NumVectors_);
  }

  virtual double& LHS(const int i, const int Vector = 0);
  
  virtual double& RHS(const int i, const int Vector = 0);

  virtual int SetMatrixElement(const int row, const int col,
			       const double value);

  virtual int GetMatrixPointer(void** Matrix);

  virtual int SetInversePointer(void* Inverse);

  virtual int Compute();

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
  Epetra_Map* Map_;
  Epetra_CrsMatrix* Matrix_;
  Epetra_Operator* Inverse_;
  Epetra_MultiVector* LHS_;
  Epetra_MultiVector* RHS_;

  Epetra_LinearProblem* Solver_;

  bool IsProblemShaped_;
  bool IsProblemComputed_;
  Epetra_Comm* SerialComm_;

};

#endif // AMESOS_CONTAINEREPETRACRS_H

