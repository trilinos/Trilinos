#ifndef AMESOS_CONTAINER_H
#define AMESOS_CONTAINER_H

class Epetra_BlockMap;
#include "Amesos_InverseFactory.h"

class Amesos_Container {

public:

  virtual ~Amesos_Container() {};

  virtual int Shape(const int NumRows, const int NumVectors = 1) = 0;

  virtual int Shape(const Epetra_BlockMap& Map, const int NumVectors = 1) = 0;

  virtual int Reshape(const int NumRows, const int NumVectors = 1) = 0;

  virtual int Reshape(const Epetra_BlockMap& Map, const int NumVectors = 1) = 0;

  virtual int NumRows() const = 0;

  virtual int NumVectors() const = 0;

  virtual double& LHS(const int i, const int Vector = 0) = 0;
  
  virtual double& RHS(const int i, const int Vector = 0) = 0;

  virtual int& GID(const int i) = 0;

  virtual int SetMatrixElement(const int row, const int col,
			       const double value) = 0;

  virtual int Compute() = 0;

  virtual int ComputeInverse(char* Type,
			     Amesos_InverseFactory& Factory,
			     Teuchos::ParameterList& List) = 0;
  
  virtual bool IsProblemShaped() const = 0;

  virtual bool IsProblemComputed() const = 0;

  virtual int Apply() = 0;

  virtual int ApplyInverse() = 0;

  virtual int Destroy() = 0;

};

#endif // AMESOS_CONTAINER_H
