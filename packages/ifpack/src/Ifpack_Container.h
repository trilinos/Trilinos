#ifndef IFPACK_CONTAINER_H
#define IFPACK_CONTAINER_H

class Epetra_RowMatrix;
class Ifpack_Partitioner;
namespace Teuchos {
  class ParameterList;
}

//! Ifpack_Container: a pure virtual class for creating and solving local linear systems.

/*!
  Class Ifpack_Containers provides the abstract interfaces for 
  containers. A "container" is an object that hosts all it is necessary
  to create, populate, and solve local linear system.

  A container should be used in the following manner:
  1- create an container object.
  2- set the number of rows using Shape(). All rows are local.
  3- set matrix elements using SetMatrixElement(),
     and LHS and/or RHS elements using LHS() and  RHS().
  4- if necessary, set parameters for the solution using
     SetParameters().
  5- prepare the linear system solver using Compute().
  6- Solve the linear system using ApplyInverse().

  The number of vectors can be set using SetNumVectors().
  
  Containers are used in block preconditioners (Ifpack_BlockJacobi
  and Ifpack_BlockGaussSeidel). For these preconditioners, a container
  represents a block to be solved. Users can take advantage of
  method ID(), which can be used to set an integer value to each row 
  in the container (typically, ID can be set as restriction operator
  from the matrix to be preconditioned, to the sub-block in the container).

  Two concrete implementations are provided in classes
  SparseContainer (that stores matrices in sparse the format
  Epetra_CrsMatrix) and DenseContainer
  (for relatively small matrices, as matrices are stored as
  Epetra_SerialDenseMatrix's).
  
  \date Sep-04
  
*/

class Ifpack_Container {

public:

  virtual ~Ifpack_Container() {};

  virtual int Shape(const int NumRows, const int NumVectors = 1) = 0;

  virtual int Reshape(const int NumRows, const int NumVectors = 1) = 0;

  virtual int NumRows() const = 0;

  virtual int NumVectors() const = 0;

  virtual int SetNumVectors(const int i) = 0;

  virtual double& LHS(const int i, const int Vector = 0) = 0;
  
  virtual double& RHS(const int i, const int Vector = 0) = 0;

  virtual int& ID(const int i) = 0;

  virtual int Extract(const Epetra_RowMatrix* Matrix) = 0;
  
  virtual int SetMatrixElement(const int row, const int col,
			       const double value) = 0;

  virtual int Compute() = 0;

  virtual int SetParameters(Teuchos::ParameterList& List) = 0;

  virtual bool IsShaped() const = 0;

  virtual bool IsComputed() const = 0;

  virtual int Apply() = 0;

  virtual int ApplyInverse() = 0;

  virtual int Destroy() = 0;

};


#endif // IFPACK_CONTAINER_H
