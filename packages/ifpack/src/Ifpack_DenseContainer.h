#ifndef IFPACK_DENSECONTAINER_H
#define IFPACK_DENSECONTAINER_H

#include "Ifpack_Container.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
class Epetra_RowMatrix;

//! Ifpack_DenseContainer: a class to define containers for dense matrices.
/*!
Class Ifpack_DenseContainer enables the use of containers based on dense
matrices. In block methods like Ifpack_BlockJacobi and Ifpack_GaussSeidel,
if the blocks are small, it can be convenient to store then as dense
matrices, as efficient solvers are available through LAPACK.

A typical use of a container is as follows:
\code
#include "Ifpack_Container.h"
#include "Ifpack_DenseContainer.h"
...

Ifpack_Container* Container = new
  Ifpack_DenseContainer;
  
// local matrix of (5,5), with two vectors for solution and rhs.
Container.Shape(5,2);
// assign local rows 1, 5, 12, 13, 16 to this container
Container(0) = 1;
Container(1) = 5;
Container(2) = 12;
Container(3) = 13;
Container(4) = 16;

// Now extract the submatrix corresponding to rows and columns
// 1, 5, 12, 13, 16:
Container.Extract();

// Compute all we need (the LU factors to apply the inverse)
Container.Compute();

// We can set the RHS as follows:
Container.RHS(0) = 1.0;
Container.RHS(1) = 2.0;
Container.RHS(2) = 3.0;
Container.RHS(3) = 4.0;
Container.RHS(4) = 5.0;

// Solve as follows:
Container.ApplyInverse().
\endcode

Note that method Apply() cannot be used, as the LU factorization as
implemented in Epetra_SerialDenseMatrix (for performance reasons) 
overwrites the matrix itself.

\date Sep-04
*/

class Ifpack_DenseContainer : public Ifpack_Container {

public:

  Ifpack_DenseContainer() :
    NumRows_(-1),
    NumVectors_(1),
    IsComputed_(false),
    IsShaped_(false)
  {}

  //! Destructor.
  virtual ~Ifpack_DenseContainer()
  {
    Destroy();
  }

  //! Shapes the container, allocating space for \c NumRows and \c NumVectors.
  virtual int Shape(const int NumRows, const int NumVectors = 1);

  //! Reshapes a container.
  virtual int Reshape(const int NumRows, const int NumVectors = 1);

  //! Returns the number of rows of the matrix and LHS/RHS.
  virtual int NumRows() const;

  //! Returns the number of vectors in LHS/RHS.
  virtual int NumVectors() const
  {
    return(NumVectors_);
  }

  //! Sets the number of vectors for LHS/RHS.
  virtual int SetNumVectors(const int NumVectors)
  {
    NumVectors_ = NumVectors;
    return(0);
  }

  //! Returns the i-th component of the vector Vector of LHS.
  virtual double& LHS(const int i, const int Vector = 0);
  
  //! Returns the i-th component of the vector Vector of RHS.
  virtual double& RHS(const int i, const int Vector = 0);

  //! Returns the ID associated to local row i. 
  /*!
   * The set of (local) rows assigned to this container is defined
   * by calling ID(i) = j, where i (from 0 to NumRows()) indicates
   * the container-row, and j indicates the local row in the calling
   * process.
   *
   * This is usually used to recorder the local row ID (on calling process)
   * of the i-th row in the container.
   */
  virtual int& ID(const int i);

  //! Extract the submatrices identified by the ID set int ID().
  virtual int Extract(const Epetra_RowMatrix* Matrix);

  //! Set the matrix element (row,col) to \c value.
  virtual int SetMatrixElement(const int row, const int col,
			       const double value);

  //! Sets all necessary parameters.
  virtual int SetParameters(Teuchos::ParameterList& List)
  {
    return(0);
  }

  //! Finalizes the linear system matrix and prepares for the application of the inverse.
  virtual int Compute()
  {
    if (IsComputed() == true) {
      IFPACK_CHK_ERR(-6); // already computed
    }

    // factorize the matrix using LAPACK
    IFPACK_CHK_ERR(Solver_.Factor());

    Label_ = "Ifpack_DenseContainer";
    IsComputed_ = true;
    return(0);
  }

  //! Returns \c true is the containers is shaped.
  virtual bool IsShaped() const
  {
    return(IsShaped_);
  }

  //! Returns \c true is the container has successfully called \c Compute().
  virtual bool IsComputed() const
  {
    return(IsComputed_);
  }

  //! Apply the matrix to RHS, results are stored in LHS.
  virtual int Apply()
  {
    IFPACK_CHK_ERR(-1);
  }

  //! Apply the inverse of the matrix to RHS, results are stored in LHS.
  virtual int ApplyInverse();

  //! Destroys all data.
  virtual int Destroy()
  {
    IsShaped_ = false;
    IsComputed_ = false;
    return(0);
  }    

  //! Returns the label of \e this container.
  virtual const char* Label() const
  {
    return(Label_.c_str());
  }

private:
  
  //! Number of rows in the container.
  int NumRows_; 
  //! Number of vectors in the container.
  int NumVectors_;
  //! Dense matrix.
  Epetra_SerialDenseMatrix Matrix_;
  //! Dense vector representing the LHS.
  Epetra_SerialDenseMatrix LHS_;
  //! Dense vector representing the RHS.
  Epetra_SerialDenseMatrix RHS_;
  //! Dense solver (solution will be get using LAPACK).
  Epetra_SerialDenseSolver Solver_;
  //! Sets of local rows.
  Epetra_IntSerialDenseVector ID_;
  //! If \c true, the container has been shaped.
  bool IsShaped_;
  //! If \c true, the container has been computed.
  bool IsComputed_;
  //! Label for \c this object
  string Label_;

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

  // Set to -1 ID_'s
  for (int i = 0 ; i < NumRows_ ; ++i)
    ID_(i) = -1;  

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

  return(0);
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

  if (Matrix == 0)
    IFPACK_CHK_ERR(-1);

  for (int j = 0 ; j < NumRows_ ; ++j) {
    // be sure that the user has set all the ID's
    if (ID(j) == -1)
      IFPACK_CHK_ERR(-1);
    // be sure that all are local indices
    if (ID(j) > Matrix->NumMyRows())
      IFPACK_CHK_ERR(-2);
  }

  // allocate storage to extract matrix rows.
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
