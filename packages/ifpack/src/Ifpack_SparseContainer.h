#ifndef IFPACK_SPARSECONTAINER_H
#define IFPACK_SPARSECONTAINER_H

#include "Ifpack_Container.h"
#ifdef HAVE_IFPACK_TEUCHOS
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Teuchos_ParameterList.hpp"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

//! Ifpack_SparseContianer: a class for storing and solving linear systems.

/*!
  Containers are templated with a type T, which represent the 
  class to use in the application of the inverse. (T is not
  used in Ifpack_DenseContainer). In SparseContainer, T must be
  an Ifpack_Preconditioner derived class. The container will allocate
  a \c T object, use SetParameters() and Compute(), then
  use \c T every time the linear system as to be solved (using the
  ApplyInverse() method of \c T).
*/

template<typename T>
class Ifpack_SparseContainer : public Ifpack_Container {

public:

  //! Constructor.
  Ifpack_SparseContainer();

  //! Destructor.
  virtual ~Ifpack_SparseContainer();

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

  //! Finalizes the linear system matrix and prepares for the application of the inverse.
  virtual int Compute();

  //! Returns \c true is the container has successfully called \c Compute().
  virtual bool IsComputed() const
  {
    return(IsComputed_);
  }

  //! Sets all necessary parameters.
  virtual int SetParameters(Teuchos::ParameterList& List);

  virtual bool IsShaped() const
  {
    return(IsShaped_);
  }

  //! Apply the matrix to RHS, results are stored in LHS.
  virtual int Apply();

  //! Apply the inverse of the matrix to RHS, results are stored in LHS.
  virtual int ApplyInverse();

  //! Destroys all data.
  virtual int Destroy();

  //! Returns the label of \e this container.
  virtual const char* Label() const
  {
    return(Label_.c_str());
  }
  
private:
  
  int NumRows_; 
  int NumVectors_; 
  Epetra_Map* Map_;
  Epetra_CrsMatrix* Matrix_;
  Epetra_MultiVector* LHS_;
  Epetra_MultiVector* RHS_;
  Epetra_IntSerialDenseVector GID_;

  Epetra_LinearProblem* Problem_;

  bool IsShaped_;
  bool IsComputed_;
  Epetra_Comm* SerialComm_;
  T* Inverse_;
  //! Label for \c this object
  string Label_;

};

//==============================================================================
template<typename T>
Ifpack_SparseContainer<T>::Ifpack_SparseContainer() :
  NumRows_(-1),
  NumVectors_(-1),
  Map_(0),
  Matrix_(0),
  Inverse_(0),
  LHS_(0),
  RHS_(0),
  Problem_(0),
  IsShaped_(false),
  IsComputed_(false),
  SerialComm_(0)
{

#ifdef HAVE_MPI
  SerialComm_ = new Epetra_MpiComm(MPI_COMM_SELF);
#else
  SerialComm_ = new Epetra_SerialComm;
#endif

}

//==============================================================================
template<typename T>
Ifpack_SparseContainer<T>::~Ifpack_SparseContainer()
{
  Destroy();
}

//==============================================================================
template<typename T>
int Ifpack_SparseContainer<T>::NumRows() const
{
  if (IsShaped() == false)
    return(0);
  else
    return(NumRows_);
}

//==============================================================================
template<typename T>
int Ifpack_SparseContainer<T>::Shape(const int NumRows, const int NumVectors) 
{
  
  if (IsShaped() == true)
    IFPACK_CHK_ERR(-1); // already computed

  IsShaped_ = true;

  NumRows_ = NumRows;
  NumVectors_ = NumVectors;

  Map_ = new Epetra_Map(NumRows_,0,*SerialComm_);

  LHS_ = new Epetra_MultiVector(*Map_,NumVectors);
  RHS_ = new Epetra_MultiVector(*Map_,NumVectors);
  GID_.Reshape(NumRows,1);

  // FIXME: try View??
  Matrix_ = new Epetra_CrsMatrix(Copy,*Map_,0);
  Problem_ = new Epetra_LinearProblem;
  
  Problem_->SetOperator(Matrix_);
  Problem_->SetLHS(LHS_);
  Problem_->SetRHS(RHS_);

  // create the inverse
  Inverse_ = new T(Matrix_);

  if (Inverse_ == 0)
    IFPACK_CHK_ERR(-10);

  return(0);
  
}

//==============================================================================
template<typename T>
int Ifpack_SparseContainer<T>::Reshape(const int NumRows, const int NumVector) 
{

  if (IsShaped() == true)
    Destroy();

  Shape(NumRows);
  
}

//==============================================================================
template<typename T>
double& Ifpack_SparseContainer<T>::LHS(const int i, const int Vector)
{
  return(((*LHS_)(Vector))->Values()[i]);
}
  
//==============================================================================
template<typename T>
double& Ifpack_SparseContainer<T>::RHS(const int i, const int Vector)
{
  return(((*RHS_)(Vector))->Values()[i]);
}

//==============================================================================
template<typename T>
int Ifpack_SparseContainer<T>::
SetMatrixElement(const int row, const int col, double value)
{
  if (IsShaped() == false)
    IFPACK_CHK_ERR(-5); // problem not shaped yet

  if ((row < 0) || (row >= NumRows())) {
    IFPACK_CHK_ERR(-2); // not in range
  }

  if ((col < 0) || (col >= NumRows())) {
    IFPACK_CHK_ERR(-2); // not in range
  }

  int ierr = Matrix_->InsertGlobalValues((int)row,1,(double*)&value,(int*)&col);
  if (ierr < 0) {
    ierr = Matrix_->SumIntoGlobalValues((int)row,1,(double*)&value,(int*)&col);
    if (ierr < 0)
      IFPACK_CHK_ERR(-1);
  }

  return(0);

}

//==============================================================================
template<typename T>
int Ifpack_SparseContainer<T>::Compute()
{
  if (IsComputed() == true) {
    IFPACK_CHK_ERR(-6); // already computed
  }

  // complete the construction of Matrix
  IFPACK_CHK_ERR(Matrix_->FillComplete());
  // compute the inverse operator
  IFPACK_CHK_ERR(Inverse_->Compute());

  Label_ = "Ifpack_SparseContainer";
  
  IsComputed_ = true;

  return(0);
}

//==============================================================================
template<typename T>
int Ifpack_SparseContainer<T>::Apply()
{
  if (IsShaped() == false) {
    IFPACK_CHK_ERR(-5); // not yet shaped
  }

  if (IsComputed() == false) {
    IFPACK_CHK_ERR(-6); // not yet computed
  }
  
  IFPACK_CHK_ERR(Matrix_->Apply(*RHS_, *LHS_));

  return(0);
}

//==============================================================================
template<typename T>
int Ifpack_SparseContainer<T>::ApplyInverse()
{
  if (IsShaped() == false) {
    IFPACK_CHK_ERR(-5); // not yet shaped
  }

  if (IsComputed() == false) {
    IFPACK_CHK_ERR(-6); // not yet computed
  }
  
  if (Inverse_ == 0) {
    IFPACK_CHK_ERR(-7);
  }

  IFPACK_CHK_ERR(Inverse_->ApplyInverse(*RHS_, *LHS_));

  return(0);
}
 

//==============================================================================
template<typename T>
int Ifpack_SparseContainer<T>::Destroy()
{
  if (Map_)
    delete Map_;

  if (Matrix_)
    delete Matrix_;

  if (LHS_)
    delete LHS_;

  if (RHS_)
    delete RHS_;

  if (Inverse_)
    delete Inverse_;

  if (Problem_)
    delete Problem_;

  Map_ = 0;
  Matrix_ = 0;
  Problem_ = 0;
  Inverse_ = 0;
  LHS_ = 0;
  RHS_ = 0;

  IsShaped_ = false;
  IsComputed_ = false;
  return(0);
}

//==============================================================================
template<typename T>
int& Ifpack_SparseContainer<T>::ID(const int i)
{
  return(GID_[i]);
}

//==============================================================================
template<typename T>
int Ifpack_SparseContainer<T>::
SetParameters(Teuchos::ParameterList& List)
{
  Inverse_->SetParameters(List);
}

//==============================================================================
// FIXME: optimize performances of this guy...
template<typename T>
int Ifpack_SparseContainer<T>::Extract(const Epetra_RowMatrix* Matrix)
{

  for (int j = 0 ; j < NumRows_ ; ++j) {
    // be sure that the user has set all the ID's
    if (ID(j) == -1)
      IFPACK_CHK_ERR(-1);
    // be sure that all are local indices
    if (ID(j) > Matrix->NumMyRows())
      IFPACK_CHK_ERR(-2);
  }

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
#endif // HAVE_IFPACK_TEUCHOS
#endif // IFPACK_SPARSECONTAINER_H
