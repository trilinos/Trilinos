#ifndef IFPACK_DENSECONTAINER_H
#define IFPACK_DENSECONTAINER_H

#include "Ifpack_Container.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_RowMatrix.h"

//! Ifpack_DenseContainer: a class to define containers for dense matrices.
/*!

<P>To understand what an IFPACK container is, please refer to the documentation 
of the pure virtual class Ifpack_Container. Currently, containers are
used by class Ifpack_BlockRelaxation.

<P>Using block methods, one needs to store all diagonal blocks and
to be also to apply the inverse of each diagonal block. Using
class Ifpack_DenseContainer, one can store the blocks as dense
matrices, which can be advantageous when the 
blocks are small. Otherwise,
class Ifpack_SparseContainer is probably more appropriate.

<P>A typical use of a container is as follows:
\code
#include "Ifpack_DenseContainer.h"
...

// local matrix of (5,5), with two vectors for solution and rhs.
Ifpack_Container* Container = new
  Ifpack_DenseContainer(5,5);
  
// assign local rows 1, 5, 12, 13, 16 to this container
Container(0) = 1;
Container(1) = 5;
Container(2) = 12;
Container(3) = 13;
Container(4) = 16;

// Now extract the submatrix corresponding to rows and columns:
// 1. initialize the container.
Container.Initialize();
// 2. extract matrix values from an Epetra_RowMatrix A,
// and compute LU factors of the submatrix identified by rows
// and columns 1, 5, 12, 13 and 16 using LAPACK
Container.Compute(A);

// We can set the RHS as follows:
Container.RHS(0) = 1.0;
Container.RHS(1) = 2.0;
Container.RHS(2) = 3.0;
Container.RHS(3) = 4.0;
Container.RHS(4) = 5.0;

// The linear system with the submatrix is solved as follows:
Container.ApplyInverse().
\endcode

A call to Compute() computes the LU factorization of the
linear system matrix, using LAPACK (more precisely, by calling
the corresponding routines in Epetra_SerialDenseSolver). 
The default behavior is 
to store the matrix factors by overwriting the linear system matrix
itself. This way, method Apply() fails, as the original matrix
does no longer exists. An alternative is to call
\c KeepNonFactoredMatrix(true), which forces Ifpack_DenseContainer to
maintain in memory a copy of the non-factored matrix.

\author Marzio Sala, SNL 9214.

\date Last update Nov-04.
*/

class Ifpack_DenseContainer : public Ifpack_Container {

public:

  //@{ Constructors/Destructors

  //! Default constructor
  Ifpack_DenseContainer(const int NumRows, const int NumVectors = 1) :
    NumRows_(NumRows),
    NumVectors_(NumVectors),
    KeepNonFactoredMatrix_(false),
    IsInitialized_(false),
    IsComputed_(false),
    ComputeFlops_(0),
    ApplyFlops_(0),
    ApplyInverseFlops_(0)
  {}

  //! Copy constructor
  Ifpack_DenseContainer(const Ifpack_DenseContainer& rhs) :
    NumRows_(rhs.NumRows()),
    NumVectors_(rhs.NumVectors()),
    KeepNonFactoredMatrix_(rhs.KeepNonFactoredMatrix()),
    IsInitialized_(rhs.IsInitialized()),
    IsComputed_(rhs.IsComputed())
  {
    Matrix_ = rhs.Matrix();
    if (KeepNonFactoredMatrix_)
      NonFactoredMatrix_ = rhs.NonFactoredMatrix();
    LHS_ = rhs.LHS();
    RHS_ = rhs.RHS();
    ID_ = rhs.ID();
  }
    
  //! Destructor.
  virtual ~Ifpack_DenseContainer()
  {}
  //@}

  //@{ Overloaded operators.

  //! Operator=
  Ifpack_DenseContainer& operator=(const Ifpack_DenseContainer& rhs)
  {
    if (&rhs == this)
      return(*this);

    NumRows_ = rhs.NumRows();
    NumVectors_ = rhs.NumVectors();
    IsComputed_ = rhs.IsComputed();
    KeepNonFactoredMatrix_ = rhs.KeepNonFactoredMatrix();
    Matrix_ = rhs.Matrix();
    if (KeepNonFactoredMatrix_)
      NonFactoredMatrix_ = rhs.NonFactoredMatrix();
    LHS_ = rhs.LHS();
    RHS_ = rhs.RHS();
    ID_ = rhs.ID();

    return(*this);
  }

  //@}

  //@{ Get/Set methods.

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
    if (NumVectors_ == NumVectors) 
      return(0);

    NumVectors_ = NumVectors;
    IFPACK_CHK_ERR(RHS_.Reshape(NumRows_,NumVectors_));
    IFPACK_CHK_ERR(RHS_.Reshape(NumRows_,NumVectors_));
    // zero out vector elements
    for (int i = 0 ; i < NumRows_ ; ++i)
      for (int j = 0 ; j < NumVectors_ ; ++j) {
	LHS_(i,j) = 0.0;
	RHS_(i,j) = 0.0;
      }

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

  //! Set the matrix element (row,col) to \c value.
  virtual int SetMatrixElement(const int row, const int col,
			       const double value);

  //! Sets all necessary parameters.
  virtual int SetParameters(Teuchos::ParameterList& List)
  {
    return(0);
  }

  //! Returns \c true is the container has been successfully initialized.
  virtual bool IsInitialized() const
  {
    return(IsInitialized_);
  }

  //! Returns \c true is the container has been successfully computed.
  virtual bool IsComputed() const
  {
    return(IsComputed_);
  }

  //! Returns the label of \e this container.
  virtual const char* Label() const
  {
    return(Label_.c_str());
  }

  //! If \c flag is \c true, keeps a copy of the non-factored matrix.
  virtual int SetKeepNonFactoredMatrix(const bool flag)
  {
    KeepNonFactoredMatrix_ = flag;
    return(0);
  }

  //! Returns KeepNonFactoredMatrix_.
  virtual bool KeepNonFactoredMatrix() const
  {
    return(KeepNonFactoredMatrix_);
  }

  //! Returns the dense vector containing the LHS.
  virtual const Epetra_SerialDenseMatrix& LHS() const
  {
    return(LHS_);
  }

  //! Returns the dense vector containing the RHS.
  virtual const Epetra_SerialDenseMatrix& RHS() const
  {
    return(RHS_);
  }

  //! Returns the dense matrix or its factors.
  virtual const Epetra_SerialDenseMatrix& Matrix() const
  {
    return(Matrix_);
  }

  //! Returns the non-factored dense matrix (only if stored).
  virtual const Epetra_SerialDenseMatrix& NonFactoredMatrix() const
  {
    return(NonFactoredMatrix_);
  }

  //! Returns the integer dense vector of IDs.
  virtual const Epetra_IntSerialDenseVector& ID() const
  {
    return(ID_);
  }

  //@}

  //@{ Mathematical methods.
  //! Initialize the container.
  virtual int Initialize();

  //! Finalizes the linear system matrix and prepares for the application of the inverse.
  virtual int Compute(const Epetra_RowMatrix& Matrix);

  //! Apply the matrix to RHS, results are stored in LHS.
  virtual int Apply();

  //! Apply the inverse of the matrix to RHS, results are stored in LHS.
  virtual int ApplyInverse();

  //@}

  virtual double ComputeFlops() const
  {
    return(ComputeFlops_);
  }

  virtual double ApplyFlops() const
  {
    return(ApplyFlops_);
  }

  virtual double ApplyInverseFlops() const
  {
    return(ApplyInverseFlops_);
  }

  //! Prints basic information on iostream. This function is used by operator<<.
  virtual ostream& Print(std::ostream& os) const;

private:
  
  //! Extract the submatrices identified by the ID set int ID().
  virtual int Extract(const Epetra_RowMatrix& Matrix);

  //! Number of rows in the container.
  int NumRows_; 
  //! Number of vectors in the container.
  int NumVectors_;
  //! Dense matrix, that contains the non-factored matrix.
  Epetra_SerialDenseMatrix NonFactoredMatrix_;
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
  //! If \c true, keeps a copy of the non-factored matrix.
  bool KeepNonFactoredMatrix_;
  //! If \c true, the container has been successfully initialized.
  bool IsInitialized_;
  //! If \c true, the container has been successfully computed.
  bool IsComputed_;
  //! Label for \c this object
  string Label_;

  //! Flops in Compute().
  double ComputeFlops_;
  //! Flops in Apply().
  double ApplyFlops_;
  //! Flops in ApplyInverse().
  double ApplyInverseFlops_;
};

//==============================================================================
int Ifpack_DenseContainer::NumRows() const
{
  return(NumRows_);
}

//==============================================================================
int Ifpack_DenseContainer::Initialize()
{
  
  IsInitialized_ = false;

  if (Matrix_.N() == 0) {
    IFPACK_CHK_ERR(LHS_.Shape(NumRows_,NumVectors_));
    IFPACK_CHK_ERR(RHS_.Shape(NumRows_,NumVectors_));
    IFPACK_CHK_ERR(ID_.Shape(NumRows_,NumVectors_));
    IFPACK_CHK_ERR(Matrix_.Shape(NumRows_,NumRows_));
  }
  else {
    IFPACK_CHK_ERR(LHS_.Reshape(NumRows_,NumVectors_));
    IFPACK_CHK_ERR(RHS_.Reshape(NumRows_,NumVectors_));
    IFPACK_CHK_ERR(ID_.Reshape(NumRows_,NumVectors_));
    IFPACK_CHK_ERR(Matrix_.Reshape(NumRows_,NumRows_));
  }

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

  IsInitialized_ = true;
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
  if (IsInitialized() == false)
    IFPACK_CHK_ERR(Initialize());

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

  if (IsComputed() == false) {
    IFPACK_CHK_ERR(-3); // not yet computed
  }
  
  IFPACK_CHK_ERR(Solver_.Solve());

  ApplyInverseFlops_ += 2.0 * NumVectors_ * NumRows_ * NumRows_;
  return(0);
}

//==============================================================================
int& Ifpack_DenseContainer::ID(const int i)
{
  return(ID_[i]);
}

//==============================================================================
// FIXME: optimize performances of this guy...
int Ifpack_DenseContainer::Extract(const Epetra_RowMatrix& Matrix)
{

  for (int j = 0 ; j < NumRows_ ; ++j) {
    // be sure that the user has set all the ID's
    if (ID(j) == -1)
      IFPACK_CHK_ERR(-2);
    // be sure that all are local indices
    if (ID(j) > Matrix.NumMyRows())
      IFPACK_CHK_ERR(-2);
  }

  // allocate storage to extract matrix rows.
  int Length = Matrix.MaxNumEntries();
  vector<double> Values;
  Values.resize(Length);
  vector<int> Indices;
  Indices.resize(Length);

  for (int j = 0 ; j < NumRows_ ; ++j) {

    int LRID = ID(j);

    int NumEntries;

    int ierr = 
      Matrix.ExtractMyRowCopy(LRID, Length, NumEntries, 
			      &Values[0], &Indices[0]);
    IFPACK_CHK_ERR(ierr);

    for (int k = 0 ; k < NumEntries ; ++k) {

      int LCID = Indices[k];

      // skip off-processor elements
      if (LCID >= Matrix.NumMyRows()) 
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

//==============================================================================
int Ifpack_DenseContainer::Compute(const Epetra_RowMatrix& Matrix)
{
  IsComputed_ = false;
  if (IsInitialized() == false) {
    IFPACK_CHK_ERR(Initialize());
  }

  if (KeepNonFactoredMatrix_)
    NonFactoredMatrix_ = Matrix_;

  // extract local rows and columns
  IFPACK_CHK_ERR(Extract(Matrix));

  if (KeepNonFactoredMatrix_)
    NonFactoredMatrix_ = Matrix_;

  // factorize the matrix using LAPACK
  IFPACK_CHK_ERR(Solver_.Factor());

  Label_ = "Ifpack_DenseContainer";

  // not sure of count
  ComputeFlops_ += 2.0 * NumRows_ * NumRows_ * NumRows_ / 3;
  IsComputed_ = true;

  return(0);
}

//==============================================================================
int Ifpack_DenseContainer::Apply()
{
  if (IsComputed() == false)
    IFPACK_CHK_ERR(-3);

  if (KeepNonFactoredMatrix_) {
    IFPACK_CHK_ERR(RHS_.Multiply('N','N', 1.0,NonFactoredMatrix_,LHS_,0.0));
  }
  else
    IFPACK_CHK_ERR(RHS_.Multiply('N','N', 1.0,Matrix_,LHS_,0.0));

  ApplyFlops_ += 2 * NumRows_ * NumRows_;
  return(0);
}

//==============================================================================
ostream& Ifpack_DenseContainer::Print(ostream & os) const
{
    os << "================================================================================" << endl;
  os << "Ifpack_DenseContainer" << endl;
  os << "Number of rows          = " << NumRows_ << endl;
  os << "Number of vectors       = " << NumVectors_ << endl;
  os << "IsInitialized()         = " << IsInitialized_ << endl;
  os << "IsComputed()            = " << IsComputed_ << endl;
  os << "Flops in Compute()      = " << ComputeFlops_ << endl; 
  os << "Flops in ApplyInverse() = " << ApplyInverseFlops_ << endl; 
  os << "================================================================================" << endl;
  os << endl;

  return(os);
}
#endif
