/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK_SPARSECONTAINER_H
#define IFPACK_SPARSECONTAINER_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include "Ifpack_Container.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

/*!
\brief Ifpack_SparseContainer: a class for storing and solving linear systems
using sparse matrices.

<P>To understand what an IFPACK container is, please refer to the documentation
of the pure virtual class Ifpack_Container. Currently, containers are
used by class Ifpack_BlockRelaxation.

<P>Using block methods, one needs to store all diagonal blocks and
to be also to apply the inverse of each diagonal block. Using
class Ifpack_DenseContainer, one can store the blocks as sparse
matrices (Epetra_CrsMatrix), which can be advantageous when the
blocks are large. Otherwise,
class Ifpack_DenseContainer is probably more appropriate.

<P>Sparse containers are templated with a type T, which represent the
class to use in the application of the inverse. (T is not
used in Ifpack_DenseContainer). In SparseContainer, T must be
an Ifpack_Preconditioner derived class. The container will allocate
a \c T object, use SetParameters() and Compute(), then
use \c T every time the linear system as to be solved (using the
ApplyInverse() method of \c T).

\author Marzio Sala, SNL 9214.

\date Last modified on Nov-04.

*/

template<typename T>
class Ifpack_SparseContainer : public Ifpack_Container {

public:

  //@{ Constructors/Destructors.
  //! Constructor.
  Ifpack_SparseContainer(const int NumRows, const int NumVectors = 1);

  //! Copy constructor.
  Ifpack_SparseContainer(const Ifpack_SparseContainer<T>& rhs);

  //! Destructor.
  virtual ~Ifpack_SparseContainer();
  //@}

  //@{ Overloaded operators.

  //! Operator =
  Ifpack_SparseContainer& operator=(const Ifpack_SparseContainer<T>& rhs);
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
   virtual int SetNumVectors(const int NumVectors_in)
   {
     if (NumVectors_ != NumVectors_in)
       {
       NumVectors_=NumVectors_in;
       LHS_=Teuchos::rcp(new Epetra_MultiVector(*Map_,NumVectors_));
       RHS_=Teuchos::rcp(new Epetra_MultiVector(*Map_,NumVectors_));
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

  //! Sets all necessary parameters.
  virtual int SetParameters(Teuchos::ParameterList& List);

  //! Returns the label of \e this container.
  virtual const char* Label() const
  {
    return(Label_.c_str());
  }

  //! Returns a pointer to the internally stored map.
  Teuchos::RCP<const Epetra_Map> Map() const
  {
    return(Map_);
  }

  //! Returns a pointer to the internally stored solution multi-vector.
  Teuchos::RCP<const Epetra_MultiVector> LHS() const
  {
    return(LHS_);
  }

  //! Returns a pointer to the internally stored rhs multi-vector.
  Teuchos::RCP<const Epetra_MultiVector> RHS() const
  {
    return(RHS_);
  }

  //! Returns a pointer to the internally stored matrix.
  Teuchos::RCP<const Epetra_CrsMatrix> Matrix() const
  {
    return(Matrix_);
  }

  //! Returns a pointer to the internally stored ID's.
  const Epetra_IntSerialDenseVector& ID() const
  {
    return(GID_);
  }

  //! Returns a pointer to the internally stored inverse operator.
  Teuchos::RCP<const T> Inverse() const
  {
    return(Inverse_);
  }
  //@}

  //@{ Mathematical functions.
  /*!
   * \brief Initializes the container, by completing all the operations based
   * on matrix structure.
   *
   * \note After a call to Initialize(), no new matrix entries can be
   * added.
   */
  virtual int Initialize();
  //! Finalizes the linear system matrix and prepares for the application of the inverse.
  virtual int Compute(const Epetra_RowMatrix& Matrix_in);
  //! Apply the matrix to RHS, result is stored in LHS.
  virtual int Apply();

  //! Apply the inverse of the matrix to RHS, result is stored in LHS.
  virtual int ApplyInverse();

  //@}

  //@{ Miscellaneous methods
  //! Destroys all data.
  virtual int Destroy();
  //@}

  //! Returns the flops in Compute().
  virtual double InitializeFlops() const
  {
    if (Inverse_ == Teuchos::null)
      return (0.0);
    else
      return(Inverse_->InitializeFlops());
  }

  //! Returns the flops in Compute().
  virtual double ComputeFlops() const
  {
    if (Inverse_ == Teuchos::null)
      return (0.0);
    else
      return(Inverse_->ComputeFlops());
  }

  //! Returns the flops in Apply().
  virtual double ApplyFlops() const
  {
    return(ApplyFlops_);
  }

  //! Returns the flops in ApplyInverse().
  virtual double ApplyInverseFlops() const
  {
    if (Inverse_ == Teuchos::null)
      return (0.0);
    else
      return(Inverse_->ApplyInverseFlops());
  }

  //! Prints basic information on iostream. This function is used by operator<<.
  virtual std::ostream& Print(std::ostream& os) const;

private:

  //! Extract the submatrices identified by the ID set int ID().
  virtual int Extract(const Epetra_RowMatrix& Matrix_in);

  //! Number of rows in the local matrix.
  int NumRows_;
  //! Number of vectors in the local linear system.
  int NumVectors_;
  //! Linear map on which the local matrix is based.
  Teuchos::RefCountPtr<Epetra_Map> Map_;
  //! Pointer to the local matrix.
  Teuchos::RefCountPtr<Epetra_CrsMatrix> Matrix_;
  //! Solution vector.
  Teuchos::RefCountPtr<Epetra_MultiVector> LHS_;
  //! right-hand side for local problems.
  Teuchos::RefCountPtr<Epetra_MultiVector> RHS_;
  //! Contains the subrows/subcols of A that will be inserted in Matrix_.
  Epetra_IntSerialDenseVector GID_;
  //! If \c true, the container has been successfully initialized.
  bool IsInitialized_;
  //! If \c true, the container has been successfully computed.
  bool IsComputed_;
  //! Serial communicator (containing only MPI_COMM_SELF if MPI is used).
  Teuchos::RefCountPtr<Epetra_Comm> SerialComm_;
  //! Pointer to an Ifpack_Preconditioner object whose ApplyInverse() defined the action of the inverse of the local matrix.
  Teuchos::RefCountPtr<T> Inverse_;
  //! Label for \c this object
  std::string Label_;
  Teuchos::ParameterList List_;
  double ApplyFlops_;

};

//==============================================================================
template<typename T>
Ifpack_SparseContainer<T>::
Ifpack_SparseContainer(const int NumRows_in, const int NumVectors_in) :
  NumRows_(NumRows_in),
  NumVectors_(NumVectors_in),
  IsInitialized_(false),
  IsComputed_(false),
  ApplyFlops_(0.0)
{

#ifdef HAVE_MPI
  SerialComm_ = Teuchos::rcp( new Epetra_MpiComm(MPI_COMM_SELF) );
#else
  SerialComm_ = Teuchos::rcp( new Epetra_SerialComm );
#endif

}

//==============================================================================
template<typename T>
Ifpack_SparseContainer<T>::
Ifpack_SparseContainer(const Ifpack_SparseContainer<T>& rhs) :
  NumRows_(rhs.NumRows()),
  NumVectors_(rhs.NumVectors()),
  IsInitialized_(rhs.IsInitialized()),
  IsComputed_(rhs.IsComputed())
{

#ifdef HAVE_MPI
  SerialComm_ = Teuchos::rcp( new Epetra_MpiComm(MPI_COMM_SELF) );
#else
  SerialComm_ = Teuchos::rcp( new Epetra_SerialComm );
#endif

  if (!rhs.Map().is_null())
    Map_ = Teuchos::rcp( new Epetra_Map(*rhs.Map()) );

  if (!rhs.Matrix().is_null())
    Matrix_ = Teuchos::rcp( new Epetra_CrsMatrix(*rhs.Matrix()) );

  if (!rhs.LHS().is_null())
    LHS_ = Teuchos::rcp( new Epetra_MultiVector(*rhs.LHS()) );

  if (!rhs.RHS().is_null())
    RHS_ = Teuchos::rcp( new Epetra_MultiVector(*rhs.RHS()) );

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
  if (IsInitialized() == false)
    return(0);
  else
    return(NumRows_);
}

//==============================================================================
template<typename T>
int Ifpack_SparseContainer<T>::Initialize()
{

  if (IsInitialized_ == true)
    Destroy();

  IsInitialized_ = false;

#if !defined(EPETRA_NO_32BIT_GLOBAL_INDICES) || !defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
  Map_ = Teuchos::rcp( new Epetra_Map(NumRows_,0,*SerialComm_) );
#endif

  LHS_ = Teuchos::rcp( new Epetra_MultiVector(*Map_,NumVectors_) );
  RHS_ = Teuchos::rcp( new Epetra_MultiVector(*Map_,NumVectors_) );
  GID_.Reshape(NumRows_,1);

#if defined(HAVE_TEUCHOSCORE_CXX11)
  Matrix_ = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy,*Map_,0) );
#else
  Matrix_ = Teuchos::rcp( new Epetra_CrsMatrix(::Copy,*Map_,0) );
#endif

  // create the inverse
  Inverse_ = Teuchos::rcp( new T(Matrix_.get()) );

  if (Inverse_ == Teuchos::null)
    IFPACK_CHK_ERR(-5);

  IFPACK_CHK_ERR(Inverse_->SetParameters(List_));

  // Call Inverse_->Initialize() in Compute(). This saves
  // some time, because I can extract the diagonal blocks faster,
  // and only once.

  Label_ = "Ifpack_SparseContainer";

  IsInitialized_ = true;
  return(0);

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
SetMatrixElement(const int row, const int col, const double value)
{
  if (!IsInitialized())
    IFPACK_CHK_ERR(-3); // problem not shaped yet

  if ((row < 0) || (row >= NumRows())) {
    IFPACK_CHK_ERR(-2); // not in range
  }

  if ((col < 0) || (col >= NumRows())) {
    IFPACK_CHK_ERR(-2); // not in range
  }

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
   if(Matrix_->RowMatrixRowMap().GlobalIndicesInt()) {
     int ierr = Matrix_->InsertGlobalValues((int)row,1,(double*)&value,(int*)&col);
     if (ierr < 0) {
       ierr = Matrix_->SumIntoGlobalValues((int)row,1,(double*)&value,(int*)&col);
       if (ierr < 0)
         IFPACK_CHK_ERR(-1);
     }
   }
   else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
   if(Matrix_->RowMatrixRowMap().GlobalIndicesLongLong()) {
     long long col_LL = col;
     int ierr = Matrix_->InsertGlobalValues(row,1,(double*)&value,&col_LL);
     if (ierr < 0) {
       ierr = Matrix_->SumIntoGlobalValues(row,1,(double*)&value,&col_LL);
       if (ierr < 0)
         IFPACK_CHK_ERR(-1);
     }
   }
   else
#endif
     throw "Ifpack_SparseContainer<T>::SetMatrixElement: GlobalIndices type unknown";

  return(0);

}

//==============================================================================
template<typename T>
int Ifpack_SparseContainer<T>::Compute(const Epetra_RowMatrix& Matrix_in)
{

  IsComputed_ = false;
  if (!IsInitialized()) {
    IFPACK_CHK_ERR(Initialize());
  }

  // extract the submatrices
  IFPACK_CHK_ERR(Extract(Matrix_in));

  // initialize the inverse operator
  IFPACK_CHK_ERR(Inverse_->Initialize());

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
  if (IsComputed() == false) {
    IFPACK_CHK_ERR(-3); // not yet computed
  }

  IFPACK_CHK_ERR(Matrix_->Apply(*RHS_, *LHS_));

  ApplyFlops_ += 2 * Matrix_->NumGlobalNonzeros64();
  return(0);
}

//==============================================================================
template<typename T>
int Ifpack_SparseContainer<T>::ApplyInverse()
{
  if (!IsComputed())
    IFPACK_CHK_ERR(-1);

  IFPACK_CHK_ERR(Inverse_->ApplyInverse(*RHS_, *LHS_));

  return(0);
}


//==============================================================================
template<typename T>
int Ifpack_SparseContainer<T>::Destroy()
{
  IsInitialized_ = false;
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
  List_ = List;
  return(0);
}

//==============================================================================
// FIXME: optimize performances of this guy...
template<typename T>
int Ifpack_SparseContainer<T>::Extract(const Epetra_RowMatrix& Matrix_in)
{

  for (int j = 0 ; j < NumRows_ ; ++j) {
    // be sure that the user has set all the ID's
    if (ID(j) == -1)
      IFPACK_CHK_ERR(-1);
    // be sure that all are local indices
    if (ID(j) > Matrix_in.NumMyRows())
      IFPACK_CHK_ERR(-1);
  }

  int Length = Matrix_in.MaxNumEntries();
  std::vector<double> Values;
  Values.resize(Length);
  std::vector<int> Indices;
  Indices.resize(Length);

  for (int j = 0 ; j < NumRows_ ; ++j) {

    int LRID = ID(j);

    int NumEntries;

    int ierr =
      Matrix_in.ExtractMyRowCopy(LRID, Length, NumEntries,
                               &Values[0], &Indices[0]);
    IFPACK_CHK_ERR(ierr);

    for (int k = 0 ; k < NumEntries ; ++k) {

      int LCID = Indices[k];

      // skip off-processor elements
      if (LCID >= Matrix_in.NumMyRows())
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

  IFPACK_CHK_ERR(Matrix_->FillComplete());

  return(0);
}

//==============================================================================
template<typename T>
std::ostream& Ifpack_SparseContainer<T>::Print(std::ostream & os) const
{
  using std::endl;

  os << "================================================================================" << endl;
  os << "Ifpack_SparseContainer" << endl;
  os << "Number of rows          = " << NumRows() << endl;
  os << "Number of vectors       = " << NumVectors() << endl;
  os << "IsInitialized()         = " << IsInitialized() << endl;
  os << "IsComputed()            = " << IsComputed() << endl;
  os << "Flops in Initialize()   = " << InitializeFlops() << endl;
  os << "Flops in Compute()      = " << ComputeFlops() << endl;
  os << "Flops in ApplyInverse() = " << ApplyInverseFlops() << endl;
  os << "================================================================================" << endl;
  os << endl;

  return(os);
}
#endif // IFPACK_SPARSECONTAINER_H
