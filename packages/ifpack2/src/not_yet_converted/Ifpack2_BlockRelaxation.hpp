
//@HEADER
// ***********************************************************************
// 
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#ifndef IFPACK2_BLOCKPRECONDITIONER_HPP
#define IFPACK2_BLOCKPRECONDITIONER_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Preconditioner.hpp" 
#include "Ifpack2_Partitioner.hpp"
#include "Ifpack2_LinearPartitioner.hpp"
#include "Ifpack2_GreedyPartitioner.hpp"
#include "Ifpack2_METISPartitioner.hpp"
#include "Ifpack2_EquationPartitioner.hpp"
#include "Ifpack2_UserPartitioner.hpp"
#include "Ifpack2_Graph_Tpetra_RowMatrix.hpp"
#include "Ifpack2_DenseContainer.hpp" 
#include "Ifpack2_Utils.hpp" 
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Time.hpp"
#include "Tpetra_Import.hpp"

static const int IFPACK2_JACOBI = 0;
static const int IFPACK2_GS = 1;
static const int IFPACK2_SGS = 2;


//! Ifpack2_BlockRelaxation: a class to define block relaxation preconditioners of Tpetra_RowMatrix's.

/*! The Ifpack2_BlockRelaxation class enables the construction of
  block relaxation
  preconditioners of an Tpetra_RowMatrix. Ifpack2_PointRelaxation 
  is derived from 
  the Ifpack2_Preconditioner class, which is derived from Tpetra_Operator.
  Therefore this object can be used as preconditioner everywhere an
  ApplyInverse() method is required in the preconditioning step.
 
  The class currently support:
  - block Jacobi;
  - block Gauss-Seidel;
  - symmetric block Gauss-Seidel.
  
  The idea of block relaxation method is to extend their point relaxation
  counterpart (implemented in Ifpack2_PointRelaxation), by working on a
  group of equation simulteneously. Generally, larger blocks result
  in better convergence and increased robusteness.

  The user can decide:
  - the number of blocks (say, NumBlocks). If NumBlocks is equal to the
    number of rows, then the resulting scheme is equivalent to
    a point relaxation scheme;
  - how to apply the inverse of each diagonal block, by choosing a dense
    container or a sparse container. The implementation of
    block relaxation schemes requires the application of the
    inverse of each diagonal block. This can be done using LAPACK (dense 
    container), or any Ifpack2_Preconditioner derived class (sparse
    container);
  - blocks can be defined using a linear decomposition, by a simple greedy
    algorithm, or by resorting to METIS.

The following is an example of usage of this preconditioner with dense
containers. First, we include the header files:
\code
#include "Ifpack2_AdditiveSchwarz.hpp"
#include "Ifpack2_BlockPreconditioner.hpp"
#include "Ifpack2_DenseContainer.hpp"
\endcode

Then, we declare the preconditioner. Note that this is done through
the class Ifpack2_AdditiveSchwarz (see note below in this section).
\code
// A is an Tpetra_RowMatrix
// List is a Teuchos::ParameterList
Ifpack2_AdditiveSchwarz<Ifpack2_BlockRelaxation<Ifpack2_DenseContainer> > > Prec(A);
IFPACK2_CHK_ERR(Prec.SetParameters(List));
IFPACK2_CHK_ERR(Prec.Initialize());
IFPACK2_CHK_ERR(Prec.Compute());

// action of the preconditioner is given by ApplyInverse()
// Now use it in AztecOO, solver is an AztecOO object
solver.SetPrecOperator(&Prec);
\endcode

<P>The complete list of supported parameters is reported in page \ref ifp_params. For a presentation of basic relaxation schemes, please refer to page
\ref Ifpack2_PointRelaxation.

\author Michael Heroux, SNL 9214.

\date Last modified on 25-Jan-05.
  
*/
template<typename T>
class Ifpack2_BlockRelaxation : public Ifpack2_Preconditioner {

public:

  //@{ \name Constructors/Destructors
  //! Ifpack2_BlockRelaxation constructor with given Tpetra_RowMatrix.
  /*! Creates an Ifpack2_Preconditioner preconditioner. 
   *
   * \param In
   * Matrix - Pointer to matrix to be preconditioned.
   */
  Ifpack2_BlockRelaxation(const Tpetra_RowMatrix* Matrix);

  virtual ~Ifpack2_BlockRelaxation();

  //@}

  //@{ \name Mathematical functions.

  //! Applies the matrix to an Tpetra_MultiVector.
  /*! 
    \param In
    X - A Tpetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
    Y -A Tpetra_MultiVector of dimension NumVectors containing the result.

    \return Integer error code, set to 0 if successful.
    */
  virtual int Apply(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
		    Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  //! Applies the block Jacobi preconditioner to X, returns the result in Y.
  /*! 
    \param In
    X - A Tpetra_MultiVector of dimension NumVectors to be preconditioned.
    \param Out
    Y -A Tpetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    */
  virtual int ApplyInverse(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
			   Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  //! Returns the infinity norm of the global matrix (not implemented)
  virtual double NormInf() const
  {
    return(-1.0);
  }
  //@}

  //@{ \name Atribute access functions

  virtual int SetUseTranspose(bool UseTranspose_in)
  {
    if (UseTranspose_in)
      IFPACK2_CHK_ERR(-98); // FIXME: can I work with the transpose?
    return(0);
  }

  virtual const char* Label() const;
 
  //! Returns the current UseTranspose setting.
  virtual bool UseTranspose() const
  {
    return(false);
  }

  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  virtual bool HasNormInf() const
  {
    return(false);
  }

  //! Returns a pointer to the Tpetra_Comm communicator associated with this operator.
  virtual const Tpetra_Comm & Comm() const;

  //! Returns the Tpetra_Map object associated with the domain of this operator.
  virtual const Tpetra_Map & OperatorDomainMap() const;

  //! Returns the Tpetra_Map object associated with the range of this operator.
  virtual const Tpetra_Map & OperatorRangeMap() const;
  //@}

  //! Returns the number local blocks.
  int NumLocalBlocks() const 
  {
    return(NumLocalBlocks_);
  }

  //! Returns \c true if the preconditioner has been successfully computed.
  virtual bool IsInitialized() const
  {
    return(IsInitialized_);
  }

  //! Returns \c true if the preconditioner has been successfully computed.
  virtual bool IsComputed() const
  {
    return(IsComputed_);
  }

  //! Sets all the parameters for the preconditioner.
  virtual int SetParameters(Teuchos::ParameterList& List);

  //! Initializes the preconditioner.
  virtual int Initialize();

  //! Computes the preconditioner.
  virtual int Compute();

  virtual const Tpetra_RowMatrix& Matrix() const
  {
    return(*Matrix_);
  }

  virtual double Condest(const Ifpack2_CondestType CT = Ifpack2_Cheap,
                         const int MaxIters = 1550,
                         const double Tol = 1e-9,
			 Tpetra_RowMatrix* Matrix_in = 0)
  {
    return(-1.0);
  }

  virtual double Condest() const
  {
    return(-1.0);
  }

  std::ostream& Print(std::ostream& os) const;

  //! Returns the number of calls to Initialize().
  virtual int NumInitialize() const
  {
    return(NumInitialize_);
  }

  //! Returns the number of calls to Compute().
  virtual int NumCompute() const
  {
    return(NumCompute_);
  }

  //! Returns the number of calls to ApplyInverse().
  virtual int NumApplyInverse() const
  {
    return(NumApplyInverse_);
  }

  //! Returns the time spent in Initialize().
  virtual double InitializeTime() const
  {
    return(InitializeTime_);
  }

  //! Returns the time spent in Compute().
  virtual double ComputeTime() const
  {
    return(ComputeTime_);
  }

  //! Returns the time spent in ApplyInverse().
  virtual double ApplyInverseTime() const
  {
    return(ApplyInverseTime_);
  }

  //! Returns the number of flops in the initialization phase.
  virtual double InitializeFlops() const
  {
    if (Containers_.size() == 0)
      return(0.0);

    // the total number of flops is computed each time InitializeFlops() is
    // called. This is becase I also have to add the contribution from each
    // container.
    double total = InitializeFlops_;
    for (unsigned int i = 0 ; i < Containers_.size() ; ++i)
      total += Containers_[i]->InitializeFlops();
    return(total);
  }

  virtual double ComputeFlops() const
  {
    if (Containers_.size() == 0)
      return(0.0);
    
    double total = ComputeFlops_;
    for (unsigned int i = 0 ; i < Containers_.size() ; ++i)
      total += Containers_[i]->ComputeFlops();
    return(total);
  }

  virtual double ApplyInverseFlops() const
  {
    if (Containers_.size() == 0)
      return(0.0);

    double total = ApplyInverseFlops_;
    for (unsigned int i = 0 ; i < Containers_.size() ; ++i) {
      total += Containers_[i]->ApplyInverseFlops();
    }
    return(total);
  }

private:

  //! Copy constructor (PRIVATE, should not be used).
  Ifpack2_BlockRelaxation(const Ifpack2_BlockRelaxation& rhs);

  //! operator= (PRIVATE, should not be used).
  Ifpack2_BlockRelaxation & operator=(const Ifpack2_BlockRelaxation& rhs)
  {
    return(*this);
  }

  virtual int ApplyInverseJacobi(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
                                 Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  virtual int DoJacobi(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
                                  Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  virtual int ApplyInverseGS(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
                             Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  virtual int DoGaussSeidel(Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
                            Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  virtual int ApplyInverseSGS(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
                              Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  virtual int DoSGS(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                    Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Xtmp,
                    Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  int ExtractSubmatrices();

  // @{ Initializations, timing and flops

  //! If true, the preconditioner has been successfully initialized.
  bool IsInitialized_;
  //! If true, the preconditioner has been successfully computed.
  bool IsComputed_;
  //! Contains the number of successful calls to Initialize().
  int NumInitialize_;
  //! Contains the number of successful call to Compute().
  int NumCompute_;
  //! Contains the number of successful call to ApplyInverse().
  mutable int NumApplyInverse_;
  //! Contains the time for all successful calls to Initialize().
  double InitializeTime_;
  //! Contains the time for all successful calls to Compute().
  double ComputeTime_;
  //! Contains the time for all successful calls to ApplyInverse().
  mutable double ApplyInverseTime_;
  //! Contains the number of flops for Initialize().
  double InitializeFlops_;
  //! Contains the number of flops for Compute().
  double ComputeFlops_;
  //! Contain sthe number of flops for ApplyInverse().
  mutable double ApplyInverseFlops_;
  // @}

  // @{ Settings
  //! Number of preconditioning sweeps.
  int NumSweeps_;
  //! Damping parameter.
  double DampingFactor_;
  //! Number of local blocks
  int NumLocalBlocks_;
  //! Parameters list to be used to solve on each subblock
  Teuchos::ParameterList List_;
  // @}

  // @{ Other data
  //! Containers_[i] contains all the necessary information to solve on each subblock.
  //! Pointers to the matrix to be preconditioned.
  Teuchos::RCP< const Tpetra_RowMatrix > Matrix_;
  mutable std::vector<Teuchos::RCP<T> > Containers_;
  //! Contains information about non-overlapping partitions.
  Teuchos::RCP<Ifpack2_Partitioner> Partitioner_;
  string PartitionerType_;
  int PrecType_;
  //! Label for \c this object
  string Label_;
  //! If \c true, starting solution is the zero vector.
  bool ZeroStartingSolution_;
  Teuchos::RCP<Ifpack2_Graph> Graph_;
  //! Weights for overlapping Jacobi only.
  Teuchos::RCP<Tpetra_Vector> W_;
  // Level of overlap among blocks (for Jacobi only).
  int OverlapLevel_;
  mutable Tpetra_Time Time_;
  bool IsParallel_;
  Teuchos::RCP<Tpetra_Import> Importer_;
  // @}
  
}; // class Ifpack2_BlockRelaxation

//==============================================================================
template<typename T>
Ifpack2_BlockRelaxation<T>::
Ifpack2_BlockRelaxation(const Tpetra_RowMatrix* Matrix_in) :
  IsInitialized_(false),
  IsComputed_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyInverseTime_(0.0),
  InitializeFlops_(0.0),
  ComputeFlops_(0.0),
  ApplyInverseFlops_(0.0),
  NumSweeps_(1),
  DampingFactor_(1.0),
  NumLocalBlocks_(1),
  Matrix_(Teuchos::rcp(Matrix_in,false)),
  PartitionerType_("greedy"),
  PrecType_(IFPACK2_JACOBI),
  ZeroStartingSolution_(true),
  OverlapLevel_(0),
  Time_(Comm()),
  IsParallel_(false)
{
  if (Matrix_in->Comm().NumProc() != 1)
    IsParallel_ = true;
}

//==============================================================================
template<typename T>
Ifpack2_BlockRelaxation<T>::~Ifpack2_BlockRelaxation()
{
}

//==============================================================================
template<typename T>
const char* Ifpack2_BlockRelaxation<T>::Label() const
{
  return(Label_.c_str());
}

//==============================================================================
template<typename T>
int Ifpack2_BlockRelaxation<T>::
Apply(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  IFPACK2_RETURN(Matrix().Apply(X,Y));
}

//==============================================================================
template<typename T>
const Tpetra_Comm& Ifpack2_BlockRelaxation<T>::
Comm() const
{
  return(Matrix().Comm());
}

//==============================================================================
template<typename T>
const Tpetra_Map& Ifpack2_BlockRelaxation<T>::
OperatorDomainMap() const
{
  return(Matrix().OperatorDomainMap());
}

//==============================================================================
template<typename T>
const Tpetra_Map& Ifpack2_BlockRelaxation<T>::
OperatorRangeMap() const
{
  return(Matrix().OperatorRangeMap());
}

//==============================================================================
template<typename T>
int Ifpack2_BlockRelaxation<T>::ExtractSubmatrices()
{

  if (Partitioner_ == Teuchos::null)
    IFPACK2_CHK_ERR(-3);

  NumLocalBlocks_ = Partitioner_->NumLocalParts();

  Containers_.resize(NumLocalBlocks());

  for (int i = 0 ; i < NumLocalBlocks() ; ++i) {

    int rows = Partitioner_->NumRowsInPart(i);
    Containers_[i] = Teuchos::rcp( new T(rows) );
    
    //Ifpack2_DenseContainer* DC = 0;
    //DC = dynamic_cast<Ifpack2_DenseContainer*>(Containers_[i]);

    if (Containers_[i] == Teuchos::null)
      IFPACK2_CHK_ERR(-5);
    
    IFPACK2_CHK_ERR(Containers_[i]->SetParameters(List_));
    IFPACK2_CHK_ERR(Containers_[i]->Initialize());
    // flops in Initialize() will be computed on-the-fly in method InitializeFlops().

    // set "global" ID of each partitioner row
    for (int j = 0 ; j < rows ; ++j) {
      int LRID = (*Partitioner_)(i,j);
      Containers_[i]->ID(j) = LRID;
    }

    IFPACK2_CHK_ERR(Containers_[i]->Compute(*Matrix_));
    // flops in Compute() will be computed on-the-fly in method ComputeFlops().

  }

  return(0);
}

//==============================================================================
template<typename T>
int Ifpack2_BlockRelaxation<T>::Compute()
{

  if (!IsInitialized())
    IFPACK2_CHK_ERR(Initialize());

  Time_.ResetStartTime();

  IsComputed_ = false;

  if (Matrix().NumGlobalRows() != Matrix().NumGlobalCols())
    IFPACK2_CHK_ERR(-2); // only square matrices

  IFPACK2_CHK_ERR(ExtractSubmatrices());
  
  if (IsParallel_ && PrecType_ != IFPACK2_JACOBI) {
    // not needed by Jacobi (done by matvec)
    Importer_ = Teuchos::rcp( new Tpetra_Import(Matrix().RowMatrixColMap(),
                                                Matrix().RowMatrixRowMap()) );

    if (Importer_ == Teuchos::null) IFPACK2_CHK_ERR(-5);
  }
  IsComputed_ = true;
  ComputeTime_ += Time_.ElapsedTime();
  ++NumCompute_;

  return(0);

}

//==============================================================================
template<typename T>
int Ifpack2_BlockRelaxation<T>::
ApplyInverse(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  if (!IsComputed())
    IFPACK2_CHK_ERR(-3);

  if (X.NumVectors() != Y.NumVectors())
    IFPACK2_CHK_ERR(-2);

  Time_.ResetStartTime();

  // AztecOO gives X and Y pointing to the same memory location,
  // need to create an auxiliary vector, Xcopy
  Teuchos::RCP< const Tpetra_MultiVector > Xcopy;
  if (X.Pointers()[0] == Y.Pointers()[0])
    Xcopy = Teuchos::rcp( new Tpetra_MultiVector(X) );
  else
    Xcopy = Teuchos::rcp( &X, false );

  switch (PrecType_) {
  case IFPACK2_JACOBI:
    IFPACK2_CHK_ERR(ApplyInverseJacobi(*Xcopy,Y));
    break;
  case IFPACK2_GS:
    IFPACK2_CHK_ERR(ApplyInverseGS(*Xcopy,Y));
    break;
  case IFPACK2_SGS:
    IFPACK2_CHK_ERR(ApplyInverseSGS(*Xcopy,Y));
    break;
  }

  ApplyInverseTime_ += Time_.ElapsedTime();
  ++NumApplyInverse_;

  return(0);
}

//==============================================================================
// This method in general will not work with AztecOO if used
// outside Ifpack2_AdditiveSchwarz and OverlapLevel_ != 0
//
template<typename T>
int Ifpack2_BlockRelaxation<T>::
ApplyInverseJacobi(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
                   Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{

  if (ZeroStartingSolution_)
    Y.PutScalar(0.0);

  // do not compute the residual in this case
  if (NumSweeps_ == 1 && ZeroStartingSolution_) {
    IFPACK2_RETURN(DoJacobi(X,Y));
  }

  Tpetra_MultiVector AX(Y);

  for (int j = 0; j < NumSweeps_ ; j++) {
    IFPACK2_CHK_ERR(Apply(Y,AX));
    ApplyInverseFlops_ += X.NumVectors() * 2 * Matrix_->NumGlobalNonzeros();
    IFPACK2_CHK_ERR(AX.Update(1.0,X,-1.0));
    ApplyInverseFlops_ += X.NumVectors() * 2 * Matrix_->NumGlobalRows();
    IFPACK2_CHK_ERR(DoJacobi(AX,Y));
    // flops counted in DoJacobi()
  }


  return(0);
}

//==============================================================================
template<typename T>
int Ifpack2_BlockRelaxation<T>::
DoJacobi(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  int NumVectors = X.NumVectors();

  if (OverlapLevel_ == 0) {

    for (int i = 0 ; i < NumLocalBlocks() ; ++i) {
     
      // may happen that a partition is empty
      if (Containers_[i]->NumRows() == 0) 
        continue;

      int LID;

      // extract RHS from X
      for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
        LID = Containers_[i]->ID(j);
        for (int k = 0 ; k < NumVectors ; ++k) {
          Containers_[i]->RHS(j,k) = X[k][LID];
        }
      }

      // apply the inverse of each block. NOTE: flops occurred
      // in ApplyInverse() of each block are summed up in method
      // ApplyInverseFlops().
      IFPACK2_CHK_ERR(Containers_[i]->ApplyInverse());

      // copy back into solution vector Y
      for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
        LID = Containers_[i]->ID(j);
        for (int k = 0 ; k < NumVectors ; ++k) {
          Y[k][LID] += DampingFactor_ * Containers_[i]->LHS(j,k);
        }
      }

    }
    // NOTE: flops for ApplyInverse() of each block are summed up
    // in method ApplyInverseFlops()
    ApplyInverseFlops_ += NumVectors * 2 * Matrix_->NumGlobalRows();

  }
  else {

    for (int i = 0 ; i < NumLocalBlocks() ; ++i) {

      // may happen that a partition is empty
      if (Containers_[i]->NumRows() == 0) 
        continue;

      int LID;

      // extract RHS from X
      for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
        LID = Containers_[i]->ID(j);
        for (int k = 0 ; k < NumVectors ; ++k) {
          Containers_[i]->RHS(j,k) = (*W_)[LID] * X[k][LID];
        }
      }

      // apply the inverse of each block
      IFPACK2_CHK_ERR(Containers_[i]->ApplyInverse());

      // copy back into solution vector Y
      for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
        LID = Containers_[i]->ID(j);
        for (int k = 0 ; k < NumVectors ; ++k) {
          Y[k][LID] += DampingFactor_ * (*W_)[LID] * Containers_[i]->LHS(j,k);
        }
      }

    }
    // NOTE: flops for ApplyInverse() of each block are summed up
    // in method ApplyInverseFlops()
    // NOTE: do not count for simplicity the flops due to overlapping rows
    ApplyInverseFlops_ += NumVectors * 4 * Matrix_->NumGlobalRows();
  }

  return(0);
}

//==============================================================================
template<typename T>
int Ifpack2_BlockRelaxation<T>::
ApplyInverseGS(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
               Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{

  if (ZeroStartingSolution_)
    Y.PutScalar(0.0);

  Tpetra_MultiVector Xcopy(X);
  for (int j = 0; j < NumSweeps_ ; j++) {
    IFPACK2_CHK_ERR(DoGaussSeidel(Xcopy,Y));
    if (j != NumSweeps_ - 1)
      Xcopy = X;
  }

  return(0);

}

//==============================================================================
template<typename T>
int Ifpack2_BlockRelaxation<T>::
DoGaussSeidel(Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{

  // cycle over all local subdomains

  int Length = Matrix().MaxNumEntries();
  std::vector<int> Indices(Length);
  std::vector<double> Values(Length);

  int NumMyRows = Matrix().NumMyRows();
  int NumVectors = X.NumVectors();

  // an additonal vector is needed by parallel computations
  // (note that applications through Ifpack2_AdditiveSchwarz
  // are always seen are serial)
  Teuchos::RCP< Tpetra_MultiVector > Y2;
  if (IsParallel_)
    Y2 = Teuchos::rcp( new Tpetra_MultiVector(Importer_->TargetMap(), NumVectors) );
  else
    Y2 = Teuchos::rcp( &Y, false );

  double** y_ptr;
  double** y2_ptr;
  Y.ExtractView(&y_ptr);
  Y2->ExtractView(&y2_ptr);

  // data exchange is here, once per sweep
  if (IsParallel_)
    IFPACK2_CHK_ERR(Y2->Import(Y,*Importer_,Insert));

  for (int i = 0 ; i < NumLocalBlocks() ; ++i) {

    // may happen that a partition is empty
    if (Containers_[i]->NumRows() == 0) 
      continue;

    int LID;

    // update from previous block

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);

      int NumEntries;
      IFPACK2_CHK_ERR(Matrix().ExtractMyRowCopy(LID, Length,NumEntries,
                                               &Values[0], &Indices[0]));

      for (int k = 0 ; k < NumEntries ; ++k) {
        int col = Indices[k];

          for (int kk = 0 ; kk < NumVectors ; ++kk) {
            X[kk][LID] -= Values[k] * y2_ptr[kk][col];
          }
      }
    }

    // solve with this block

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);
      for (int k = 0 ; k < NumVectors ; ++k) {
        Containers_[i]->RHS(j,k) = X[k][LID];
      }
    }

    IFPACK2_CHK_ERR(Containers_[i]->ApplyInverse());
    ApplyInverseFlops_ += Containers_[i]->ApplyInverseFlops();

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);
      for (int k = 0 ; k < NumVectors ; ++k) {
        y2_ptr[k][LID] += DampingFactor_ * Containers_[i]->LHS(j,k);
      }
    }

  }

  // operations for all getrow()'s
  // NOTE: flops for ApplyInverse() of each block are summed up
  // in method ApplyInverseFlops()
  ApplyInverseFlops_ += NumVectors * 2 * Matrix_->NumGlobalNonzeros();
  ApplyInverseFlops_ += NumVectors * 2 * Matrix_->NumGlobalRows();

  // Attention: this is delicate... Not all combinations
  // of Y2 and Y will always work (tough for ML it should be ok)
  if (IsParallel_)
    for (int m = 0 ; m < NumVectors ; ++m) 
      for (int i = 0 ; i < NumMyRows ; ++i)
        y_ptr[m][i] = y2_ptr[m][i];

  return(0);
}

//==============================================================================
template<typename T>
int Ifpack2_BlockRelaxation<T>::
ApplyInverseSGS(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{

  if (ZeroStartingSolution_)
    Y.PutScalar(0.0);

  Tpetra_MultiVector Xcopy(X);
  for (int j = 0; j < NumSweeps_ ; j++) {
    IFPACK2_CHK_ERR(DoSGS(X,Xcopy,Y));
    if (j != NumSweeps_ - 1)
      Xcopy = X;
  }
  return(0);
}

//==============================================================================
template<typename T>
int Ifpack2_BlockRelaxation<T>::
DoSGS(const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Xcopy, 
      Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{

  int NumMyRows = Matrix().NumMyRows();
  int NumVectors = X.NumVectors();

  int Length = Matrix().MaxNumEntries();
  std::vector<int> Indices;
  std::vector<double> Values;
  Indices.resize(Length);
  Values.resize(Length);

  // an additonal vector is needed by parallel computations
  // (note that applications through Ifpack2_AdditiveSchwarz
  // are always seen are serial)
  Teuchos::RCP< Tpetra_MultiVector > Y2;
  if (IsParallel_)
    Y2 = Teuchos::rcp( new Tpetra_MultiVector(Importer_->TargetMap(), NumVectors) );
  else
    Y2 = Teuchos::rcp( &Y, false );

  double** y_ptr;
  double** y2_ptr;
  Y.ExtractView(&y_ptr);
  Y2->ExtractView(&y2_ptr);

  // data exchange is here, once per sweep
  if (IsParallel_)
    IFPACK2_CHK_ERR(Y2->Import(Y,*Importer_,Insert));

  for (int i = 0 ; i < NumLocalBlocks() ; ++i) {

    // may happen that a partition is empty
    if (Containers_[i]->NumRows() == 0) 
      continue;

    int LID;

    // update from previous block

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);

      int NumEntries;
      IFPACK2_CHK_ERR(Matrix().ExtractMyRowCopy(LID, Length,NumEntries,
                                               &Values[0], &Indices[0]));

      for (int k = 0 ; k < NumEntries ; ++k) {
        int col = Indices[k];

        for (int kk = 0 ; kk < NumVectors ; ++kk) {
          Xcopy[kk][LID] -= Values[k] * y2_ptr[kk][col];
        }
      }
    }

    // solve with this block

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);
      for (int k = 0 ; k < NumVectors ; ++k) {
        Containers_[i]->RHS(j,k) = Xcopy[k][LID];
      }
    }

    IFPACK2_CHK_ERR(Containers_[i]->ApplyInverse());
    ApplyInverseFlops_ += Containers_[i]->ApplyInverseFlops();

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);
      for (int k = 0 ; k < NumVectors ; ++k) {
        y2_ptr[k][LID] += DampingFactor_ * Containers_[i]->LHS(j,k);
      }
    }
  }

  // operations for all getrow()'s
  ApplyInverseFlops_ += NumVectors * 2 * Matrix_->NumGlobalNonzeros();
  ApplyInverseFlops_ += NumVectors * 2 * Matrix_->NumGlobalRows();

  Xcopy = X;

  for (int i = NumLocalBlocks() - 1; i >=0 ; --i) {

    if (Containers_[i]->NumRows() == 0) 
      continue;

    int LID;

    // update from previous block

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);

      int NumEntries;
      IFPACK2_CHK_ERR(Matrix().ExtractMyRowCopy(LID, Length,NumEntries,
                                               &Values[0], &Indices[0]));

      for (int k = 0 ; k < NumEntries ; ++k) {
        int col = Indices[k];

          for (int kk = 0 ; kk < NumVectors ; ++kk) {
            Xcopy[kk][LID] -= Values[k] * y2_ptr[kk][col];
          }
      }
    }

    // solve with this block

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);
      for (int k = 0 ; k < NumVectors ; ++k) {
        Containers_[i]->RHS(j,k) = Xcopy[k][LID];
      }
    }

    IFPACK2_CHK_ERR(Containers_[i]->ApplyInverse());
    ApplyInverseFlops_ += Containers_[i]->ApplyInverseFlops();

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);
      for (int k = 0 ; k < NumVectors ; ++k) {
        y2_ptr[k][LID] += DampingFactor_ * Containers_[i]->LHS(j,k);
      }
    }
  }

  // operations for all getrow()'s
  ApplyInverseFlops_ += NumVectors * 2 * Matrix_->NumGlobalNonzeros();
  ApplyInverseFlops_ += NumVectors * 2 * Matrix_->NumGlobalRows();

  // Attention: this is delicate... Not all combinations
  // of Y2 and Y will always work (tough for ML it should be ok)
  if (IsParallel_)
    for (int m = 0 ; m < NumVectors ; ++m) 
      for (int i = 0 ; i < NumMyRows ; ++i)
        y_ptr[m][i] = y2_ptr[m][i];

  return(0);
}

//==============================================================================
template<typename T>
ostream& Ifpack2_BlockRelaxation<T>::Print(ostream & os) const
{

  string PT;
  if (PrecType_ == IFPACK2_JACOBI)
    PT = "Jacobi";
  else if (PrecType_ == IFPACK2_GS)
    PT = "Gauss-Seidel";
  else if (PrecType_ == IFPACK2_SGS)
    PT = "symmetric Gauss-Seidel";

  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack2_BlockRelaxation, " << PT << endl;
    os << "Sweeps = " << NumSweeps_ << endl;
    os << "Damping factor = " << DampingFactor_;
    if (ZeroStartingSolution_) 
      os << ", using zero starting solution" << endl;
    else
      os << ", using input starting solution" << endl;
    os << "Number of local blocks = " << Partitioner_->NumLocalParts() << endl;
    //os << "Condition number estimate = " << Condest_ << endl;
    os << "Global number of rows            = " << Matrix_->NumGlobalRows() << endl;
    os << endl;
    os << "Phase           # calls   Total Time (s)       Total MFlops     MFlops/s" << endl;
    os << "-----           -------   --------------       ------------     --------" << endl;
    os << "Initialize()    "   << std::setw(5) << NumInitialize() 
       << "  " << std::setw(15) << InitializeTime() 
       << "  " << std::setw(15) << 1.0e-6 * InitializeFlops();
    if (InitializeTime() != 0.0)
      os << "  " << std::setw(15) << 1.0e-6 * InitializeFlops() / InitializeTime() << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "Compute()       "   << std::setw(5) << NumCompute() 
       << "  " << std::setw(15) << ComputeTime()
       << "  " << std::setw(15) << 1.0e-6 * ComputeFlops();
    if (ComputeTime() != 0.0) 
      os << "  " << std::setw(15) << 1.0e-6 * ComputeFlops() / ComputeTime() << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "ApplyInverse()  "   << std::setw(5) << NumApplyInverse() 
       << "  " << std::setw(15) << ApplyInverseTime()
       << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops();
    if (ApplyInverseTime() != 0.0) 
      os << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops() / ApplyInverseTime() << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "================================================================================" << endl;
    os << endl;
  }

  return(os);
}

//==============================================================================
template<typename T>
int Ifpack2_BlockRelaxation<T>::SetParameters(Teuchos::ParameterList& List)
{

  string PT;
  if (PrecType_ == IFPACK2_JACOBI)
    PT = "Jacobi";
  else if (PrecType_ == IFPACK2_GS)
    PT = "Gauss-Seidel";
  else if (PrecType_ == IFPACK2_SGS)
    PT = "symmetric Gauss-Seidel";

  PT = List.get("relaxation: type", PT);

  if (PT == "Jacobi") {
    PrecType_ = IFPACK2_JACOBI;
  }
  else if (PT == "Gauss-Seidel") {
    PrecType_ = IFPACK2_GS;
  }
  else if (PT == "symmetric Gauss-Seidel") {
    PrecType_ = IFPACK2_SGS;
  } else {
    cerr << "Option `relaxation: type' has an incorrect value ("
      << PT << ")'" << endl;
    cerr << "(file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
    exit(EXIT_FAILURE);
  }

  NumSweeps_            = List.get("relaxation: sweeps", NumSweeps_);
  DampingFactor_        = List.get("relaxation: damping factor", 
                                   DampingFactor_);
  ZeroStartingSolution_ = List.get("relaxation: zero starting solution", 
                                   ZeroStartingSolution_);
  PartitionerType_      = List.get("partitioner: type", 
                                   PartitionerType_);
  NumLocalBlocks_       = List.get("partitioner: local parts", 
                                   NumLocalBlocks_);
  // only Jacobi can work with overlap among local domains,
  OverlapLevel_         = List.get("partitioner: overlap", 
                                   OverlapLevel_);

  // check parameters
  if (PrecType_ != IFPACK2_JACOBI)
    OverlapLevel_ = 0;
  if (NumLocalBlocks_ < 0)
    NumLocalBlocks_ = Matrix().NumMyRows() / (-NumLocalBlocks_);
  // other checks are performed in Partitioner_
  
  // copy the list as each subblock's constructor will
  // require it later
  List_ = List;

  // set the label
  string PT2;
  if (PrecType_ == IFPACK2_JACOBI)
    PT2 = "BJ";
  else if (PrecType_ == IFPACK2_GS)
    PT2 = "BGS";
  else if (PrecType_ == IFPACK2_SGS)
    PT2 = "BSGS";
  Label_ = "TIFPACK (" + PT2 + ", sweeps=" 
    + Ifpack2_toString(NumSweeps_) + ", damping="
    + Ifpack2_toString(DampingFactor_) + ", blocks="
    + Ifpack2_toString(NumLocalBlocks()) + ")";

  return(0);
}

//==============================================================================
template<typename T>
int Ifpack2_BlockRelaxation<T>::Initialize()
{
  IsInitialized_ = false;
  Time_.ResetStartTime();

  Graph_ = Teuchos::rcp( new Ifpack2_Graph_Tpetra_RowMatrix(Teuchos::rcp(&Matrix(),false)) );
  if (Graph_ == Teuchos::null) IFPACK2_CHK_ERR(-5);

  if (PartitionerType_ == "linear")
    Partitioner_ = Teuchos::rcp( new Ifpack2_LinearPartitioner(&*Graph_) );
  else if (PartitionerType_ == "greedy")
    Partitioner_ = Teuchos::rcp( new Ifpack2_GreedyPartitioner(&*Graph_) );
  else if (PartitionerType_ == "metis")
    Partitioner_ = Teuchos::rcp( new Ifpack2_METISPartitioner(&*Graph_) );
  else if (PartitionerType_ == "equation")
    Partitioner_ = Teuchos::rcp( new Ifpack2_EquationPartitioner(&*Graph_) );
  else if (PartitionerType_ == "user")
    Partitioner_ = Teuchos::rcp( new Ifpack2_UserPartitioner(&*Graph_) );
  else
    IFPACK2_CHK_ERR(-2);

  if (Partitioner_ == Teuchos::null) IFPACK2_CHK_ERR(-5);

  // need to partition the graph of A
  IFPACK2_CHK_ERR(Partitioner_->SetParameters(List_));
  IFPACK2_CHK_ERR(Partitioner_->Compute());

  // get actual number of partitions
  NumLocalBlocks_ = Partitioner_->NumLocalParts();

  // weight of each vertex
  W_ = Teuchos::rcp( new Tpetra_Vector(Matrix().RowMatrixRowMap()) );
  W_->PutScalar(0.0);

  for (int i = 0 ; i < NumLocalBlocks() ; ++i) {

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      int LID = (*Partitioner_)(i,j);
      (*W_)[LID]++;
    }
  }
  W_->Reciprocal(*W_);

  InitializeTime_ += Time_.ElapsedTime();
  IsInitialized_ = true;
  ++NumInitialize_;

  return(0);
}

//==============================================================================
#endif // IFPACK2_BLOCKPRECONDITIONER_HPP
