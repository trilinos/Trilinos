#ifndef IFPACK_BLOCKPRECONDITIONER_H
#define IFPACK_BLOCKPRECONDITIONER_H

#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_IFPACK_TEUCHOS
#include "Ifpack_Preconditioner.h" 
#include "Ifpack_Partitioner.h"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_RowMatrix.h"
#include "Epetra_MultiVector.h"

//! Ifpack_BlockPreconditioner: a class to define block preconditioners preconditioners of Epetra_RowMatrix's.

/*! The Ifpack_BlockPreconditioner class enables the construction of
  block-preconditioners
  preconditioners of an Epetra_RowMatrix. The blocks can have arbitrary
  sizes. Ifpack_Partitioner is used to create the blocks (that is, to
  partition the graph of the local matrix). Hence, valid partitioning schemes 
  are:
  - a linear decomposition;
  - a simple greedy algorithm;
  - METIS.

 \date Sep-04
  
 */
template<typename T>
class Ifpack_BlockPreconditioner : public Ifpack_Preconditioner {

public:

  //@{ \name Constructors/Destructors
  //! Ifpack_BlockPreconditioner constructor with given Epetra_RowMatrix.
  /*! Creates an Ifpack_Preconditioner preconditioner. 
   *
   * \param In
   * Matrix - Pointer to matrix to be preconditioned.
   */
  Ifpack_BlockPreconditioner(Epetra_RowMatrix* Matrix);

  virtual ~Ifpack_BlockPreconditioner();

  //@}

  //@{ \name Mathematical functions.

  //! Applies the matrix to an Epetra_MultiVector.
  /*! 
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectors containing the result.

    \return Integer error code, set to 0 if successful.
    */
  virtual int Apply(const Epetra_MultiVector& X, 
		    Epetra_MultiVector& Y) const;

  //! Applies the block Jacobi preconditioner to X, returns the result in Y.
  /*! 
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to be preconditioned.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning In order to work with AztecOO, any implementation of this method 
    must support the case where X and Y are the same object.
    */
  virtual int ApplyInverse(const Epetra_MultiVector& X, 
			   Epetra_MultiVector& Y) const = 0;

  //! Returns the infinity norm of the global matrix (not implemented)
  virtual double NormInf() const
  {
    return(-1.0);
  }

  //@}

  //@{ \name Atribute access functions

  virtual int SetUseTranspose(bool UseTranspose)
  {
    IFPACK_CHK_ERR(-1); // FIXME: can I work with the transpose?
  }

  virtual char * Label() const;
 
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

  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  virtual const Epetra_Comm & Comm() const;

  //! Returns the Epetra_Map object associated with the domain of this operator.
  virtual const Epetra_Map & OperatorDomainMap() const;

  //! Returns the Epetra_Map object associated with the range of this operator.
  virtual const Epetra_Map & OperatorRangeMap() const;
  //@}

  //! Returns the number local blocks.
  int NumLocalBlocks() const 
  {
    return(NumLocalBlocks_);
  }

  //! Sets the number of local blocks.
  int SetNumLocalBlocks(const int NumLocalBlocks) const 
  {
    NumLocalBlocks_ = NumLocalBlocks;
    return(0);
  }
  //! Returns \c true if the preconditioner has been successfully computed.
  virtual bool IsComputed() const
  {
    return(IsComputed_);
  }

  //! Sets all the parameters for the preconditioner.
  virtual int SetParameters(Teuchos::ParameterList& List)
  {
    
    // copy the list as each subblock's constructor will
    // require it later
    List_ = List;

    SetNumSweeps(List.get("sweeps",NumSweeps()));
    SetDampingFactor(List.get("omega", DampingFactor()));
    PrintLevel_ = List.get("print level", PrintLevel());
    Partitioner_ = List.get("partitioner object", (Ifpack_Partitioner*)0);

    // derived class will set the appropriate label
    IFPACK_CHK_ERR(SetLabel());
    return(0);
  }

  //! Computes the preconditioners.
  /*! Computes the preconditioners: 
   *
   * \param In
   * List - list specifying the parameters for Jacobi. See above.
   *    
   * \return
   * 0 if successful, 1 if a zero element has been found on the
   * diagonal. In this latter case, the preconditioner can still be
   * applied, as the inverse of zero elements is replaced by 1.0.
   */
  virtual int Compute();

  //! Sets the number of sweeps.
  int SetNumSweeps(const int NumSweeps)
  {
    NumSweeps_ = NumSweeps;
    return(0);
  }

  //! Gets the number of sweeps.
  int NumSweeps() const
  {
    return(NumSweeps_);
  }
 
  //! Sets the damping factor.
  int SetDampingFactor(const double DampingFactor)
  {
    DampingFactor_ = DampingFactor;
    return(0);
  }

  //! Gets the damping parameter.
  double DampingFactor() const
  {
    return(DampingFactor_);
  }

  virtual Epetra_RowMatrix* Matrix() const
  {
    return(Matrix_);
  }

  virtual int PrintLevel() const
  {
    return(PrintLevel_);
  }

  virtual int SetPrintLevel(const int PrintLevel)
  {
    PrintLevel_ = PrintLevel;
    return(0);
  }
  //@}

protected:

  virtual int SetLabel() = 0;

  int ExtractSubmatrices();

  //! Containers_[i] contains all the necessary information to solve on each subblock.
  mutable vector<T*> Containers_;
  //! Contains information about non-overlapping partitions.
  Ifpack_Partitioner* Partitioner_;
  //! Label for \c this object
  string Label_;

private:

  bool IsSymmetric() const 
  {
    return(IsSymmetric_);
  }

  //! Pointers to the matrix to be preconditioned.
  Epetra_RowMatrix* Matrix_;
  // FIXME
  bool IsSymmetric_;
  //! Solver type.
  string SubmatricesType_;
  //! Number of preconditioning sweeps.
  int NumSweeps_;
  //! Damping parameter.
  double DampingFactor_;
  //! Number of local blocks
  int NumLocalBlocks_;
  //! Parameters list to be used to solve on each subblock
  Teuchos::ParameterList List_;
  //! If true, the preconditioner has been successfully computed.
  bool IsComputed_;
  //! Print level, from 0 (silent) to 10 (verbose)
  int PrintLevel_;
  
}; // class Ifpack_BlockPreconditioner

//==============================================================================
template<typename T>
Ifpack_BlockPreconditioner<T>::Ifpack_BlockPreconditioner(Epetra_RowMatrix* Matrix) :
  Matrix_(Matrix),
  IsSymmetric_(true),
  NumSweeps_(1),
  DampingFactor_(1.0),
  NumLocalBlocks_(1),
  Partitioner_(0),
  PrintLevel_(0)
{
}

//==============================================================================
template<typename T>
Ifpack_BlockPreconditioner<T>::~Ifpack_BlockPreconditioner()
{
  // FIXME: move in delete??
  for (int i = 0 ; i < NumLocalBlocks() ; ++i)
    delete Containers_[i];
}

//==============================================================================
template<typename T>
char* Ifpack_BlockPreconditioner<T>::Label() const
{
  return(const_cast<char*>(Label_.c_str()));
}

//==============================================================================
template<typename T>
int Ifpack_BlockPreconditioner<T>::
Apply(const Epetra_MultiVector& X, 
      Epetra_MultiVector& Y) const
{
  IFPACK_RETURN(Matrix_->Apply(X,Y));
}

//==============================================================================
template<typename T>
const Epetra_Comm& Ifpack_BlockPreconditioner<T>::
Comm() const
{
  return(Matrix_->Comm());
}

//==============================================================================
template<typename T>
const Epetra_Map& Ifpack_BlockPreconditioner<T>::
OperatorDomainMap() const
{
  return(Matrix_->OperatorDomainMap());
}

//==============================================================================
template<typename T>
const Epetra_Map& Ifpack_BlockPreconditioner<T>::
OperatorRangeMap() const
{
  return(Matrix_->OperatorRangeMap());
}

//==============================================================================
template<typename T>
int Ifpack_BlockPreconditioner<T>::
ExtractSubmatrices()
{

  if (Partitioner_ == 0)
    IFPACK_CHK_ERR(-1);

  NumLocalBlocks_ = Partitioner_->NumLocalParts();

  Containers_.resize(NumLocalBlocks());

  for (int i = 0 ; i < NumLocalBlocks() ; ++i) {

    Containers_[i] = new T;
    
    if (Containers_[i] == 0)
      IFPACK_CHK_ERR(-10);
    
    // set "global" ID of each partitioner row
    int rows = Partitioner_->NumRowsInPart(i);
    IFPACK_CHK_ERR(Containers_[i]->Shape(rows));

    for (int j = 0 ; j < rows ; ++j) {
      int LRID = (*Partitioner_)(i,j);
      Containers_[i]->ID(j) = LRID;
    }

    IFPACK_CHK_ERR(Containers_[i]->Extract(Matrix()));
    IFPACK_CHK_ERR(Containers_[i]->SetParameters(List_));
    IFPACK_CHK_ERR(Containers_[i]->Compute());

  }

  return(0);
}

//==============================================================================
template<typename T>
int Ifpack_BlockPreconditioner<T>::Compute()
{

  if (Matrix()->NumGlobalRows() != Matrix()->NumGlobalCols())
    IFPACK_CHK_ERR(-1); // only square matrices

  IFPACK_CHK_ERR(ExtractSubmatrices());
  
  IsComputed_ = true;

  return(0);

}

#endif // HAVE_IFPACK_TEUCHOS
#endif // IFPACK_BLOCKPRECONDITIONER_H
