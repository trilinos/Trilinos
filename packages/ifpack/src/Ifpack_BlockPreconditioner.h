#ifndef IFPACK_BLOCKPRECONDITIONER_H
#define IFPACK_BLOCKPRECONDITIONER_H

#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_IFPACK_TEUCHOS
#include "Epetra_RowMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Ifpack_Preconditioner.h" 
#include "Ifpack_Partitioner.h"
#include "Ifpack_LinearPartitioner.h"
#include "Ifpack_GreedyPartitioner.h"
#include "Ifpack_METISPartitioner.h"
#include "Ifpack_Graph_Epetra_RowMatrix.h"
#include "Ifpack_DenseContainer.h" 
#include "Teuchos_ParameterList.hpp"

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
  virtual int SetParameters(Teuchos::ParameterList& List)
  {
    
    SetNumSweeps(List.get("block: sweeps",NumSweeps()));
    SetDampingFactor(List.get("block: damping factor", DampingFactor()));
    PrintFrequency_ = List.get("block: print frequency", PrintFrequency());
//    Partitioner_ = List.get("partitioner: object", (Ifpack_Partitioner*)0);
    ZeroStartingSolution_ = List.get("block: zero starting solution", 
				     ZeroStartingSolution_);
    PartitionerType_ = List.get("partitioner: type", PartitionerType_);
    NumLocalBlocks_ = List.get("partitioner: local parts", NumLocalBlocks_);
    // only Jacobi can work with overlap among local domains,
    // I wanna check this out later
    OverlapLevel_ = List.get("partitioner: overlap", OverlapLevel_);

    // copy the list as each subblock's constructor will
    // require it later
    List_ = List;

    // derived class will set the appropriate label
    IFPACK_CHK_ERR(SetLabel());
    return(0);
  }

  //! Not implemented.
  virtual int SetParameter(const string Name, const int value)
  {
    return(-1);
  }

  //! Not implemented.
  virtual int SetParameter(const string Name, const double value)
  {
    return(-1);
  }

  //! Initializes the preconditioner.
  virtual int Initialize()
  {
    IsInitialized_ = false;

    if (Partitioner_)
      delete Partitioner_;
    if (Graph_)
      delete Graph_;

    Graph_ = new Ifpack_Graph_Epetra_RowMatrix(&Matrix());
    assert (Graph_ != 0);

    if (PartitionerType_ == "Linear")
      Partitioner_ = new Ifpack_LinearPartitioner(Graph_);
    else if (PartitionerType_ == "Greedy")
      Partitioner_ = new Ifpack_GreedyPartitioner(Graph_);
    else if (PartitionerType_ == "METIS")
      Partitioner_ = new Ifpack_METISPartitioner(Graph_);
    else
      IFPACK_CHK_ERR(-1);

    assert (Partitioner_ != 0);

    // need to partition the graph of A
    IFPACK_CHK_ERR(Partitioner_->SetParameters(List_));
    IFPACK_CHK_ERR(Partitioner_->Compute());

    // get actual number of partitions
    NumLocalBlocks_ = Partitioner_->NumLocalParts();
    
    // weight of each vertex
    if (W_)
      delete W_;
    W_ = new Epetra_Vector(Matrix().RowMatrixRowMap());
    W_->PutScalar(0.0);

    for (int i = 0 ; i < NumLocalBlocks() ; ++i) {

      for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
	int LID = (*Partitioner_)(i,j);
        (*W_)[LID]++;
      }
    }
    W_->Reciprocal(*W_);

    IsInitialized_ = true;

    return(0);
  }

  //! Computes the preconditioner.
  virtual int Compute();

  //! Sets the number of sweeps.
  inline int SetNumSweeps(const int NumSweeps)
  {
    NumSweeps_ = NumSweeps;
    return(0);
  }

  //! Gets the number of sweeps.
  inline int NumSweeps() const
  {
    return(NumSweeps_);
  }
 
  //! Sets the damping factor.
  inline int SetDampingFactor(const double DampingFactor)
  {
    DampingFactor_ = DampingFactor;
    return(0);
  }

  //! Gets the damping parameter.
  inline double DampingFactor() const
  {
    return(DampingFactor_);
  }

  virtual const Epetra_RowMatrix& Matrix() const
  {
    return(*Matrix_);
  }

  virtual Epetra_RowMatrix& Matrix()
  {
    return(*Matrix_);
  }

  inline int SetPrintFrequency(const int PrintFrequency)
  {
    PrintFrequency_ = PrintFrequency;
    return(0);
  }

  inline int PrintFrequency() const
  {
    return(PrintFrequency_);
  }

  virtual double Condest(const Ifpack_CondestType CT = Ifpack_Cheap,
			 Epetra_RowMatrix* Matrix = 0)
  {
    return(-1.0);
  }

  std::ostream& Print(std::ostream& os) const
  {
    if (Matrix().Comm().MyPID())
      return(os);
    os << "*** " << Label() << endl;
    os << "*** Container label: " << Containers_[0]->Label() << endl;
    return(os);
  }

  //@}

protected:

  virtual int SetLabel() = 0;

  int ExtractSubmatrices();

  //! Containers_[i] contains all the necessary information to solve on each subblock.
  mutable vector<T*> Containers_;
  //! Contains information about non-overlapping partitions.
  Ifpack_Partitioner* Partitioner_;
  string PartitionerType_;
  //! Label for \c this object
  string Label_;
  //! If \c true, starting solution is the zero vector.
  bool ZeroStartingSolution_;
  Ifpack_Graph* Graph_;
  Epetra_Vector* W_;

  bool KeepNonFactoredMatrix_;
  int OverlapLevel_;

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
  //! If true, the preconditioner has been successfully initialized.
  bool IsInitialized_;
  //! If true, the preconditioner has been successfully computed.
  bool IsComputed_;
  //! Print level, from 0 (silent) to 10 (verbose)
  int PrintFrequency_;
  
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
  Graph_(0),
  W_(0),
  PrintFrequency_(0),
  ZeroStartingSolution_(true),
  PartitionerType_("Greedy"),
  KeepNonFactoredMatrix_(false),
  OverlapLevel_(0)
{
}

//==============================================================================
template<typename T>
Ifpack_BlockPreconditioner<T>::~Ifpack_BlockPreconditioner()
{
  // FIXME: move in delete??
  for (int i = 0 ; i < NumLocalBlocks() ; ++i)
    delete Containers_[i];
  if (Partitioner_)
    delete Partitioner_;
  if (Graph_)
    delete Graph_;
  if (W_)
    delete W_;
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
  IFPACK_RETURN(Matrix().Apply(X,Y));
}

//==============================================================================
template<typename T>
const Epetra_Comm& Ifpack_BlockPreconditioner<T>::
Comm() const
{
  return(Matrix().Comm());
}

//==============================================================================
template<typename T>
const Epetra_Map& Ifpack_BlockPreconditioner<T>::
OperatorDomainMap() const
{
  return(Matrix().OperatorDomainMap());
}

//==============================================================================
template<typename T>
const Epetra_Map& Ifpack_BlockPreconditioner<T>::
OperatorRangeMap() const
{
  return(Matrix().OperatorRangeMap());
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
    
    Ifpack_DenseContainer* DC = 0;
    DC = dynamic_cast<Ifpack_DenseContainer*>(Containers_[i]);
    if (DC != 0)
      DC->KeepNonFactoredMatrix(KeepNonFactoredMatrix_);

    if (Containers_[i] == 0)
      IFPACK_CHK_ERR(-10);
    
    // set "global" ID of each partitioner row
    int rows = Partitioner_->NumRowsInPart(i);
    IFPACK_CHK_ERR(Containers_[i]->Shape(rows));

    for (int j = 0 ; j < rows ; ++j) {
      int LRID = (*Partitioner_)(i,j);
      Containers_[i]->ID(j) = LRID;
    }

    IFPACK_CHK_ERR(Containers_[i]->Extract(&Matrix()));
    IFPACK_CHK_ERR(Containers_[i]->SetParameters(List_));
    IFPACK_CHK_ERR(Containers_[i]->Compute());

  }

  return(0);
}

//==============================================================================
template<typename T>
int Ifpack_BlockPreconditioner<T>::Compute()
{

  if (!IsInitialized())
    IFPACK_CHK_ERR(Initialize());
  IsComputed_ = false;

  if (Matrix().NumGlobalRows() != Matrix().NumGlobalCols())
    IFPACK_CHK_ERR(-1); // only square matrices

  IFPACK_CHK_ERR(ExtractSubmatrices());
  
  IsComputed_ = true;

  return(0);

}

#endif // HAVE_IFPACK_TEUCHOS
#endif // IFPACK_BLOCKPRECONDITIONER_H
