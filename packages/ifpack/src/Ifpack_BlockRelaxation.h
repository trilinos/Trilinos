#ifndef IFPACK_BLOCKPRECONDITIONER_H
#define IFPACK_BLOCKPRECONDITIONER_H

#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_IFPACK_TEUCHOS
#include "Ifpack_Preconditioner.h" 
#include "Ifpack_Partitioner.h"
#include "Ifpack_LinearPartitioner.h"
#include "Ifpack_GreedyPartitioner.h"
#include "Ifpack_METISPartitioner.h"
#include "Ifpack_Graph_Epetra_RowMatrix.h"
#include "Ifpack_DenseContainer.h" 
#include "Teuchos_ParameterList.hpp"
#include "Epetra_RowMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Time.h"

static const int IFPACK_JACOBI = 0;
static const int IFPACK_GS = 1;
static const int IFPACK_SGS = 2;


//! Ifpack_BlockRelaxation: a class to define block preconditioners preconditioners of Epetra_RowMatrix's.

/*! The Ifpack_BlockRelaxation class enables the construction of
  block-preconditioners
  preconditioners of an Epetra_RowMatrix. The blocks can have arbitrary
  sizes. Ifpack_Partitioner is used to create the blocks (that is, to
  partition the graph of the local matrix). 

  \author Marzio Sala, SNL 9214.

 \date Last modified: Oct-04.
  
 */
template<typename T>
class Ifpack_BlockRelaxation : public Ifpack_Preconditioner {

public:

  //@{ \name Constructors/Destructors
  //! Ifpack_BlockRelaxation constructor with given Epetra_RowMatrix.
  /*! Creates an Ifpack_Preconditioner preconditioner. 
   *
   * \param In
   * Matrix - Pointer to matrix to be preconditioned.
   */
  Ifpack_BlockRelaxation(const Epetra_RowMatrix* Matrix);

  //! Copy constructor.
  Ifpack_BlockRelaxation(const Ifpack_BlockRelaxation& rhs);

#ifdef FIXME
  //! operator=
  Ifpack_BlockRelaxation & operator=(const Ifpack_BlockRelaxation& rhs);
#endif

  virtual ~Ifpack_BlockRelaxation();

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
			   Epetra_MultiVector& Y) const;

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

    string PT;
    if (PrecType_ == IFPACK_JACOBI)
      PT = "Jacobi";
    else if (PrecType_ == IFPACK_GS)
      PT = "Gauss-Seidel";
    else if (PrecType_ == IFPACK_SGS)
      PT = "symmetric Gauss-Seidel";

    PT = List.get("block: type", PT);

    if (PT == "Jacobi") {
      PrecType_ = IFPACK_JACOBI;
      KeepNonFactoredMatrix_ = false;
    }
    else if (PT == "Gauss-Seidel") {
      PrecType_ = IFPACK_GS;
      KeepNonFactoredMatrix_ = false;
    }
    else if (PT == "symmetric Gauss-Seidel") {
      KeepNonFactoredMatrix_ = true;
      PrecType_ = IFPACK_SGS;
    }

    SetNumSweeps(List.get("block: sweeps",NumSweeps()));
    SetDampingFactor(List.get("block: damping factor", DampingFactor()));
    PrintFrequency_ = List.get("block: print frequency", PrintFrequency());
    ZeroStartingSolution_ = List.get("block: zero starting solution", 
                                     ZeroStartingSolution_);
    PartitionerType_ = List.get("partitioner: type", PartitionerType_);
    NumLocalBlocks_ = List.get("partitioner: local parts", NumLocalBlocks_);
    // only Jacobi can work with overlap among local domains,
    OverlapLevel_ = List.get("partitioner: overlap", OverlapLevel_);
    if (PrecType_ != IFPACK_JACOBI)
      OverlapLevel_ = 0;

    // copy the list as each subblock's constructor will
    // require it later
    List_ = List;

    Label_ = "IFPACK block " + PT + ", # blocks = "
      + Ifpack_toString(NumLocalBlocks());
    return(0);
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

    if (PartitionerType_ == "linear")
      Partitioner_ = new Ifpack_LinearPartitioner(Graph_);
    else if (PartitionerType_ == "greedy")
      Partitioner_ = new Ifpack_GreedyPartitioner(Graph_);
    else if (PartitionerType_ == "metis")
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
                         const int MaxIters = 1550,
                         const double Tol = 1e-9,
			 Epetra_RowMatrix* Matrix = 0)
  {
    return(-1.0);
  }

  virtual double Condest() const
  {
    return(-1.0);
  }

  std::ostream& Print(std::ostream& os) const
  {
    os << "*** " << Label() << endl << endl;
    os << "Container label: " << Containers_[0]->Label() << endl;
    os << endl;
    os << "Number of rows         = " << Matrix().NumMyRows() << endl;
    os << "Number of sweeps       = " << NumSweeps_ << endl;
    os << "Damping Factor         = " << DampingFactor_ << endl;
    os << "Number of local blocks = " << Partitioner_->NumLocalParts() << endl;
    os << "Print frequency        = " << PrintFrequency_ << endl;
    os << "IsInitialized()        = " << IsInitialized_ << endl;
    os << "IsComputed()           = " << IsComputed_ << endl;
    os << endl;
    return(os);
  }

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

  virtual long int ComputeFlops() const
  {
    return(ComputeFlops_);
  }

  virtual long int ApplyInverseFlops() const
  {
    return(ApplyInverseFlops_);
  }

  virtual int OverlapLevel() const
  {
    return(OverlapLevel_);
  }

  virtual int PrecType() const
  {
    return(PrecType_);
  }

  virtual bool ZeroStartingSolution() const
  {
    return(ZeroStartingSolution_);
  }

  virtual string PartitionerType() const
  {
    return(PartitionerType_);
  }

protected:

  virtual int ApplyInverseJacobi(const Epetra_MultiVector& X, 
                                 Epetra_MultiVector& Y) const;

  virtual int ApplyInverseJacobi2(const Epetra_MultiVector& X, 
                                  Epetra_MultiVector& Y) const;

  virtual int ApplyInverseGS(const Epetra_MultiVector& X, 
                             Epetra_MultiVector& Y) const;

  virtual int ApplyInverseGS2(const Epetra_MultiVector& X, 
                              Epetra_MultiVector& Y) const;

  virtual int ApplyInverseSGS(const Epetra_MultiVector& X, 
                              Epetra_MultiVector& Y) const;

  virtual int ApplyInverseSGS2(Epetra_MultiVector& X, 
                             Epetra_MultiVector& Y) const;

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

  //! Pointers to the matrix to be preconditioned.
  const Epetra_RowMatrix* Matrix_;
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
  Epetra_Time* Time_;

  //! Contains the number of flops for Compute().
  long int ComputeFlops_;
  //! Contain sthe number of flops for ApplyInverse().
  long int ApplyInverseFlops_;
  int PrecType_;

}; // class Ifpack_BlockRelaxation

//==============================================================================
template<typename T>
Ifpack_BlockRelaxation<T>::
Ifpack_BlockRelaxation(const Epetra_RowMatrix* Matrix) :
  Matrix_(Matrix),
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
  OverlapLevel_(0),
  PrecType_(IFPACK_JACOBI),
  Time_(0)
{
}

//==============================================================================
template<typename T>
Ifpack_BlockRelaxation<T>::
Ifpack_BlockRelaxation(const Ifpack_BlockRelaxation& rhs) :
  Matrix_(&rhs.Matrix()),
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
  OverlapLevel_(0),
  PrecType_(IFPACK_JACOBI),
  Time_(0)
{
  if (rhs.IsInitialized())
    Initialize();

  if (rhs.IsComputed())
    Compute();
}

//==============================================================================
template<typename T>
Ifpack_BlockRelaxation<T>::~Ifpack_BlockRelaxation()
{
  for (int i = 0 ; i < NumLocalBlocks() ; ++i)
    delete Containers_[i];
  if (Partitioner_)
    delete Partitioner_;
  if (Graph_)
    delete Graph_;
  if (W_)
    delete W_;
  if (Time_)
    delete Time_;
}

//==============================================================================
template<typename T>
char* Ifpack_BlockRelaxation<T>::Label() const
{
  return(const_cast<char*>(Label_.c_str()));
}

//==============================================================================
template<typename T>
int Ifpack_BlockRelaxation<T>::
Apply(const Epetra_MultiVector& X, 
      Epetra_MultiVector& Y) const
{
  IFPACK_RETURN(Matrix().Apply(X,Y));
}

//==============================================================================
template<typename T>
const Epetra_Comm& Ifpack_BlockRelaxation<T>::
Comm() const
{
  return(Matrix().Comm());
}

//==============================================================================
template<typename T>
const Epetra_Map& Ifpack_BlockRelaxation<T>::
OperatorDomainMap() const
{
  return(Matrix().OperatorDomainMap());
}

//==============================================================================
template<typename T>
const Epetra_Map& Ifpack_BlockRelaxation<T>::
OperatorRangeMap() const
{
  return(Matrix().OperatorRangeMap());
}

//==============================================================================
template<typename T>
int Ifpack_BlockRelaxation<T>::
ExtractSubmatrices()
{

  if (Partitioner_ == 0)
    IFPACK_CHK_ERR(-1);

  NumLocalBlocks_ = Partitioner_->NumLocalParts();

  Containers_.resize(NumLocalBlocks());

  for (int i = 0 ; i < NumLocalBlocks() ; ++i) {

    int rows = Partitioner_->NumRowsInPart(i);
    Containers_[i] = new T(rows);
    
    Ifpack_DenseContainer* DC = 0;
    DC = dynamic_cast<Ifpack_DenseContainer*>(Containers_[i]);
    if (DC != 0)
      DC->SetKeepNonFactoredMatrix(KeepNonFactoredMatrix_);

    if (Containers_[i] == 0)
      IFPACK_CHK_ERR(-10);
    
    IFPACK_CHK_ERR(Containers_[i]->SetParameters(List_));
    IFPACK_CHK_ERR(Containers_[i]->Initialize());

    // set "global" ID of each partitioner row
    for (int j = 0 ; j < rows ; ++j) {
      int LRID = (*Partitioner_)(i,j);
      Containers_[i]->ID(j) = LRID;
    }

    IFPACK_CHK_ERR(Containers_[i]->Compute(*Matrix_));

  }

  return(0);
}

//==============================================================================
template<typename T>
int Ifpack_BlockRelaxation<T>::Compute()
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

//==============================================================================
template<typename T>
int Ifpack_BlockRelaxation<T>::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (IsComputed() == false)
    IFPACK_CHK_ERR(-4);

  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-3);

  switch (PrecType_) {
  case IFPACK_JACOBI:
    IFPACK_CHK_ERR(ApplyInverseJacobi(X,Y));
    break;
  case IFPACK_GS:
    IFPACK_CHK_ERR(ApplyInverseGS(X,Y));
    break;
  case IFPACK_SGS:
    IFPACK_CHK_ERR(ApplyInverseSGS(X,Y));
    break;
  }
  return(0);
}

//==============================================================================
template<typename T>
int Ifpack_BlockRelaxation<T>::
ApplyInverseJacobi(const Epetra_MultiVector& X, 
                   Epetra_MultiVector& Y) const
{

  // ------------ //
  // single sweep //
  // ------------ //

  if (NumSweeps_ == 1 && ZeroStartingSolution_
      && (PrintFrequency_ == 0)) {
    IFPACK_CHK_ERR(ApplyInverseJacobi2(X,Y));
    return(0);
  }

  // --------------------- //
  // general case (solver) //
  // --------------------- //

  if (PrintFrequency())
    Ifpack_PrintResidual(Label(),Matrix(),Y,X);

  // starting solution
  Epetra_MultiVector AX(Y);

  if (ZeroStartingSolution_)
    Y.PutScalar(0.0);

  for (int j = 0; j < NumSweeps() ; j++) {

    // compute the residual
    IFPACK_CHK_ERR(Apply(Y,AX));

    AX.Update(1.0,X,-1.0);

    // apply the block diagonal of A and update
    // the residual
    ApplyInverseJacobi2(AX,Y);

    if (PrintFrequency() && (j != 0) && (j % PrintFrequency() == 0))
      Ifpack_PrintResidual(j,Matrix(),Y,X);

  }

  if (PrintFrequency())
    Ifpack_PrintResidual(NumSweeps(),Matrix(),Y,X);

  return(0);
}

//==============================================================================
template<typename T>
int Ifpack_BlockRelaxation<T>::
ApplyInverseJacobi2(const Epetra_MultiVector& X, 
                    Epetra_MultiVector& Y) const
{
  // cycle over all local subdomains
  for (int i = 0 ; i < NumLocalBlocks() ; ++i) {

    int LID, GID;

    // extract RHS from X
    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);
      for (int k = 0 ; k < X.NumVectors() ; ++k) {
        Containers_[i]->RHS(j,k) = X[k][LID];
      }
    }

    // apply the inverse of each block
    IFPACK_CHK_ERR(Containers_[i]->ApplyInverse());

    // copy back into solution vector Y
    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);
      for (int k = 0 ; k < Y.NumVectors() ; ++k) {
        Y[k][LID] = Y[k][LID] + DampingFactor() * (*W_)[LID] * Containers_[i]->LHS(j,k);
      }
    }
  }
  return(0);
}

//==============================================================================
template<typename T>
int Ifpack_BlockRelaxation<T>::
ApplyInverseGS(const Epetra_MultiVector& X, 
               Epetra_MultiVector& Y) const
{

  // ------------ //
  // single sweep //
  // ------------ //

  if ((NumSweeps_ == 1) && ZeroStartingSolution_
      && (PrintFrequency_ == 0)) {
    IFPACK_CHK_ERR(ApplyInverseGS2(X,Y));
    return(0);
  }

  // --------------------- //
  // general case (solver) //
  // --------------------- //

  if (PrintFrequency())
    Ifpack_PrintResidual(Label(),Matrix(),Y,X);

  // starting solution
  Epetra_MultiVector AX(X);
  Epetra_MultiVector Ynew(Y);

  if (ZeroStartingSolution_)
    Y.PutScalar(0.0);

  for (int j = 0; j < NumSweeps() ; j++) {

    // compute the residual, I can skip first iteration
    // if starting solution is zero
    if (j || !ZeroStartingSolution_) {
      IFPACK_CHK_ERR(Apply(Y,AX));
      AX.Update(1.0,X,-1.0);
    }

    // apply the block diagonal of A and update the residual
    ApplyInverseGS2(AX,Ynew);

    Y.Update(DampingFactor(),Ynew,1.0);

    if (PrintFrequency() && (j != 0) && (j % PrintFrequency() == 0))
      Ifpack_PrintResidual(j,Matrix(),Y,X);

  }

  if (PrintFrequency())
    Ifpack_PrintResidual(NumSweeps(),Matrix(),Y,X);

  return(0);

}

//==============================================================================
template<typename T>
int Ifpack_BlockRelaxation<T>::
ApplyInverseGS2(const Epetra_MultiVector& X, 
                Epetra_MultiVector& Y) const
{

  // cycle over all local subdomains

  int Length = Matrix().MaxNumEntries();
  vector<int> Indices;
  vector<double> Values;
  Indices.resize(Length);
  Values.resize(Length);

  int NumMyRows = Matrix().NumMyRows();

// DELETE  Y.PutScalar(0.0);

  for (int i = 0 ; i < NumLocalBlocks() ; ++i) {

    int LID, GID;

    // update from previous block

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);

      int NumEntries;
      IFPACK_CHK_ERR(Matrix().ExtractMyRowCopy(LID, Length,NumEntries,
                                               &Values[0], &Indices[0]));

      for (int k = 0 ; k < NumEntries ; ++k) {
        int col = Indices[k];
        if (col >= NumMyRows) 
          continue;

        if ((*Partitioner_)(col) < i) {
          for (int kk = 0 ; kk < Y.NumVectors() ; ++kk) {
            X[kk][LID] -= Values[k] * Y[kk][col];
          }
        }
      }
    }

    // solve with this block

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);
      for (int k = 0 ; k < Y.NumVectors() ; ++k) {
        Containers_[i]->RHS(j,k) = X[k][LID];
      }
    }

    IFPACK_CHK_ERR(Containers_[i]->ApplyInverse());

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);
      for (int k = 0 ; k < Y.NumVectors() ; ++k) {
        Y[k][LID] += Containers_[i]->LHS(j,k);
      }
    }
  }

  return(0);
}

//==============================================================================
template<typename T>
int Ifpack_BlockRelaxation<T>::
ApplyInverseSGS(const Epetra_MultiVector& X, 
                Epetra_MultiVector& Y) const
{

  // ------------ //
  // single sweep //
  // ------------ //

  if (NumSweeps_ == 1 && ZeroStartingSolution_
      && (PrintFrequency_ == 0)) {
    Epetra_MultiVector Xtmp(X);
    Y.PutScalar(0.0);
    IFPACK_CHK_ERR(ApplyInverseSGS2(Xtmp,Y));
    return(0);
  }

  // --------------------- //
  // general case (solver) //
  // --------------------- //

  Epetra_MultiVector Xtmp(X);
  if (PrintFrequency())
    Ifpack_PrintResidual(Label(),Matrix(),Y,Xtmp);

  // starting solution
  Epetra_MultiVector AX(Xtmp);
  Epetra_MultiVector Ynew(Y);

  if (ZeroStartingSolution_)
    Y.PutScalar(0.0);

  for (int j = 0; j < NumSweeps() ; j++) {

    // compute the residual, I can skip first iteration
    // if starting solution is zero
    if (j || !ZeroStartingSolution_) {
      IFPACK_CHK_ERR(Apply(Y,AX));
      AX.Update(1.0,Xtmp,-1.0);
    }

    // apply the block diagonal of A and update the residual
    ApplyInverseSGS2(AX,Ynew);

    Y.Update(DampingFactor(),Ynew,1.0);

    if (PrintFrequency() && (j != 0) && (j % PrintFrequency() == 0))
      Ifpack_PrintResidual(j,Matrix(),Y,Xtmp);

  }

  if (PrintFrequency())
    Ifpack_PrintResidual(NumSweeps(),Matrix(),Y,Xtmp);

  return(0);

}

//==============================================================================
template<typename T>
int Ifpack_BlockRelaxation<T>::
ApplyInverseSGS2(Epetra_MultiVector& X, 
                 Epetra_MultiVector& Y) const
{

  // Y : in input, previous solution
  //     in output, update solution
  // X : tmp vector containing in input the rhs. It will be
  //     overwritten by temporary data

  // cycle over all local subdomains
  int Length = Matrix().MaxNumEntries();
  vector<int> Indices;
  vector<double> Values;
  Indices.resize(Length);
  Values.resize(Length);

  int NumMyRows = Matrix().NumMyRows();

  // ================== //
  // apply (D - E)^{-1} //
  // ================== //

  for (int i = 0 ; i < NumLocalBlocks() ; ++i) {

    int LID, GID;

    // update from previous block

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);

      int NumEntries;
      IFPACK_CHK_ERR(Matrix().ExtractMyRowCopy(LID, Length,NumEntries,
                                               &Values[0], &Indices[0]));

      for (int k = 0 ; k < NumEntries ; ++k) {
        int col = Indices[k];
        if (col >= NumMyRows ) 
          continue;

        if ((*Partitioner_)(col) < i) {
          for (int kk = 0 ; kk < Y.NumVectors() ; ++kk) {
            X[kk][LID] -= Values[k] * Y[kk][col];
          }
        }
      }
    }

    // solve with this block

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);
      for (int k = 0 ; k < Y.NumVectors() ; ++k) {
        Containers_[i]->RHS(j,k) = X[k][LID];
      }
    }

    IFPACK_CHK_ERR(Containers_[i]->ApplyInverse());

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);
      for (int k = 0 ; k < Y.NumVectors() ; ++k) {
        Y[k][LID] = Containers_[i]->LHS(j,k);
      }
    }
  }

  // ============ //
  // apply D to Y //
  // ============ //

  for (int i = 0 ; i < NumLocalBlocks() ; ++i) {

    int LID, GID;

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);
      for (int k = 0 ; k < Y.NumVectors() ; ++k) {
        Containers_[i]->RHS(j,k) = Y[k][LID];
      }
    }

    IFPACK_CHK_ERR(Containers_[i]->Apply());

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);
      for (int k = 0 ; k < Y.NumVectors() ; ++k) {
        Y[k][LID] = Containers_[i]->LHS(j,k);
      }
    }
  }

  X = Y;

  // ================== //
  // apply (D - F)^{-1} //
  // ================== //

  for (int i = NumLocalBlocks() - 1; i >=0 ; --i) {

    int LID, GID;

    // update from previous block

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);

      int NumEntries;
      IFPACK_CHK_ERR(Matrix().ExtractMyRowCopy(LID, Length,NumEntries,
                                               &Values[0], &Indices[0]));

      for (int k = 0 ; k < NumEntries ; ++k) {
        int col = Indices[k];
        if (col >= NumMyRows ) 
          continue;

        if ((*Partitioner_)(col) > i) {
          for (int kk = 0 ; kk < Y.NumVectors() ; ++kk) {
            X[kk][LID] -= Values[k] * Y[kk][col];
          }
        }
      }
    }

    // solve with this block

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);
      for (int k = 0 ; k < Y.NumVectors() ; ++k) {
        Containers_[i]->RHS(j,k) = X[k][LID];
      }
    }

    IFPACK_CHK_ERR(Containers_[i]->ApplyInverse());

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);
      for (int k = 0 ; k < Y.NumVectors() ; ++k) {
        Y[k][LID] = Containers_[i]->LHS(j,k);
      }
    }
  }
  return(0);
}

//==============================================================================
#endif // HAVE_IFPACK_TEUCHOS
#endif // IFPACK_BLOCKPRECONDITIONER_H
