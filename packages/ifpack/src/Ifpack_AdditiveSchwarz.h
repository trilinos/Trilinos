#ifndef IFPACK_ADDITIVESCHWARZ_H
#define IFPACK_ADDITIVESCHWARZ_H

#if defined(HAVE_IFPACK_TEUCHOS)

#include "Ifpack_Preconditioner.h"
#include "Teuchos_ParameterList.hpp"
#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_Reordering.h"
#include "Ifpack_RCMReordering.h"
#include "Ifpack_METISReordering.h"
#include "Ifpack_LocalFilter.h"
#include "Ifpack_ReorderFilter.h"
#include "Ifpack_DropFilter.h"
#include "Ifpack_SparsityFilter.h"
#include "Ifpack_SingletonFilter.h"
#include "Ifpack_Utils.h"
#include "Ifpack_OverlappingRowMatrix.h"
#include "Epetra_CombineMode.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Epetra_Time.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_ParameterList.hpp"

//! Ifpack_AdditiveSchwarz: a class to define Additive Schwarz preconditioners.

/*!
  Class  Ifpack_AdditiveSchwarz enables the construction of Additive
  Schwarz (one-level overlapping domain decomposition) preconditioners, 
  for a given Epetra_RowMatrix.
  Ifpack_AdditiveSchwarz is derived from Ifpack_Preconditioner,
  itself derived from Epetra_Operator. An application of
  the Additive Schwarz preconditioner can be obtained 
  by calling method ApplyInverse().

  One-level overlapping domain decomposition preconditioners use 
  local solvers, of Dirichlet type. This means that the inverse of
  the local matrix (with minimal or wider overlap) is applied to
  the residual to be preconditioned.

  The preconditioner can be written as:
  \f[
  P_{AS}^{-1} = \sum_{i=1}^M R_i^T A_i^{-1} R_i ,
  \f]
  where \f$M\f$ is the number of subdomains (that is, the number of 
  processors in
  the computation), \f$R_i\f$ is an operator that restricts the global
  vector to the vector lying on subdomain \f$i\f$, and
  \f[
  A_i = R_i A R_i^T.
  \f]

  The construction of Schwarz preconditioners is mainly composed by
  two steps:
  - definition of the restriction and prolongation operator
  \f$R_i\f$ and \f$R_i^T\f$. If minimal overlap is chosen, their
  implementation is trivial, \f$R_i\f$ will return all the local
  components. For wider overlaps, instead, Epetra_Import and
  Epetra_Export will be used to import/export data. The user
  must provide both the matrix to be preconditioned (which is suppose
  to have minimal-overlap) and the matrix with wider overlap.
  - definition of a technique to apply the inverse of \f$A_i\f$.
  To solve on each subdomain, the user can adopt any class, derived
  from Ifpack_Preconditioner. This can be easily accomplished, as
  Ifpack_AdditiveSchwarz is templated with the solver for each subdomain.

*/

template<typename T>
class Ifpack_AdditiveSchwarz : public virtual Ifpack_Preconditioner {
      
public:

  //@{ \name Constructors/Destructors
  //! Ifpack_AdditiveSchwarz constructor with given Epetra_RowMatrix.
  /*! Creates an Ifpack_AdditiveSchwarz preconditioner with overlap.
   * To use minimal-overlap, OverlappingMatrix is omitted
   * (as defaulted to 0).
   *
   * \param
   * Matrix - (In) Pointer to matrix to be preconditioned
   *
   * \param
   * OverlappingMatrix - (In) Pointer to the matrix extended with the
   *                     desired level of overlap.
   */
  Ifpack_AdditiveSchwarz(Epetra_RowMatrix* Matrix,
			 int OverlapLevel = 0);
  
  Ifpack_AdditiveSchwarz(const Ifpack_AdditiveSchwarz& RHS);

  //! Destructor
  virtual ~Ifpack_AdditiveSchwarz();
  //@}

  //@{ \name Atribute set methods.

    //! If set true, transpose of this operator will be applied (not implemented).
    /*! This flag allows the transpose of the given operator to be used 
     * implicitly.  
      
    \param 
	   UseTranspose - (In) If true, multiply by the transpose of operator, 
	   otherwise just use operator.

    \return Integer error code, set to 0 if successful.  Set to -1 if this implementation does not support transpose.
  */
    virtual int SetUseTranspose(bool UseTranspose);
  //@}
  
  //@{ \name Mathematical functions.

  //! Applies the matrix to an Epetra_MultiVector.
  /*! 
    \param
    X - (In) A Epetra_MultiVector of dimension NumVectors 
       to multiply with matrix.
    \param
    Y -(Out) A Epetra_MultiVector of dimension NumVectors 
       containing the result.

    \return Integer error code, set to 0 if successful.
    */
    virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //! Applies the preconditioner to X, returns the result in Y.
  /*! 
    \param
    X - (In) A Epetra_MultiVector of dimension NumVectors to be preconditioned.
    \param
    Y -(Out) A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning In order to work with AztecOO, any implementation of this method 
    must support the case where X and Y are the same object.
    */
    virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //! Returns the infinity norm of the global matrix (not implemented)
    virtual double NormInf() const;
  //@}
  
  //@{ \name Atribute access functions

    //! Returns a character string describing the operator
    virtual char * Label() const;

    //! Returns the current UseTranspose setting.
    virtual bool UseTranspose() const;

    //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
    virtual bool HasNormInf() const;

    //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
    virtual const Epetra_Comm & Comm() const;

    //! Returns the Epetra_Map object associated with the domain of this operator.
    virtual const Epetra_Map & OperatorDomainMap() const;

    //! Returns the Epetra_Map object associated with the range of this operator.
    virtual const Epetra_Map & OperatorRangeMap() const;
  //@}

  //! Returns \c true if the preconditioner has been successfully initialized.
  virtual bool IsInitialized() const
  {
    return(IsInitialized_);
  }

  //! Returns \c true if the preconditioner has been successfully computed.
  virtual bool IsComputed() const
  {
    return(IsComputed_);
  }

  //! Sets the parameters.
  /*! Sets the parameter for the additive Schwarz preconditioner,
   * as well as for all the preconditioners that may need to
   * be defined on each subblock.
   * Parameters accepted by List are:
   * - \c "schwarz: combine mode" : It must be an Epetra_CombineMode.
   *     Default: \c Zero.
   *     It Can be assume of the following values:
   *   - Add: Components on the receiving processor will be added together;
   *   - Zero: Off-processor components will be ignored;
   *   - Insert: Off-processor components will be inserted into locations on 
   *     receiving processor replacing existing values.
   *   - Average: Off-processor components will be averaged with existing;
   *   - AbsMax: Magnitudes of Off-processor components will be 
   *     maxed with magnitudes of existing components on the receiving 
   *     processor.
   * - \c "schwarz: compute condest" : if \c true, \c Compute() will
   *    estimate the condition number of the preconditioner. 
   *    Default: \c false.
   */
  virtual int SetParameters(Teuchos::ParameterList& List);

  //! Initialized the preconditioner.
  virtual int Initialize();

  //! Computes the preconditioner.
  virtual int Compute();

  //! Returns the estimated condition number (if computed).
  virtual double Condest(const Ifpack_CondestType CT = Ifpack_Cheap,
                         const int MaxIters = 1550,
                         const double Tol = 1e-9,
			 Epetra_RowMatrix* Matrix = 0);

  //! Returns the estimated condition number, never computes it.
  virtual double Condest() const
  {
    return(Condest_);
  }

  //! Returns a refernence to the internally stored matrix.
  virtual const Epetra_RowMatrix& Matrix() const
  {
    return(*Matrix_);
  }

  //! Returns \c true is an overlapping matrix is present.
  virtual bool IsOverlapping() const
  {
    return(IsOverlapping_);
  }

  //! Prints major information about this preconditioner.
  virtual std::ostream& Print(std::ostream&) const;
  
  virtual const T* Inverse() const
  {
    return(Inverse_);
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

  virtual const Teuchos::ParameterList& List() const
  {
    return(List_);
  }

protected:

  //! Sets up the localized matrix and the singleton filter.
  int Setup();
  
  //! Destroys all allocated data
  void Destroy();

  //! Pointers to the matrix to be preconditioned.
  const Epetra_RowMatrix* Matrix_;
  //! Pointers to the overlapping matrix.
  Ifpack_OverlappingRowMatrix* OverlappingMatrix_;
  //! Localized version of Matrix_ or OverlappingMatrix_.
  Ifpack_LocalFilter* LocalizedMatrix_;
  //! Contains the label of \c this object.
  string Label_;
  //! If true, the preconditioner has been successfully initialized.
  bool IsInitialized_;
  //! If true, the preconditioner has been successfully computed.
  bool IsComputed_;
  //! Pointer to the local solver.
  T* Inverse_;
  //! If true, overlapping is used
  bool IsOverlapping_;
  int OverlapLevel_;
  //! Stores a copy of the list given in SetParameters()
  Teuchos::ParameterList List_;
  //! Combine mode for off-process elements (only if overlap is used)
  Epetra_CombineMode CombineMode_;
  //! Contains the estimated condition number.
  double Condest_;
  
  bool UseReordering_;
  string ReorderingType_;
  Ifpack_Reordering* Reordering_;
  Ifpack_ReorderFilter* ReorderedLocalizedMatrix_;
  bool UseFilter_;
  bool FilterSingletons_;
  Epetra_RowMatrix* FilteredMatrix_;
  Ifpack_SingletonFilter* SingletonFilter_;
 
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

  //! Contains the number of flops for Compute().
  long int ComputeFlops_;
  //! Contain sthe number of flops for ApplyInverse().
  mutable long int ApplyInverseFlops_;

  Epetra_Time* Time_;
};

//==============================================================================
template<typename T>
Ifpack_AdditiveSchwarz<T>::
Ifpack_AdditiveSchwarz(Epetra_RowMatrix* Matrix,
		       int OverlapLevel) :
  Matrix_(Matrix),
  OverlappingMatrix_(0),
  LocalizedMatrix_(0),
  IsInitialized_(false),
  IsComputed_(false),
  Inverse_(0),
  OverlapLevel_(OverlapLevel),
  IsOverlapping_(false),
  CombineMode_(Zero),
  Condest_(-1.0),
  UseReordering_(false),
  ReorderedLocalizedMatrix_(0),
  Reordering_(0),
  UseFilter_(false),
  FilterSingletons_(false),
  FilteredMatrix_(0),
  SingletonFilter_(0),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyInverseTime_(0.0),
  ComputeFlops_(0),
  ApplyInverseFlops_(0),
  Time_(0)
{
  if (Matrix_->Comm().NumProc() == 1)
    OverlapLevel_ = 0;

  if ((OverlapLevel_ != 0) && (Matrix_->Comm().NumProc() > 1))
    IsOverlapping_ = true;
  // Sets parameters to default values
  Teuchos::ParameterList List;
  SetParameters(List);
}

//==============================================================================
template<typename T>
Ifpack_AdditiveSchwarz<T>::
Ifpack_AdditiveSchwarz(const Ifpack_AdditiveSchwarz& RHS) :
  Matrix_(&RHS.Matrix()),
  OverlappingMatrix_(0),
  LocalizedMatrix_(0),
  IsInitialized_(false),
  IsComputed_(false),
  Inverse_(0),
  OverlapLevel_(RHS.OverlapLevel()),
  IsOverlapping_(RHS.IsOverlapping()),
  CombineMode_(Zero),
  Condest_(-1.0),
  UseReordering_(false),
  ReorderedLocalizedMatrix_(0),
  Reordering_(0),
  UseFilter_(false),
  FilterSingletons_(false),
  FilteredMatrix_(0),
  SingletonFilter_(0),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyInverseTime_(0.0),
  ComputeFlops_(0),
  ApplyInverseFlops_(0),
  Time_(0)
{
  // Sets parameters to default values
  Teuchos::ParameterList List(RHS.List());
  SetParameters(List);

  // as Inverse_ is a generic (without copy constructor) 
  // Ifpack_Preconditioner, I simply recall Initialize() and
  // Construct() is RHS has been initialized and/or constructed
  if (RHS.IsInitialized())
    Initialize();

  if (RHS.IsComputed())
    Compute();
}

//==============================================================================
template<typename T>
Ifpack_AdditiveSchwarz<T>::~Ifpack_AdditiveSchwarz()
{
  Destroy();
}

//==============================================================================
template<typename T>
void Ifpack_AdditiveSchwarz<T>::Destroy() 
{
  if (OverlappingMatrix_)
    delete OverlappingMatrix_;
  OverlappingMatrix_ = 0;

  if (Inverse_)
    delete Inverse_;
  Inverse_ = 0;

  if (LocalizedMatrix_)
    delete LocalizedMatrix_;
  LocalizedMatrix_ = 0;

  if (ReorderedLocalizedMatrix_)
    delete ReorderedLocalizedMatrix_;
  ReorderedLocalizedMatrix_ = 0;

  if (FilteredMatrix_)
    delete FilteredMatrix_;
  FilteredMatrix_ = 0;

  if (SingletonFilter_)
    delete SingletonFilter_;
  SingletonFilter_ = 0;

  if (Reordering_)
    delete Reordering_;
  Reordering_ = 0;

  if (Time_)
    delete Time_;
  Time_ = 0;
}

//==============================================================================
template<typename T>
int Ifpack_AdditiveSchwarz<T>::Setup()
{

  double AddToDiag = List_.get("filter: add to diagonal", 0.0);
  Epetra_RowMatrix* MatrixPtr;

  if (OverlappingMatrix_)
    LocalizedMatrix_ = new Ifpack_LocalFilter(OverlappingMatrix_,
					      AddToDiag);
  else
    LocalizedMatrix_ = new Ifpack_LocalFilter(Matrix_,
					      AddToDiag);

  if (LocalizedMatrix_ == 0)
    IFPACK_CHK_ERR(-1);

  // users may want to skip singleton check
  if (FilterSingletons_) {
    SingletonFilter_ = new Ifpack_SingletonFilter(LocalizedMatrix_);
    MatrixPtr = SingletonFilter_;
  }
  else
    MatrixPtr = LocalizedMatrix_;

  if (UseReordering_) {

    // create reordeing and compute it
    if (ReorderingType_ == "rcm")
      Reordering_ = new Ifpack_RCMReordering();
    else if (ReorderingType_ == "metis")
      Reordering_ = new Ifpack_METISReordering();
    else {
      cerr << "reordering type not correct (" << ReorderingType_ << ")" << endl;
      exit(EXIT_FAILURE);
    }
    IFPACK_CHK_ERR(Reordering_->SetParameters(List_));
    IFPACK_CHK_ERR(Reordering_->Compute(*MatrixPtr));

    // now create reordered localized matrix
    ReorderedLocalizedMatrix_ = 
      new Ifpack_ReorderFilter(MatrixPtr,Reordering_);
    assert(ReorderedLocalizedMatrix_ != 0);
    MatrixPtr = ReorderedLocalizedMatrix_;
  }

  if (UseFilter_) {
    string DropScheme = List_.get("filter: type", "by-value");

    if (DropScheme == "by-value") {
      double DropValue = List_.get("filter: drop value", 1e-9);
      FilteredMatrix_ = new Ifpack_DropFilter(MatrixPtr,DropValue);
    }
    else if (DropScheme == "by-sparsity") {
      int AllowedEntries = List_.get("filter: allowed entries", 1);
      int AllowedBandwidth = List_.get("filter: allowed bandwidth", 
				   MatrixPtr->NumMyRows());
      FilteredMatrix_ = new Ifpack_SparsityFilter(MatrixPtr,AllowedEntries,
						  AllowedBandwidth);
    }
    else {
      cerr << "Option `filter: type' not recognized ("
	   << DropScheme << ")." << endl;
      exit(EXIT_FAILURE);
    }
      
    assert (FilteredMatrix_ != 0);
    MatrixPtr = FilteredMatrix_;
  }

  Inverse_ = new T(MatrixPtr);

  if (Inverse_ == 0)
    IFPACK_CHK_ERR(-1);

  return(0);
}

//==============================================================================
template<typename T>
int Ifpack_AdditiveSchwarz<T>::SetParameters(Teuchos::ParameterList& List)
{
 
  CombineMode_ = List.get("schwarz: combine mode", CombineMode_);
  UseReordering_ = List.get("schwarz: use reordering", UseReordering_);
  ReorderingType_ = List.get("schwarz: reordering type", "rcm");
  UseFilter_ = List.get("schwarz: use filter", UseFilter_);
  FilterSingletons_ = List.get("schwarz: filter singletons", FilterSingletons_);

  // This copy may be needed by Amesos or other preconditioners.
  List_ = List;

  return(0);
}

//==============================================================================
template<typename T>
int Ifpack_AdditiveSchwarz<T>::Initialize()
{

  IsInitialized_ = false;
  IsComputed_ = false; // values required
  Condest_ = -1.0; // zero-out condest

  Destroy();

  if (Time_ == 0)
    Time_ = new Epetra_Time(Comm());

  Time_->ResetStartTime();

  // compute the overlapping matrix if necessary
  if (IsOverlapping_) {
    OverlappingMatrix_ = 
      new Ifpack_OverlappingRowMatrix(Matrix_, OverlapLevel_);
    if (OverlappingMatrix_ == 0)
      IFPACK_CHK_ERR(-1);
  }

  IFPACK_CHK_ERR(Setup());

  if (Inverse_ == 0)
    IFPACK_CHK_ERR(-1);

  if (LocalizedMatrix_ == 0)
    IFPACK_CHK_ERR(-1);

  IFPACK_CHK_ERR(Inverse_->SetParameters(List_));
  IFPACK_CHK_ERR(Inverse_->Initialize());

  // Label is for Aztec-like solvers
  Label_ = "Ifpack_AdditiveSchwarz, ov = " + Ifpack_toString(OverlapLevel_)
    + ", local solver = \n\t\t***** `" + string(Inverse_->Label()) + "'";

  IsInitialized_ = true;
  ++NumInitialize_;
  InitializeTime_ += Time_->ElapsedTime();

  return(0);
}
//==============================================================================
template<typename T>
int Ifpack_AdditiveSchwarz<T>::Compute()
{

  if (IsInitialized() == false)
    IFPACK_CHK_ERR(Initialize());

  Time_->ResetStartTime();
  IsComputed_ = false;
  Condest_ = -1.0;
  
  IFPACK_CHK_ERR(Inverse_->Compute());

  string R = "";
  if (UseReordering_)
    R = ReorderingType_ + " reord, ";

  // A call to Condest() will update the condition number estimate
  Condest();
  
  // reset lavel with condest()
  Label_ = "Ifpack_AdditiveSchwarz, ov = " + Ifpack_toString(OverlapLevel_)
    + ", local solver = \n\t\t***** `" + string(Inverse_->Label()) + "'"
    + "\n\t\t***** " + R + "Condition number estimate = "
    + Ifpack_toString(Condest());

  IsComputed_ = true;
  ++NumCompute_;
  ComputeTime_ += Time_->ElapsedTime();

  return(0);
}

//==============================================================================
template<typename T>
int Ifpack_AdditiveSchwarz<T>::SetUseTranspose(bool UseTranspose)
{
  IFPACK_CHK_ERR(-99); // not implemented
}

//==============================================================================
template<typename T>
int Ifpack_AdditiveSchwarz<T>::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_CHK_ERR(Matrix_->Apply(X,Y));
  return(0);
}

//==============================================================================
template<typename T>
double Ifpack_AdditiveSchwarz<T>::NormInf() const
{
  return(-1.0);
}

//==============================================================================
template<typename T>
char * Ifpack_AdditiveSchwarz<T>::Label() const
{
  return((char*)Label_.c_str());
}

//==============================================================================
template<typename T>
bool Ifpack_AdditiveSchwarz<T>::UseTranspose() const
{
  return(false);
}

//==============================================================================
template<typename T>
bool Ifpack_AdditiveSchwarz<T>::HasNormInf() const
{
  return(false);
}

//==============================================================================
template<typename T>
const Epetra_Comm & Ifpack_AdditiveSchwarz<T>::Comm() const
{
  return(Matrix_->Comm());
}

//==============================================================================
template<typename T>
const Epetra_Map & Ifpack_AdditiveSchwarz<T>::OperatorDomainMap() const
{
  return(Matrix_->OperatorDomainMap());
}

//==============================================================================
template<typename T>
const Epetra_Map & Ifpack_AdditiveSchwarz<T>::OperatorRangeMap() const
{
  return(Matrix_->OperatorRangeMap());
}

//==============================================================================
template<typename T>
int Ifpack_AdditiveSchwarz<T>::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int NumVectors = X.NumVectors();

  if (NumVectors != Y.NumVectors())
    IFPACK_CHK_ERR(-1); // wrong input

  // FIXME
  if (!X.Map().SameAs(Matrix().RowMatrixRowMap()))
    IFPACK_CHK_ERR(-99);
  if (!Y.Map().SameAs(Matrix().RowMatrixRowMap()))
    IFPACK_CHK_ERR(-99);

  Time_->ResetStartTime();

  Epetra_MultiVector* OverlappingX;
  Epetra_MultiVector* OverlappingY;

  // process overlap, may need to create vectors and import data
  if (IsOverlapping()) {
    OverlappingX = new Epetra_MultiVector(OverlappingMatrix_->RowMatrixRowMap(),
					  X.NumVectors());
    OverlappingY = new Epetra_MultiVector(OverlappingMatrix_->RowMatrixRowMap(),
                                          Y.NumVectors());
    assert (OverlappingY != 0);

    OverlappingY->PutScalar(0.0);
    OverlappingX->PutScalar(0.0);
    IFPACK_CHK_ERR(OverlappingMatrix_->ImportMultiVector(X,*OverlappingX,Insert));
  }
  else {
    OverlappingX = (Epetra_MultiVector*)&X;
    OverlappingY = &Y;
  }

  if (FilterSingletons_) {
    // process singleton filter
    Epetra_MultiVector ReducedX(SingletonFilter_->Map(),NumVectors);
    Epetra_MultiVector ReducedY(SingletonFilter_->Map(),NumVectors);
    IFPACK_CHK_ERR(SingletonFilter_->SolveSingletons(*OverlappingX,*OverlappingY));
    IFPACK_CHK_ERR(SingletonFilter_->CreateReducedRHS(*OverlappingY,*OverlappingX,ReducedX));

    // process reordering
    if (!UseReordering_) {
      IFPACK_CHK_ERR(Inverse_->ApplyInverse(ReducedX,ReducedY));
    }
    else {
      Epetra_MultiVector ReorderedX(ReducedX);
      Epetra_MultiVector ReorderedY(ReducedY);
      IFPACK_CHK_ERR(Reordering_->P(ReducedX,ReorderedX));
      IFPACK_CHK_ERR(Inverse_->ApplyInverse(ReorderedX,ReorderedY));
      IFPACK_CHK_ERR(Reordering_->Pinv(ReorderedY,ReducedY));
    }

    // finish up with singletons
    IFPACK_CHK_ERR(SingletonFilter_->UpdateLHS(ReducedY,*OverlappingY));
  }
  else {
    // process reordering
    if (!UseReordering_) {
      IFPACK_CHK_ERR(Inverse_->ApplyInverse(*OverlappingX,*OverlappingY));
    }
    else {
      Epetra_MultiVector ReorderedX(*OverlappingX);
      Epetra_MultiVector ReorderedY(*OverlappingY);
      IFPACK_CHK_ERR(Reordering_->P(*OverlappingX,ReorderedX));
      IFPACK_CHK_ERR(Inverse_->ApplyInverse(ReorderedX,ReorderedY));
      IFPACK_CHK_ERR(Reordering_->Pinv(ReorderedY,*OverlappingY));
    }
  }

  if (IsOverlapping()) {
    IFPACK_CHK_ERR(OverlappingMatrix_->ExportMultiVector(*OverlappingY,Y,
                                                        CombineMode_));

    delete OverlappingX;
    delete OverlappingY;
  }

  ++NumApplyInverse_;
  ApplyInverseTime_ += Time_->ElapsedTime();

  return(0);
 
}

//==============================================================================
template<typename T>
std::ostream& Ifpack_AdditiveSchwarz<T>::
Print(std::ostream& os) const
{
  if( Matrix().Comm().MyPID())
    return(os);

  os << "*** Ifpack_AdditiveSchwarz" << endl << endl;
  os << "Cond number estimate = " << Condest_ << endl;
  os << "Overlap level        = " << OverlapLevel_ << endl;
  if (CombineMode_ == Insert)
    os << "Combine mode         = Insert" << endl;
  else if (CombineMode_ == Add)
    os << "Combine mode         = Add" << endl;
  else if (CombineMode_ == Zero)
    os << "Combine mode         = Zero" << endl;
  else if (CombineMode_ == Average)
    os << "Combine mode         = Average" << endl;
  else if (CombineMode_ == AbsMax)
    os << "Combine mode         = AbsMax" << endl;

  os << endl;
  os << "Number of initialization phases = " << NumInitialize_ << endl;
  os << "Number of computation phases    = " << NumCompute_ << endl;
  os << "Number of applications          = " << NumApplyInverse_ << endl;
  os << endl;
  os << "Total time for Initialize()     = " << InitializeTime_ << " (s)\n";
  os << "Total time for Compute()        = " << ComputeTime_ << " (s)\n";
  os << "Total time for ApplyInverse()   = " << ApplyInverseTime_ << " (s)\n";
  os << endl;

  if (Comm().MyPID() == 0) {
    os << endl << "Printing Inverse() on process 0..." << endl;
    os << *Inverse_;
  }
  os << endl;

  return(os);
}

#include "Ifpack_Condest.h"
//==============================================================================
template<typename T>
double Ifpack_AdditiveSchwarz<T>::
Condest(const Ifpack_CondestType CT, const int MaxIters, 
        const double Tol, Epetra_RowMatrix* Matrix)
{
  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  if (Condest_ == -1.0)
    Condest_ = Ifpack_Condest(*this, CT, MaxIters, Tol, Matrix);

  return(Condest_);
}

#endif // HAVE_IFPACK_TEUCHOS
#endif // IFPACK_ADDITIVESCHWARZ_H
