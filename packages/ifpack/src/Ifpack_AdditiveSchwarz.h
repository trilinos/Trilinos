#ifndef IFPACK_ADDITIVESCHWARZ_H
#define IFPACK_ADDITIVESCHWARZ_H

#include "Ifpack_Preconditioner.h"
#include "Teuchos_ParameterList.hpp"
#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_LocalRowMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_RowMatrix.h"
#include "Teuchos_ParameterList.hpp"

//! Ifpack_AdditiveSchwarz: a class to define Additive Schwarz preconditioners.

/*!
  Class  Ifpack_AdditiveSchwarz enables the construction of Additive
  Schwarz preconditioners, for a given Epetra_RowMatrix.
  Ifpack_AdditiveSchwarz is derived from Ifpack_Preconditioner,
  itself derived from Epetra_Operator. An application of
  the Additive Schwarz preconditioner can be obtained 
  by calling method ApplyInverse().

  To solve on each subdomain, the user can adopt any class, derived
  from Ifpack_Preconditioner. This can be easily accomplished, as
  Ifpack_AdditiveSchwarz is templated with the solver for each subdomain.

  The following example shows how to define an additive Schwarz
  preconditioner, that use Jacobi to solve the local subproblems.

  \code
  #include "Ifpack_AdditiveSchwarz.h"
  #include "Ifpack_Jacobi.h"

  Teuchos::ParameterList List;
  List.set("omega", 0.67);
  List.set("local blocks", 4);
  List.set("overlap level", OverlapLevel);

  Ifpack_Preconditioner* Prec = new
    Ifpack_AdditiveSchwarz<Ifpack_Jacobi>(A);

  assert(Prec != 0);

  IFPACK_CHK_ERR(Prec->SetParameters(List));
  IFPACK_CHK_ERR(Prec->Compute());
  \endcode
  
  At this point, \c Prec can be used, for instance, by an AztecOO solver.
*/

template<typename T>
class Ifpack_AdditiveSchwarz : public Ifpack_Preconditioner {
      
public:

  //@{ \name Constructors/Destructors
  //! Ifpack_AdditiveSchwarz constructor with given Epetra_RowMatrix.
  /*! Creates an Ifpack_AdditiveSchwarz preconditioner with overlap.
   *
   * \param In
   * Matrix - Pointer to matrix to be preconditioned
   *
   * \param In
   * OverlappingMatrix - Pointer to the matrix extended with the
   *                     desired level of overlap.
   */
  Ifpack_AdditiveSchwarz(Epetra_RowMatrix* Matrix,
			 Epetra_RowMatrix* OverlappingMatrix = 0);
  
  //@{ \name Destructor.
  //! Destructor
  virtual ~Ifpack_AdditiveSchwarz();
  //@}

  //@{ \name Atribute set methods.

    //! If set true, transpose of this operator will be applied (not implemented).
    /*! This flag allows the transpose of the given operator to be used 
     * implicitly.  
      
    \param In
	   UseTranspose -If true, multiply by the transpose of operator, 
	   otherwise just use operator.

    \return Integer error code, set to 0 if successful.  Set to -1 if this implementation does not support transpose.
  */
    virtual int SetUseTranspose(bool UseTranspose);
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
    virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //! Applies the preconditioner to X, returns the result in Y.
  /*! 
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to be preconditioned.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectors containing result.

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

  //! Returns \c true if the preconditioner has been successfully computed.
  virtual bool IsComputed() const
  {
    return(IsComputed_);
  }

  //! Sets the parameters.
  virtual int SetParameters(Teuchos::ParameterList& List);

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

  virtual const Epetra_RowMatrix* Matrix() const
  {
    return(Matrix_);
  }

  virtual bool IsOverlapping() const
  {
    return(IsOverlapping_);
  }

protected:

  //! Sets up LocalizedMatrix_, Importer and Exporter_.
  int SetUp();

  //! Pointers to the matrix to be preconditioned.
  Epetra_RowMatrix* Matrix_;
  //! Pointers to the overlapping matrix.
  Epetra_RowMatrix* OverlappingMatrix_;
  //! Localized version of Matrix_ or OverlappingMatrix_.
  Ifpack_LocalRowMatrix* LocalizedMatrix_;
  //! Contains the label of \c this object.
  string Label_;
  //! If true, the preconditioner has been successfully computed.
  bool IsComputed_;
  //! Pointer to the local solver.
  T* Inverse_;
  //! If true, overlapping is used
  bool IsOverlapping_;
  Epetra_Import* Importer_;
  Epetra_Export* Exporter_;
  //! Stores a copy of the list given in SetParameters()
  Teuchos::ParameterList List_;

};

//==============================================================================
template<typename T>
Ifpack_AdditiveSchwarz<T>::
Ifpack_AdditiveSchwarz(Epetra_RowMatrix* Matrix,
		       Epetra_RowMatrix* OverlappingMatrix) :
  Matrix_(Matrix),
  OverlappingMatrix_(OverlappingMatrix),
  LocalizedMatrix_(0),
  IsComputed_(false),
  Inverse_(0),
  IsOverlapping_(false),
  Importer_(0),
  Exporter_(0)
{
  if ((OverlappingMatrix_ != 0) && (Matrix_->Comm().NumProc() > 1))
    IsOverlapping_ = true;
}

//==============================================================================
template<typename T>
Ifpack_AdditiveSchwarz<T>::~Ifpack_AdditiveSchwarz()
{
  if (Inverse_)
    delete Inverse_;

  if (LocalizedMatrix_)
    delete LocalizedMatrix_;

  if (Importer_)
    delete Importer_;

  if (Exporter_)
    delete Exporter_;

}

//==============================================================================
template<typename T>
int Ifpack_AdditiveSchwarz<T>::SetUp()
{

  if (OverlappingMatrix_)
    LocalizedMatrix_ = new Ifpack_LocalRowMatrix(OverlappingMatrix_);
  else
    LocalizedMatrix_ = new Ifpack_LocalRowMatrix(Matrix_);

  if (LocalizedMatrix_ == 0)
    IFPACK_CHK_ERR(-1);

  Inverse_ = new T(LocalizedMatrix_);
  if (Inverse_ == 0)
    IFPACK_CHK_ERR(-1);

  if (IsOverlapping()) {
    Importer_ = new Epetra_Import(OverlappingMatrix_->RowMatrixRowMap(), 
				  Matrix_->RowMatrixRowMap());
    Exporter_ = new Epetra_Export(Matrix_->RowMatrixRowMap(),
				  OverlappingMatrix_->RowMatrixRowMap());
    if ((Importer_ == 0) || (Exporter_ == 0))
      IFPACK_CHK_ERR(-1);
  }

  return(0);
}

//==============================================================================
template<typename T>
int Ifpack_AdditiveSchwarz<T>::SetParameters(Teuchos::ParameterList& List)
{
 
  List_ = List;

  return(0);
}

//==============================================================================
template<typename T>
int Ifpack_AdditiveSchwarz<T>::Compute()
{

  IFPACK_CHK_ERR(SetUp());

  if (Inverse_ == 0)
    IFPACK_CHK_ERR(-1);

  if (LocalizedMatrix_ == 0)
    IFPACK_CHK_ERR(-1);

  // FIXME: add overlap
  Label_ = "Ifpack Additive Schwarz";

  IFPACK_CHK_ERR(Inverse_->SetParameters(List_));

  IFPACK_CHK_ERR(Inverse_->Compute());

  IsComputed_ = true;

  return(0);
}

//==============================================================================
template<typename T>
int Ifpack_AdditiveSchwarz<T>::SetUseTranspose(bool UseTranspose)
{
  IFPACK_CHK_ERR(-99); // not implemented
  return(-99);
}

//==============================================================================
template<typename T>
int Ifpack_AdditiveSchwarz<T>::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_CHK_ERR(Matrix_->Apply(X,Y));
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

  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-1); // wrong input
  // check for maps??

  if (IsOverlapping() == false) {
    // fairly easy without overlap
    Inverse_->ApplyInverse(X,Y);
  }
  else {
    // a bit more with overlap
    Epetra_MultiVector OverlappingX(OverlappingMatrix_->RowMatrixRowMap(),
				     X.NumVectors());
    Epetra_MultiVector OverlappingY(OverlappingMatrix_->RowMatrixRowMap(),
				     Y.NumVectors());

    IFPACK_CHK_ERR(OverlappingX.Export(X,*Exporter_,Insert));
    Epetra_MultiVector LocalX(View,LocalizedMatrix_->RowMatrixRowMap(),
			      OverlappingX.Pointers(),X.NumVectors());
    Epetra_MultiVector LocalY(View,LocalizedMatrix_->RowMatrixRowMap(),
			      OverlappingY.Pointers(),Y.NumVectors());

    IFPACK_CHK_ERR(Inverse_->ApplyInverse(LocalX,LocalY));

    // FIXME: add more combine mode
    IFPACK_CHK_ERR(Y.Export(OverlappingY,*Importer_,Add));
  }

  return(0);
 
}

#endif // IFPACK_ADDITIVESCHWARZ_H
