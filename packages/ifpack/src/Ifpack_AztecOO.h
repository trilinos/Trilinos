#ifndef IFPACK_AZTECOO_H
#define IFPACK_AZTECOO_H

#include "Ifpack_ConfigDefs.h"
#if defined(HAVE_IFPACK_TEUCHOS) && defined(HAVE_IFPACK_AZTECOO)
#include "Ifpack_Preconditioner.h"
#include "Epetra_Operator.h"
#include "Ifpack_AztecOO.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "AztecOO.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_RowMatrix.h"
#include "Teuchos_ParameterList.hpp"

template<class T>
class Ifpack_AztecOO : public Ifpack_Preconditioner {
      
public:

  Ifpack_AztecOO<T>(Epetra_RowMatrix* Matrix);

  //@{ \name Destructor.
  //! Destructor
  virtual ~Ifpack_AztecOO();
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

  //! Computes the preconditioners.
  /*! Computes the preconditioners: 
   *
   * \param In
   * List - list specifying the parameters for Jacobi. See above.
   *    
   * \return
   * 0 if successful, 1 problems occurred.
   */
  virtual int Compute();

  //! Sets all the parameters for the preconditioner.
  virtual int SetParameters(Teuchos::ParameterList& List);

  //! Returns a pointer to the internally stored matrix.
  virtual const Epetra_RowMatrix* Matrix() const
  {
    return(Matrix_);
  }
private:

  //! Pointers to the matrix to be preconditioned.
  Epetra_RowMatrix* Matrix_;
  //! AztecOO solver, use to apply the inverse of the local matrix.
  AztecOO* Solver_;
  //! Linear problem required by Solver_.
  Epetra_LinearProblem* Problem_;
  //! Contains the label of \c this object.
  string Label_;
  //! If true, the preconditioner has been successfully computed.
  bool IsComputed_;
  //! Pointer to the AztecOO preconditioner
  Ifpack_Preconditioner* Prec_;

};

//==============================================================================
template<class T>
Ifpack_AztecOO<T>::Ifpack_AztecOO(Epetra_RowMatrix* Matrix) :
  Matrix_(Matrix),
  Solver_(0),
  Problem_(0),
  IsComputed_(false),
  Prec_(0)
{
}

//==============================================================================
template<class T>
Ifpack_AztecOO<T>::~Ifpack_AztecOO()
{
  if (Problem_)
    delete Problem_;

  if (Solver_)
    delete Solver_;

  if (Prec_)
    delete Prec_;
}

#include "Ifpack_LocalRowMatrix.h"
//==============================================================================
template<class T>
int Ifpack_AztecOO<T>::SetParameters(Teuchos::ParameterList& List)
{

  Problem_ = new Epetra_LinearProblem;
  //FIXME
  Ifpack_LocalRowMatrix* LocalMatrix_ = new
    Ifpack_LocalRowMatrix(Matrix_);
  Problem_->SetOperator(LocalMatrix_);
  Solver_ = new AztecOO(*Problem_);

  Ifpack_Preconditioner* Prec_ = new T(Matrix_);

  IFPACK_CHK_ERR(Prec_->SetParameters(List));
  IFPACK_CHK_ERR(Prec_->Compute());

  Solver_->SetPrecOperator(Prec_);
  return(0);
}

//==============================================================================
template<class T>
int Ifpack_AztecOO<T>::Compute()
{

  if (Solver_ == 0)
    IFPACK_CHK_ERR(-1);

  Solver_->SetAztecOption(AZ_solver,AZ_gmres);
  Solver_->SetAztecOption(AZ_output,1);

  Solver_->SetAztecOption(AZ_precond, AZ_dom_decomp);
  Solver_->SetAztecOption(AZ_subdomain_solve, AZ_ilu);

  IsComputed_ = true;

  return(0);
}

//==============================================================================
template<class T>
int Ifpack_AztecOO<T>::SetUseTranspose(bool UseTranspose)
{
  IFPACK_CHK_ERR(-99); // not implemented
  return(-99);
}

//==============================================================================
template<class T>
int Ifpack_AztecOO<T>::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // check for maps ???
  IFPACK_CHK_ERR(Matrix_->Apply(X,Y));
}

//==============================================================================
template<class T>
int Ifpack_AztecOO<T>::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

  if (IsComputed() == false)
    IFPACK_CHK_ERR(-1);

  // check for maps??

  int NumVectors = X.NumVectors();

  if (NumVectors != Y.NumVectors())
    IFPACK_CHK_ERR(-1); // wrong input
  
  Epetra_MultiVector Xtmp(X);
  Solver_->SetLHS(&Y);
  Solver_->SetRHS((Epetra_MultiVector*)&Xtmp);
  IFPACK_CHK_ERR(Solver_->Iterate(1550, 1e-9));

  return(0);
}

//=============================================================================
template<class T>
double Ifpack_AztecOO<T>::NormInf() const
{
  return(-1.0);
}

//==============================================================================
template<class T>
char * Ifpack_AztecOO<T>::Label() const
{
  return((char*)Label_.c_str());
}

//==============================================================================
template<class T>
bool Ifpack_AztecOO<T>::UseTranspose() const
{
  return(false);
}

//==============================================================================
template<class T>
bool Ifpack_AztecOO<T>::HasNormInf() const
{
  return(false);
}

//==============================================================================
template<class T>
const Epetra_Comm & Ifpack_AztecOO<T>::Comm() const
{
  return(Matrix_->Comm());
}

//==============================================================================
template<class T>
const Epetra_Map & Ifpack_AztecOO<T>::OperatorDomainMap() const
{
  return(Matrix_->OperatorDomainMap());
}

//==============================================================================
template<class T>
const Epetra_Map & Ifpack_AztecOO<T>::OperatorRangeMap() const
{
  return(Matrix_->OperatorRangeMap());
}
#endif // HAVE_IFPACK_TEUCHOS && HAVE_IFPACK_AZTECOO
#endif // IFPACK_AZTECOO_H
