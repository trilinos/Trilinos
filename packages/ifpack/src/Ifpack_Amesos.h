#ifndef IFPACK_AMESOS_H
#define IFPACK_AMESOS_H

#include "Ifpack_ConfigDefs.h"
#if defined(HAVE_IFPACK_AMESOS) && defined(HAVE_IFPACK_TEUCHOS)

#include "Ifpack_Preconditioner.h"
#include "Epetra_Operator.h"
namespace Teuchos {
  class ParameterList;
}

class Epetra_Map;
class Epetra_Comm;
class Amesos_BaseSolver;
class Epetra_LinearProblem;
class Epetra_RowMatrix;

//! Ifpack_Amesos: a class to use Amesos's LU solvers as preconditioners.

class Ifpack_Amesos : public Ifpack_Preconditioner {
      
public:

  Ifpack_Amesos(Epetra_RowMatrix* Matrix);

  //@{ \name Destructor.
  //! Destructor
  virtual ~Ifpack_Amesos();
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

  virtual const Epetra_RowMatrix* Matrix() const
  {
    return(Matrix_);
  }
  
private:

  //! Pointers to the matrix to be preconditioned.
  Epetra_RowMatrix* Matrix_;
  //! Amesos solver, use to apply the inverse of the local matrix.
  Amesos_BaseSolver* Solver_;
  //! Linear problem required by Solver_.
  Epetra_LinearProblem* Problem_;
  //! Contains the label of \c this object.
  string Label_;
  //! If true, the preconditioner has been successfully computed.
  bool IsComputed_;

};

#endif // HAVE_IFPACK_AMESOS && HAVE_IFPAC_TEUCHOS
#endif // IFPACK_AMESOS_H
