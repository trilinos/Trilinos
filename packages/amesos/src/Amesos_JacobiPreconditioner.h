#ifndef AMESOS_JACOBIPRECONDITIONER_H
#define AMESOS_JACOBIPRECONDITIONER_H

#include "Amesos_ConfigDefs.h"
#include "Epetra_Operator.h"
class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_Map;
class Epetra_Comm;
class Epetra_Vector;
class Epetra_RowMatrix;
#include "Teuchos_ParameterList.hpp"

class Amesos_JacobiPreconditioner : public Epetra_Operator {

public:

  //@{ \name Constructor.
  Amesos_JacobiPreconditioner(const Epetra_RowMatrix* Matrix,
			      Teuchos::ParameterList& List);

  //@}
 
  //@{ \name Destructor.
  ~Amesos_JacobiPreconditioner();

  //@{ \name Atribute set methods.

  //! If set true, transpose of this operator will be applied.
  /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
    affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface 
    does not support transpose use, this method should return a value of -1.

    \param In
    UseTranspose -If true, multiply by the transpose of operator, otherwise just use operator.

    \return Integer error code, set to 0 if successful.  Set to -1 if this implementation does not support transpose.
    */
  virtual int SetUseTranspose(bool UseTranspose)
  {
    UseTranspose_ = UseTranspose;
  }

  //@}

  //@{ \name Mathematical functions.

  //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
  /*! 
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
    */
  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
  /*! 
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning In order to work with AztecOO, any implementation of this method must 
    support the case where X and Y are the same object.
    */
  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the infinity norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_\infty\f$ such that
     \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].

     \warning This method must not be called unless HasNormInf() returns true.
     */ 
  virtual double NormInf() const
  {
    return(-1.0);
  }
  //@}

  //@{ \name Atribute access functions

  //! Returns a character string describing the operator
  virtual char * Label() const
  {
    return((char*)Label_);
  }

  //! Returns the current UseTranspose setting.
  virtual bool UseTranspose() const
  {
    return(UseTranspose_);
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

  virtual bool IsComputed() const
  {
    return(IsComputed_);
  }

private:

  int NumApplications_;
  double Omega_;
  bool DebugSmoother_;
  int NumMyRows_;

  const Epetra_RowMatrix* Matrix_;
  mutable Epetra_MultiVector* AX_;
  Epetra_Vector* InvDiagonal_;

  bool UseTranspose_;
  bool IsComputed_;
  char Label_[80];

};

#endif // AMESOS_JACOBIPRECONDITIONER_H
