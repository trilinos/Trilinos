// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_EPETRA_SCHUROPERATOR_H
#define NOX_EPETRA_SCHUROPERATOR_H

#include "NOX_Epetra.H"
#include "NOX_Epetra_SchurInterface.H"

//-----------------------------------------------------------------------------

// This is a new class that may eventually get moved into NOX::Epetra.  For now,
// this is simply used as a testbed for driving new coupling algorithms - RWH 10/23/2006

namespace NOX {

namespace Epetra {

class SchurOp : public Epetra_Operator
{
public:

  SchurOp( int probId, int depId, SchurInterface & );

  virtual ~SchurOp();


  // Inherited methods from Epetra_Operator

  //! @name Attribute set methods
  //@{

    //! If set true, transpose of this operator will be applied.
    /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
        affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface
    does not support transpose use, this method should return a value of -1.

    \param In
       UseTranspose -If true, multiply by the transpose of operator, otherwise just use operator.

    \return Integer error code, set to 0 if successful.  Set to -1 if this implementation does not support transpose.
  */
    virtual int SetUseTranspose(bool UseTranspose);
  //@}

  //! @name Mathematical functions
  //@{

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
    virtual double NormInf() const;
  //@}

  //! @name Attribute access functions
  //@{

    //! Returns a character std::string describing the operator
    virtual const char * Label() const;

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

  //! @name Schur-based modifying methods
  //@{

    //! Modifies the RHS vector using block operations
    virtual void modifyRHS( Teuchos::RCP<Epetra_Vector> &, Teuchos::RCP<Epetra_Vector> & );
  //@}

protected:

  // Methods

  // Data

  //! The interface to the action of matrix blocks
  SchurInterface & schurInterface ;

  //! Matrix row block identifier
  int probId;

  //! Matrix column block identifier
  int depId;

  //! Label for operator
  std::string label;

};
} // Epetra
} // NOX

#endif
