
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef _NOX_EPETRA_OPERATOR_H_
#define _NOX_EPETRA_OPERATOR_H_

#include "Epetra_Comm.h"
#include "Epetra_Operator.h"	// base class

#include "NOX_Common.H"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"

// Forward Declarations
class Epetra_Comm;
class Epetra_Map;

namespace NOX {

namespace Epetra {

/*! \brief Pure virtual class for creating an Operator for NOX.
  This class derives from the Epetra_Operator class with a single
  additional method: compute().    This allows the user to supply their own
  preconditioner or jacobian operators.
 */
class Operator : public Epetra_Operator {
      
 public:

  //! Constructor 
  Operator(){};

  //! Pure virtual destructor
  virtual ~Operator(){};

  //! If set true, transpose of this operator will be applied.
  /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
    affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface 
    does not support transpose use, this method should return a value of -1.      
    \param  UseTranspose -If true, multiply by the transpose of operator, otherwise just use operator.
    
    \return Integer error code, set to 0 if successful.  Set to -1 if this implementation does not support transpose.
  */
  virtual int SetUseTranspose(bool UseTranspose) = 0;

  //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
  /*! 
    \param X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Y - A Epetra_MultiVector of dimension NumVectors containing result.
    
    \return Integer error code, set to 0 if successful.
  */
  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;

  //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
  /*! 
    \param X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Y -A Epetra_MultiVector of dimension NumVectors containing result.
    
    \return Integer error code, set to 0 if successful.
    
    \warning In order to work with AztecOO, any implementation of this method must 
    support the case where X and Y are the same object.
  */
  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;
  
  //! Returns the infinity norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_\infty\f$ such that
     \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].
     
     \warning This method must not be called unless HasNormInf() returns true.    */ 
  virtual double NormInf() const = 0;
  
  //! Returns a character string describing the operator
  virtual char* Label () const = 0;
  
  //! Returns the current UseTranspose setting.
  virtual bool UseTranspose() const = 0;

  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  virtual bool HasNormInf() const = 0;

  //! Returns a reference to the Epetra_Comm communicator associated with this operator.
  virtual const Epetra_Comm & Comm() const = 0;

  //! Returns the Epetra_BlockMap object associated with the domain of this matrix operator.
  virtual const Epetra_Map& OperatorDomainMap () const = 0;

  //! Returns the Epetra_BlockMap object associated with the range of this matrix operator.
  virtual const Epetra_Map& OperatorRangeMap () const = 0;

  //! Compute the Operator given the specified input vector, x. In some cases the compute function may need additional information from an auxiliary operator, so we allow for one to be supplied.  For example, if a preconditioner is implemented that relies on a separate Jacobian matrix, it could take the Jacobian as the second argument and extract the required information.
  virtual bool compute(const Epetra_Vector& x, const Epetra_Operator* J) = 0;

};
} // namespace Epetra
} // namespace NOX

#endif /* _NOX_EPETRA_OPERATOR_H_ */
