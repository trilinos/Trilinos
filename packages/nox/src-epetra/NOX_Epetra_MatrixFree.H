
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

#ifndef _NOX_EPETRA_MATRIXFREE_H_
#define _NOX_EPETRA_MATRIXFREE_H_

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
class Interface;
}
}

namespace NOX {

namespace Epetra {

/*! \brief Concrete implementation for creating an Epetra_Operator Jacobian based on the Matrix-Free Newton-Krylov method.

Matrix-Free Newton-Krylov is a method that takes advantage of the fact the Newton Krylov solvers do not require an explicit Jacobian matrix.  Newton-Krylov solvers only require the matrix-vector product \f$Jy\f$ in the iteration sequence.  This product can approximated by the following:

\f[ Jy = \frac{F(x + \delta y) - F(x)}{\delta} \f]

where \f$J\f$ is the Jacobian, \f$F\f$ is the function evaluation, \f$x\f$ is the solution vector, \f$y\f$ is the vector to be operated on, and \f$\delta\f$ is a scalar perturbation calculated by:

\f[ \delta = \lambda * (\lambda + \frac{\| x\|}{\| y\|} ) \f]

where \f$ \lambda = 1.0e-6 \f$.
 */
class MatrixFree : public Epetra_Operator {
      
 public:

  //! Constructor 
  MatrixFree(Interface& i, const Epetra_Vector& x);

  //! Pure virtual destructor
  virtual ~MatrixFree();

  //! If set true, transpose of this operator will be applied.
  /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
    affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface 
    does not support transpose use, this method should return a value of -1.      
    \param     UseTranspose -If true, multiply by the transpose of operator, otherwise just use operator.
    
    \return Integer error code, set to 0 if successful.  Set to -1 if this implementation does not support transpose.
  */
  virtual int SetUseTranspose(bool UseTranspose);

  //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
  /*! 
    \param     X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param     Y -A Epetra_MultiVector of dimension NumVectors containing result.
    
    \return Integer error code, set to 0 if successful.
  */
  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
  /*! 
    \param     X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param     Y -A Epetra_MultiVector of dimension NumVectors containing result.
    
    \return Integer error code, set to 0 if successful.
    
    \warning In order to work with AztecOO, any implementation of this method must 
    support the case where X and Y are the same object.
  */
  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  
  //! Returns the infinity norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_\infty\f$ such that
     \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].
     
     \warning This method must not be called unless HasNormInf() returns true.    */ 
  virtual double NormInf() const;
  
  //! Returns a character string describing the operator
  virtual char* Label () const;
  
  //! Returns the current UseTranspose setting.
  virtual bool UseTranspose() const;

  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  virtual bool HasNormInf() const;

  //! Returns a reference to the Epetra_Comm communicator associated with this operator.
  virtual const Epetra_Comm & Comm() const;

  //! Returns the Epetra_BlockMap object associated with the domain of this matrix operator.
  virtual const Epetra_Map& OperatorDomainMap () const;

  //! Returns the Epetra_BlockMap object associated with the range of this matrix operator.
  virtual const Epetra_Map& OperatorRangeMap () const;

  //! Compute Jacobian given the specified input vector, x.  Returns true if computation was successful.
  virtual bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);
 
protected:

  //! Label for matrix
  string label;
  
  //! User provided interface function
  Interface& interface;

  //! The current solution vector 
  Epetra_Vector currentX;

  //! Perturbed solution vector 
  mutable Epetra_Vector perturbX;

  //! Function evaluation at currentX 
  mutable Epetra_Vector fo;

  //! Function evaluation at perturbX 
  mutable Epetra_Vector fp;

  //! Epetra_Map object used in the returns of the Epetra_Operator derived methods.
  /*! If the user is using Epetra_BlockMaps, then NOX::Epetra::MatrixFree must create an equivalent Epetra_Map from the Epetra_BlockMap that can be used as the return object of the OperatorDomainMap() and OperatorRangeMap() methods.
   */ 
  const Epetra_Map* epetraMap;

  //! Set to true if the epetraMap is created by NOX::Epetra::MatrixFree.
  /*! If the user is using Epetra_BlockMaps, then NOX::Epetra::MatrixFree must create an equivalent Epetra_Map from the Epetra_BlockMap that can be used as the return object of the OperatorDomainMap() and OperatorRangeMap() methods.  This flag tells the destructor whether or not to destroy the Epetra_Map.
   */ 
  bool ownsMap;

};
} // namespace Epetra
} // namespace NOX

#endif /* _EPETRA_ROWMATRIX_H_ */
