//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER
                                                  
#ifndef NOX_EPETRA_MATRIXFREE_H
#define NOX_EPETRA_MATRIXFREE_H

#include "Epetra_Comm.h"

#include "Epetra_Operator.h"	           // base class
#include "NOX_EpetraNew_Interface_Jacobian.H" // base class

#include "NOX_Common.H"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"

// Forward Declarations
class Epetra_Comm;
class Epetra_Map;

namespace NOX {
  namespace Abstract {
    class Group;
  }
  namespace EpetraNew {
    namespace Interface {
      class Required;
    }
  }
}

namespace NOX {

namespace EpetraNew {

/*! \brief Concrete implementation for creating an Epetra_Operator Jacobian based on the Matrix-Free Newton-Krylov method.

Matrix-Free Newton-Krylov is a method that takes advantage of the fact the Newton Krylov solvers do not require an explicit Jacobian matrix.  Newton-Krylov solvers only require the matrix-vector product \f$Jy\f$ in the iteration sequence.  This product can approximated by the following:

\f[ Jy = \frac{F(x + \delta y) - F(x)}{\delta} \f]

where \f$J\f$ is the Jacobian, \f$F\f$ is the function evaluation, \f$x\f$ is the solution vector, \f$y\f$ is the vector to be operated on, and \f$\delta\f$ is a scalar perturbation calculated by:

\f[ \delta = \lambda * (\lambda + \frac{\| x\|}{\| y\|} ) \f]

where \f$ \lambda = 1.0e-6 \f$.
 */
class MatrixFree : public Epetra_Operator, 
		   public virtual NOX::EpetraNew::Interface::Jacobian {
      
 public:

  //! Define types for use of the perturbation parameter \f$ \delta\f$.
  enum DifferenceType {Forward, Backward, Centered};

  //! Constructor 
  MatrixFree(NOX::EpetraNew::Interface::Required& i, const Epetra_Vector& x);

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
  virtual bool computeJacobian(const Epetra_Vector& x);

  //! Set the type of perturbation method used (default is Forward)
  virtual void setDifferenceMethod( DifferenceType type );
 
  //! Allows the user to change the value of \f$ \lambda \f$ in the perturbation calculation.
  void setLambda(double lambda_);

  //! Flag that toggles whether MatrixFree should compute the perturbation parameter \f$ \eta \f$ or use a value supplied by the user through setPerturbation().
  void setComputePerturbation(bool bVal);

  //! Set the perturbation parameter \f$ \eta \f$.
  void setPerturbation(double eta_);
  
  //! Returns the most recently used value of the perturbation parameter \f$ \eta \f$.
  double getPerturbation() const;

  //! Register a NOX::Abstract::Group derived object and use the computeF() method of that group for the perturbation instead of the NOX::Epetra::Interface::Required::computeF() method.  This is required for LOCA to get the operators correct during homotopy.
  void setGroupForComputeF(NOX::Abstract::Group& group);

protected:

  //! Label for matrix
  string label;
  
  //! User provided interface function
  NOX::EpetraNew::Interface::Required& interface;

  //! The current solution vector 
  Epetra_Vector currentX;

  //! Perturbed solution vector 
  mutable Epetra_Vector perturbX;

  //! Function evaluation at currentX 
  mutable Epetra_Vector fo;

  //! Function evaluation at perturbX 
  mutable Epetra_Vector fp;

  //! Optional pointer to function evaluation at -perturbX - needed only for centered finite differencing
  mutable Epetra_Vector* fmPtr;

  //! Epetra_Map object used in the returns of the Epetra_Operator derived methods.
  /*! If the user is using Epetra_BlockMaps, then NOX::EpetraNew::MatrixFree must create an equivalent Epetra_Map from the Epetra_BlockMap that can be used as the return object of the OperatorDomainMap() and OperatorRangeMap() methods.
   */ 
  const Epetra_Map* epetraMap;

  //! Set to true if the epetraMap is created by NOX::EpetraNew::MatrixFree.
  /*! If the user is using Epetra_BlockMaps, then NOX::EpetraNew::MatrixFree must create an equivalent Epetra_Map from the Epetra_BlockMap that can be used as the return object of the OperatorDomainMap() and OperatorRangeMap() methods.  This flag tells the destructor whether or not to destroy the Epetra_Map.
   */ 
  bool ownsMap;

  //! Define types for use of the perturbation parameter \f$ \delta\f$.
  DifferenceType diffType;

  //! Scale factor for eta calculation
  double lambda;

  //! Perturbation value to use in the directional derivative
  mutable double eta;

  //! User specified perturbation value to use in the directional derivative.  Set by setPerturbation().
  double userEta;

  //! Flag that determines if we should calculate eta or use a value set by the user.
  bool computeEta;

  //! Flag to enables the use of a group instead of the interface for the computeF() calls in the directional difference calculation.
  bool useGroupForComputeF;
  
  //! Pointer to the group for computeF() calls.
  NOX::Abstract::Group* groupPtr;

};
} // namespace Epetra
} // namespace NOX

#endif /* NOX_EPETRA_MATRIXFREE_H */
