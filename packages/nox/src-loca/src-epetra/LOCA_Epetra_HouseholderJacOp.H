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
                                                                                
#ifndef LOCA_EPETRA_HOUSEHOLDERJACOP_H
#define LOCA_EPETRA_HOUSEHOLDERJACOP_H

#include "Epetra_Operator.h"	// base class
#include "NOX_Common.H"         // for string data member

// Forward declarations
class Epetra_Vector;
class Epetra_Map;

namespace LOCA {

  namespace Epetra {

    /*!
     * \brief %Epetra operator representing a \f$n+1\f$ bordered matrix.
     */
    /*!
     * The %LOCA::Epetra::HouseholderJacOp is an Epetra_Operator representing 
     * the matrix
     * \f[
     *     \begin{bmatrix}
     *        J & \frac{df}{dp} 
     *     \end{bmatrix} Q
     * \f]
     * where \f$J\f$ is an Epetra_Operator representing an \f$n\times n\f$
     * Jacobian matrix \f$\frac{df}{dx}\f$, and \f$Q\in\Re^{n+1\times n}\f$ is
     * orthogonal and whose columns span \f$v^{\perp}\f$ with 
     * \f$v\in\Re^{n+1}\f$
     * a supplied vector.  The action of \f$Q\f$ is defined for any 
     * \f$z\in\Re^n\f$ by
     * \f[
     *    Q z = P \begin{bmatrix} z \\ 0 \end{bmatrix}.
     * \f]
     * Here \f$P = I - \beta u u^T\in\Re^{n+1\times n+1}\f$ is a Householder
     * transformation such that \f$P v = \|v\|_2 e_{n+1}\f$ and 
     * \f$u_{n+1} = 1\f$ where \f$e_{n+1}\f$ is the \f$n+1\f$ column of the 
     * identity.
     *
     * The Apply() method performs the matrix-vector multiplication using 
     * the Householder transformation \f$P\f$ and the matrix operator \f$J\f$.
     * ApplyInverse() is not implemented.
     *
     * This operator is used to compute solutions to the system of equations
     * \f[
     *    \begin{bmatrix}
     *      J & \frac{df}{dp}
     *    \end{bmatrix} 
     *    \begin{bmatrix}
     *      \Delta x \\
     *      \Delta p
     *    \end{bmatrix} = f \quad \text{s.t.}\quad v^T\begin{bmatrix}\Delta x \\ \Delta p\end{bmatrix} = 0
     * \f]
     * in an iterative solver context for pseudo-arclength continuation 
     * where \f$v\f$ is the predictor vector.
     */
    class HouseholderJacOp : public Epetra_Operator {
      
    public:

      //! Constructor 
      /*!
       * Builds the operator using the supplied operator \em jac 
       * and Epetra_Vector's \em dfdpVec and \em houseVec_x.  Here
       * \em houseVec_x is the x (solution) component of the Householder
       * vector, \em houseVec_p is the parameter component, and \em b is
       * the scalar in the Householder transformation.
       */
      HouseholderJacOp(Epetra_Operator& jac, 
		       const Epetra_MultiVector& dfdpVec, 
		       const Epetra_MultiVector& houseVec_x,
		       double houseVec_p, double b);

      //! Destructor
      virtual ~HouseholderJacOp();

      //! This operator does not support a transpose
      /*! 
       * Setting this to true throws an error.
       */
      virtual int SetUseTranspose(bool UseTranspose);

      /*! 
       * \brief Returns the result of a Epetra_Operator applied to a 
       * Epetra_MultiVector Input in Result as described above.
       */
      virtual int Apply(const Epetra_MultiVector& Input, 
			Epetra_MultiVector& Result) const;

      /*!
       * \brief This method does nothing.
       */
      virtual int ApplyInverse(const Epetra_MultiVector& X, 
			       Epetra_MultiVector& Y) const;
  
      //! Returns an approximate infinity norm of the Householder jac matrix.
      /*!
       * This is defined only if %NormInf() of the underlying operator \f$J\f$
       * is defined and is given by 
       * \f$\|J\|_\infty+\|\frac{df}{dp}\|_\infty\f$.
       */
      virtual double NormInf() const;
    
      //! Returns a character string describing the operator
      virtual char* Label () const;
  
      //! Returns the current UseTranspose setting.  Always returns false.
      virtual bool UseTranspose() const;

      /*! 
       * \brief Returns true if the \e this object can provide an 
       * approximate Inf-norm, false otherwise.
       */
      virtual bool HasNormInf() const;

      /*! 
       * \brief Returns a reference to the Epetra_Comm communicator 
       * associated with this operator.
       */
      virtual const Epetra_Comm & Comm() const;

      /*!
       * \brief Returns the Epetra_Map object associated with the 
       * domain of this matrix operator.
       */
      virtual const Epetra_Map& OperatorDomainMap () const;

      /*! 
       * \brief Returns the Epetra_Map object associated with the 
       * range of this matrix operator.
       */
      virtual const Epetra_Map& OperatorRangeMap () const;

      /*!
       * \brief Initialize operator.  Call this before starting a linear
       * solve.  The Epetra_MultiVector argument \em x must be of the
       * same size and distribution as arguments to Apply().
       */
      void init(const Epetra_MultiVector& x);

      /*!
       * \brief Finish up solve.  Call this after a linear solve is
       * finished to inform the operator that the solve is completed.
       */
      void finish();

      /*!
       * \brief This method applies the Housholder transformation from
       * \f$\Re^{n+1}\f$ to \f$\Re^{n+1}\f$ and overwrites its arguments.
       */
      void applyHouse(Epetra_MultiVector& x, double& p) const; 

      /*!
       * \brief This version of applyHouse does not overwrite its arguments.
       */
      void applyHouse(const Epetra_MultiVector& x, double p,
		      Epetra_MultiVector& result_x, double& result_p) const;

    protected:

      //! Label for operator
      string label;
    
      //! Stores operator representing \f$J\f$
      Epetra_Operator& jacOperator;
  
      //! Stores reference to \f$df/dp\f$ vector
      const Epetra_MultiVector& dfdp;

      //! Stores reference to x component of Householder vector
      const Epetra_MultiVector& house_x;

      //! Stores parameter component of Householder vector
      double house_p;

      //! Stores constant in Householder transformation
      double beta;

      //! Temporary multi-vec used in Apply()
      Epetra_MultiVector *tmp;

    };
  } // namespace Epetra
} // namespace LOCA

#endif 
