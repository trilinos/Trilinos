// $Id$
// $Source$

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

#ifndef NOX_BELOS_GROUP_H
#define NOX_BELOS_GROUP_H

#include "NOX_Abstract_Group.H"         // base class and data member
#include "NOX_Utils.H"

namespace NOX { 

  //! %NOX Belos support
  namespace Belos { 

    //! %NOX %Group interface to Belos.
    /*!
     * This group is an implementation of the NOX::Abstract::Group to use
     * Belos for Jacobian solves.  Belos is a block, generic iterative
     * linear solver package in Trilinos.
     *
     * This group is designed to be a wrapper for another group tied to a 
     * specific application code.  It forwards all group methods to this
     * group except applyJacobianInverse, which it implements using Belos.
     * The iterative solve is implemented using the supplied group's 
     * applyJacobian and applyRightPreconditioning methods.
     */
    class Group : public NOX::Abstract::Group {

      public:

      //! Constructor
      Group(NOX::Abstract::Group& g, NOX::Parameter::List& printParams);

      //! Copy constructor.
      Group(const Group& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~Group();

      //! Assignment operator
      virtual Group& operator=(const Group& source);

      /*! 
       * @name Implementation of NOX::Abstract::Group virtual methods 
       */
      //@{

      //! Assignment operator
      virtual NOX::Abstract::Group& 
      operator=(const NOX::Abstract::Group& source);
     
      //! Cloning function
      virtual NOX::Abstract::Group* 
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Set the solution vector, x, to y.
      virtual void setX(const NOX::Abstract::Vector& y);

      //! Compute this.x = grp.x + step * d.
      virtual void computeX(const NOX::Abstract::Group& g, 
			    const NOX::Abstract::Vector& d,
			    double step);

      //! Compute F
      virtual NOX::Abstract::Group::ReturnType computeF();

      //! Compute Jacobian
      virtual NOX::Abstract::Group::ReturnType computeJacobian();

      //! Compute gradient
      virtual NOX::Abstract::Group::ReturnType computeGradient();

      //! Compute %Newton direction using applyJacobianInverse().
      virtual NOX::Abstract::Group::ReturnType 
      computeNewton(NOX::Parameter::List& params);
  
      //! Computes the Jacobian vector product
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobian(const NOX::Abstract::Vector& input, 
		    NOX::Abstract::Vector& result) const;

      //! Computes the Jacobian-transpose vector product
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianTranspose(const NOX::Abstract::Vector& input, 
			     NOX::Abstract::Vector& result) const;

      //! Applies the inverse of the Jacobian matrix using Belos.
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianInverse(NOX::Parameter::List& params, 
			   const NOX::Abstract::Vector& input, 
			   NOX::Abstract::Vector& result) const;

      //! Applies the right preconditioner.
      virtual NOX::Abstract::Group::ReturnType 
      applyRightPreconditioning(bool useTranspose,
				NOX::Parameter::List& params,
				const NOX::Abstract::Vector& input, 
				NOX::Abstract::Vector& result) const;

      //! applyJacobian for multiple right-hand sides 
      virtual NOX::Abstract::Group::ReturnType
      applyJacobianMultiVector(const NOX::Abstract::MultiVector& input, 
			       NOX::Abstract::MultiVector& result) const;

      //! applyJacobianTranspose for multiple right-hand sides 
      virtual NOX::Abstract::Group::ReturnType
      applyJacobianTransposeMultiVector(
				    const NOX::Abstract::MultiVector& input, 
				    NOX::Abstract::MultiVector& result) const;

      //! applyJacobianInverse for multiple right-hand sides
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianInverseMultiVector(
				    NOX::Parameter::List& params, 
				    const NOX::Abstract::MultiVector& input, 
				    NOX::Abstract::MultiVector& result) const;

      //! applyRightPreconditioning for multiple right-hand sides
      virtual NOX::Abstract::Group::ReturnType 
      applyRightPreconditioningMultiVector(
				   bool useTranspose,
				   NOX::Parameter::List& params,
				   const NOX::Abstract::MultiVector& input, 
				   NOX::Abstract::MultiVector& result) const;
  
      //! Return \c true if the residual \f$F\f$ is valid.
      virtual bool isF() const;

      //! Return \c true if the Jacobian is valid.
      virtual bool isJacobian() const;
      
      //! Return \c true if the gradient is valid.
      virtual bool isGradient() const;
      
      //! Return \c true if the Newton direction is valid.
      virtual bool isNewton() const;
      
      //! Return solution vector \f$x\f$.  
      virtual const NOX::Abstract::Vector& getX() const;

      //! Return residual \f$F\f$
      virtual const NOX::Abstract::Vector& getF() const;

      //! Return 2-norm of \f$F(x)\f$. 
      virtual double getNormF() const;

      //! Return gradient
      virtual const NOX::Abstract::Vector& getGradient() const;

      //! Return Newton direction.
      virtual const NOX::Abstract::Vector& getNewton() const;

      //! Return the norm of the %Newton solve residual.
      virtual double getNormNewtonSolveResidual() const;

      //@}

    protected:

      //! Rest all is-valid flags
      void resetIsValid();

    protected:

      //! Stores the underlying group
      NOX::Abstract::Group* grpPtr;  

      //! Flag indicating whether we own the underlying group
      bool ownsGroup;

      //! Vector to store the Newton vector
      NOX::Abstract::Vector *newtonVecPtr;

      //! Flag indicating whether Newton vector is valid
      bool isValidNewton;

      //! %Utils object
      NOX::Utils utils;

      //! ID of processor
      int myPID;
  
    };
  } // namespace Belos
} // namespace LOCA

#endif
