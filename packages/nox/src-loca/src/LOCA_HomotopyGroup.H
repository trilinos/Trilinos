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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_HOMOTOPYGROUP_H
#define LOCA_HOMOTOPYGROUP_H

#include "LOCA_Abstract_Group.H"    // base class

#include "LOCA_Parameter_Vector.H"  // class data member
#include "LOCA_ErrorCheck.H"        // class data member

namespace NOX {
  namespace Abstract {
    class Vector;
  }
}

namespace LOCA {
  class ParameterVector;
}

namespace LOCA { 
  
  /*! 
    \brief %LOCA's Homotopy Group.  Implements a noninvasive homotopy
    algorithm to solve a set of nonlinear equations.  
  */
  /*!
    The %HomotopyGroup is a concrete implementation of the 
    LOCA::Abstract::Group that modifies the set of nonlinear equations to be solved to allow for a noninvasive Homotopy algorithm to be applied to the system.  This algorithm solves a system of nonlinear equations supplied by the user (\f$ F(x) \f$) through continuation.  An artificial parameter \f$ \lambda \f$ is used to control the continuation.  The idea is to solve a simple equation starting at \f$ \lambda \f$ = 0 and, using the solution from the previous step, solve systems of equations that gets progressively closer to the true system of interest ( at \f$ \lambda \f$ = 1.0 we recover the original equations \f$ F(x) \f$).  
 
    The modified system of equations, \f$ g(x, \lambda) \f$, supplied by the HomotopyGroup is defined as:

    \f[ g(x, \lambda) = \lambda F(x) + (1.0 - \lambda)(x - a) \f]
    
    where \f$x\f$ is the solution vector, \f$ \lambda \f$ is an artificial 
    parameter, \f$ F(x) \f$ is the set of nonlinear equations the user 
    supplies, \f$ g(x) \f$ is the corresponding set of 
    homotopy equations that LOCA will solve, and \f$ a \f$ is a random vector.

    This group requires the loca Stepper for continuation
    from \f$ \lambda \f$ = 0.0 (a simple set of equations to solve) to 
    \f$ \lambda \f$ = 1.0 (the set of equations requested by the user, 
    \f$ F(x) \f$).  The HomotopyGroup will generate the Stepper parameter sublist in the parameter list that is passed in to the constructor.  The user is free to modify this list (it sets default values) before passing it into the stepper object but should NOT change the starting and stopping values for the continuation parameter.
   
    This class implements all of the NOX::Abstract::Group and 
    LOCA::Abstract::Group methods for this extended set of equations.
  */

    class HomotopyGroup : public virtual LOCA::Abstract::Group {

    public:
  
      //! Constructor to set the base group and generate the "%Stepper" sublist for homotopy continuation.
      /*! 
	The locaSublist variable is the "LOCA" sublist that will be used in loca continuatioin runs.

      IMPORTANT: For homotopy to work correctly you should not change the starting and stopping parameter values (0.0 and 1.0 respectively) set in the "%Stepper" sublist.
       */
      HomotopyGroup(NOX::Parameter::List& locaSublist,
		    const LOCA::Abstract::Group& g);

      //! Copy constructor.
      HomotopyGroup(const HomotopyGroup& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~HomotopyGroup();
  
      //! Assignment operator
      virtual LOCA::Abstract::Group& operator=(const LOCA::Abstract::Group& source);
      //! Assignment operator
      virtual NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);
      //! Assignment operator
      virtual HomotopyGroup& operator=(const HomotopyGroup& source);
     
      //! Cloning function
      virtual NOX::Abstract::Group* clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Set the parameter list in the group to p (pVector = p).  
      virtual void setParams(const ParameterVector& p);

      //! Return a const reference to the ParameterVector owned by the group. 
      virtual const ParameterVector& getParams() const;

      //! Set parameter indexed by paramID
      virtual void setParam(int paramID, double val);

      //! Return copy of parameter indexed by paramID
      virtual double getParam(int paramID) const;

      //! Set parameter indexed by paramID
      virtual void setParam(string paramID, double val);

      //! Return copy of parameter indexed by paramID
      virtual double getParam(string paramID) const;

      //@{ \name "Compute" functions.

      //! Set the solution vector, x, to y.
      virtual void setX(const NOX::Abstract::Vector& y);

      //! Compute and return solution vector, x, where this.x = grp.x + step * d.
      virtual void computeX(const NOX::Abstract::Group& g, 
			    const NOX::Abstract::Vector& d,
			    double step);

      //! Compute and return solution vector, x, where this.x = grp.x + step * d.
      virtual void computeX(const LOCA::HomotopyGroup& g, 
			    const NOX::Abstract::Vector& d,
			    double step);

      //! Compute the turning point equation residual $G$ 
      virtual NOX::Abstract::Group::ReturnType computeF();

      //! Compute the blocks of the Jacobian derivative of $G$
      /*!
	Certain second derivatives such as (Jn)_x are not calculated since
	only their action on vectors are needed for later solves.
      */
      virtual NOX::Abstract::Group::ReturnType computeJacobian();

      //! Gradient computation is not defined for this group
      virtual NOX::Abstract::Group::ReturnType computeGradient();

      //! Compute %Newton direction using bordering algorithm
      virtual NOX::Abstract::Group::ReturnType 
      computeNewton(NOX::Parameter::List& params);

      //@}

      /** @name Jacobian operations. */

      //@{
  
      //! Computes the extended Jacobian vector product
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobian(const NOX::Abstract::Vector& input, 
		    NOX::Abstract::Vector& result) const;

      //! Jacobian transpose product is not defined by this group
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianTranspose(const NOX::Abstract::Vector& input, 
			     NOX::Abstract::Vector& result) const;

      //! Applies the inverse of the extended Jacobian matrix using the bordering algorithm
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianInverse(NOX::Parameter::List& params, 
			   const NOX::Abstract::Vector& input, 
			   NOX::Abstract::Vector& result) const;

      virtual NOX::Abstract::Group::ReturnType 
      applyRightPreconditioning(NOX::Parameter::List& params,
				const NOX::Abstract::Vector& input, 
				NOX::Abstract::Vector& result) const;
      //@}

      /** @name "Is" functions.
   
      Checks to see if various objects have been computed. Returns \c true
      if the corresponding "compute" function has been called since the
      last update to the solution vector (via instantiation or
      computeX).
      */

      //@{
  
      //! Return \c true if the extended residual \f$G\f$ is valid.
      virtual bool isF() const;

      //! Return \c true if the extended Jacobian is valid.
      virtual bool isJacobian() const;

      //! Always returns false.
      virtual bool isGradient() const;

      //! Return \c true if the extended Newton direction is valid.
      virtual bool isNewton() const;
      //@}

      /** @name "Get" functions.
   
      Note that these function do not check whether or not the vectors
      are valid. Must use the "Is" functions for that purpose.
      */
      //@{ 

      //! Return extended solution vector \f$z\f$.  
      virtual const NOX::Abstract::Vector& getX() const;

      //! Return extended equation residual \f$G(z)\f$
      virtual const NOX::Abstract::Vector& getF() const;

      //! Return 2-norm of \f$G(z)\f$. 
      virtual double getNormF() const;

      //! Vector returned is not valid.
      virtual const NOX::Abstract::Vector& getGradient() const;

      //! Return extended Newton direction.
      virtual const NOX::Abstract::Vector& getNewton() const;

      //@}

      //! Calls groups print function for solution vector null vec, and params 
      virtual void printSolution(const double conParam) const;

      //! Set scale vector
      virtual void setScaleVec(const NOX::Abstract::Vector& s);

      //! Get scale vector
      virtual const NOX::Abstract::Vector& getScaleVec() const;

    protected:
      
      //! Reset the isValid flags to false.  
      /*! This is called when the solution vector or parameter vector 
	is changed.
      */
      void resetIsValidFlags();

      //! Creates and sets the "Stepper" parameter sublist
      void setStepperParameters(NOX::Parameter::List& params);

    protected:

      //! Stores the underlying group, which stores the solution vector,
      //! jacobian, and solution residual
      LOCA::Abstract::Group* grpPtr;   
  
      //! Stores the homotopy residual vector, \f$ g \f$.
      NOX::Abstract::Vector* gVecPtr;

      //! Stores the random Vector, \f$ a \f$.
      NOX::Abstract::Vector* randomVecPtr;

      //! Stores the homotopy Newton vector, \f$ \frac{\patrial r}{\partial x} \f$.
      NOX::Abstract::Vector* newtonVecPtr;

      //! Stores the homotopy gradient vector if needed, \f$ \frac{\partial r}{\partial x} \f$.
      NOX::Abstract::Vector* gradVecPtr;

      //! Is residual vector valid
      bool isValidF;

      //! Is Jacobian matrix valid
      bool isValidJacobian;

      //! Is Newton vector valid
      bool isValidNewton;

      //! Is gradient vector valid
      bool isValidGradient;
      
      //! Copy of the ParameterVector for the underlying grpPtr.
      /*! We copy this and then add the homotopy parameter to the list. */
      LOCA::ParameterVector paramVec;
      
      //! Continuation parameter.  For homotopy this parameter will go
      //! from 0.0 (easy solution) to 1.0 (solution to the system of 
      //! interest).
      double conParam;

      //! Continuatioin parameter ID number from the ParameterVector.
      int conParamID;
      
      //! Contains the string used to identify the homotopy parameter in the ParameterVector object.
      const string conParamLabel;

      //! Utility to check return types, print warnings and throw errors.
      LOCA::ErrorCheck utils;

    };
} // namespace LOCA

#endif
