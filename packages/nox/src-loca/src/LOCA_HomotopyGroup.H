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
    \brief %LOCA's Homotopy Algorithm.  
  */
  /*!
    The %HomotopyGroup is a concrete implementation of the 
    LOCA::Abstract::Group that modifies the set of nonlinear equations to be solved to allow for Homotopy to be applied to the system.  This object should be used in conjunction with the LOCA::Stepper object to drive the continuation.  This algorithm solves a system of nonlinear equations supplied by the user (\f$ F(x) \f$) through continuation.  An artificial parameter \f$ \lambda \f$ is used to control the continuation.  The idea is to solve a simple equation starting at \f$ \lambda \f$ = 0 and, using the solution from the previous step, solve systems of equations that gets progressively closer to the true system of interest ( at \f$ \lambda \f$ = 1.0 we recover the original equations \f$ F(x) \f$).  By constraining the definition of \f$ g(x, \lambda) \f$ and using artificial parameter contiuation, the continuation branch should be free of multiplicity and bifurcation phenomena.
 
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
   
    References:

     - ALGORITHM 652 HOMPACK: A Suite of Codes for Globally Convergent Homotopy Algorithms, Watson, L.T., Billups, S.C, and Morgan, A.P., ACM Transactions on Mathematical Software, Vol. 13, No. 3, September 1987, pp281-310.

  */

    class HomotopyGroup : public virtual LOCA::Abstract::Group {

    public:
  
      //! Constructor to set the base group and generate the "%Stepper" sublist for homotopy continuation.
      /*! 
	The locaSublist variable is the "LOCA" sublist (of type NOX::Parameter::List) that will be used in loca continuation runs.

	The variables scalarRandomVector and scalarInitialGuess are used to give some control over the generation of the random vector.  In certain instances we have seen the random vector force the solution to a set of variables that are unphysical and could break the function evaluations (cause them to return nan).  For example, in heat transfer problems, the temperature could be the dependent variable.  If the solution vector has an unphysical temperature (the random vector could force the temperature to negative or near zero values for the solution at \f$ \lambda = 0\f$) then property evaluations could break.  The random vector can be modified to keep the values near the initial guess based on values supplied to the constructor of the HomotopyGroup: 

\f[
    a = abs(r) * \mbox{scalarRandom} + x_o * \mbox{scalarInitialGuess}
\f]

where \f$ r \f$ is the random vector generated by a call to NOX::Abstract::Vector::random(), \f$ \mbox{scalarRandom} \f$ is a scalar value, \f$ x_o \f$ is the initial guess to the solution vector, and \f$ \mbox{scalarInitialGuess} \f$ is a scalar value.  The defualt values force the random vector to be calculated as:

\f[
    a = abs(r)
\f]

      IMPORTANT: For homotopy to work correctly you should not change the starting and stopping parameter values (0.0 and 1.0 respectively) set in the "%Stepper" sublist.
       */
      HomotopyGroup(NOX::Parameter::List& locaSublist,
		    const LOCA::Abstract::Group& g,
		    double scaleRandom = 1.0,
		    double scaleInitialGuess = 0.0);

      //! Constructor with a user supplied random vector.
      HomotopyGroup(NOX::Parameter::List& locaSublist,
		    const LOCA::Abstract::Group& g,
		    const NOX::Abstract::Vector& randomVector);

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

      //! Stores the underlying loca group.
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
      
      //! Value of the homotopy continuation parameter.
      /*! Ranges from 0.0 (easy solution) to 1.0 (solution to the system of 
	  interest).
      */
      double conParam;

      //! Continuatioin parameter ID number from the ParameterVector.
      int conParamID;
      
      //! Contains the string used to identify the homotopy parameter in the ParameterVector object.
      const string conParamLabel;

      //! Utility to check return types, print warnings and throw errors.
      LOCA::ErrorCheck utils;

      //! Tells whether the underlying group was created by the HomotopyGroup.  If so, the destructor should destroy the group.  
      bool ownsGroup;

    };
} // namespace LOCA

#endif
