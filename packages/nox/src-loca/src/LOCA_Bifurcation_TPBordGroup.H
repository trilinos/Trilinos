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

#ifndef LOCA_BIFURCATION_TPBORDGROUP_H
#define LOCA_BIFURCATION_TPBORDGROUP_H

#include "LOCA_Abstract_Group.H"
#include "LOCA_Bifurcation_TPBordVector.H"
#include "LOCA_DerivUtils.H"

namespace LOCA {
  class ParameterVector;
}
namespace LOCA { 

  //! %LOCA bifurcation vectors and groups
  namespace Bifurcation { 

    /*! 
      \brief %LOCA's generic turning point group.  Implements the bordering
      algorithm to solve turning point equations.  
    */
    /*!
      The %TPBordGroup is a concrete implementation of the 
      LOCA::Abstract::Group that defines the following extended set of 
      equations:
      \f[
       G(z) = \left[
              \begin{array}{c}
	        F(x,p) \\
		Jn \\
		l^Tn-1
	      \end{array}
              \right] = 0
    \f]
    where \f$z = [x, n, p]\in\Re^{2n+1}\f$, \f$x\f$ is the solution vector,
    \f$n\f$ is the null vector, \f$l\f$ is the length normalization vector
    and \f$J\f$ is the Jacobian of F.

    This class implements all of the NOX::Abstract::Group and 
    LOCA::Abstract::Group methods for this extended set of equations.
    %Newton solves are computed by the bordering algorithm.
    */

    class TPBordGroup : public virtual LOCA::Abstract::Group {

    public:
  
      //! Constructor to set the base group, bifurcation paramter, length vector, and derivative objects.
      TPBordGroup(const Abstract::Group& g, 
		  const NOX::Abstract::Vector& lenVec, int paramId = 0,
		  const DerivUtils& d = DerivUtils());

      //! Copy constructor.
      TPBordGroup(const TPBordGroup& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~TPBordGroup();
  
      //! Assignment operator
      virtual LOCA::Abstract::Group& operator=(const LOCA::Abstract::Group& source);
      //! Assignment operator
      virtual NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);
      //! Assignment operator
      virtual TPBordGroup& operator=(const TPBordGroup& source);
     
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

      //! Set the solution vector, x, to y.
      virtual void setX(const TPBordVector& y);

      //! Compute and return solution vector, x, where this.x = grp.x + step * d.
      virtual void computeX(const NOX::Abstract::Group& g, 
			    const NOX::Abstract::Vector& d,
			    double step);

      //! Compute and return solution vector, x, where this.x = grp.x + step * d.
      virtual void computeX(const TPBordGroup& g, const TPBordVector& d,
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

      //! Applies the inverse of the extended Jacobian matrix using the bordering algorithm for multiple right-hand sides
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianInverseMulti(NOX::Parameter::List& params, 
			   const NOX::Abstract::Vector* const* inputs, 
			   NOX::Abstract::Vector** results, int nVecs) const;

      /*! 
	Applies Jacobian inverse of underlying group on multiple right hand sides
      */
      virtual NOX::Abstract::Group::ReturnType 
      underlyingGroupApplyJacobianInverseManager(NOX::Parameter::List& params, 
			  const NOX::Abstract::Vector& input1, 
			  const NOX::Abstract::Vector& input2, 
			  const NOX::Abstract::Vector& approxNullVec,
			  const NOX::Abstract::Vector& jacTimesApproxNullVec,
			  NOX::Abstract::Vector& result1,
			  NOX::Abstract::Vector& result2) const;

      /*! 
	Applies Jacobian inverse of underlying group on multiple right hand sides
      */
      virtual NOX::Abstract::Group::ReturnType 
      underlyingGroupApplyJacobianInverseManagerMulti(
                           NOX::Parameter::List& params, 
			   const NOX::Abstract::Vector*const* inputs, 
			   const NOX::Abstract::Vector& approxNullVec,
			   const NOX::Abstract::Vector& jacTimesApproxNullVec,
			   NOX::Abstract::Vector** results, int nVecs) const;

      //! Apply right preconditiong to the given input vector
      /*!
	This uses the same bordering algorithm as applyJacobianInverse with
	calls of the base group's applyJacobianInverse replaced by 
	applyRightPreconditioning.
      */
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

      //! Return the norm of the Newton solve residual.
      virtual double getNormNewtonSolveResidual() const;

      //@}

      //! Calls groups print function for solution vector null vec, and params 
      virtual void printSolution(const double conParam) const;

      //! Set scale vector
      virtual void setScaleVec(const NOX::Abstract::Vector& s);

      //! Set scale vector
      void setScaleVec(const LOCA::Bifurcation::TPBordVector& s);

      //! Get scale vector
      virtual const NOX::Abstract::Vector& getScaleVec() const;

    protected:

      //! Get bifurcation parameter
      double getBifParam() const;

      //! Set bifurcation parameter
      void setBifParam(double param);

    protected:

      //! Stores the underlying group, which stores the solution vector,
      //! jacobian, and solution residual
      LOCA::Abstract::Group* grpPtr;   
  
      //! Stores the turning point solution vector
      TPBordVector tpXVec;

      //! Stores the turning point residual vector
      TPBordVector tpFVec;

      //! Stores the turning point Newton vector
      TPBordVector tpNewtonVec;

      //! Stores the constant length vector
      NOX::Abstract::Vector* lengthVecPtr;

      //! Stores the scale vector
      TPBordVector tpScaleVec;

      //! Stores the bifurcation parameter index
      int bifParamId;

      //! Stores the derivative of the solution residual w.r.t. bifparam
      NOX::Abstract::Vector* derivResidualParamPtr;

      //! Stores the derivative of the null vector residual w.r.t. bifparam
      NOX::Abstract::Vector* derivNullResidualParamPtr;

      //! Stores the directional derivative object used to compute derivatives
      //! of residual equations
      DerivUtils* derivPtr;

      //! Is residual vector valid
      bool isValidF;

      //! Is Jacobian matrix valid
      bool isValidJacobian;

      //! Is Newton vector valid
      bool isValidNewton;

      //! Is Tangent vector valid
      bool isValidTangent;
    };
  } // namespace Bifurcation
} // namespace LOCA

#endif
