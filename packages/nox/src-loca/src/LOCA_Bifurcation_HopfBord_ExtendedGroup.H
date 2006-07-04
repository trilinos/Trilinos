// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
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

#ifndef LOCA_BIFURCATION_HOPFBORD_EXTENDEDGROUP_H
#define LOCA_BIFURCATION_HOPFBORD_EXTENDEDGROUP_H

#include "LOCA_Continuation_AbstractGroup.H"           // base class
#include "LOCA_Extended_AbstractGroup.H"               // base class
#include "LOCA_Bifurcation_HopfBord_ExtendedVector.H"  // class data element

namespace LOCA {
  class ParameterVector;
  
  namespace Bifurcation {
    namespace HopfBord {
      class AbstractGroup;
    }
  }
}

namespace LOCA { 

  namespace Bifurcation { 

    namespace HopfBord {

      /*! 
	\brief %LOCA's Hopf point group.  Implements the bordering
	algorithm to solve turning point equations.  
      */
      /*!
	The %HopfBord::Group is a concrete implementation of the 
	LOCA::Continuation::AbstractGroup that defines the following extended 
	set of equations:
	\f[
	G(z) = \left[
	\begin{array}{c}
	F(x,p) \\
	Jy-\omega Bz \\
	Jz+\omega By \\
	l^Ty-1 \\
	l^Tz
	\end{array}
	\right] = 0
	\f]
	where \f$z = [x, y, z, \omega, p]\in\Re^{3n+2}\f$, \f$x\f$ is the 
	solution vector, \f$y+i\omega z\f$ is the complex eigenvector of \f$J\f$
	with corresponding eigenvalues \f$\pm i\omega\f$, \f$l\f$ is the length 
	normalization vector and \f$J\f$ is the Jacobian of F w.r.t \f$x\f$.

	This class implements all of the NOX::Abstract::Group and 
	LOCA::Continuation::AbstractGroup methods for this extended set of 
	equations. %Newton solves are computed by the bordering algorithm.
      */

      class ExtendedGroup : 
	public LOCA::Continuation::AbstractGroup,
	public LOCA::Extended::AbstractGroup {

      public:

	/*! 
	 * \brief Constructor with initial data passed through parameter lists.
	 */
	ExtendedGroup(LOCA::Bifurcation::HopfBord::AbstractGroup& g, 
		      Teuchos::ParameterList& bifParamList);
  
	//! Constructor to set the base group, bifurcation paramter, length vector.
	ExtendedGroup(LOCA::Bifurcation::HopfBord::AbstractGroup& g,
		      const NOX::Abstract::Vector& realEigenVec,
		      const NOX::Abstract::Vector& imaginaryEigenVec,
		      NOX::Abstract::Vector& lenVec, double frequency, 
		      int paramId = 0);

	//! Constructor to set the base group, bifurcation paramter, length vector.
	ExtendedGroup(const LOCA::Bifurcation::HopfBord::AbstractGroup& g,
		      const NOX::Abstract::Vector& realEigenVec,
		      const NOX::Abstract::Vector& imaginaryEigenVec,
		      NOX::Abstract::Vector& lenVec, double frequency, 
		      int paramId = 0);

	//! Copy constructor.
	ExtendedGroup(const ExtendedGroup& source, 
		      NOX::CopyType type = NOX::DeepCopy);

	//! Destructor.
	virtual ~ExtendedGroup();
  
	//! Assignment operator
	virtual LOCA::Continuation::AbstractGroup& 
	operator=(const LOCA::Continuation::AbstractGroup& source);

	//! Assignment operator
	virtual LOCA::Extended::AbstractGroup& 
	operator=(const LOCA::Extended::AbstractGroup& source);

	//! Assignment operator
	virtual NOX::Abstract::Group& 
	operator=(const NOX::Abstract::Group& source);

	//! Assignment operator
	virtual ExtendedGroup& operator=(const ExtendedGroup& source);
     
	//! Cloning function
	virtual NOX::Abstract::Group* 
	clone(NOX::CopyType type = NOX::DeepCopy) const;

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

	//! Compute \f$dF/dp\f$ for the extended system
	virtual NOX::Abstract::Group::ReturnType
	computeDfDp(int paramID, NOX::Abstract::Vector& result);

	//@{ \name "Compute" functions.

	//! Set the solution vector, x, to y.
	virtual void setX(const NOX::Abstract::Vector& y);

	//! Set the solution vector, x, to y.
	virtual void setX(const LOCA::Bifurcation::HopfBord::ExtendedVector& y);

	//! Compute and return solution vector, x, where this.x = grp.x + step * d.
	virtual void computeX(const NOX::Abstract::Group& g, 
			      const NOX::Abstract::Vector& d,
			      double step);

	//! Compute and return solution vector, x, where this.x = grp.x + step * d.
	virtual void computeX(const ExtendedGroup& g, 
			      const LOCA::Bifurcation::HopfBord::ExtendedVector& d,
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
	computeNewton(Teuchos::ParameterList& params);

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
	applyJacobianInverse(Teuchos::ParameterList& params, 
			     const NOX::Abstract::Vector& input, 
			     NOX::Abstract::Vector& result) const;

	//! Applies the inverse of the extended Jacobian matrix using the bordering algorithm for multiple right-hand sides
	virtual NOX::Abstract::Group::ReturnType 
	applyJacobianInverseMulti(Teuchos::ParameterList& params, 
				  const NOX::Abstract::Vector* const* inputs, 
				  NOX::Abstract::Vector** results, int nVecs) const;
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

	//! Calls groups print function for solution vector null vec, and params 
	virtual void printSolution(const NOX::Abstract::Vector& x_,
				   const double conParam) const;

	//! Return underlying group
	virtual const LOCA::Continuation::AbstractGroup& 
	getUnderlyingGroup() const;
      
	//! Return underlying group
	virtual LOCA::Continuation::AbstractGroup& 
	getUnderlyingGroup();

      protected:

	//! Initializes group
	void init(bool perturbSoln = false, double perturbSize = 0.0);

	//! Get bifurcation parameter
	double getBifParam() const;

	//! Set bifurcation parameter
	void setBifParam(double param);

	//! Defines null vector normalization \f$l^Tn\f$
	double lTransNorm(const NOX::Abstract::Vector& n) const;

      protected:

	//! Stores the underlying group, which stores the solution vector,
	//! jacobian, and solution residual
	LOCA::Bifurcation::HopfBord::AbstractGroup* grpPtr;   
  
	//! Stores the Hopf point solution vector
	LOCA::Bifurcation::HopfBord::ExtendedVector hopfXVec;

	//! Stores the Hopf point residual vector
	LOCA::Bifurcation::HopfBord::ExtendedVector hopfFVec;

	//! Stores the Hopf point Newton vector
	LOCA::Bifurcation::HopfBord::ExtendedVector hopfNewtonVec;

	//! Stores the constant length vector
	NOX::Abstract::Vector* lengthVecPtr;

	//! Stores the bifurcation parameter index
	int bifParamId;

	//! Stores the derivative of the solution residual w.r.t. bifparam
	NOX::Abstract::Vector* derivResidualParamPtr;

	//! Stores the real component of derivative of the eigenvector residual 
	//! w.r.t. bifparam
	NOX::Abstract::Vector* derivRealEigenResidualParamPtr;

	//! Stores the imaginary component of derivative of the eigenvector 
	//! residual w.r.t. bifparam
	NOX::Abstract::Vector* derivImagEigenResidualParamPtr;

	//! Stores the mass matrix times the real component of 
	//! the eigenvector, i.e., B*y
	NOX::Abstract::Vector* massTimesYPtr;

	//! Stores the negative of the mass matrix times the imaginary component
	//! of the eigenvector, i.e., -B*z
	NOX::Abstract::Vector* minusMassTimesZPtr;

	//! Flag indicating whether we have our own copy of the underlying group
	bool ownsGroup;

	//! Is residual vector valid
	bool isValidF;

	//! Is Jacobian matrix valid
	bool isValidJacobian;

	//! Is Newton vector valid
	bool isValidNewton;
      };
    } // namespace HopfBord
  } // namespace Bifurcation
} // namespace LOCA

#endif
