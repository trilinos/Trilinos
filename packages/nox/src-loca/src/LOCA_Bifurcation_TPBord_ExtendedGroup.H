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

#ifndef LOCA_BIFURCATION_TPBORD_EXTENDEDGROUP_H
#define LOCA_BIFURCATION_TPBORD_EXTENDEDGROUP_H

#include "LOCA_Continuation_AbstractGroup.H"         // base class
#include "LOCA_Extended_AbstractGroup.H"             // base class
#include "LOCA_Bifurcation_TPBord_ExtendedVector.H"  // class data element

namespace LOCA {
  class ParameterVector;

  namespace Bifurcation {
    namespace TPBord {
      class AbstractGroup;
    }
  }
}
namespace LOCA { 

  namespace Bifurcation { 
    
    namespace TPBord {

      /*! 
       * \brief %LOCA's generic turning point group.  Implements the bordering
       * algorithm to solve turning point equations.  
       */
      /*!
       The LOCA::Bifurcation::TPBord::Group is a concrete implementation of 
       the NOX::Abstract::Group, LOCA::Continuation::AbstractGroup and 
       LOCA::Extended::AbstractGroup that defines the following extended 
       set of equations that are regular at a generic turning point:
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

       The group stores an underlying group of type 
       LOCA::Bifurcation::TPBord::AbstractGroup to represent the equations
       \f$F(x,p) = 0\f$ and to manipulate the underlying Jacobian \f$J\f$.
       Note that the entire extended Jacobian \f$D_z G\f$ is not stored in
       memory, rather a block-elimination algorithm (bordering algorithm) is
       used to compute linear solves of the extended Jacobian (see 
       LOCA::Bifurcation::TPBord::ExtendedGroup::applyJacobianInverse(),
       LOCA::Bifurcation::TPBord::ExtendedGroup::applyJacobianInverseMulti()
       for more details).

       This class implements all of the NOX::Abstract::Group, 
       LOCA::Continuation::AbstractGroup, and  LOCA::Extended::AbstractGroup 
       methods for this extended set of equations and therefore is a complete
       group which can be passed to most %NOX solvers to locate a single 
       turning point or to the LOCA::Stepper to compute a family of turning 
       points in a second parameter.

       However, Jacobian-tranpose operations and gradient calculations cannot 
       be implemented efficiently and therefore gradient-base nonlinear solvers
       such as steepest descent and Trust region methods cannot be used to 
       solve the extended turning point equations.
       */

      class ExtendedGroup 
	: public LOCA::Continuation::AbstractGroup,
	  public LOCA::Extended::AbstractGroup {

      public:

	/*! 
	 * \brief Constructor with initial data passed through parameter lists.
	 */
	ExtendedGroup(LOCA::Bifurcation::TPBord::AbstractGroup& g, 
		      NOX::Parameter::List& bifParamList);
  
	/*! 
	 * \brief Constructor to set the base group, length normalization 
	 * vector \f$l\f$, initial guess for the null vector \f$n\f$, and 
	 * bifurcation paramter id (integer).
	 */
	/*!
	 * The initial guess to the null vector is re-scaled so that the 
	 * normalization condition \f$l^Tn = 1\f$ is satisfied.
	 *
	 * In this version of the constructor, a non-const reference to the
	 * underlying group \em g is passed and therefore the group will
	 * not be copied.  
	 */
	ExtendedGroup(LOCA::Bifurcation::TPBord::AbstractGroup& g, 
		      const NOX::Abstract::Vector& lenVec, 
		      const NOX::Abstract::Vector& nullVec, int paramId);

	/*! 
	 * \brief Constructor to set the base group, length normalization 
	 * vector \f$l\f$, initial guess for the null vector \f$n\f$, and 
	 * bifurcation paramter id (integer).
	 */
	/*!
	 * The initial guess to the null vector is re-scaled so that the 
	 * normalization condition \f$l^Tn = 1\f$ is satisfied.
	 *
	 * In this version of the constructor, a const reference to the
	 * underlying group \em g is passed and therefore the group will
	 * be copied.  
	 */
	ExtendedGroup(const LOCA::Bifurcation::TPBord::AbstractGroup& g, 
		      const NOX::Abstract::Vector& lenVec, 
		      const NOX::Abstract::Vector& nullVec, int paramId);

	//! Copy constructor.
	ExtendedGroup(const ExtendedGroup& source, 
		      NOX::CopyType type = NOX::DeepCopy);

	//! Destructor.
	virtual ~ExtendedGroup();

	//! Assignment operator
	virtual ExtendedGroup& operator=(const ExtendedGroup& source);

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

		//! Set the solution vector, x, to y.
	virtual void setX(const LOCA::Bifurcation::TPBord::ExtendedVector& y);

	//! Compute this.x = grp.x + step * d.
	virtual void computeX(const NOX::Abstract::Group& g, 
			      const NOX::Abstract::Vector& d,
			      double step);

	//! Compute this.x = grp.x + step * d.
	virtual void computeX(
			  const LOCA::Bifurcation::TPBord::ExtendedGroup& g, 
			  const LOCA::Bifurcation::TPBord::ExtendedVector& d,
			  double step);

	//! Compute the turning point equation residual \f$G\f$.
	/*!
	  This method fills the extended residual
	  \f[
	      G(z) = \left[
	      \begin{array}{c}
	         F(x,p) \\
	         Jn \\
	         l^Tn-1
	      \end{array}
	      \right].
	  \f]
	  The solution component residual \f$F(x,p)\f$ and the null-vector 
	  residual \f$Jn\f$ are calculated via the computeF and applyJacobian 
	  methods of the underlying group.
	*/
	virtual NOX::Abstract::Group::ReturnType computeF();

	//! Compute the blocks of the Jacobian derivative of \f$G\f$
	/*!
	  This method computes the \f$J\f$, \f$\partial F/\partial p\f$, and
	  \f$\partial Jn/\partial p\f$ blocks of the extended Jacobian:
	  \f[
	      D_z G(z) = 
	      \begin{bmatrix}
	         J & 0 & \frac{\partial F}{\partial p} \\
		 \frac{\partial Jn}{\partial x} & J & \frac{\partial Jn}{\partial p} \\
		 0 & l^T & 0
	      \end{bmatrix}
	  \f]
	  by calling the computeJacobian, computeDfDp, and computeDJnDp
	  methods of the underlying group.  The second derivative matrix
	  \f$\partial Jn/\partial x\f$ is not calculated since only its
	  action on vectors is needed for linear solves using the bordering
	  algorithm.
	*/
	virtual NOX::Abstract::Group::ReturnType computeJacobian();

	//! Gradient computation is not defined for this group
	virtual NOX::Abstract::Group::ReturnType computeGradient();

	//! Compute %Newton direction using applyJacobianInverse().
	virtual NOX::Abstract::Group::ReturnType 
	computeNewton(NOX::Parameter::List& params);
  
	//! Computes the extended Jacobian vector product
	/*! 
	  This method computes the extended Jacobian vector product
	  \f[
	      \begin{bmatrix}
	         J & 0 & \frac{\partial F}{\partial p} \\
		 \frac{\partial Jn}{\partial x} & J & \frac{\partial Jn}{\partial p} \\
		 0 & l^T & 0
	      \end{bmatrix}
	      \begin{bmatrix}
	         a \\
		 b \\
		 c
	      \end{bmatrix}
	      = 
	      \begin{bmatrix}
	         Ja + \frac{\partial F}{\partial p}c \\
		 \frac{\partial Jn}{\partial x}a + Jb + \frac{\partial Jn}{\partial p}c \\
		 l^T b
	      \end{bmatrix}
	  \f]
	  using the applyJacobian and computeDJnDxa methods of the underlying
	  group where \f$a\f$, \f$b\f$, and \f$c\f$ are the solution, 
	  null-vector, and paramter components of the given vector \em input.
	  Vectors \em input and \em result must be of type 
	  LOCA::Bifurcation::TPBord::ExtendedVector, otherwise an error is 
	  thrown.
	*/
	virtual NOX::Abstract::Group::ReturnType 
	applyJacobian(const NOX::Abstract::Vector& input, 
		      NOX::Abstract::Vector& result) const;

	//! Jacobian transpose product is not defined by this group
	virtual NOX::Abstract::Group::ReturnType 
	applyJacobianTranspose(const NOX::Abstract::Vector& input, 
			       NOX::Abstract::Vector& result) const;

	/*! 
	 * \brief Applies the inverse of the extended Jacobian matrix using 
	 * the bordering algorithm.
	 */
	/*!
	 * This method is a special case of applyJacobianInverseMulti() for
	 * a single right-hand-side.
	 */
	virtual NOX::Abstract::Group::ReturnType 
	applyJacobianInverse(NOX::Parameter::List& params, 
			     const NOX::Abstract::Vector& input, 
			     NOX::Abstract::Vector& result) const;

	//! Applies the extended right preconditioner.
	/*!
	 * This method is implemented using a similar bordering algorithm
	 * to applyJacobianInverse(), except calls to the underlying group's
	 * applyJacobianInvese method are replaced by 
	 * applyRightPreconditioning.
	 */
	virtual NOX::Abstract::Group::ReturnType 
	applyRightPreconditioning(bool useTranspose,
				  NOX::Parameter::List& params,
				  const NOX::Abstract::Vector& input, 
				  NOX::Abstract::Vector& result) const;
  
	//! Return \c true if the extended residual \f$G\f$ is valid.
	virtual bool isF() const;

	//! Return \c true if the extended Jacobian is valid.
	virtual bool isJacobian() const;

	//! Always returns false.
	virtual bool isGradient() const;

	//! Return \c true if the extended Newton direction is valid.
	virtual bool isNewton() const;

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

	//! Return the norm of the %Newton solve residual.
	virtual double getNormNewtonSolveResidual() const;

	//@}

	/*! 
	 * @name Implementation of LOCA::Continuation::AbstractGroup virtual methods 
	 */
	//@{

	//! Assignment operator
	virtual LOCA::Continuation::AbstractGroup& 
	operator=(const LOCA::Continuation::AbstractGroup& source);

	//! Set the parameter vector in the group to p.
	virtual void setParams(const ParameterVector& p);

	//! Return a const reference to the paramter vector owned by the group.
	virtual const ParameterVector& getParams() const;

	//! Set parameter indexed by paramID
	virtual void setParam(int paramID, double val);

	//! Return copy of parameter indexed by paramID
	virtual double getParam(int paramID) const;

	//! Set parameter indexed by paramID
	virtual void setParam(string paramID, double val);

	//! Return copy of parameter indexed by paramID
	virtual double getParam(string paramID) const;

	/*!
	 * \brief Compute parameter derivative \f$dG/d\lambda\f$ for the 
	 * extended system.
	 */
	/*!
	   For a given parameter \f$\lambda\f$, this method computes the 
	   derivative
	   \f[
	      \frac{\partial G(z,\lambda)}{\partial\lambda} = \left[
	      \begin{array}{c}
	         \frac{\partial F(x,p,\lambda)}{\partial\lambda} \\
	         \frac{\partial Jn}{\partial\lambda} \\
	         0
	      \end{array}
	      \right].
	  \f]
	*/
	virtual NOX::Abstract::Group::ReturnType
	computeDfDp(int paramID, NOX::Abstract::Vector& result);

	//! applyJacobian for multiple right-hand sides 
	/*! 
	 * The default implementation here calls applyJacobian() for
	 * each right hand side serially but should be overloaded if a 
	 * block solver is available.
	 */
// 	virtual NOX::Abstract::Group::ReturnType
// 	applyJacobianMultiVector(const NOX::Abstract::MultiVector& input, 
// 				 NOX::Abstract::MultiVector& result) const;

	/*! 
	 * \brief Applies the inverse of the extended Jacobian matrix using 
	 * the bordering algorithm for multiple right-hand sides.
	 */
	/*!
	   This method performs the following multiple right-hand-side linear 
	   solve of the extended Jacobian matrix:
	   \f[
	      \begin{bmatrix}
	         J & 0 & \frac{\partial F}{\partial p} \\
		 \frac{\partial Jn}{\partial x} & J & \frac{\partial Jn}{\partial p} \\
		 0 & l^T & 0
	      \end{bmatrix}
	      \begin{bmatrix}
	         U \\
		 V \\
		 w^T
	      \end{bmatrix}
	      = 
	      \begin{bmatrix}
	         X \\
		 Y \\
		 z^T
	      \end{bmatrix}
	  \f]
	  where \f$U\f$, \f$V\f$, \f$X\f$, and \f$Y\f$ each have \em nVecs 
	  columns and \f$w\f$, and \f$z\f$ each have \em nVecs components
	  (i.e., each column of \f$X\f$, \f$Y\f$, and \f$z^T\f$ have the 
	  solution, null-vector, and parameter components of one 
	  right-hand-side).  The solution components \f$U\f$, \f$V\f$, and 
	  \f$w\f$ are computed using a multiple right-hand-side block
	  elimination (bordering) algorithm:
	  \f[
	       \begin{split}
	          J\begin{bmatrix}
		      A & b
		   \end{bmatrix} &= 
		   \begin{bmatrix} 
		      X & \frac{\partial F}{\partial p}
		   \end{bmatrix} \\
		  J\begin{bmatrix}
		      C & d
		   \end{bmatrix} &= 
		   \begin{bmatrix} 
		      Y & \frac{\partial Jn}{\partial p}
		   \end{bmatrix} - 
		  \frac{\partial Jn}{\partial x} 
		   \begin{bmatrix} 
		      A & b
		   \end{bmatrix} \\
		  w &= \frac{C^Tl - z}{l^Td} \\
		  U &= A - b w^T \\
		  V &= C - d w^T
	       \end{split}
	  \f]
	  which requires \f$ 2(\text{nVecs} + 1) \f$ solves of the underlying 
	  Jacobian \f$J\f$ calculated by 2 calls to the 
	  applyJacobianInverseMulti method of the underlying group.

	  Each vector in the arrays \em inputs and \em outputs must be of
	  type LOCA::Bifurcation::TPBord::ExtendedVector, otherwise an 
	  error is thrown.
	*/
	virtual NOX::Abstract::Group::ReturnType 
	applyJacobianInverseMulti(NOX::Parameter::List& params, 
				  const NOX::Abstract::Vector* const* inputs, 
				  NOX::Abstract::Vector** results, int nVecs) const;

	/*! 
	 * \brief Function to print out extended solution and continuation 
	 * parameter after successful continuation step.
	 */
	/*!
	 * This method prints the solution, null-vector, and parameter 
	 * components of the extended solution vector using the printSolution
	 * method of the underlying group.
	 */
	virtual void printSolution(const double conParam) const;

	//@}

	/*! 
	 * @name Implementation of LOCA::Extended::AbstractGroup virtual methods 
	 */
	//@{

	//! Assignment operator
	virtual LOCA::Extended::AbstractGroup& 
	operator=(const LOCA::Extended::AbstractGroup& source);

	//! Return underlying group
	virtual const LOCA::Continuation::AbstractGroup& 
	getUnderlyingGroup() const;
      
	//! Return underlying group
	virtual LOCA::Continuation::AbstractGroup& 
	getUnderlyingGroup();
	
	//@}

	//! Get bifurcation parameter
	double getBifParam() const;

      protected:

	//! Initializes group
	void init(bool perturbSoln = false, double perturbSize = 0.0);

	//! Set bifurcation parameter
	void setBifParam(double param);

	//! Defines null vector normalization \f$l^Tn\f$
	double lTransNorm(const NOX::Abstract::Vector& n) const;

      protected:

	/*! 
	 * \brief Stores the underlying group, which stores the solution 
	 * vector, jacobian, and solution residual
	 */
	LOCA::Bifurcation::TPBord::AbstractGroup* grpPtr;   
  
	//! Stores the turning point solution vector
	LOCA::Bifurcation::TPBord::ExtendedVector tpXVec;

	//! Stores the turning point residual vector
	LOCA::Bifurcation::TPBord::ExtendedVector tpFVec;

	//! Stores the turning point Newton vector
	LOCA::Bifurcation::TPBord::ExtendedVector tpNewtonVec;

	//! Stores the constant length vector
	NOX::Abstract::Vector* lengthVecPtr;

	//! Stores the bifurcation parameter index
	int bifParamId;

	//! Stores the derivative of the solution residual w.r.t. bifparam
	NOX::Abstract::Vector* derivResidualParamPtr;

	//! Stores the derivative of the null vector residual w.r.t. bifparam
	NOX::Abstract::Vector* derivNullResidualParamPtr;

	//! Flag indicating whether we have our own copy of the underlying grp
	bool ownsGroup;

	//! Is residual vector valid
	bool isValidF;

	//! Is Jacobian matrix valid
	bool isValidJacobian;

	//! Is Newton vector valid
	bool isValidNewton;
      };
    } // namespace TPBord
  } // namespace Bifurcation
} // namespace LOCA

#endif
