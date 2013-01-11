// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef LOCA_TURNINGPOINT_MOORESPENCE_EXTENDEDGROUP_H
#define LOCA_TURNINGPOINT_MOORESPENCE_EXTENDEDGROUP_H

#include "Teuchos_RCP.hpp"

#include "LOCA_MultiContinuation_AbstractGroup.H"         // base class
#include "LOCA_Extended_MultiAbstractGroup.H"             // base class
#include "LOCA_Abstract_TransposeSolveGroup.H"            // base class
#include "LOCA_TurningPoint_MooreSpence_ExtendedVector.H" // class data element
#include "LOCA_TurningPoint_MooreSpence_ExtendedMultiVector.H" // class data element

namespace LOCA {
  class ParameterVector;
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }

  namespace TurningPoint {
    namespace MooreSpence {
      class AbstractGroup;
      class SolverStrategy;
    }
  }

  namespace TimeDependent {
    class AbstractGroup;
  }
}
namespace LOCA { 

  namespace TurningPoint { 
    
    namespace MooreSpence {

      /*! 
       * \brief A group representing the Moore-Spence turning point equations. 
       */
      /*!
       The LOCA::TurningPoint::MooreSpence::ExtendedGroup is a concrete 
       implementation of 
       the NOX::Abstract::Group, LOCA::MultiContinuation::AbstractGroup and 
       LOCA::Extended::MultiAbstractGroup that defines the following extended 
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
       LOCA::TurningPoint::MooreSpence AbstractGroup to represent the equations
       \f$F(x,p) = 0\f$ and to manipulate the underlying Jacobian \f$J\f$.
       Note that the entire extended Jacobian \f$D_z G\f$ is not stored in
       memory, rather a block-elimination algorithm (bordering algorithm) is
       used to compute linear solves of the extended Jacobian (see 
       LOCA::TurningPoint::MooreSpence::SolverFactory for more details).

       This class implements all of the NOX::Abstract::Group, 
       LOCA::MultiContinuation::AbstractGroup, and 
       LOCA::Extended::MultiAbstractGroup 
       methods for this extended set of equations and therefore is a complete
       group which can be passed to most %NOX solvers to locate a single 
       turning point or to the LOCA::Stepper to compute a family of turning 
       points in a second parameter.

       However, Jacobian-tranpose operations and gradient calculations cannot 
       be implemented efficiently and therefore gradient-base nonlinear solvers
       such as steepest descent and Trust region methods cannot be used to 
       solve the extended turning point equations.

       The class is intialized via the \c tpParams parameter list argument to 
       the constructor.  The parameters this class recognizes are:
       <ul>
       <li> "Bifurcation Parameter" -- [string] (Must be supplied) Name of the bifurcation parameter \f$p\f$
       <li> "Length Normalization Vector" -- 
            [Teuchos::RCP<NOX::Abstract::Vector>] (Must be supplied) 
            Vector storing length normalization vector \f$l\f$
       <li> "Initial Null Vector"  -- 
            [Teuchos::RCP<NOX::Abstract::Vector>] (Must be supplied) - 
            Vector storing initial guess for the null vector \f$n\f$
       <li> "Perturb Initial Solution" -- [bool] (default: false) Flag 
            indicating whether to perturb the initial solution
       <li> "Relative Perturbation Size" -- [double] (default: 1.0e-3) 
            Relative perturbation size if perturbing the initial solution
       <li> "Null Vector Scaling" -- [string] (default:  "Order N") - 
             Method to scale \f$l\f$.  This determines the norm of \f$l\f$
	     (and hence, \f$n\f$).  Valid choices are:
             <ul>
             <li> "None" -- Use initial scaling
             <li> "Order 1" -- Scale to unit norm
             <li> "Order N" -- Use vector length scaling
             </ul>
       <li> "Update Null Vectors Every Continuation Step"  -- 
           [bool] (default: false) - 
           Flag indicating whether to update length normalization vector 
	   \f$l\f$ via \f$l = n\f$ every continuation step
       <li> "Multiply Null Vectors by Mass Matrix"  -- 
           [bool] (default: false) - 
           Flag indicating whether to multiply length scaling vector \f$l\f$ by
           the mass matrix \f$M = \partial f/\partial\dot{x}\f$ at the start
           of a turning point calculation, and each time it is
           updated.  This can improve the scaling of this vector, and
           may orthogonalize it against structural null spaces (i.e., 
           pressure null space for incompressible Navier-Stokes).
       </ul>
       */

      class ExtendedGroup 
	: public virtual LOCA::Extended::MultiAbstractGroup,
	  public virtual LOCA::MultiContinuation::AbstractGroup,
	  public virtual LOCA::Abstract::TransposeSolveGroup {

      public:

	/*! 
	 * \brief Constructor with initial data passed through parameter lists.
	 */
	ExtendedGroup(
	 const Teuchos::RCP<LOCA::GlobalData>& global_data,
         const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& tpParams,
	 const Teuchos::RCP<LOCA::TurningPoint::MooreSpence::AbstractGroup>& g);

	//! Copy constructor.
	ExtendedGroup(const ExtendedGroup& source, 
		      NOX::CopyType type = NOX::DeepCopy);

	//! Destructor.
	virtual ~ExtendedGroup();

	/*! 
	 * @name Implementation of NOX::Abstract::Group virtual methods 
	 */
	//@{

	//! Assignment operator
	virtual NOX::Abstract::Group& 
	operator=(const NOX::Abstract::Group& source);
     
	//! Cloning function
	virtual Teuchos::RCP<NOX::Abstract::Group>
	clone(NOX::CopyType type = NOX::DeepCopy) const;

	//! Set the solution vector, x, to y.
	virtual void setX(const NOX::Abstract::Vector& y);

	//! Compute this.x = grp.x + step * d.
	virtual void computeX(const NOX::Abstract::Group& g, 
			      const NOX::Abstract::Vector& d,
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
	computeNewton(Teuchos::ParameterList& params);
  
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
	  LOCA::TurningPoint::MooreSpence::ExtendedVector, otherwise an error 
	  is thrown.
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
	 * This method is a special case of applyJacobianInverseMultiVector() 
	 * for a single right-hand-side.
	 */
	virtual NOX::Abstract::Group::ReturnType 
	applyJacobianInverse(Teuchos::ParameterList& params, 
			     const NOX::Abstract::Vector& input, 
			     NOX::Abstract::Vector& result) const;

	//! Applies Jacobian for extended system
	virtual NOX::Abstract::Group::ReturnType 
	applyJacobianMultiVector(const NOX::Abstract::MultiVector& input, 
				 NOX::Abstract::MultiVector& result) const;

	//! Jacobian transpose not defined for this system
	virtual NOX::Abstract::Group::ReturnType 
	applyJacobianTransposeMultiVector(
				     const NOX::Abstract::MultiVector& input, 
				     NOX::Abstract::MultiVector& result) const;

	//! Applies Jacobian inverse for extended system
	/*!
	 * Uses a LOCA::TurningPoint::MooreSpence::SolverStrategy instantiated
	 * by the LOCA::TurningPoint::MooreSpence::SolverFactory to 
	 * implement the solve.
	 */
	virtual NOX::Abstract::Group::ReturnType 
	applyJacobianInverseMultiVector(
				    Teuchos::ParameterList& params, 
				    const NOX::Abstract::MultiVector& input, 
				    NOX::Abstract::MultiVector& result) const;
  
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

	//! Return extended solution vector \f$z\f$.  
	virtual Teuchos::RCP< const NOX::Abstract::Vector > getXPtr() const;

	//! Return extended equation residual \f$G(z)\f$
	virtual Teuchos::RCP< const NOX::Abstract::Vector > getFPtr() const;

	//! Vector returned is not valid.
	virtual Teuchos::RCP< const NOX::Abstract::Vector > getGradientPtr() const;

	//! Return extended Newton direction.
	virtual Teuchos::RCP< const NOX::Abstract::Vector > getNewtonPtr() const;

	//! Return the norm of the %Newton solve residual.
	virtual double getNormNewtonSolveResidual() const;

	//@}

	/*! 
	 * @name Implementation of LOCA::Extended::MultiAbstractGroup 
	 * virtual methods 
	 */
	//@{
	
	//! Return underlying group
	virtual 
	Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup> 
	getUnderlyingGroup() const;
	
	//! Return underlying group
	virtual 
	Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
	getUnderlyingGroup();

	//@}

	/*! 
	 * @name Implementation of LOCA::MultiContinuation::AbstractGroup
	 * virtual methods 
	 */
	//@{

	//! Assignment operator
	virtual void copy(const NOX::Abstract::Group& source);

	//! Set parameters indexed by (integer) paramIDs
	virtual void setParamsMulti(
			  const std::vector<int>& paramIDs, 
			  const NOX::Abstract::MultiVector::DenseMatrix& vals);

	//! Set the parameter vector in the group to p.
	virtual void setParams(const ParameterVector& p);

	//! Set parameter indexed by paramID
	virtual void setParam(int paramID, double val);

	//! Set parameter indexed by paramID
	virtual void setParam(std::string paramID, double val);

	//! Return a const reference to the paramter vector owned by the group.
	virtual const ParameterVector& getParams() const;

	//! Return copy of parameter indexed by paramID
	virtual double getParam(int paramID) const;

	//! Return copy of parameter indexed by paramID
	virtual double getParam(std::string paramID) const;

	/*! 
	 * Compute \f$\partial F/\partial p\f$ for each parameter \f$p\f$ 
	 * indexed by paramIDs.  The first column of \em dfdp holds F,
	 * which is valid if \em isValidF is true.  Otherwise F must be
	 * computed.
	 */
	virtual NOX::Abstract::Group::ReturnType
	computeDfDpMulti(const std::vector<int>& paramIDs, 
			 NOX::Abstract::MultiVector& dfdp, 
			 bool isValidF);

	//! Perform any preprocessing before a continuation step starts.
	/*!
	 * The \c stepStatus argument indicates whether the previous step was
	 * successful.
	 */
	virtual void 
	preProcessContinuationStep(
			 LOCA::Abstract::Iterator::StepStatus stepStatus);
	
	//! Perform any postprocessing after a continuation step finishes.
	/*!
	 * The \c stepStatus argument indicates whether the step was
	 * successful.
	 */
	virtual void 
	postProcessContinuationStep(
			 LOCA::Abstract::Iterator::StepStatus stepStatus);

	//! Projects solution to a few scalars for multiparameter continuation
	virtual void projectToDraw(const NOX::Abstract::Vector& x,
				   double *px) const;

	//! Returns the dimension of the project to draw array
	virtual int projectToDrawDimension() const;

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

	/*! 
	 * \brief Function to print out extended solution and continuation 
	 * parameter after successful continuation step.
	 */
	/*!
	 * This method prints the solution, null-vector, and parameter 
	 * components of the extended solution vector using the printSolution
	 * method of the underlying group.
	 */
	 virtual void printSolution(const NOX::Abstract::Vector& x_,
				    const double conParam) const;

	//@}

	/*! 
	 * @name Implementation of LOCA::Abstract::TransposeSolveGroup
	 * virtual methods 
	 */
	//@{

	//! Solve Jacobian-tranpose system
	virtual NOX::Abstract::Group::ReturnType 
	applyJacobianTransposeInverse(Teuchos::ParameterList& params, 
				      const NOX::Abstract::Vector& input, 
				      NOX::Abstract::Vector& result) const;

	//! Solve Jacobian-tranpose system with multiple right-hand sides
	virtual NOX::Abstract::Group::ReturnType 
	applyJacobianTransposeInverseMultiVector(
				   Teuchos::ParameterList& params, 
				   const NOX::Abstract::MultiVector& input, 
				   NOX::Abstract::MultiVector& result) const;

	//@}

	//! Get bifurcation parameter
	double getBifParam() const;

	//! Defines null vector normalization \f$l^Tn\f$ equation
	double lTransNorm(const NOX::Abstract::Vector& n) const;

	//! null vector normalization for multivectors
	/*!
	 * Note:  result should have 1 row and n.numVectors() columns.
	 */
	void lTransNorm(const NOX::Abstract::MultiVector& n,
			NOX::Abstract::MultiVector::DenseMatrix& result) const;

	//! Return length normalization vector
	Teuchos::RCP<NOX::Abstract::Vector>
	getLengthVector() const;

      protected:

	//! Enumerated type determining type of scaling
	enum NullVectorScaling {
	  NVS_None,      //! No scaling
	  NVS_OrderOne,  //! Scale to O(1)
	  NVS_OrderN     //! Scale to O(N) when N is the vector length
	};

	//! Set bifurcation parameter
	void setBifParam(double param);

	//! Sets up multivector views
	void setupViews();

	//! Initializes group
	void init(bool perturbSoln = false, double perturbSize = 0.0);

	//! Scale null vector
	void scaleNullVector(NOX::Abstract::Vector& a);

      private:

	//! Prohibit generation and use of operator=()
	ExtendedGroup& operator=(const ExtendedGroup&);

      protected:

	//! Pointer LOCA global data object
	Teuchos::RCP<LOCA::GlobalData> globalData;

	//! Parsed top-level parameters
	Teuchos::RCP<LOCA::Parameter::SublistParser> parsedParams;

	//! Bifurcation parameter list
	Teuchos::RCP<Teuchos::ParameterList> turningPointParams;

	//! Pointer to base group that defines \f$F\f$
	Teuchos::RCP<LOCA::TurningPoint::MooreSpence::AbstractGroup> grpPtr;

	//! Stores the extended solution vector
	LOCA::TurningPoint::MooreSpence::ExtendedMultiVector xMultiVec;
	
	//! Stores the extended residual vector and df/dp
	LOCA::TurningPoint::MooreSpence::ExtendedMultiVector fMultiVec;
	
	//! Stores the extended Newton vector
	LOCA::TurningPoint::MooreSpence::ExtendedMultiVector newtonMultiVec;

	//! Stores the length normalization vector
	Teuchos::RCP<NOX::Abstract::MultiVector> lengthMultiVec;
	
	//! Stores view of first column of xMultiVec
	Teuchos::RCP<LOCA::TurningPoint::MooreSpence::ExtendedVector> xVec;
	
	//! Stores view of first column of fMultiVec
	Teuchos::RCP<LOCA::TurningPoint::MooreSpence::ExtendedVector> fVec;
	
	//! Stores view of first column of fMultiVec as a multivec
	Teuchos::RCP<LOCA::TurningPoint::MooreSpence::ExtendedMultiVector> ffMultiVec;

	//! Stores view of df/dp columns of fMultiVec
	Teuchos::RCP<LOCA::TurningPoint::MooreSpence::ExtendedMultiVector> dfdpMultiVec;

	//! Stores view of first column of newtonMultiVec
	Teuchos::RCP<LOCA::TurningPoint::MooreSpence::ExtendedVector> newtonVec;

	//! Stores view of first column of lengthMultiVec
	Teuchos::RCP<NOX::Abstract::Vector> lengthVec;

	//! Stores bordering strategy
	Teuchos::RCP<LOCA::TurningPoint::MooreSpence::SolverStrategy> solverStrategy;

	//! Stores indices for getting f part of fMultiVec
	std::vector<int> index_f;

	//! Stores indices for getting df/dp part of fMultiVec
	std::vector<int> index_dfdp;

	//! Stores the bifurcation parameter index
	std::vector<int> bifParamID;

	//! Is residual vector valid
	bool isValidF;

	//! Is Jacobian matrix valid
	bool isValidJacobian;

	//! Is Newton vector valid
	bool isValidNewton;

	/*! 
	 * \brief Flag indicating whether to update \f$n\f$
	 * every continuation step.
	 */
	bool updateVectorsEveryContinuationStep;

	//! Null vector scaling method
	NullVectorScaling nullVecScaling;

	//! Multiply null vectors by mass matrix
	bool multiplyMass;

	//! Time dependent group interface for multiplying by mass matrix
	Teuchos::RCP<LOCA::TimeDependent::AbstractGroup> tdGrp;

	//! Temporary vector for mulplying null vectors by mass matrix
	Teuchos::RCP<NOX::Abstract::Vector> tmp_mass;

      }; // class ExtendedGroup

    } // namespace MooreSpence

  } // namespace TurningPoint

} // namespace LOCA

#endif
