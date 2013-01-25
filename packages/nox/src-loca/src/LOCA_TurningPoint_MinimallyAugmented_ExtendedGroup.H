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

#ifndef LOCA_TURNINGPOINT_MINIMALLYAUGMENTED_EXTENDEDGROUP_H
#define LOCA_TURNINGPOINT_MINIMALLYAUGMENTED_EXTENDEDGROUP_H

#include "Teuchos_RCP.hpp"

#include "LOCA_MultiContinuation_AbstractGroup.H"         // base class
#include "LOCA_Extended_MultiAbstractGroup.H"             // base class
#include "LOCA_BorderedSystem_AbstractGroup.H"          // base class
#include "LOCA_Abstract_TransposeSolveGroup.H"          // base class

namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
  namespace TurningPoint {
    namespace MinimallyAugmented {
      class AbstractGroup;
      class Constraint;
    }
  }
  namespace MultiContinuation {
    class ConstrainedGroup;
  }
}
namespace LOCA { 

  namespace TurningPoint { 
    
    namespace MinimallyAugmented {

      /*! 
       * \brief A group representing the minimally augemented turning point 
       * equations. 
       */
      /*!
       * The LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup is a 
       * concrete implementation of the NOX::Abstract::Group, 
       * LOCA::MultiContinuation::AbstractGroup and 
       * LOCA::Extended::MultiAbstractGroup that defines the following 
       * extended set of equations that are regular at a generic turning point:
       * \f[
       *     G(z) = \left[
       *      \begin{array}{c}
       *         F(x,p) \\
       *         \sigma
       *      \end{array}
       *      \right] = 0
       *   \f]
       * where \f$z = [x, p]\in\Re^{n+1}\f$, \f$x\f$ is the solution vector,
       * \f$p\f$ is the bifurcation parameter, \f$J\f$ is the Jacobian of 
       * \f$F\f$, and \f$\sigma\in\Re\f$ is a measure of the singularity of 
       * \f$F\f$ and is defined via
       * \f[
       *     \begin{bmatrix}
       *       J   & a \\
       *       b^T & 0 
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *       v \\
       *       \sigma_1
       *     \end{bmatrix} = 
       *     \begin{bmatrix}
       *       0 \\
       *       n
       *     \end{bmatrix},
       * \f]
       * \f[
       *     \begin{bmatrix}
       *       J^T & b \\
       *       a^T & 0 
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *       w \\
       *       \sigma_2
       *     \end{bmatrix} = 
       *     \begin{bmatrix}
       *       0 \\
       *       n
       *     \end{bmatrix},
       * \f]
       * \f[
       *     \sigma = w^T J v/n
       * \f]
       * for any vectors \f$a\f$ and \f$b\f$ in \f$\Re^n\f$.  Using these
       * relationships, it is easy to show
       * \f[
       *     \begin{split}
       *        \sigma_x &= (w^T J v)_x/n = w^T J_x v/n \\
       *        \sigma_p &= (w^T J v)_p/n = w^T J_p v/n
       *     \end{split}
       * \f]
       *
       * The group stores an underlying group of type 
       * LOCA::TurningPoint::MinimallyAugmented::AbstractGroup to represent 
       * the equations \f$F(x,p) = 0\f$ and to manipulate the underlying 
       * Jacobian \f$J\f$.  This interface defines methods for computing
       * the derivatives \f$(w^T J v)_x\f$ and \f$(w^T J v)_p\f$ as well.
       *
       * This class implements all of the NOX::Abstract::Group, 
       * LOCA::MultiContinuation::AbstractGroup, and 
       * LOCA::Extended::MultiAbstractGroup 
       * methods for this extended set of equations and therefore is a complete
       * group which can be passed to most %NOX solvers to locate a single 
       * turning point or to the LOCA::Stepper to compute a family of 
       * turning points in a second parameter.
       *
       * The class is intialized via the \c tpParams parameter list argument 
       * to the constructor.  The parameters this class recognizes are:
       * <ul>
       * <li> "Bifurcation Parameter" -- [string] (Must be supplied) - 
       *    Name of the bifurcation parameter \f$p\f$
       * <li> "Bordered Solver Method" -- [string] (default "Bordering") 
       *      Method for solving bordered systems of equations.  See 
       *      LOCA::BorderedSolver::Factory for a description.
       * <li> "Symmetric Jacobian"  -- [bool] (default: false) - 
       *    Flag indicating whether Jacobian matrix \f$J\f$ is symmetric,
       *    in which case we force \f$a = b\f$ and therefore the second
       *    tranpose solve for \f$w\f$ is unnecessary
       * <li> "Constraint Method" -- [string] (default: "Default") -
       *    Type of constraint method to use.  Valid choices are
       *    <ul>
       *    <li> "Default" (LOCA::TurningPoint::MinimallyAugmented::Constraint)
       *         Default method described above.
       *    <li> "Modified" 
       *         (LOCA::TurningPoint::MinimallyAugmented::ModifedConstraint)
       *          A modified method that computes updates to the null vectors
       *          every nonlinear interation, instead of directly solving
       *          for them
       *    </ul>
       * <li> "Initial Null Vector Compuation" -- [string] 
       *      (default: "User Provided") - Method to compute initial \f$a\f$ 
       *      and \f$b\f$ vectors.  Valid choices are:
       *      <ul>
       *      <li> "User Provided" - Initial vectors are provided in the
       *           parameter list, in which case the following parameters
       *           are relevant:
       *           <ul>
       *           <li> "Initial A Vector"  -- 
       *                [Teuchos::RCP<NOX::Abstract::Vector>] 
       *                (Must be supplied) - 
       *                 Vector storing initial value for \f$a\f$ vector
       *           <li> "Initial B Vector"  -- 
       *                 [Teuchos::RCP<NOX::Abstract::Vector>] 
       *                 (Must be supplied for nonsymmetric Jacobians) - 
       *                 Vector storing initial value for \f$b\f$ vector
       *           </ul>
       *      <li> "Solve df/dp" - Compute \f$a = J^{-T}df/dp\f$ and 
       *           \f$b = J^{-1} df/dp\f$ where \f$p\f$ is the bifurcation
       *           parameter.
       *      <li> "Constant" - Entries of $a$ and $b$ are set to 1.0
       *      </ul>
       * <li> "Null Vector Scaling Method" -- [string] (default:  "Order N") - 
       *      Method to scale \f$a\f$ and \f$b\f$.  This determines the norm of
       *      these vectors and the scaling of \f$\sigma\f$.  Valid choices are:
       *      <ul>
       *      <li> "None" -- Use initial scaling
       *      <li> "Order 1" -- Scale to unit norm
       *      <li> "Order N" -- Use vector length scaling
       *      </ul>
       * <li> "Update Null Vectors Every Continuation Step"  -- 
       *    [bool] (default: true) - 
       *    Flag indicating whether to update \f$a\f$ and \f$b\f$ vectors via
       *    \f$a = w\f$ and \f$b = v\f$ every continuation step
       * <li> "Update Null Vectors Every Nonlinear Iteration"  -- 
       *    [bool] (default: false) - 
       *    Flag indicating whether to update \f$a\f$ and \f$b\f$ vectors via
       *    \f$a = w\f$ and \f$b = v\f$ every nonlinear iteration
       * <li> "Multiply Null Vectors by Mass Matrix"  -- 
       *    [bool] (default: false) - 
       *    Flag indicating whether to multiply \f$a\f$ and \f$b\f$ vectors by
       *    the mass matrix \f$M = \partial f/\partial\dot{x}\f$ at the strart
       *    of a turning point calculation, and each time \f$a\f$ and \f$b\f$
       *    are updated.  This can improve the scaling of these vectors, and
       *    may orthogonalize them against structural null spaces (i.e., 
       *    pressure null space for incompressible Navier-Stokes).
       * </ul>
       */
      class ExtendedGroup 
	: public virtual LOCA::Extended::MultiAbstractGroup,
	  public virtual LOCA::MultiContinuation::AbstractGroup,
	  public virtual LOCA::BorderedSystem::AbstractGroup,
	  public virtual LOCA::Abstract::TransposeSolveGroup {

      public:

	/*! 
	 * \brief Constructor with initial data passed through parameter lists.
	 */
	ExtendedGroup(
	 const Teuchos::RCP<LOCA::GlobalData>& global_data,
         const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& tpParams,
	 const Teuchos::RCP<LOCA::TurningPoint::MinimallyAugmented::AbstractGroup>& g);

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
	 * This method fills the extended residual
	 * \f[
	 *     G(z) = \left[
	 *     \begin{array}{c}
	 *        F(x,p) \\
	 *       \sigma
	 *    \end{array}
	 *     \right].
	 * \f]
	 */
	virtual NOX::Abstract::Group::ReturnType computeF();

	//! Compute the blocks of the Jacobian derivative of \f$G\f$
	/*!
	 * This method computes the \f$J\f$, \f$F_p\f$, \f$\sigma_x\f$ and 
	 * \f$\sigma_p\f$ blocks of the extended Jacobian:
	 * \f[
	 *    D_z G(z) = 
	 *    \begin{bmatrix}
	 *       J        & F_p \\
	 *       \sigma_x & \sigma_p
	 *    \end{bmatrix}
	 * \f]
	 */
	virtual NOX::Abstract::Group::ReturnType computeJacobian();

	//! Gradient computation is not defined for this group
	virtual NOX::Abstract::Group::ReturnType computeGradient();

	//! Compute %Newton direction using applyJacobianInverse().
	virtual NOX::Abstract::Group::ReturnType 
	computeNewton(Teuchos::ParameterList& params);
  
	//! Computes the extended Jacobian vector product
	virtual NOX::Abstract::Group::ReturnType 
	applyJacobian(const NOX::Abstract::Vector& input, 
		      NOX::Abstract::Vector& result) const;

	//! Computes the extended Jacobian transpose vector product
	virtual NOX::Abstract::Group::ReturnType 
	applyJacobianTranspose(const NOX::Abstract::Vector& input, 
			       NOX::Abstract::Vector& result) const;

	//! Applies the inverse of the extended Jacobian matrix 
	virtual NOX::Abstract::Group::ReturnType 
	applyJacobianInverse(Teuchos::ParameterList& params, 
			     const NOX::Abstract::Vector& input, 
			     NOX::Abstract::Vector& result) const;

	//! Applies Jacobian for extended system
	virtual NOX::Abstract::Group::ReturnType 
	applyJacobianMultiVector(const NOX::Abstract::MultiVector& input, 
				 NOX::Abstract::MultiVector& result) const;

	//! Jacobian transpose for extended system
	virtual NOX::Abstract::Group::ReturnType 
	applyJacobianTransposeMultiVector(
				     const NOX::Abstract::MultiVector& input, 
				     NOX::Abstract::MultiVector& result) const;

	//! Applies Jacobian inverse for extended system
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

	//! Return RCP to extended solution vector \f$z\f$.  
	virtual Teuchos::RCP< const NOX::Abstract::Vector > getXPtr() const;

	//! Return RCP to extended equation residual \f$G(z)\f$
	virtual Teuchos::RCP< const NOX::Abstract::Vector > getFPtr() const;

	//! Vector returned is not valid.
	virtual Teuchos::RCP< const NOX::Abstract::Vector > getGradientPtr() const;

	//! Return RCP to extended Newton direction.
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
	virtual void setParam(std::string paramID, double val);

	//! Set parameter indexed by paramID
	virtual void setParam(int paramID, double val);

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
	 * @name Implementation of 
	 * LOCA::BorderedSystem::AbstractGroup virtual methods 
	 */
	//@{
	
	//! Return the total width of the bordered rows/columns
	virtual int getBorderedWidth() const;
	
	//! Get bottom-level unbordered group
	virtual Teuchos::RCP<const NOX::Abstract::Group>
	getUnborderedGroup() const;
	
	//! Indicates whether combined A block is zero
	virtual bool isCombinedAZero() const;
	
	//! Indicates whether combined B block is zero
	virtual bool isCombinedBZero() const;
	
	//! Indicates whether combined C block is zero
	virtual bool isCombinedCZero() const;
	
	/*!
	 * Given the vector \c v, extract the underlying solution component
	 * corresponding to the unbordered group.
	 */
	virtual void
	extractSolutionComponent(const NOX::Abstract::MultiVector& v,
				 NOX::Abstract::MultiVector& v_x) const;

	/*!
	 * Given the vector \c v, extract the parameter components of all
	 * of the nested subvectors in \c v down to the solution component
	 * for the unbordered group.
	 */
	virtual void 
	extractParameterComponent(
			   bool use_transpose,
                           const NOX::Abstract::MultiVector& v,
                           NOX::Abstract::MultiVector::DenseMatrix& v_p) const;

	/*!
	 * Given the solution component \c v_x and combined parameter
	 * components \c v_p, distribute these components through the nested
	 * sub-vectors in \c v.
	 */
	virtual void
	loadNestedComponents(
			  const NOX::Abstract::MultiVector& v_x,
			  const NOX::Abstract::MultiVector::DenseMatrix& v_p,
			  NOX::Abstract::MultiVector& v) const;

	//! Fill the combined A block as described above
	virtual void fillA(NOX::Abstract::MultiVector& A) const;
	
	//! Fill the combined B block as described above
	virtual void fillB(NOX::Abstract::MultiVector& B) const;
	
	//! Fill the combined C block as described above
	virtual void fillC(NOX::Abstract::MultiVector::DenseMatrix& C) const;
	
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

	//! Returns left null vector
	Teuchos::RCP<const NOX::Abstract::Vector>
	getLeftNullVec() const;

	//! Returns right null vector v
	Teuchos::RCP<const NOX::Abstract::Vector>
	getRightNullVec() const;

	//! Returns "A" vector
	Teuchos::RCP<const NOX::Abstract::Vector>
	getAVec() const;

	//! Returns "B" 
	Teuchos::RCP<const NOX::Abstract::Vector>
	getBVec() const;

      protected:

	//! Set bifurcation parameter
	void setBifParam(double param);

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
	Teuchos::RCP<LOCA::TurningPoint::MinimallyAugmented::AbstractGroup> grpPtr;

	//! Pointer to the constraint equation
	Teuchos::RCP<LOCA::TurningPoint::MinimallyAugmented::Constraint> constraint;

	//! Pointer to constrained group implementation
	Teuchos::RCP<LOCA::MultiContinuation::ConstrainedGroup> conGroup;

	//! Stores the bifurcation parameter index
	int bifParamID;

      }; // class ExtendedGroup

    } // namespace MinimallyAugmented

  } // namespace TurningPoint

} // namespace LOCA

#endif
