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

#ifndef LOCA_HOPF_MINIMALLYAUGMENTED_EXTENDEDGROUP_H
#define LOCA_HOPF_MINIMALLYAUGMENTED_EXTENDEDGROUP_H

#include "Teuchos_RCP.hpp"

#include "LOCA_Extended_MultiAbstractGroup.H"           // base class
#include "LOCA_MultiContinuation_AbstractGroup.H"       // base class
#include "LOCA_BorderedSystem_AbstractGroup.H"          // base class
#include "LOCA_MultiContinuation_ExtendedVector.H"      // class data element
#include "LOCA_MultiContinuation_ExtendedMultiVector.H" // class data element

// forward declarations
namespace Teuchos {
  class ParameterList;
}
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
  namespace Hopf {
    namespace MinimallyAugmented {
      class AbstractGroup;
      class Constraint;
    }
  }
  namespace BorderedSolver {
    class AbstractStrategy;
    class JacobianOperator;
  }
}

namespace LOCA {

  namespace Hopf {
    
    namespace MinimallyAugmented {

      /*! 
       * \brief A group representing the minimally augemented Hopf
       * equations. 
       */
      /*!
       * The LOCA::Hopf::MinimallyAugmented::ExtendedGroup is a 
       * concrete implementation of the NOX::Abstract::Group, 
       * LOCA::MultiContinuation::AbstractGroup and 
       * LOCA::Extended::MultiAbstractGroup that defines the following 
       * extended set of equations that are regular at a generic Hopf:
       * \f[
       *     G(z) = \left[
       *      \begin{array}{c}
       *         F(x,p)\\
       *         \sigma 
       *      \end{array}
       *      \right] = 0
       *   \f]
       * where \f$z = [x, p, \omega]\in\Re^{n+2}\f$, \f$x\f$ is the 
       * solution vector, \f$p\f$ is the bifurcation parameter, \f$\omega\f$
       * is the Hopf frequency and \f$\sigma\in\Re\f$ is a measure of the 
       * singularity of \f$J+i\omega M\f$ and is defined via
       * \f[
       *     \begin{bmatrix}
       *       J+i\omega M & a \\
       *       b^H         & 0 
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
       *       J^T-i\omega M^T & b \\
       *       a^H             & 0 
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
       *     \sigma = w^H J+i\omega M v/n
       * \f]
       * for any vectors \f$a\f$ and \f$b\f$ in \f$\mathbb{C}^n\f$.  Using 
       * these relationships, it is easy to show
       * \f[
       *     \begin{split}
       *        \sigma_x &= (w^H(J+i\omega M)v)_x/n = w^H(J+i\omega M)_x v/n \\
       *        \sigma_p &= (w^H(J+i\omega M)v)_p/n = w^H(J+i\omega M)_p v/n
       *     \end{split}
       * \f]
       *
       * The group stores an underlying group of type 
       * LOCA::Hopf::MinimallyAugmented::AbstractGroup to represent 
       * the equations \f$F(x,p) = 0\f$ and to manipulate the underlying 
       * complex matrix \f$C = J+i\omega M\f$.  This interface defines methods 
       * for computing the derivatives \f$(w^H C v)_x\f$ and \f$(w^H C v)_p\f$.
       * Since %LOCA is not able to deal with complex vectors and matrices
       * directly, the real-equivalent formulation is used for all complex
       * calculations.
       *
       * This class implements all of the NOX::Abstract::Group, 
       * LOCA::MultiContinuation::AbstractGroup, and 
       * LOCA::Extended::MultiAbstractGroup 
       * methods for this extended set of equations and therefore is a complete
       * group which can be passed to most %NOX solvers to locate a single 
       * pitchfork or to the LOCA::Stepper to compute a family of 
       * pitchforks in a second parameter.
       *
       * The class is intialized via the \c hpfParams parameter list argument 
       * to the constructor.  The parameters this class recognizes are:
       * <ul>
       * <li> "Bifurcation Parameter" -- [string] (Must be supplied) - 
       *    Name of the bifurcation parameter \f$p\f$
       * <li> "Initial Frequency" -- [double] (Must be supplied)
       *      Initial guess for the Hopf frequency \f$\omega\f$.
       * <li> "Bordered Solver Method" -- [string] (default "Bordering") 
       *      Method for solving bordered systems of equations.  See 
       *      LOCA::BorderedSolver::Factory for a description.
       * <li> "Symmetric Jacobian"  -- [bool] (default: false) - 
       *    Flag indicating whether Jacobian matrix \f$J\f$ is symmetric,
       *    in which case we force \f$a = b\f$ and therefore the second
       *    tranpose solve for \f$w\f$ is unnecessary
       * <li> "Initial Null Vector Compuation" -- [string] 
       *      (default: "User Provided") - Method to compute initial \f$a\f$ 
       *      and \f$b\f$ vectors.  Valid choices are:
       *      <ul>
       *      <li> "User Provided" - Initial vectors are provided in the
       *           parameter list, in which case the following parameters
       *           are relevant:
       *           <ul>
       *           <li> "Initial Real A Vector"  -- 
       *                [Teuchos::RCP<NOX::Abstract::Vector>] 
       *                (Must be supplied) - 
       *                 Vector storing initial value for the real component
       *                 of the \f$a\f$ vector
       *           <li> "Initial Imaginary A Vector"  -- 
       *                [Teuchos::RCP<NOX::Abstract::Vector>] 
       *                (Must be supplied) - 
       *                 Vector storing initial value for the imaginary 
       *                 component of the \f$a\f$ vector
       *           <li> "Initial Real B Vector"  -- 
       *                 [Teuchos::RCP<NOX::Abstract::Vector>] 
       *                 (Must be supplied for nonsymmetric Jacobians) - 
       *                 Vector storing initial value for the real component
       *                 of the \f$b\f$ vector
       *           <li> "Initial Imaginary B Vector"  -- 
       *                [Teuchos::RCP<NOX::Abstract::Vector>] 
       *                (Must be supplied) - 
       *                 Vector storing initial value for the imaginary 
       *                 component of the \f$b\f$ vector
       *           </ul>
       *      <li> "Solve df/dp" - Compute \f$a = J^{-T}df/dp\f$ and 
       *           \f$b = J^{-1} df/dp\f$ where \f$p\f$ is the bifurcation
       *           parameter.
       *      </ul>
       * <li> "Update Null Vectors Every Continuation Step"  -- 
       *    [bool] (default: true) - 
       *    Flag indicating whether to update \f$a\f$ and \f$b\f$ vectors via
       *    \f$a = w\f$ and \f$b = v\f$ every continuation step
       * <li> "Update Null Vectors Every Nonlinear Iteration"  -- 
       *    [bool] (default: false) - 
       *    Flag indicating whether to update \f$a\f$ and \f$b\f$ vectors via
       *    \f$a = w\f$ and \f$b = v\f$ every nonlinear iteration
       * </ul>
       */
      class ExtendedGroup : 
	public virtual LOCA::Extended::MultiAbstractGroup,
	public virtual LOCA::MultiContinuation::AbstractGroup,
	public virtual LOCA::BorderedSystem::AbstractGroup {

      public:

	//! Constructor 
	/*!
	 * \param global_data [in] Global data object
	 * \param topParams [in] Parsed top-level parameter list.
	 * \param hpfParams [in] Parameter list determining the
	 * bordered solver method.
	 * \param grp [in] Group representing \f$f\f$.
	 */
      ExtendedGroup(
       const Teuchos::RCP<LOCA::GlobalData>& global_data,
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& hpfParams,
       const Teuchos::RCP<LOCA::Hopf::MinimallyAugmented::AbstractGroup>& grp);

	//! Copy constructor
	ExtendedGroup(const ExtendedGroup& source, 
		      NOX::CopyType type = NOX::DeepCopy);

	//! Destructor.
	virtual ~ExtendedGroup();

	//! Get bifurcation parameter
	double getBifParam() const;

	/*! 
	 * @name Implementation of NOX::Abstract::Group virtual methods 
	 */
	//@{

	//! Assignment operator
	virtual NOX::Abstract::Group& 
	operator=(const NOX::Abstract::Group& source);

	//! Clone function
	virtual Teuchos::RCP<NOX::Abstract::Group>
	clone(NOX::CopyType type = NOX::DeepCopy) const;

	//! Set the solution vector to y.
	virtual void setX(const NOX::Abstract::Vector& y);
	
	/*! 
	 * \brief Compute and return solution vector, x, where 
	 * this.x = grp.x + step * d.
	 */
	virtual void computeX(const NOX::Abstract::Group& g, 
			      const NOX::Abstract::Vector& d,
			      double step);

	//! Compute extended continuation equations 
	virtual NOX::Abstract::Group::ReturnType computeF();
	
	//! Compute extended continuation jacobian
	virtual NOX::Abstract::Group::ReturnType computeJacobian();

	//! Gradient is not defined for this system
	virtual NOX::Abstract::Group::ReturnType computeGradient();

	//! Compute Newton direction for extended continuation system
	virtual NOX::Abstract::Group::ReturnType 
	computeNewton(Teuchos::ParameterList& params);

	//! Applies Jacobian for extended system
	virtual NOX::Abstract::Group::ReturnType 
	applyJacobian(const NOX::Abstract::Vector& input, 
		      NOX::Abstract::Vector& result) const;

	//! Jacobian transpose not defined for this system
	virtual NOX::Abstract::Group::ReturnType 
	applyJacobianTranspose(const NOX::Abstract::Vector& input, 
			       NOX::Abstract::Vector& result) const;

	//! Applies Jacobian inverse for extended system
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
	virtual NOX::Abstract::Group::ReturnType 
	applyJacobianInverseMultiVector(
				    Teuchos::ParameterList& params, 
				    const NOX::Abstract::MultiVector& input, 
				    NOX::Abstract::MultiVector& result) const;

	//! Return \c true if extended residual is valid.
	virtual bool isF() const;

	//! Return \c true if the extended Jacobian is valid.
	virtual bool isJacobian() const;
      
	//! Always returns false
	virtual bool isGradient() const;
      
	//! Return \c true if the extended Newton direction is valid
	virtual bool isNewton() const;

	//! Return extended solution vector.  
	virtual const NOX::Abstract::Vector& getX() const;

	//! Return extended residual
	virtual const NOX::Abstract::Vector& getF() const;

	//! Return 2-norm of extended residual. 
	virtual double getNormF() const;

	//! Gradient is never valid
	virtual const NOX::Abstract::Vector& getGradient() const;

	//! Return extended Newton direction.
	virtual const NOX::Abstract::Vector& getNewton() const;

	//! Return RCP to extended solution vector.  
	virtual Teuchos::RCP< const NOX::Abstract::Vector > getXPtr() const;

	//! Return RCP to extended residual
	virtual Teuchos::RCP< const NOX::Abstract::Vector > getFPtr() const;

	//! Gradient is never valid
	virtual Teuchos::RCP< const NOX::Abstract::Vector > getGradientPtr() const;

	//! Return RCP to extended Newton direction.
	virtual Teuchos::RCP< const NOX::Abstract::Vector > getNewtonPtr() const;
	
	//! Returns 2-norm of extended Newton solve residual
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

	//! Set the parameter vector in the group to p (pVector = p).  
	virtual void setParams(const ParameterVector& p);

	//! Set parameter indexed by (integer) paramID
	virtual void setParam(int paramID, double val);

	//! Set parameter indexed by (std::string) paramID
	virtual void setParam(std::string paramID, double val);

	//! Return a const reference to the ParameterVector owned by the group.
	virtual const ParameterVector& getParams() const;

	//! Return copy of parameter indexed by (integer) paramID
	virtual double getParam(int paramID) const;

	//! Return copy of parameter indexed by (std::string) paramID
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

	//! Compute a scaled dot product
	virtual double
	computeScaledDotProduct(const NOX::Abstract::Vector& a,
				const NOX::Abstract::Vector& b) const;

	//! Function to print out solution and parameter after successful step
	virtual void printSolution(const double conParam) const;

	//! Function to print out a vector and parameter after successful step
	virtual void printSolution(const NOX::Abstract::Vector& x,
				   const double conParam) const;

	//! Scales a vector using scaling vector
	virtual void
	scaleVector(NOX::Abstract::Vector& x) const;

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

      protected:
	
	//! Resets all isValid flags to false
	virtual void resetIsValid();

	//! Sets up multivector views
	virtual void setupViews();

	//! Set bifurcation parameter
	void setBifParam(double param);

	//! Computes initial "a" and "b" vectors
	void getInitialVectors(
		    Teuchos::RCP<NOX::Abstract::Vector>& aVecPtr_real,
		    Teuchos::RCP<NOX::Abstract::Vector>& aVecPtr_imag,
		    Teuchos::RCP<NOX::Abstract::Vector>& bVecPtr_real,
		    Teuchos::RCP<NOX::Abstract::Vector>& bVecPtr_imag,
		    bool isSymmetric);

      private:
	
	//! Prevent generation and use of operator=()
	ExtendedGroup& operator=(const ExtendedGroup& source);

      protected:

	//! Pointer LOCA global data object
	Teuchos::RCP<LOCA::GlobalData> globalData;

	//! Parsed top-level parameters
	Teuchos::RCP<LOCA::Parameter::SublistParser> parsedParams;

	//! Hopf parameter list
	Teuchos::RCP<Teuchos::ParameterList> hopfParams;

	//! Pointer to base group that defines \f$F\f$
	Teuchos::RCP<LOCA::Hopf::MinimallyAugmented::AbstractGroup> grpPtr;

	//! Pointer to base group as a bordered group
	Teuchos::RCP<LOCA::BorderedSystem::AbstractGroup> bordered_grp;

	//! Pointer to constraint object
	Teuchos::RCP<LOCA::Hopf::MinimallyAugmented::Constraint> constraintsPtr;

	//! Stores the extended solution vector
	LOCA::MultiContinuation::ExtendedMultiVector xMultiVec;

	//! Stores the extended residual vector and df/dp
	LOCA::MultiContinuation::ExtendedMultiVector fMultiVec;

	//! Stores the extended Newton vector
	LOCA::MultiContinuation::ExtendedMultiVector newtonMultiVec;

	//! Stores the extended gradient vector
	LOCA::MultiContinuation::ExtendedMultiVector gradientMultiVec;

	//! Stores view of first column of xMultiVec
	Teuchos::RCP<LOCA::MultiContinuation::ExtendedVector> xVec;

	//! Stores view of first column of fMultiVec
	Teuchos::RCP<LOCA::MultiContinuation::ExtendedVector> fVec;
      
	//! Stores view of first column of fMultiVec as a multivec
	Teuchos::RCP<LOCA::MultiContinuation::ExtendedMultiVector> ffMultiVec;

	//! Stores view of df/dp columns of fMultiVec
	Teuchos::RCP<LOCA::MultiContinuation::ExtendedMultiVector> dfdpMultiVec;

	//! Stores view of f and first column of df/dp
	Teuchos::RCP<LOCA::MultiContinuation::ExtendedMultiVector> fBifMultiVec;

	//! Stores view of first column of newtonMultiVec
	Teuchos::RCP<LOCA::MultiContinuation::ExtendedVector> newtonVec;

	//! Stores view of first column of gradientMultiVec
	Teuchos::RCP<LOCA::MultiContinuation::ExtendedVector> gradientVec;

	// Stores operator for bordered solver
	Teuchos::RCP<LOCA::BorderedSolver::JacobianOperator> jacOp;

	//! Stores bordered solver strategy
	Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy> borderedSolver;

	//! Stores indices for getting f part of fMultiVec
	std::vector<int> index_f;

	//! Stores indices for getting df/dp part of fMultiVec
	std::vector<int> index_dfdp;

	//! Stores the bifurcation parameter index
	int bifParamID;

	//! Is residual vector valid
	bool isValidF;

	//! Is Jacobian matrix valid
	bool isValidJacobian;

	//! Is Newton vector valid
	bool isValidNewton;

	//! Is Gradient vector valid
	bool isValidGradient;

	//! Flag that indicates whether underlying group is a bordered group
	bool isBordered;

      }; // Class ExtendedGroup

    } // Namespace MinimallyAugmented

  } // Namespace Hopf

} // Namespace LOCA

#endif // LOCA_HOPF_MINIMALLYAUGMENTED_EXTENDEDGROUP_H
