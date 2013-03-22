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

#ifndef LOCA_HOPF_MOORESPENCE_EXTENDEDGROUP_H
#define LOCA_HOPF_MOORESPENCE_EXTENDEDGROUP_H

#include "Teuchos_RCP.hpp"

#include "LOCA_MultiContinuation_AbstractGroup.H"         // base class
#include "LOCA_Extended_MultiAbstractGroup.H"             // base class
#include "LOCA_Hopf_MooreSpence_ExtendedVector.H" // class data element
#include "LOCA_Hopf_MooreSpence_ExtendedMultiVector.H" // class data element

namespace LOCA {
  class ParameterVector;
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }

  namespace Hopf {
    namespace MooreSpence {
      class AbstractGroup;
      class SolverStrategy;
    }
  }
}
namespace LOCA { 

  namespace Hopf { 
    
    namespace MooreSpence {

      /*! 
       * \brief A group representing the Moore-Spence Hopf equations. 
       */
      /*!
       * The LOCA::Hopf::MooreSpence::ExtendedGroup is a concrete 
       * implementation of the NOX::Abstract::Group, 
       * LOCA::MultiContinuation::AbstractGroup and 
       * LOCA::Extended::MultiAbstractGroup that defines the following 
       * extended set of equations that are regular at a generic Hopf point:
       * \f[
       *     G(z) = \left[
       *     \begin{array}{c}
       *         F(x,p) \\
       *         Jy-\omega Bz \\
       *         Jz+\omega By \\
       *         l^Ty-1 \\
       *         l^Tz
       *      \end{array}
       *      \right] = 0
       * \f]
       * where \f$z = [x, y, z, \omega, p]\in\Re^{3n+2}\f$, \f$x\f$ is the 
       * solution vector, \f$y+i\omega z\f$ is the complex eigenvector of 
       * \f$J\f$ with corresponding eigenvalues \f$\pm i\omega\f$, \f$l\f$ 
       * is the length normalization vector and \f$J\f$ is the Jacobian of F 
       * w.r.t \f$x\f$.
       * 
       * The group stores an underlying group of type 
       * LOCA::Hopf::MooreSpence AbstractGroup to represent the equations
       * \f$F(x,p) = 0\f$ and to manipulate the underlying Jacobian \f$J\f$.
       * Note that the entire extended Jacobian \f$D_z G\f$ is not stored in
       * memory, rather a block-elimination algorithm (bordering algorithm) is
       * used to compute linear solves of the extended Jacobian (see 
       * LOCA::Hopf::MooreSpence::SolverFactory() for more details).
       *
       * This class implements all of the NOX::Abstract::Group, 
       * LOCA::MultiContinuation::AbstractGroup, and 
       * LOCA::Extended::MultiAbstractGroup 
       * methods for this extended set of equations and therefore is a complete
       * group which can be passed to most %NOX solvers to locate a single 
       * Hopf point or to the LOCA::Stepper to compute a family of Hopf 
       * points in a second parameter.
       *
       * However, Jacobian-tranpose operations and gradient calculations 
       * cannot be implemented efficiently and therefore gradient-base 
       * nonlinear solvers such as steepest descent and Trust region methods 
       * cannot be used to solve the extended Hopf point equations.
       *
       * The class is intialized via the \c hopfParams parameter list argument 
       * to the constructor.  The parameters this class recognizes are:
       * <ul>
       * <li> "Bifurcation Parameter" -- [string] (Must be supplied) Name of 
       *      the bifurcation parameter \f$p\f$
       * <li> "Length Normalization Vector" -- 
       *      [Teuchos::RCP<NOX::Abstract::Vector>] (Must be supplied) 
       *     Vector storing length normalization vector \f$l\f$
       * <li> "Initial Real Eigenvector"  -- 
       *      [Teuchos::RCP<NOX::Abstract::Vector>] (Must be supplied)
       *     Vector storing initial guess for the real component of the  
       *     eigenvector \f$y\f$
       * <li> "Initial Imaginary Eigenvector"  -- 
       *      [Teuchos::RCP<NOX::Abstract::Vector>] (Must be supplied)
       *     Vector storing initial guess for the imaginary component of the  
       *     eigenvector \f$z\f$
       * <li> "Initial Frequency" -- [double] (Must be supplied)
       *      Initial guess for the Hopf frequency \f$\omega\f$.
       * <li> "Perturb Initial Solution" -- [bool] (default: false) Flag 
       *      indicating whether to perturb the initial solution
       * <li> "Relative Perturbation Size" -- [double] (default: 1.0e-3) 
       *      Relative perturbation size if perturbing the initial solution
       * </ul>
       */
      class ExtendedGroup 
	: public virtual LOCA::Extended::MultiAbstractGroup,
	  public virtual LOCA::MultiContinuation::AbstractGroup {

      public:

	/*! 
	 * \brief Constructor with initial data passed through parameter lists.
	 */
	ExtendedGroup(
	 const Teuchos::RCP<LOCA::GlobalData>& global_data,
         const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& hopfParams,
	 const Teuchos::RCP<LOCA::Hopf::MooreSpence::AbstractGroup>& g);

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

	//! Compute the Hopf point equation residual \f$G\f$.
	virtual NOX::Abstract::Group::ReturnType computeF();

	//! Compute the blocks of the Jacobian derivative of \f$G\f$
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
	 * Uses a LOCA::Hopf::MooreSpence::SolverStrategy instantiated
	 * by the LOCA::Hopf::MooreSpence::SolverFactory to 
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

	//! Get bifurcation parameter
	double getBifParam() const;

	//! Get bifurcation frequency
	double getFrequency() const;

	//! Defines null vector normalization \f$l^Tz\f$ equation
	double lTransNorm(const NOX::Abstract::Vector& z) const;

	//! null vector normalization for multivectors
	/*!
	 * Note:  result should have 1 row and z.numVectors() columns.
	 */
	void lTransNorm(const NOX::Abstract::MultiVector& z,
			NOX::Abstract::MultiVector::DenseMatrix& result) const;

      protected:

	//! Set bifurcation parameter
	void setBifParam(double param);

	//! Sets up multivector views
	void setupViews();

	//! Initializes group
	void init(bool perturbSoln = false, double perturbSize = 0.0);

      private:

	//! Prohibit generation and use of operator=()
	ExtendedGroup& operator=(const ExtendedGroup&);

      protected:

	//! Pointer LOCA global data object
	Teuchos::RCP<LOCA::GlobalData> globalData;

	//! Parsed top-level parameters
	Teuchos::RCP<LOCA::Parameter::SublistParser> parsedParams;

	//! Bifurcation parameter list
	Teuchos::RCP<Teuchos::ParameterList> hopfParams;

	//! Pointer to base group that defines \f$F\f$
	Teuchos::RCP<LOCA::Hopf::MooreSpence::AbstractGroup> grpPtr;

	//! Stores the extended solution vector
	LOCA::Hopf::MooreSpence::ExtendedMultiVector xMultiVec;
	
	//! Stores the extended residual vector and df/dp
	LOCA::Hopf::MooreSpence::ExtendedMultiVector fMultiVec;
	
	//! Stores the extended Newton vector
	LOCA::Hopf::MooreSpence::ExtendedMultiVector newtonMultiVec;

	//! Stores the length normalization vector
	Teuchos::RCP<NOX::Abstract::MultiVector> lengthMultiVec;
	
	//! Stores view of first column of xMultiVec
	Teuchos::RCP<LOCA::Hopf::MooreSpence::ExtendedVector> xVec;
	
	//! Stores view of first column of fMultiVec
	Teuchos::RCP<LOCA::Hopf::MooreSpence::ExtendedVector> fVec;
	
	//! Stores view of first column of fMultiVec as a multivec
	Teuchos::RCP<LOCA::Hopf::MooreSpence::ExtendedMultiVector> ffMultiVec;

	//! Stores view of df/dp columns of fMultiVec
	Teuchos::RCP<LOCA::Hopf::MooreSpence::ExtendedMultiVector> dfdpMultiVec;

	//! Stores view of first column of newtonMultiVec
	Teuchos::RCP<LOCA::Hopf::MooreSpence::ExtendedVector> newtonVec;

	//! Stores view of first column of lengthMultiVec
	Teuchos::RCP<NOX::Abstract::Vector> lengthVec;

	//! Stores mass matrix times real component of eigenvector
	Teuchos::RCP<NOX::Abstract::MultiVector> massTimesY;

	/*!
	 * \brief Stores negative of mass matrix times imaginary component of 
	 eigenvector
	*/
	Teuchos::RCP<NOX::Abstract::MultiVector> minusMassTimesZ;

	//! Stores bordering strategy
	Teuchos::RCP<LOCA::Hopf::MooreSpence::SolverStrategy> solverStrategy;

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

      }; // class ExtendedGroup

    } // namespace MooreSpence

  } // namespace Hopf

} // namespace LOCA

#endif
