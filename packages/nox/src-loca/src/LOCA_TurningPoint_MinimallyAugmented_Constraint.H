// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
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

#ifndef LOCA_TURNINGPOINT_MINIMALLYAUGMENTED_CONSTRAINT_H
#define LOCA_TURNINGPOINT_MINIMALLYAUGMENTED_CONSTRAINT_H

#include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H" // base class
#include "NOX_Abstract_MultiVector.H"                       // for dense matrix

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
  namespace TurningPoint {
    namespace MinimallyAugmented {
      class AbstractGroup;
    }
  }
  namespace BorderedSolver {
    class AbstractStrategy;
  }
}

namespace LOCA { 

  namespace TurningPoint {

    namespace MinimallyAugmented {

      /*! 
       * \brief Implementation of 
       * LOCA::MultiContinuation::ConstraintInterfaceMVDX for computing
       * turning points for the minimally augmented turning point formulation.
       */
      /*!
       * This class implements the turning point constraint equation 
       * \f$\sigma(x,p) = 0\f$ for the minimally augmented turning point 
       * formulation where \f$\sigma\f$ is defined via
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
       *     \sigma = -w^T J v/n
       * \f]
       * for any vectors \f$a\f$ and \f$b\f$ in \f$\Re^n\f$.  Using these
       * relationships, it is easy to show
       * \f[
       *     \begin{split}
       *        \sigma_x &= -(w^T J v)_x/n = -w^T J_x v/n \\
       *        \sigma_p &= -(w^T J v)_p/n = -w^T J_p v/n
       *     \end{split}
       * \f]
       *
       * The class is intialized via the \c tpParams parameter list argument 
       * to the constructor.  The parameters this class recognizes are:
       * <ul>
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
      class Constraint : 
	public LOCA::MultiContinuation::ConstraintInterfaceMVDX {

      public:

	//! Constructor
	Constraint(
	 const Teuchos::RCP<LOCA::GlobalData>& global_data,
	 const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& tpParams,
	 const Teuchos::RCP<LOCA::TurningPoint::MinimallyAugmented::AbstractGroup>& g,
	 bool is_symmetric,
	 const NOX::Abstract::Vector& a,
	 const NOX::Abstract::Vector* b,
	 int bif_param);

	//! Copy constructor
	Constraint(const Constraint& source, 
		   NOX::CopyType type = NOX::DeepCopy);

	//! Destructor
	virtual ~Constraint();

	//! Set the group pointer
	/*!
	 * This method should be called when ever the constrained group
	 * is copied, since we don't explicitly copy the underlying group
	 * here.
	 */
	virtual void setGroup(const Teuchos::RCP<LOCA::TurningPoint::MinimallyAugmented::AbstractGroup>& g);

	//! Returns left null vector w
	virtual Teuchos::RCP<const NOX::Abstract::Vector>
	getLeftNullVec() const;

	//! Returns right null vector v
	virtual Teuchos::RCP<const NOX::Abstract::Vector>
	getRightNullVec() const;

	//! Returns a vector
	virtual Teuchos::RCP<const NOX::Abstract::Vector>
	getAVec() const;

	//! Returns b vector
	virtual Teuchos::RCP<const NOX::Abstract::Vector>
	getBVec() const;

	//! Returns sigma
	virtual double getSigma() const;

	/*! 
	 * @name Implementation of LOCA::MultiContinuation::ConstraintInterface
	 * virtual methods 
	 */
	//@{

	//! Copy
	virtual void 
	copy(const LOCA::MultiContinuation::ConstraintInterface& source);

	//! Cloning function
	virtual 
	Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
	clone(NOX::CopyType type = NOX::DeepCopy) const;

	//! Return number of constraints
	virtual int numConstraints() const;

	//! Set the solution vector to y.
	virtual void setX(const NOX::Abstract::Vector& y);

	//! Sets parameter indexed by paramID
	virtual void setParam(int paramID, double val);

	//! Sets parameters indexed by paramIDs
	virtual void setParams(
			  const vector<int>& paramIDs, 
			  const NOX::Abstract::MultiVector::DenseMatrix& vals);

	//! Compute continuation constraint equations
	virtual NOX::Abstract::Group::ReturnType
	computeConstraints();

	//! Compute derivative of constraints w.r.t. solution vector x
	virtual NOX::Abstract::Group::ReturnType
	computeDX();

	//! Compute derivative of constraints w.r.t. supplied parameters.
	/*!
	 * The first column of \c dgdp should be filled with the constraint
	 * residuals \f$g\f$ if \c isValidG is \c false.  If \c isValidG is
	 * \c true, then the \c dgdp contains \f$g\f$ on input.
	 */
	virtual NOX::Abstract::Group::ReturnType
	computeDP(const vector<int>& paramIDs, 
		  NOX::Abstract::MultiVector::DenseMatrix& dgdp, 
		  bool isValidG);

	//! Return \c true if constraint residuals are valid
	virtual bool isConstraints() const;

	//! Return \c true if derivatives of constraints w.r.t. x are valid
	virtual bool isDX() const;

	//! Return constraint residuals
	virtual const NOX::Abstract::MultiVector::DenseMatrix&
	getConstraints() const;

	//! Return solution component of constraint derivatives
	virtual const NOX::Abstract::MultiVector*
	getDX() const;

	/*! 
	 * \brief Return \c true if solution component of constraint 
	 * derivatives is zero
	 */
	virtual bool isDXZero() const;

	//! Perform any postprocessing after a continuation step finishes.
	/*!
	 * The \c stepStatus argument indicates whether the step was
	 * successful.  Here we update the \f$a\f$ and \f$b\f$ vectors to 
	 * \f$w\f$ and \f$v\f$ respectively if requested.
	 */
	virtual void 
	postProcessContinuationStep(
			      LOCA::Abstract::Iterator::StepStatus stepStatus);

	//@}

      private:

	//! Prohibit generation and use of operator=()
	Constraint& operator=(const Constraint& source);

      protected:

	//! Pointer LOCA global data object
	Teuchos::RCP<LOCA::GlobalData> globalData;

	//! Parsed top-level parameters
	Teuchos::RCP<LOCA::Parameter::SublistParser> parsedParams;

	//! Bifurcation parameter list
	Teuchos::RCP<Teuchos::ParameterList> turningPointParams;

	//! Pointer to base group that defines \f$F\f$
	Teuchos::RCP<LOCA::TurningPoint::MinimallyAugmented::AbstractGroup> grpPtr;

	//! Vector for \f$a\f$
	Teuchos::RCP<NOX::Abstract::MultiVector> a_vector;

	//! Vector for \f$b\f$
	Teuchos::RCP<NOX::Abstract::MultiVector> b_vector;

	//! Stores left null vector
	Teuchos::RCP<NOX::Abstract::MultiVector> w_vector;

	//! Stores right null vector
	Teuchos::RCP<NOX::Abstract::MultiVector> v_vector;

	//! Stores J*v
	Teuchos::RCP<NOX::Abstract::MultiVector> Jv_vector;

	//! Stores sigma_x
	Teuchos::RCP<NOX::Abstract::MultiVector> sigma_x;

	//! Constraint values
	NOX::Abstract::MultiVector::DenseMatrix constraints;

	//! Stores bordered solver strategy
	Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy> borderedSolver;

	//! Stores vector length as a double
	double dn;

	//! Stores scale factor on sigma
	double sigma_scale;

	/*!
	 * \brief Flag indicating whether Jacobian is symmetric
	 */
	bool isSymmetric;

	//! Flag indicating whether constraints are valid
	bool isValidConstraints;

	//! Flag indicating whether sigma_x is valid
	bool isValidDX;

	//! Stores the bifurcation parameter index
	vector<int> bifParamID;

	/*! 
	 * \brief Flag indicating whether to update \f$a\f$ and \f$b\f$
	 * every continuation step.
	 */
	bool updateVectorsEveryContinuationStep;

	/*! 
	 * \brief Flag indicating whether to update \f$a\f$ and \f$b\f$
	 * every nonlinear iteration.
	 */
	bool updateVectorsEveryIteration;

      }; // Class Constraint

    } // namespace MinimallyAugmented

  } // namespace TurningPoint

} // namespace LOCA

#endif // LOCA_MULTICONTINUATION_ARCLENGTHCONSTRAINT_H
