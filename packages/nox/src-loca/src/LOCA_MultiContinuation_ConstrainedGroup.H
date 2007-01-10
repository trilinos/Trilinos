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

#ifndef LOCA_MULTICONTINUATION_CONSTRAINEDGROUP_H
#define LOCA_MULTICONTINUATION_CONSTRAINEDGROUP_H

#include "Teuchos_RefCountPtr.hpp"

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
  namespace MultiContinuation {
    class ConstraintInterface;
  }
  namespace BorderedSolver {
    class AbstractStrategy;
  }
}

namespace LOCA {

  namespace MultiContinuation {

    /*!
     * \brief Extended group representing a constrained nonlinear problem.
     */
    /*!
     * This class represents a constrained system of nonlinear equations:
     * \f[
     *     \begin{split}
     *         f(x,p) &= 0 \\
     *         g(x,p) &= 0
     *     \end{split}
     * \f]
     * where \f$x\in\Re^n\f$ is the solution vector, \f$p\in\Re^m\f$ is a 
     * set of constraint parameters, \f$f(x,p)\in\Re^n\f$ is represented by 
     * some LOCA::MultiContinuation::AbstractGroup, and \f$g(x,p)\in\Re^m\f$ 
     * is a constraint represented by a 
     * LOCA::MultiContinuation::ConstraintInterface object.  Newton steps
     * for this system are computed via some 
     * LOCA::BorderedSolver::AbstractStrategy which is specified via   
     * the \c constraintParams argument to the constructor.
     */
    class ConstrainedGroup : 
      public virtual LOCA::Extended::MultiAbstractGroup,
      public virtual LOCA::MultiContinuation::AbstractGroup,
      public virtual LOCA::BorderedSystem::AbstractGroup {

    public:

      //! Constructor 
      /*!
       * \param global_data [in] Global data object
       * \param topParams [in] Parsed top-level parameter list.
       * \param constraintParams [in] Parameter list determining the
       * bordered solver method.
       * \param grp [in] Group representing \f$f\f$.
       * \param constraints [in] Constraint object representing \f$g\f$.
       * \param paramIDs [in] Parameter IDs of the constraint parameters.
       */
      ConstrainedGroup(
       const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
       const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RefCountPtr<Teuchos::ParameterList>& constraintParams,
       const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& grp,
       const Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstraintInterface>& constraints,
       const vector<int>& paramIDs);

      //! Copy constructor
      ConstrainedGroup(const ConstrainedGroup& source, 
		       NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~ConstrainedGroup();

      //! Set constraint parameter \c i to value \c val
      virtual void setConstraintParameter(int i, double val);

      //! Get constraint parameter \c i 
      virtual double getConstraintParameter(int i) const;

      //! Get constraint parameter IDs
      virtual const vector<int>&
      getConstraintParamIDs() const;

      /*! 
       * @name Implementation of NOX::Abstract::Group virtual methods 
       */
      //@{

      //! Assignment operator
      virtual NOX::Abstract::Group& 
      operator=(const NOX::Abstract::Group& source);

      //! Clone function
      virtual Teuchos::RefCountPtr<NOX::Abstract::Group>
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
      Teuchos::RefCountPtr<const LOCA::MultiContinuation::AbstractGroup> 
      getUnderlyingGroup() const;
      
      //! Return underlying group
      virtual 
      Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>
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
		     const vector<int>& paramIDs, 
		     const NOX::Abstract::MultiVector::DenseMatrix& vals);

      //! Set the parameter vector in the group to p (pVector = p).  
      virtual void setParams(const ParameterVector& p);

       //! Set parameter indexed by (integer) paramID
      virtual void setParam(int paramID, double val);

      //! Set parameter indexed by (string) paramID
      virtual void setParam(string paramID, double val);

      //! Return a const reference to the ParameterVector owned by the group. 
      virtual const ParameterVector& getParams() const;

      //! Return copy of parameter indexed by (integer) paramID
      virtual double getParam(int paramID) const;

      //! Return copy of parameter indexed by (string) paramID
      virtual double getParam(string paramID) const;

      /*! 
       * Compute \f$\partial F/\partial p\f$ for each parameter \f$p\f$ 
       * indexed by paramIDs.  The first column of \em dfdp holds F,
       * which is valid if \em isValidF is true.  Otherwise F must be
       * computed.
       */
      virtual NOX::Abstract::Group::ReturnType
      computeDfDpMulti(const vector<int>& paramIDs, 
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
      virtual Teuchos::RefCountPtr<const NOX::Abstract::Group>
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
      loadNestedComponents(const NOX::Abstract::MultiVector& v_x,
			   const NOX::Abstract::MultiVector::DenseMatrix& v_p,
			   NOX::Abstract::MultiVector& v) const;

      //! Fill the combined A block as described above
      virtual void fillA(NOX::Abstract::MultiVector& A) const;

      //! Fill the combined B block as described above
      virtual void fillB(NOX::Abstract::MultiVector& B) const;

      //! Fill the combined C block as described above
      virtual void fillC(NOX::Abstract::MultiVector::DenseMatrix& C) const;

      //@}

      //! Get group
      virtual Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>
      getGroup();

      //! Get constraints
      virtual 
      Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstraintInterface>
      getConstraints();

    protected:

      //! Resets all isValid flags to false
      virtual void resetIsValid();

      //! Sets up multivector views
      virtual void setupViews();

    private:

      //! Prevent generation and use of operator=()
      ConstrainedGroup& operator=(const ConstrainedGroup& source);

    protected:

      //! Pointer LOCA global data object
      Teuchos::RefCountPtr<LOCA::GlobalData> globalData;

      //! Parsed top-level parameters
      Teuchos::RefCountPtr<LOCA::Parameter::SublistParser> parsedParams;

      //! Constraint parameter list
      Teuchos::RefCountPtr<Teuchos::ParameterList> constraintParams;

      //! Pointer to base group that defines \f$F\f$
      Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup> grpPtr;

      //! Pointer to base group as a bordered group
      Teuchos::RefCountPtr<LOCA::BorderedSystem::AbstractGroup> bordered_grp;

      //! Pointer to constraint object
      Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstraintInterface> constraintsPtr;

      //! Number of parameters
      int numParams;

      //! Stores the extended solution vector
      LOCA::MultiContinuation::ExtendedMultiVector xMultiVec;

      //! Stores the extended residual vector and df/dp
      LOCA::MultiContinuation::ExtendedMultiVector fMultiVec;

      //! Stores the extended Newton vector
      LOCA::MultiContinuation::ExtendedMultiVector newtonMultiVec;

      //! Stores the extended gradient vector
      LOCA::MultiContinuation::ExtendedMultiVector gradientMultiVec;

      //! Stores view of first column of xMultiVec
      Teuchos::RefCountPtr<LOCA::MultiContinuation::ExtendedVector> xVec;

      //! Stores view of first column of fMultiVec
      Teuchos::RefCountPtr<LOCA::MultiContinuation::ExtendedVector> fVec;
      
      //! Stores view of first column of fMultiVec as a multivec
      Teuchos::RefCountPtr<LOCA::MultiContinuation::ExtendedMultiVector> ffMultiVec;

      //! Stores view of df/dp columns of fMultiVec
      Teuchos::RefCountPtr<LOCA::MultiContinuation::ExtendedMultiVector> dfdpMultiVec;

      //! Stores view of first column of newtonMultiVec
      Teuchos::RefCountPtr<LOCA::MultiContinuation::ExtendedVector> newtonVec;

      //! Stores view of first column of gradientMultiVec
      Teuchos::RefCountPtr<LOCA::MultiContinuation::ExtendedVector> gradientVec;

      //! Stores bordered solver strategy
      Teuchos::RefCountPtr<LOCA::BorderedSolver::AbstractStrategy> borderedSolver;

      //! Stores indices for getting f part of fMultiVec
      vector<int> index_f;

      //! Stores indices for getting df/dp part of fMultiVec
      vector<int> index_dfdp;

      //! integer id of constraint parameters
      vector<int> constraintParamIDs;

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

      //! Flag indicating whether to skip df/dp computations
      bool skipDfDp;

    }; // Class ConstrainedGroup

  } // Namespace MultiContinuation

} // Namespace LOCA

#endif // LOCA_CONSTRAINEDSYSTEM_CONSTRAINEDGROUP_H
