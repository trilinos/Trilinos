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

#ifndef LOCA_MULTICONTINUATION_EXTENDEDGROUP_H
#define LOCA_MULTICONTINUATION_EXTENDEDGROUP_H

#include "Teuchos_RefCountPtr.hpp"

#include "LOCA_MultiContinuation_AbstractStrategy.H"    // base class
#include "LOCA_MultiContinuation_ExtendedVector.H"      // class data element
#include "LOCA_MultiContinuation_ExtendedMultiVector.H" // class data element

// forward declarations
namespace NOX {
  namespace Parameter {
    class List;
  }
}
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
  namespace MultiContinuation {
    class AbstractGroup;
    class ConstraintInterface;
    class ConstrainedGroup;
  }
  namespace MultiPredictor {
    class AbstractStrategy;
  }
}

namespace LOCA { 

  namespace MultiContinuation { 

    /*!
     * \brief Base class for all continuation groups. 
     */
    /*!
     * %Continuation is defined as computing some curve 
     * \f$(x(s),p(s))\in\Re^{n+1}\f$ such that \f$F(x(s),p(s))=0\f$ for some
     * parameterization \f$s\f$.  Given some point \f$(x_0,p_0)\f$ on the
     * curve, another nearby point on the curve is calculated by first 
     * computing a predictor direction \f$v\in\Re^{n+1}\f$ and the 
     * approximate point \f$(x^\ast,p^\ast) = (x_0,p_0) + v\Delta s\f$ 
     * where \f$\Delta s\f$ is the step size.  Then the next point on the 
     * curve is computed by solving the extended set of equations
     *  \f[
     *      \begin{array}{cc}
     *         F(x,p) &= 0 \\
     *         g(x,p,x_0,p_0,x^\ast,p^\ast,v,\Delta s) &= 0
     *      \end{array}
     *  \f]
     *  for \f$(x,p)\f$.  The equation 
     *  \f$g(x,p,x_0,p_0,x^\ast,p^\ast,v,\Delta s)=0\f$ is called the 
     *  continuation equation
     *  and different choices of \f$g\f$ yield different continuation methods.
     *
     * Mathematically, this computation amounts to repeatedly computing 
     * solutions to a constrained nonlinear system.  This class provides 
     * a common implementation for all continuation groups in terms of the 
     * LOCA::MultiContinuation::ConstrainedGroup using a supplied group
     * to represent \f$F\f$ and an implementation of
     * LOCA::MultiContinuation::ConstraintInterface to represent \f$g\f$.
     *
     * Note that this class has no public constructor other than the 
     * copy constructor since it is intended to only provide an implemenation
     * of much of the continuation work.  Each derived class that 
     * implements a specific continuation strategy should provide its own
     * public constructor.
     */
    class ExtendedGroup : 
      public virtual LOCA::MultiContinuation::AbstractStrategy {

    public:

      //! Copy constructor
      ExtendedGroup(const ExtendedGroup& source, 
		    NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~ExtendedGroup();

       //! Assignment operator
      virtual ExtendedGroup& 
      operator=(const ExtendedGroup& source);

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
      computeNewton(NOX::Parameter::List& params);

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
      applyJacobianInverse(NOX::Parameter::List& params, 
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
				    NOX::Parameter::List& params, 
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

      //! Assignment operator
      virtual LOCA::Extended::MultiAbstractGroup& 
      operator=(const LOCA::Extended::MultiAbstractGroup& source);

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
       * @name Implementation of LOCA::MultiContinuation::AbstractStrategy 
       * virtual methods 
       */
      //@{

      //! Assignment operator
      virtual LOCA::MultiContinuation::AbstractStrategy&
      operator=(const MultiContinuation::AbstractStrategy& source);

      //! Returns number of parameters
      virtual int getNumParams() const;

      //! Notifies group that continuation step is completed
      virtual void notifyCompletedStep();

      //! Compute predictor directions
      virtual NOX::Abstract::Group::ReturnType 
      computePredictor();

      //! Is Predictor valid
      virtual bool isPredictor() const;

      //! Scales tangent to predictor
      virtual void scaleTangent();

      //! Sets tangent to predictor
      /*!
       * This is required by MF which takes the tangent space, 
       * orthogonalizes it, and then sets it back in the group.
       */
      virtual void
      setPredictorTangentDirection(
			    const LOCA::MultiContinuation::ExtendedVector& v,
			    int i);

      //! Returns tangent to predictor
      virtual const LOCA::MultiContinuation::ExtendedMultiVector&
      getPredictorTangent() const;

      //! Returns scaled tangent to predictor
      virtual const LOCA::MultiContinuation::ExtendedMultiVector&
      getScaledPredictorTangent() const;

      //! Set the previous solution vector y
      virtual void setPrevX(const NOX::Abstract::Vector& y);

      //! Gets the previous solution vector
      virtual const LOCA::MultiContinuation::ExtendedVector& 
      getPrevX() const ;

      //! Set step size for continuation constraint equation \em i
      virtual void setStepSize(double deltaS, int i = 0);

      //! Get step size for continuation constraint equation \em i
      virtual double getStepSize(int i = 0) const;

      //! Sets the value for continuation parameter \em i
      virtual void setContinuationParameter(double val, int i = 0);

      //! Returns the value for continuation parameter \em i
      virtual double getContinuationParameter(int i = 0) const;

      //! Get the continuation parameter id for parameter \em i
      virtual int getContinuationParameterID(int i = 0) const;

      //! Get the continuation parameter ids
      virtual const vector<int>& getContinuationParameterIDs() const;

      //! Get the continuation parameter id for parameter \em i
      virtual string getContinuationParameterName(int i = 0) const;

      //! Returns step size scale factor for constraint equation \em i
      virtual double getStepSizeScaleFactor(int i = 0) const;

      //! Prints the group
      virtual void printSolution() const;

      //! Computes a scaled dot product between two continuation vectors
      virtual double computeScaledDotProduct(
			       const NOX::Abstract::Vector& x, 
			       const NOX::Abstract::Vector& y) const;

      //! Returns dimension of project to draw array
      virtual int projectToDrawDimension() const;

      //! Fills the project to draw array
      virtual void projectToDraw(
			    const LOCA::MultiContinuation::ExtendedVector& x, 
			    double *px) const;

      //@}

    protected:

      //! Constructor used by derived classes
      ExtendedGroup(
      const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
      const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RefCountPtr<NOX::Parameter::List>& continuationParams,
      const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& grp,
      const Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy>& pred,
      const vector<int>& paramIDs);

      //! Set constraint object
      /*!
       * This allows the constraint object to be set after the group
       * is constructed using the above constructor.
       */
      virtual void setConstraints(const Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstraintInterface>& constraints);

    protected:

      //! Pointer LOCA global data object
      Teuchos::RefCountPtr<LOCA::GlobalData> globalData;

      //! Parsed top-level parameters
      Teuchos::RefCountPtr<LOCA::Parameter::SublistParser> parsedParams;

      //! Continuation parameter list
      Teuchos::RefCountPtr<NOX::Parameter::List> continuationParams;

      //! Pointer to underlying group
      Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup> grpPtr;

      //! Pointer to predictor object
      Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy> predictor;

      //! Pointer to constrained group implementation
      Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstrainedGroup> conGroup;

      //! Number of parameters
      int numParams;

      //! Stores the tangent to the predictor
      LOCA::MultiContinuation::ExtendedMultiVector tangentMultiVec;

      //! Stores the scaled tangent to the predictor
      LOCA::MultiContinuation::ExtendedMultiVector scaledTangentMultiVec;


      //! Stores the previous extended solution vector
      LOCA::MultiContinuation::ExtendedVector prevXVec;

      //! integer id of continuation parameters
      vector<int> conParamIDs;

      //! continuation step size
      vector<double> stepSize;

      //! step size scale factors
      vector<double> stepSizeScaleFactor;

      //! Is Predictor vector valid
      bool isValidPredictor;

      //! Flag indicating whether to base predictor direction on secant
      bool baseOnSecant;

    };

  } // namespace MultiContinuation
} // namespace LOCA

#endif
