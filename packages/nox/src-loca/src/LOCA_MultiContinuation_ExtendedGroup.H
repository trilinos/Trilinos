// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_MULTICONTINUATION_EXTENDEDGROUP_H
#define LOCA_MULTICONTINUATION_EXTENDEDGROUP_H

#include "Teuchos_RCP.hpp"

#include "LOCA_MultiContinuation_AbstractStrategy.H"    // base class
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
      public virtual LOCA::MultiContinuation::AbstractStrategy,
      public virtual LOCA::BorderedSystem::AbstractGroup {

    public:

      //! Copy constructor
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

      //! Return extended solution vector.
      virtual Teuchos::RCP< const NOX::Abstract::Vector > getXPtr() const;

      //! Return extended residual
      virtual Teuchos::RCP< const NOX::Abstract::Vector > getFPtr() const;

      //! Gradient is never valid
      virtual Teuchos::RCP< const NOX::Abstract::Vector > getGradientPtr() const;

      //! Return extended Newton direction.
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
       * @name Implementation of LOCA::MultiContinuation::AbstractStrategy
       * virtual methods
       */
      //@{

      //! Assignment operator
      virtual void copy(const NOX::Abstract::Group& source);

      //! Returns number of parameters
      virtual int getNumParams() const;

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
      virtual const std::vector<int>& getContinuationParameterIDs() const;

      //! Get the continuation parameter id for parameter \em i
      virtual std::string getContinuationParameterName(int i = 0) const;

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

    protected:

      //! Constructor used by derived classes
      ExtendedGroup(
      const Teuchos::RCP<LOCA::GlobalData>& global_data,
      const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RCP<Teuchos::ParameterList>& continuationParams,
      const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& grp,
      const Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>& pred,
      const std::vector<int>& paramIDs);

      //! Set constraint object
      /*!
       * This allows the constraint object to be set after the group
       * is constructed using the above constructor.
       */
      virtual void setConstraints(const Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>& constraints, bool skip_dfdp);

    private:

      //! Prohibit generation and use of operator=()
      ExtendedGroup& operator=(const ExtendedGroup& source);

    protected:

      //! Pointer LOCA global data object
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Parsed top-level parameters
      Teuchos::RCP<LOCA::Parameter::SublistParser> parsedParams;

      //! Continuation parameter list
      Teuchos::RCP<Teuchos::ParameterList> continuationParams;

      //! Pointer to underlying group
      Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> grpPtr;

      //! Pointer to predictor object
      Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy> predictor;

      //! Pointer to constrained group implementation
      Teuchos::RCP<LOCA::MultiContinuation::ConstrainedGroup> conGroup;

      //! Number of parameters
      int numParams;

      //! Stores the tangent to the predictor
      LOCA::MultiContinuation::ExtendedMultiVector tangentMultiVec;

      //! Stores the scaled tangent to the predictor
      LOCA::MultiContinuation::ExtendedMultiVector scaledTangentMultiVec;


      //! Stores the previous extended solution vector
      LOCA::MultiContinuation::ExtendedVector prevXVec;

      //! integer id of continuation parameters
      std::vector<int> conParamIDs;

      //! continuation step size
      std::vector<double> stepSize;

      //! step size scale factors
      std::vector<double> stepSizeScaleFactor;

      //! Is Predictor vector valid
      bool isValidPredictor;

      //! Flag indicating whether to base predictor direction on secant
      bool baseOnSecant;

    };

  } // namespace MultiContinuation
} // namespace LOCA

#endif
