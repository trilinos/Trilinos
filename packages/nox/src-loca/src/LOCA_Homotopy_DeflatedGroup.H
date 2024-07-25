// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_HOMOTOPY_DEFAULTEDGROUP_H
#define LOCA_HOMOTOPY_DEFAULTEDGROUP_H

#include "Teuchos_RCP.hpp"

#include "LOCA_Extended_MultiAbstractGroup.H"           // base class
#include "LOCA_MultiContinuation_AbstractGroup.H"       // base class
#include "LOCA_BorderedSystem_AbstractGroup.H"          // base class
#include "LOCA_MultiContinuation_ExtendedVector.H"      // class data element
#include "LOCA_MultiContinuation_ExtendedMultiVector.H" // class data element
#include "LOCA_Parameter_Vector.H"

// forward declarations
namespace Teuchos {
  class ParameterList;
}
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
  namespace Homotopy {
    class AbstractGroup;
  }
  namespace BorderedSolver {
    class AbstractStrategy;
    class JacobianOperator;
  }
}

namespace LOCA {

  namespace Homotopy {

    /*!
     * \brief %LOCA's Homotopy Algorithm.
     */
    /*!
     * The %HomotopyGroup is a concrete implementation of the
     * LOCA::Continuation::AbstractGroup that modifies the set of nonlinear
     * equations to be solved to allow for Homotopy to be applied to the
     * system.  This object should be used in conjunction with the
     * LOCA::Stepper object to drive the continuation.  This algorithm solves
     * a system of nonlinear equations supplied by the user (\f$ F(x) \f$)
     * through continuation.  An artificial parameter \f$ \lambda \f$ is used
     * to control the continuation.  The idea is to solve a simple equation
     * starting at \f$ \lambda \f$ = 0 and, using the solution from the
     * previous step, solve systems of equations that gets progressively
     * closer to the true system of interest ( at \f$ \lambda \f$ = 1.0 we
     * recover the original equations \f$ F(x) \f$).  By constraining the
     * definition of \f$ g(x, \lambda) \f$ and using artificial parameter
     * contiuation, the continuation branch should be free of multiplicity and
     * bifurcation phenomena.
     *
     * The modified system of equations, \f$ g(x, \lambda) \f$, supplied by
     * the HomotopyGroup is defined as:
     *
     * \f[ g(x, \lambda) = \lambda F(x) + (1.0 - \lambda)(x - a)(S) \f]
     *
     * where \f$ x\f$ is the solution vector, \f$ \lambda \f$ is an artificial
     * parameter, \f$ F(x) \f$ is the set of nonlinear equations the user
     * supplies, \f$ g(x) \f$ is the corresponding set of
     * homotopy equations that LOCA will solve, \f$ a \f$ is a random
     * vector, and  \f$ S \f$ is a scaling factor used to switch sign of
     * the last term (typically valued 1.0 or -1.0).
     *
     * This group requires the loca Stepper for continuation
     * from \f$ \lambda \f$ = 0.0 (a simple set of equations to solve) to
     * \f$ \lambda \f$ = 1.0 (the set of equations requested by the user,
     * \f$ F(x) \f$).  The Homotopy::Group will generate the Stepper parameter
     * sublist in the parameter list that is passed in to the constructor.
     * The user is free to modify this list (it sets default values) before
     * passing it into the stepper object but should NOT change the starting
     * and stopping values for the continuation parameter.
     *
     * References:
     *
     * - ALGORITHM 652 HOMPACK: A Suite of Codes for Globally Convergent
     * Homotopy Algorithms, Watson, L.T., Billups, S.C, and Morgan, A.P.,
     * ACM Transactions on Mathematical Software, Vol. 13, No. 3, September
     * 1987, pp281-310.
     */
    class DeflatedGroup :
      public virtual LOCA::Extended::MultiAbstractGroup,
      public virtual LOCA::MultiContinuation::AbstractGroup,
      public virtual LOCA::BorderedSystem::AbstractGroup {

    public:

      //! Constructor
      /*!
       * \param global_data [in] Global data object
       * \param topParams [in] Parsed top-level parameter list.
       * \param hParams [in] Homotopy parameters
       * \param grp [in] Group representing \f$ f\f$.
       */
      DeflatedGroup(
       const Teuchos::RCP<LOCA::GlobalData>& global_data,
       const Teuchos::RCP<Teuchos::ParameterList>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& hParams,
       const Teuchos::RCP<LOCA::Homotopy::AbstractGroup>& grp,
       const Teuchos::RCP<const NOX::Abstract::Vector>& start_vec,
       const std::vector< Teuchos::RCP<const NOX::Abstract::Vector> >& prev_solns,
       const double identity_sign = 1.0);

      //! Copy constructor
      DeflatedGroup(const DeflatedGroup& source,
            NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~DeflatedGroup();

      //! Get homotopy parameter
      double getHomotopyParam() const;

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
       * Compute \f$\partial F/\partial p\f$ for each parameter \f$ p\f$
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

      //! Set homotopy parameter
      void setHomotopyParam(double param);

      //! Creates and sets the "Stepper" parameter sublist
      void setStepperParameters(Teuchos::ParameterList& topParams);

    private:

      //! Prevent generation and use of operator=()
      DeflatedGroup& operator=(const DeflatedGroup& source);

    protected:

      //! Pointer LOCA global data object
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Parsed top-level parameters
      Teuchos::RCP<LOCA::Parameter::SublistParser> parsedParams;

      //! Pitchfork parameter list
      Teuchos::RCP<Teuchos::ParameterList> homotopyParams;

      //! Pointer to base group that defines \f$ F\f$
      Teuchos::RCP<LOCA::Homotopy::AbstractGroup> grpPtr;

      //! Pointer to base group as a bordered group
      Teuchos::RCP<LOCA::BorderedSystem::AbstractGroup> bordered_grp;

      //! Stores the extended solution vector
      LOCA::MultiContinuation::ExtendedMultiVector xMultiVec;

      //! Stores the extended residual vector
      LOCA::MultiContinuation::ExtendedMultiVector fMultiVec;

      //! Stores the extended Newton vector
      LOCA::MultiContinuation::ExtendedMultiVector newtonMultiVec;

      //! Stores the extended gradient vector
      LOCA::MultiContinuation::ExtendedMultiVector gradientMultiVec;

      //! Stores view of first column of xMultiVec
      Teuchos::RCP<LOCA::MultiContinuation::ExtendedVector> xVec;

      //! Stores view of first column of fMultiVec
      Teuchos::RCP<LOCA::MultiContinuation::ExtendedVector> fVec;

      //! Stores view of first column of newtonMultiVec
      Teuchos::RCP<LOCA::MultiContinuation::ExtendedVector> newtonVec;

      //! Stores view of first column of gradientMultiVec
      Teuchos::RCP<LOCA::MultiContinuation::ExtendedVector> gradientVec;

      //! Stores the starting vector, \f$ x^\ast \f$.
      Teuchos::RCP<const NOX::Abstract::Vector> startVec;

      //!
      double identitySign;

      //! Stores array of deflated solution vectors
      std::vector< Teuchos::RCP<const NOX::Abstract::Vector> > solns;

      //! Temporary vector for storing distance
      Teuchos::RCP<NOX::Abstract::Vector> distVec;

      //! Temporary vector for storing sum of distances
      Teuchos::RCP<NOX::Abstract::MultiVector> totalDistMultiVec;

      //! Temporary vector for storing sum of distances
      Teuchos::RCP<NOX::Abstract::Vector> totalDistVec;

      //! underlying group's F vector as a multivector
      Teuchos::RCP<NOX::Abstract::MultiVector> underlyingF;

      // Stores operator for bordered solver
      Teuchos::RCP<LOCA::BorderedSolver::JacobianOperator> jacOp;

      //! Stores bordered solver strategy
      Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy> borderedSolver;

      //! Stores -1
      Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix> minusOne;

      //! Number of previous solutions
      int numSolns;

      //! Array of distances from previous solutions
      std::vector<double> distances;

      //! Product of distances
      double distProd;

      //! Stores indices for getting f part of fMultiVec
      std::vector<int> index_f;

      //! Copy of the ParameterVector for the underlying grpPtr.
      /*! We copy this and then add the homotopy parameter to the list. */
      LOCA::ParameterVector paramVec;

      //! Value of the homotopy continuation parameter.
      /*! Ranges from 0.0 (easy solution) to 1.0 (solution to the system of
       * interest).
       */
      double conParam;

      //! Continuatioin parameter ID number from the ParameterVector.
      int conParamID;

      /*!
       * \brief Contains the std::string used to identify the homotopy parameter in
       * the ParameterVector object.
       */
      const std::string conParamLabel;

      /*!
       * \brief Tracks whether the LOCA::Homotopy::Group method
       * augmentJacobianForHomotopy is implemented.  If not, the augmentation
       * is applied during the applyJacobian assuming a matrix-free
       * implementation.
       */
      bool augmentJacForHomotopyNotImplemented;

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

    }; // Class DeflatedGroup

  } // Namespace Homotopy

} // Namespace LOCA

#endif // LOCA_HOMOTOPY_DEFLATEDGROUP_H
