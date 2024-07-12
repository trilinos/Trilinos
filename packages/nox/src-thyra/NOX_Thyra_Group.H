// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_THYRA_GROUP_H
#define NOX_THYRA_GROUP_H

#include "Teuchos_RCP.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "NOX_Abstract_Group.H"    // base class
#include "NOX_Common.H"             // class data element (std::string)
#include "NOX_Thyra_Vector.H"        // class data element
#include "NOX_SharedObjectTemplate.H"  // class data element

// Forward declares
namespace NOX {
  struct SolverStats;
  namespace Parameter {
    class List;
  }
}

namespace Thyra {
  template <class Scalar> class ModelEvaluator;
  template <class Scalar> class MultiVectorBase;
  template <class Scalar> class LinearOpBase;
  template <class Scalar> class LinearOpWithSolveBase;
  template <class Scalar> class LinearOpWithSolveFactoryBase;
  template <class Scalar> class DefaultLinearOpSource;
  template <class Scalar> class PreconditionerFactoryBase;
  template <class Scalar> class PreconditionerBase;
}

namespace NOX {

  namespace Thyra {

    /** \brief A concrete implementation of the NOX::Abstract::Group using Thyra.

        NOTE: This Group supports row sum scaling of the function
        (residual and Jacobian).  This is enabled by setting a weight
        vector on the initial guess vector in the Group constructor.
        The residual and Jacobian must be scaled before and then
        unscaled after calls to construct the preconditioner and solve
        the linear system.  This follows the nox epetra group.  This
        design should be changed in a future nox refactor, but
        requires significant changes to the Group object.
    */
    class Group : public virtual NOX::Abstract::Group {

    public:

      /** \brief The default constructor that uses the linear solver from the ModelEvaluator.

          Most users should use this constructor.  It is meant to be
          used in conjunction with a stratimikos linear solver that is
          built as part of the input model evaluator.  For finer
          grained control over the use of the preconditioner and for
          Jacobian-Free Newton-Krylov cases, the power user
          constructor should be used.

          \param[in] initialGuess Initial guess for the solution vector
          \param[in] model ModelEvaluator
          \param[in] weightVector Optional diagonal weighting vector for the model.
          \param[in] rightWeightVector Optional solution vector weighting
          \param[in] inv_rightWeightVector Optional inverse solution vector weighting
          \param[in] rightScalingFirst Optional bool to select if right scaling should be applied before left scaling
      */
      Group(const NOX::Thyra::Vector& initialGuess,
            const Teuchos::RCP<const ::Thyra::ModelEvaluator<double> >& model,
            const Teuchos::RCP<const ::Thyra::VectorBase<double> >& weightVector = Teuchos::null,
            const Teuchos::RCP<const ::Thyra::VectorBase<double> >& rightWeightVector = Teuchos::null,
            const Teuchos::RCP<::Thyra::VectorBase<double> >& inv_rightWeightVector = Teuchos::null,
            const bool rightScalingFirst = false);

      /** \brief Power user constructor that takes explicit linear solver objects to handle different combinations.

          This class allows the user to set user-defined linear
          operators and preconditioners (and corresponding
          factories). The user can set the linear_op to be a
          Jacobian-Free Newton Krylov operator (use the class
          NOX::Thyra::MatrixFreeJacobianOperator).

          \param[in] initialGuess (Required) Initial guess for the solution vector
          \param[in] model (Required) ModelEvaluator
          \param[in] linearOp (Optional) Forward operator for the Jacobian.  Must be non-null for Newton-based solvers.
          \param[in] lowsFactory (Optional) Factory for building and updating linear solver.
          \param[in] precOp (Optional) Preconditioner operator.  If set to Teuchos::null and a non-null prec_factory exists, the prec_op will be constructed using the preconditioner factory.
          \param[in] precFactory (Optional) Factory for updating the precondiitoner.  If set to Teuchos::null and there is a non-null prec_op, then the preconditioner will be updated using the model evaluator as long as the ModelEvaluator::outArgs supports W_prec.
          \param[in] weightVector (Optional) diagonal weighting vector for the model.
          \param[in] rightWeightVector Optional solution vector weighting
          \param[in] inv_rightWeightVector Optional inverse solution vector weighting
          \param[in] rightScalingFirst Optional bool to select if right scaling should be applied before left scaling
          \param[in] updatePreconditioner Optional bool to select if the Group should auotmatically update the preconditioner matrix values between Newton iterations
          \param[in] jacobianIsEvaluated Optional bool, if true this means that the input Jacobian operator (linearOp) has been evaluated externally and is consistent with the initialGuess. In this case, the isJacobian() flag is initialized to true.
      */
      Group(const NOX::Thyra::Vector& initialGuess,
            const Teuchos::RCP<const ::Thyra::ModelEvaluator<double> >& model,
            const Teuchos::RCP< ::Thyra::LinearOpBase<double> >& linearOp,
            const Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<double> >& lowsFactory,
            const Teuchos::RCP< ::Thyra::PreconditionerBase<double> >& precOp,
            const Teuchos::RCP< ::Thyra::PreconditionerFactoryBase<double> >& precFactory,
            const Teuchos::RCP<const ::Thyra::VectorBase<double> >& weightVector = Teuchos::null,
            const Teuchos::RCP<const ::Thyra::VectorBase<double> >& rightWeightVector = Teuchos::null,
            const Teuchos::RCP<::Thyra::VectorBase<double> >& inv_rightWeightVector = Teuchos::null,
            const bool rightScalingFirst = false,
            const bool updatePreconditioner = true,
            const bool jacobianIsEvaluated = false);

      //! Copy constructor
      Group(const NOX::Thyra::Group& source, NOX::CopyType type = DeepCopy);

      //! Destructor.
      ~Group();

      NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);

      NOX::Abstract::Group& operator=(const NOX::Thyra::Group& source);

      Teuchos::RCP<const ::Thyra::VectorBase<double> > get_current_x() const;

      Teuchos::RCP< ::Thyra::LinearOpBase<double> >
      getNonconstJacobianOperator();

      Teuchos::RCP<const ::Thyra::LinearOpBase<double> >
      getJacobianOperator() const;

      Teuchos::RCP<const ::Thyra::LinearOpBase<double> >
      getScaledJacobianOperator() const;

      void unscaleJacobianOperator() const;

      Teuchos::RCP< ::Thyra::LinearOpWithSolveBase<double> >
      getNonconstJacobian();

      Teuchos::RCP<const ::Thyra::LinearOpWithSolveBase<double> >
      getJacobian() const;

      Teuchos::RCP< ::Thyra::PreconditionerBase<double> >
      getNonconstPreconditioner();

      Teuchos::RCP<const ::Thyra::PreconditionerBase<double> >
      getPreconditioner() const;

      /// Dangerous power user function for LOCA Householder bordered algorithm.
      void setJacobianOperator(const Teuchos::RCP<::Thyra::LinearOpBase<double>>& op);

      /// Dangerous power user function for LOCA Householder bordered algorithm. This is the Matrix M that is used to initialize a stratimikos preconditioner. NOTE: this sets the losb_ object used to update prec_!
      void setPreconditionerMatrix(const Teuchos::RCP<const ::Thyra::DefaultLinearOpSource<double>>& op);

      /** @name "Compute" functions. */
      //@{

      void setX(const NOX::Abstract::Vector& y);
      //! See above
      void setX(const NOX::Thyra::Vector& y);

      void computeX(const NOX::Abstract::Group& grp,
            const NOX::Abstract::Vector& d,
            double step);
      //! See above.
      void computeX(const NOX::Thyra::Group& grp,
            const NOX::Thyra::Vector& d,
            double step);

      NOX::Abstract::Group::ReturnType computeF();

      NOX::Abstract::Group::ReturnType computeJacobian();

      NOX::Abstract::Group::ReturnType computeGradient();

      NOX::Abstract::Group::ReturnType
      computeNewton(Teuchos::ParameterList& params);

      //@}

      /** @name Jacobian operations.
       *
       * Operations using the Jacobian matrix. These may not be defined in
       * matrix-free scenarios. */

      //@{

      NOX::Abstract::Group::ReturnType
      applyJacobian(const NOX::Abstract::Vector& input,
            NOX::Abstract::Vector& result) const;

      NOX::Abstract::Group::ReturnType
      applyJacobian(const NOX::Thyra::Vector& input,
            NOX::Thyra::Vector& result) const;

      NOX::Abstract::Group::ReturnType
      applyJacobianMultiVector(const NOX::Abstract::MultiVector& input,
                   NOX::Abstract::MultiVector& result) const;

      NOX::Abstract::Group::ReturnType
      applyJacobianTranspose(const NOX::Abstract::Vector& input,
                 NOX::Abstract::Vector& result) const;

      NOX::Abstract::Group::ReturnType
      applyJacobianTranspose(const NOX::Thyra::Vector& input,
                 NOX::Thyra::Vector& result) const;

      NOX::Abstract::Group::ReturnType
      applyJacobianTransposeMultiVector(
                     const NOX::Abstract::MultiVector& input,
                     NOX::Abstract::MultiVector& result) const;

      NOX::Abstract::Group::ReturnType
      applyJacobianInverse(Teuchos::ParameterList& params,
               const NOX::Abstract::Vector& input,
               NOX::Abstract::Vector& result) const;

      NOX::Abstract::Group::ReturnType
      applyJacobianInverse(Teuchos::ParameterList& params,
               const NOX::Thyra::Vector& input,
               NOX::Thyra::Vector& result) const;

      NOX::Abstract::Group::ReturnType
      applyJacobianInverseMultiVector(
                    Teuchos::ParameterList& params,
                    const NOX::Abstract::MultiVector& input,
                    NOX::Abstract::MultiVector& result) const;

      NOX::Abstract::Group::ReturnType
      applyRightPreconditioning(bool useTranspose,
                Teuchos::ParameterList& params,
                const NOX::Abstract::Vector& input,
                NOX::Abstract::Vector& result) const;
      //@}

      /** @name "Is" functions
       *
       * Checks to see if various objects have been computed. Returns true
       * if the corresponding "compute" function has been called since the
       * last update to the solution vector (via instantiation or
       * computeX). */
      //@{

      bool isF() const;
      bool isJacobian() const;
      bool isGradient() const;
      bool isNewton() const;

      //@}

      /** @name "Get" functions
       *
       * Note that these function do not check whether or not the vectors
       * are valid. Must use the "Is" functions for that purpose. */
      //@{

      const NOX::Abstract::Vector& getX() const;

      const NOX::Abstract::Vector& getScaledX() const;

      const NOX::Abstract::Vector& getF() const;

      double getNormF() const;

      const NOX::Abstract::Vector& getGradient() const;

      const NOX::Abstract::Vector& getNewton() const;

      Teuchos::RCP< const NOX::Abstract::Vector > getXPtr() const;

      Teuchos::RCP< const NOX::Abstract::Vector > getFPtr() const;

      Teuchos::RCP< const NOX::Abstract::Vector > getGradientPtr() const;

      Teuchos::RCP< const NOX::Abstract::Vector > getNewtonPtr() const;

      //@}

      virtual void logLastLinearSolveStats(NOX::SolverStats& stats) const;

      virtual Teuchos::RCP<NOX::Abstract::Group>
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Print out the group
      void print() const;

      Teuchos::RCP< const ::Thyra::ModelEvaluator<double> > getModel() const;

      /** \brief Set the transient terms on the Group and use them in the underlying evalModelImpl() calls.

          \param x_dot Time derivative term. Can be set to null.
          \param alpha Model evaluator transient timer derivative multiplier
          \param beta Model evaluator transient Jacobian multiplier
          \param t Current time value
      */
      void enablePseudoTransientTerms(const Teuchos::RCP<const ::Thyra::VectorBase<double>>& x_dot,
                                      const double alpha, const double beta, const double t);

      /// Disable the pseudo traansient terms in the underlying evalModel() calls. Sets x_dot, alpha, beta and t back to steady state values.
      void disablePseudoTransientTerms();

      /// Check for whether the pseudo transient support is enabled for residual and Jacobian evaluations.
      bool usingPseudoTransientTerms() const;

      /// Set default parameters to be used with inArgs. This is a
      /// hack to support a bad design in PIRO. PIRO wraps nox as a
      /// model evaluator. The PIRO inArgs suppots all inargs, but nox
      /// only knows about nonlinear solver realted inputs. It knows
      /// nothing about parameters. PIRO should have a model evaluator
      /// wrapper that populates the extra inArgs as needed to
      /// override nominal values. Instead, we have to hard code
      /// storage of all in arg parameters in nox groups. This is only
      /// used in the ::Thyra::NOXNonlinearSolver and only when
      /// wrapped within a PIRO model evaluator that calls the
      /// setBasePoint on ::Thyra::NOXNonlinearSolver. This function
      /// should eventually be deprecated in favor of addign a wrapper
      /// ME to the piro nox solver class.
      void setBasePoint(const ::Thyra::ModelEvaluatorBase::InArgs<double>& base_point_params);

      /// Unset the base point parameters so that they are not used internally.
      void unsetBasePoint();

      /// Returns true if a base point has been set.
      bool usingBasePoint() const;

      /// Returns the thyra linear solver factory.
      Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<double>>
      getLinearOpWithSolveFactory() const;

      /// Returns the thyra preconditioner factory.
      Teuchos::RCP<::Thyra::PreconditionerFactoryBase<double>>
      getPreconditionerFactory() const;

      /// User will take control of updating an already registered preconditioner. Typically used with NOX::ObserverReusePreconditioner observer.
      void takeControlOfPreconditionerUpdates(const Teuchos::RCP< ::Thyra::PreconditionerBase<double>>& prec);

    protected:

      //! resets the isValid flags to false
      void resetIsValidFlags();

      //! Apply Jacobian inverse using Thyra objects
      NOX::Abstract::Group::ReturnType
      applyJacobianInverseMultiVector(
                  Teuchos::ParameterList& p,
                  const ::Thyra::MultiVectorBase<double>& input,
                  ::Thyra::MultiVectorBase<double>& result) const;

      ::Thyra::ESolveMeasureNormType
      getThyraNormType(const std::string& name) const;

      //! Finalizes LOWS to be a valid solver for the Jacobian
      void updateLOWS() const;

      void scaleResidualAndJacobian() const;

      void unscaleResidualAndJacobian() const;

      void computeScaledSolution();

    protected:

      //! Problem interface
      Teuchos::RCP< const ::Thyra::ModelEvaluator<double> > model_;

      /** @name IsValid flags
       *
       * True if the current solution is up-to-date with respect to the
       * currect solution vector. */
      //@{
      bool is_valid_f_;
      bool is_valid_jacobian_;
      bool is_valid_newton_dir_;
      bool is_valid_gradient_dir_;
      mutable bool is_valid_lows_;
      //@}

      //! Solution vector
      Teuchos::RCP<NOX::Thyra::Vector> x_vec_;

      //! Residual vector
      Teuchos::RCP<NOX::Thyra::Vector> f_vec_;

      //! Newton direction vector
      Teuchos::RCP<NOX::Thyra::Vector> newton_vec_;

      //! Gradient direction vector
      Teuchos::RCP<NOX::Thyra::Vector> gradient_vec_;

      //! Shared Jacobian operator with solve
      Teuchos::RCP<
        NOX::SharedObject<
          ::Thyra::LinearOpWithSolveBase<double>,
          NOX::Thyra::Group
          >
        > shared_jacobian_;

      //! Jacobian operator
      Teuchos::RCP< ::Thyra::LinearOpBase<double> > lop_;

      //! Thyra LOWS factory for building Jacobians
      Teuchos::RCP< const ::Thyra::LinearOpWithSolveFactoryBase<double> > lows_factory_;

      //! Source base needed to create preconditioner
      Teuchos::RCP< const ::Thyra::DefaultLinearOpSource<double> > losb_;

      //! Preconditioner for Jacobian
      Teuchos::RCP< ::Thyra::PreconditionerBase<double> > prec_;

      //! Preconditioner factory
      Teuchos::RCP< ::Thyra::PreconditionerFactoryBase<double> > prec_factory_;

      /** \brief Optional wieghting vector for function scaling. NOX assumes that this vector can be updated in between nonlinear iterations.

           This is pulled out of the initial guess vector
      */
      Teuchos::RCP<const ::Thyra::VectorBase<double> > weight_vec_;
      //! Optional wieghting vector for solution (right) scaling. 
      // NOX no longer assumes that this vector fixed. It is allowed to change during an entire nonlinear solve.
      // The application is in charge of updating the vector in an observer.
      Teuchos::RCP<const ::Thyra::VectorBase<double> > right_weight_vec_;

      //! Inverse of weight vector used to unscale function (left) scaling. NOX assumes that this vector can be updated in between nonlinear iterations.
      mutable Teuchos::RCP< ::Thyra::VectorBase<double> > inv_weight_vec_;
      //! Inverse of weight vector used to unscale solution (right) scaling.
      // NOX no longer assumes that this vector fixed. It is allowed to change during an entire nonlinear solve.
      // The application is in charge of updating the vector in an observer.
      Teuchos::RCP< ::Thyra::VectorBase<double> > inv_right_weight_vec_;

      //! Scaled solution vector scaled by the 
      mutable Teuchos::RCP<NOX::Thyra::Vector> scaled_x_vec_;

      //! Do right scaling before left scaling?
      bool rightScalingFirst_;

      //! If set to true, the preconditioner matrix values will be automatically updated via precFactory or ModelEvalautor. If set to false, the user must manually handle updating the preconditioner.
      bool updatePreconditioner_;

      //! The status of the last linear solve performed
      mutable NOX::Abstract::Group::ReturnType last_linear_solve_status_;

      //! Number of iterations for last linear solve performed
      mutable int last_linear_solve_num_iters_;

      //! The tolerance achieved by the last linear solver
      mutable double last_linear_solve_achieved_tol_;

      // True if pseudo transient term is enabled
      bool use_pseudo_transient_terms_;
      // Only used for pseudo transient
      Teuchos::RCP<const ::Thyra::VectorBase<double>> x_dot_;
      // Only used for pseudo transient
      double alpha_;
      // Only used for pseudo transient
      double beta_;
      // Only used for pseudo transient
      double t_;

      // True if we will use a default inArgs set on this group. This
      // is used for the thyra nonlinear solver interface that allows
      // users to set an inArgs as the base point. It is a hack to
      // work around bad PIRO design.
      bool use_base_point_;

      // Will override null parameter values on in_args with
      // these. This is used for the thyra nonlinear solver interface
      // that allows users to set an inArgs as the base point. It is a
      // hack to work around bad PIRO design.
      ::Thyra::ModelEvaluatorBase::InArgs<double> base_point_;
    };

  } // namespace LAPACK
} // namespace NOX


#endif
