// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_TURNINGPOINT_MINIMALLYAUGMENTED_MODIFIEDCONSTRAINT_H
#define LOCA_TURNINGPOINT_MINIMALLYAUGMENTED_MODIFIEDCONSTRAINT_H

#include "LOCA_TurningPoint_MinimallyAugmented_Constraint.H" // base class

namespace LOCA {

  namespace TurningPoint {

    namespace MinimallyAugmented {

      /*!
       * \brief Implementation of
       * LOCA::MultiContinuation::ConstraintInterfaceMVDX for computing
       * turning points for the minimally augmented turning point formulation.
       */
      /*!
       * This class is a modification of
       * LOCA::TurningPoint::MinimallyAugmented::Constraint where updates
       * are computed to the left and right null vectors \f$w\f$ and \f$v\f$
       * every nonlinear iteration instead of solving for them directly:
       * \f[
       *     \begin{bmatrix}
       *       J   & a \\
       *       b^T & 0
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *       \Delta v \\
       *       \Delta \sigma_1
       *     \end{bmatrix} = -
       *     \begin{bmatrix}
       *       J v + \sigma_1 a + (Jv)_x\Delta x + (Jv)_p \Delta p\\
       *       b^T a - n
       *     \end{bmatrix},
       * \f]
       * \f[
       *     \begin{bmatrix}
       *       J^T & b \\
       *       a^T & 0
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *       \Delta w \\
       *       \Delta \sigma_2
       *     \end{bmatrix} = -
       *     \begin{bmatrix}
       *       J^T w + \sigma_2 b + (J^T w)_x \Delta x + (J^T w)_p \Delta p \\
       *       a^T w - n
       *     \end{bmatrix},
       * \f]
       *
       * The class is intialized via the \c tpParams parameter list argument
       * to the constructor.  This class recognizes all paramters for
       * LOCA::TurningPoint::MinimallyAugmented::Constraint plus the following:
       * <ul>
       * <li> "Include Newton Terms"  -- [bool] (default: false) -
       *    Flag indicating whether to include the \f$\Delta x\f$ and
       *    \f$\Delta p\f$ terms above when computing the null vector updates
       * </ul>
       */
      class ModifiedConstraint : public Constraint {

      public:

    //! Constructor
    ModifiedConstraint(
     const Teuchos::RCP<LOCA::GlobalData>& global_data,
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& tpParams,
     const Teuchos::RCP<LOCA::TurningPoint::MinimallyAugmented::AbstractGroup>& g,
     int bif_param);

    //! Copy constructor
    ModifiedConstraint(const ModifiedConstraint& source,
               NOX::CopyType type = NOX::DeepCopy);

    //! Destructor
    virtual ~ModifiedConstraint();

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

    //! Compute continuation constraint equations
    virtual NOX::Abstract::Group::ReturnType
    computeConstraints();

    //! Perform any preprocessing before a continuation step starts.
    /*!
     * The \c stepStatus argument indicates whether the previous step was
     * successful.  Here we set up the constraint class to solve for
     * \f$w\f$ and \f$v\f$ for the first nonlinear iteration.
     */
    virtual void
    preProcessContinuationStep(
                  LOCA::Abstract::Iterator::StepStatus stepStatus);

    //! Perform any postprocessing after a continuation step finishes.
    /*!
     * The \c stepStatus argument indicates whether the step was
     * successful.  Here we set up the constraint class to solve for
     * \f$w\f$ and \f$v\f$ for the first nonlinear iteration.
     */
    virtual void
    postProcessContinuationStep(
                  LOCA::Abstract::Iterator::StepStatus stepStatus);

    //@}

    //! Set the newton update for x and p
    void setNewtonUpdates(const NOX::Abstract::Vector& dx, double dp,
                  double step);

      private:

    //! Prohibit generation and use of operator=()
    ModifiedConstraint& operator=(const ModifiedConstraint& source);

      protected:

    //! Stores update to left null vector
    Teuchos::RCP<NOX::Abstract::MultiVector> w_vector_update;

    //! Stores update to right null vector
    Teuchos::RCP<NOX::Abstract::MultiVector> v_vector_update;

    //! Stores left null vector residual
    Teuchos::RCP<NOX::Abstract::MultiVector> w_residual;

    //! Stores right null vector residual
    Teuchos::RCP<NOX::Abstract::MultiVector> v_residual;

    //! Stores solution update
    Teuchos::RCP<NOX::Abstract::MultiVector> deltaX;

    //! Stores sigma_1
    NOX::Abstract::MultiVector::DenseMatrix sigma1;

    //! Stores sigma_1
    NOX::Abstract::MultiVector::DenseMatrix sigma2;

    //! Stores parameter update
    double deltaP;

    /*!
     * Flag that indicates whether we're in the first solve per
     * continuation step.
     */
    bool isFirstSolve;

    //! Flag indicating whether to include the newton update terms
    bool includeNewtonTerms;


      }; // Class ModifiedConstraint

    } // namespace MinimallyAugmented

  } // namespace TurningPoint

} // namespace LOCA

#endif // LOCA_TURNINGPOINT_MINIMALLYAUGMENTED_MODIFIEDCONSTRAINT_H
