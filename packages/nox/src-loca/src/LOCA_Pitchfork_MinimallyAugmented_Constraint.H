// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_PITCHFORK_MINIMALLYAUGMENTED_CONSTRAINT_H
#define LOCA_PITCHFORK_MINIMALLYAUGMENTED_CONSTRAINT_H

#include "LOCA_TurningPoint_MinimallyAugmented_Constraint.H" // base class

// forward declarations
namespace LOCA {
  namespace Pitchfork {
    namespace MinimallyAugmented {
      class AbstractGroup;
    }
  }
}

namespace LOCA {

  namespace Pitchfork {

    namespace MinimallyAugmented {

      /*!
       * \brief Implementation of
       * LOCA::MultiContinuation::ConstraintInterfaceMVDX for computing
       * pitchforks for the minimally augmented pitchfork formulation.
       */
      /*!
       * This class implements the pitchfork constraint equations
       * \f$\sigma(x,p) = 0\f$, \f$\langle \psi,x \rangle = 0\f$ for
       * the minimally augmented pitchfork formulation where \f$\sigma\f$
       * is defined via
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
       * The class is derived from
       * LOCA::TurningPoint::MinimallyAugmented::Constraint.  See this class
       * for a description of available parameters.
       */
      class Constraint :
    public LOCA::TurningPoint::MinimallyAugmented::Constraint {

      public:

    //! Constructor
    Constraint(
     const Teuchos::RCP<LOCA::GlobalData>& global_data,
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& pfParams,
     const Teuchos::RCP<LOCA::Pitchfork::MinimallyAugmented::AbstractGroup>& g,
     const Teuchos::RCP<const NOX::Abstract::Vector>& psi,
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
    computeDP(const std::vector<int>& paramIDs,
          NOX::Abstract::MultiVector::DenseMatrix& dgdp,
          bool isValidG);

    //! Return constraint residuals
    virtual const NOX::Abstract::MultiVector::DenseMatrix&
    getConstraints() const;

    //! Return solution component of constraint derivatives
    virtual const NOX::Abstract::MultiVector*
    getDX() const;

    //@}

      private:

    //! Prohibit generation and use of operator=()
    Constraint& operator=(const Constraint& source);

      protected:

    //! Pointer to pitchfork group
    Teuchos::RCP<LOCA::Pitchfork::MinimallyAugmented::AbstractGroup> pf_grp;

    //! Vector for \f$\psi\f$
    Teuchos::RCP<const NOX::Abstract::Vector> psi_vector;

    //! Stores dg/dx
    Teuchos::RCP<NOX::Abstract::MultiVector> dgdx;

    //! Constraint values
    NOX::Abstract::MultiVector::DenseMatrix pf_constraints;

      }; // Class Constraint

    } // namespace MinimallyAugmented

  } // namespace Pitchfork

} // namespace LOCA

#endif // LOCA_PITCHFORK_MINIMALLYAUGMENTED_CONSTRAINT_H
