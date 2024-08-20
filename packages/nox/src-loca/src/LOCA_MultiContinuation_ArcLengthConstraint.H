// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_MULTICONTINUATION_ARCLENGTHCONSTRAINT_H
#define LOCA_MULTICONTINUATION_ARCLENGTHCONSTRAINT_H

#include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H" // base class
#include "NOX_Abstract_MultiVector.H"                   // for dense matrix

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace MultiContinuation {
    class ArcLengthGroup;
  }
}

namespace LOCA {

  namespace MultiContinuation {

    /*!
     * \brief Implementation of
     * LOCA::MultiContinuation::ConstraintInterfaceMVDX for arclength
     * continuation.
     */
    /*!
     * This class implements the arclength constraint equation for
     * pseudo-arclength continuation:
     * \f[
     *     g(x,p,x_0,p_0,x^\ast,p^\ast,v,\Delta s)=
     *         (x-x^\ast)^Tv_x + (p-p^\ast) v_p - \Delta s
     * \f]
     * where \f$v_x\f$, \f$v_p\f$ are the solution and  parameter components
     * of the predictor direction \f$v\f$ respectively.  Since the derivative
     * of \f$g\f$ with respect to \f$x\f$ is just \f$v\f$, the predictor
     * tangent, this class implements the MVDX version of the
     * constraint interface.
     */
    class ArcLengthConstraint :
      public LOCA::MultiContinuation::ConstraintInterfaceMVDX {

    public:

      //! Constructor
      ArcLengthConstraint(
     const Teuchos::RCP<LOCA::GlobalData>& global_data,
     const Teuchos::RCP<LOCA::MultiContinuation::ArcLengthGroup>& grp);

      //! Copy constructor
      ArcLengthConstraint(const ArcLengthConstraint& source,
              NOX::CopyType type = NOX::DeepCopy);

      //! Destructor
      ~ArcLengthConstraint();

      //! Set pointer to arclength group
      virtual void setArcLengthGroup(const Teuchos::RCP<LOCA::MultiContinuation::ArcLengthGroup>& grp);

      /*!
       * @name Implementation of LOCA::MultiContinuation::ConstraintInterface
       * virtual methods
       */
      //@{

      //! Copy
      virtual void copy(const ConstraintInterface& source);

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
             const std::vector<int>& paramIDs,
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
      computeDP(const std::vector<int>& paramIDs,
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

      //@}

    private:

      //! Prohibit generation and use of operator=()
      ArcLengthConstraint& operator=(const ArcLengthConstraint& source);

    protected:

      //! LOCA global data object
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Pointer to arc-length group
      Teuchos::RCP<
    LOCA::MultiContinuation::ArcLengthGroup> arcLengthGroup;

      //! Constraint values
      NOX::Abstract::MultiVector::DenseMatrix constraints;

      //! Flag indicating whether constraints are valid
      bool isValidConstraints;

      //! Continuation parameter IDs
      std::vector<int> conParamIDs;

    }; // Class ArcLengthConstraint

  } // namespace MultiContinuation

} // namespace LOCA

#endif // LOCA_MULTICONTINUATION_ARCLENGTHCONSTRAINT_H
