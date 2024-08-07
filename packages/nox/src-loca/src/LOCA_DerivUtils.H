// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_DERIVUTILS_H
#define LOCA_DERIVUTILS_H

#include "Teuchos_RCP.hpp"

#include "NOX_Abstract_Group.H" // For CopyType
#include "NOX_Abstract_MultiVector.H" // For DenseMatrix
#include <vector>

namespace LOCA {
  class GlobalData;
  namespace MultiContinuation {
    class AbstractGroup;
  }
  namespace Hopf {
    namespace MooreSpence {
      class AbstractGroup;
    }
     namespace MinimallyAugmented {
      class AbstractGroup;
    }
  }
}

namespace LOCA {

  /*!
   * \brief %LOCA's generic derivative computation class to compute various
   * derivatives via finite differencing.
   *
   * The %DerivUtils class provides generic methods to compute the following
   * derivatives:
   * \f[
   * \frac{\partial F}{\partial p},\quad\frac{\partial Jn}{\partial p},\quad\frac{\partial Jn}{\partial x}a
   * \f]
   * where \f$ J = \partial F/\partial x\f$ and \f$ n\f$, \f$ a\f$ are vectors.
   * These quantities are calculate by finite differencing.
   */
  class DerivUtils {

  public:

    /*!
     * \brief Default constructor.  perturb is the relative perturbation size
     * used in differencing calculations.
     */
    DerivUtils(const Teuchos::RCP<LOCA::GlobalData>& global_data,
           double perturb = 1.0e-6);

    //! Copy constructor.
    DerivUtils(const DerivUtils&);

    //! Destructor.
    virtual ~DerivUtils();

    /*!
     * \brief Cloning function.  Creates a copy of the %DerivUtils object
     * of the same type.
     */
    virtual Teuchos::RCP<DerivUtils>
    clone(NOX::CopyType type = NOX::DeepCopy) const;

    /*!
     * \brief Compute derivative of f with respect to parameter, identified
     * by param_id
     */
    virtual NOX::Abstract::Group::ReturnType
    computeDfDp(LOCA::MultiContinuation::AbstractGroup& grp,
        const std::vector<int>& param_ids,
        NOX::Abstract::MultiVector& result,
        bool isValidF) const;

    //! Compute derivative of Jn with respect to particular parameter param_id
    virtual NOX::Abstract::Group::ReturnType
    computeDJnDp(LOCA::MultiContinuation::AbstractGroup&,
         const std::vector<int>& paramIDs,
         const NOX::Abstract::Vector& nullVector,
         NOX::Abstract::MultiVector& result,
         bool isValid) const;

    //! Compute vector (Jn)_{x}a  given multi-vector a
    virtual NOX::Abstract::Group::ReturnType
    computeDJnDxa(LOCA::MultiContinuation::AbstractGroup& grp,
          const NOX::Abstract::Vector& nullVector,
          const NOX::Abstract::MultiVector& aVector,
          NOX::Abstract::MultiVector& result) const;

    //! Compute vector (Jn)_{x}a  given multi-vector a, given JnVector
    virtual NOX::Abstract::Group::ReturnType
    computeDJnDxa(LOCA::MultiContinuation::AbstractGroup& grp,
          const NOX::Abstract::Vector& nullVector,
          const NOX::Abstract::MultiVector& aVector,
          const NOX::Abstract::Vector& JnVector,
          NOX::Abstract::MultiVector& result) const;

    /*!
     * \brief Compute derivative of w^TJn with respect to particular parameter
     * param_id
     */
    virtual NOX::Abstract::Group::ReturnType
    computeDwtJnDp(LOCA::MultiContinuation::AbstractGroup& grp,
           const std::vector<int>& paramIDs,
           const NOX::Abstract::Vector& w,
           const NOX::Abstract::Vector& nullVector,
           NOX::Abstract::MultiVector::DenseMatrix& result,
           bool isValid) const;

    /*!
     * \brief Compute derivative of w^TJ with respect to particular parameter
     * param_id
     */
    virtual NOX::Abstract::Group::ReturnType
    computeDwtJDp(LOCA::MultiContinuation::AbstractGroup& grp,
          const std::vector<int>& paramIDs,
          const NOX::Abstract::Vector& w,
          NOX::Abstract::MultiVector& result,
          bool isValid) const;

    //! Compute vector (w^TJn)_{x}
    virtual NOX::Abstract::Group::ReturnType
    computeDwtJnDx(LOCA::MultiContinuation::AbstractGroup& grp,
           const NOX::Abstract::Vector& w,
           const NOX::Abstract::Vector& nullVector,
           NOX::Abstract::Vector& result) const;

        //! Compute vector (w^TJn)_{x}
    virtual NOX::Abstract::Group::ReturnType
    computeDwtJnDx(LOCA::MultiContinuation::AbstractGroup& grp,
           const NOX::Abstract::MultiVector& w,
           const NOX::Abstract::Vector& nullVector,
           NOX::Abstract::MultiVector& result) const;

    //! Compute derivative of (J+iwM)(y+iz) with respect to parameter,
    virtual NOX::Abstract::Group::ReturnType
    computeDCeDp(LOCA::Hopf::MooreSpence::AbstractGroup& grp,
         const std::vector<int>& paramIDs,
         const NOX::Abstract::Vector& yVector,
         const NOX::Abstract::Vector& zVector,
         double w,
         NOX::Abstract::MultiVector& result_real,
         NOX::Abstract::MultiVector& result_imag,
         bool isValid) const;

    //! Compute vector (J+iwM)(y+iz))_{x}a  given a
    virtual NOX::Abstract::Group::ReturnType
    computeDCeDxa(LOCA::Hopf::MooreSpence::AbstractGroup& grp,
          const NOX::Abstract::Vector& yVector,
          const NOX::Abstract::Vector& zVector,
          double w,
          const NOX::Abstract::MultiVector& aVector,
          NOX::Abstract::MultiVector& result_real,
          NOX::Abstract::MultiVector& result_imag) const;

    //! Compute vector (J+iwM)(y+iz))_{x}a  given a and (J+iwM)(y+iz) vector
    virtual NOX::Abstract::Group::ReturnType
    computeDCeDxa(LOCA::Hopf::MooreSpence::AbstractGroup& grp,
          const NOX::Abstract::Vector& yVector,
          const NOX::Abstract::Vector& zVector,
          double w,
          const NOX::Abstract::MultiVector& aVector,
          const NOX::Abstract::Vector& Ce_real,
          const NOX::Abstract::Vector& Ce_imag,
          NOX::Abstract::MultiVector& result_real,
          NOX::Abstract::MultiVector& result_imag) const;

    //! Compute derivative of (w1+iw2)^T(J+iwM)(y+iz) w.r.t. parameter p
    virtual NOX::Abstract::Group::ReturnType
    computeDwtCeDp(LOCA::Hopf::MinimallyAugmented::AbstractGroup& grp,
           const std::vector<int>& paramIDs,
           const NOX::Abstract::Vector& w1,
           const NOX::Abstract::Vector& w2,
           const NOX::Abstract::Vector& yVector,
           const NOX::Abstract::Vector& zVector,
           double omega,
           NOX::Abstract::MultiVector::DenseMatrix& result_real,
           NOX::Abstract::MultiVector::DenseMatrix& result_imag,
           bool isValid) const;

    //! Compute vector (w1+iw2)^T(J+iwM)(y+iz))_{x}
    virtual NOX::Abstract::Group::ReturnType
    computeDwtCeDx(LOCA::Hopf::MinimallyAugmented::AbstractGroup& grp,
           const NOX::Abstract::Vector& w1,
           const NOX::Abstract::Vector& w2,
           const NOX::Abstract::Vector& yVector,
           const NOX::Abstract::Vector& zVector,
           double omega,
           NOX::Abstract::Vector& result_real,
           NOX::Abstract::Vector& result_imag) const;

  protected:

    /*!
     * \brief Perturb parameter param_id in the group's parameter vector,
     * and return perturbation size
     */
    virtual double perturbParam(LOCA::MultiContinuation::AbstractGroup& grp,
                double& param,
                int param_id) const;

    /*!
     * \brief Perturb group's solution vector in direction of aVector,
     * and return perturbation size
     */
    virtual double perturbXVec(LOCA::MultiContinuation::AbstractGroup& grp,
                   const NOX::Abstract::Vector& xVector,
                   const NOX::Abstract::Vector& aVector) const;

    /*!
     * \brief Return perturbation for finite differencing around given
     * double variable
     */
    virtual double epsScalar(double p) const;

    /*!
     * \brief Return perturbation for directional derivatives around first
     * vector in direction of second vector
     */
    virtual double epsVector(const NOX::Abstract::Vector& xVector,
                 const NOX::Abstract::Vector& aVector) const;

  protected:

    //! Global data
    Teuchos::RCP<LOCA::GlobalData> globalData;

    //! Relative perturbation size using in differencing calculations.
    double perturb;

  }; // class DerivUtils
} // namespace LOCA

#endif
