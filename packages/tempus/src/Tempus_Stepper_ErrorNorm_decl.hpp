//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_Stepper_ErrorNorm_decl_hpp
#define Tempus_Stepper_ErrorNorm_decl_hpp

#include "Tempus_config.hpp"

#include "Teuchos_RCPDecl.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorSpaceFactoryBase.hpp"
namespace Tempus {

/** \brief Stepper_ErrorNorm provides error norm calcualtions for variable time
 * stepping.
 *
 */
template <class Scalar>
class Stepper_ErrorNorm {
 public:
  /// Default Constructor
  Stepper_ErrorNorm();

  /// Constructor
  Stepper_ErrorNorm(const Scalar relTol, const Scalar absTol);

  /// Destructor
  ~Stepper_ErrorNorm(){};

  /** \brief Compute the weigthed root mean square norm.
   *
   *  The WRMS norm is
   *  \f[
   *    e_{wrms} \equiv \sqrt{ \frac{1}{N} \sum_{i=1}^N \left( \frac
   *       {u^n_i}{A_{tol} + \max (|u^n_i|, |u^{n+1}_i|) R_{tol}} \right) ^2 }
   *  \f]
   *  where
   *  - \f$A_{tol}\f$ is the absolute tolerance,
   *  - \f$R_{tol}\f$ is the relative tolerance,
   *  - \f$\max\f$ is the pairwise maximum,
   *  - \f$u^n\f$ is the current solution, and
   *  - \f$u^{n+1}\f$ is the next solution.
   *  - \f$ u_i^n\f$ denotes component \f$i\f$ of time step \f$n\f$.
   *  - \f$ N\f$ denotes the number of unknowns
   */
  Scalar computeWRMSNorm(
      const Teuchos::RCP<const Thyra::VectorBase<Scalar>> &x,
      const Teuchos::RCP<const Thyra::VectorBase<Scalar>> &xNext,
      const Teuchos::RCP<const Thyra::VectorBase<Scalar>> &err);

  /** \brief Compute the error Norm.
   *
   *  The error norm is
   *  \f[
   *    e = \max_i ( u_i / ( A_{tol} + |u_i| R_{tol}))
   *  \f]
   *  where
   *  \f$A_{tol}\f$ is the absolute tolerance,
   *  \f$R_{tol}\f$ is the relative tolerance, and
   *  \f$\max_i\f$ is the maximum over all elements in \f$x\f$.
   */
  Scalar errorNorm(const Teuchos::RCP<const Thyra::VectorBase<Scalar>> &x);

  void setRelativeTolerance(const Scalar relTol) { relTol_ = relTol; }
  void setAbsoluteTolerance(const Scalar absTol) { abssTol_ = absTol; }

 protected:
  Scalar relTol_;
  Scalar abssTol_;
  Teuchos::RCP<Thyra::VectorBase<Scalar>> u_;
  Teuchos::RCP<Thyra::VectorBase<Scalar>> uNext_;
  Teuchos::RCP<Thyra::VectorBase<Scalar>> errorWeightVector_;
  Teuchos::RCP<Thyra::VectorBase<Scalar>> scratchVector_;
};

}  // namespace Tempus
#endif  // Tempus_Stepper_ErrorNorm_decl_hpp
