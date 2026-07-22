// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEU_TRUSTREGIONALGORITHM_H
#define ROL_TYPEU_TRUSTREGIONALGORITHM_H

#include "ROL_TypeU_Algorithm.hpp"
#include "ROL_TrustRegion_U_Types.hpp"
#include "ROL_TrustRegion_U.hpp"
#include "ROL_TrustRegionUtilities.hpp"
#include "ROL_Secant.hpp"

/** \class ROL::TypeU::TrustRegionAlgorithm
    \brief Provides an interface to run trust-region methods for unconstrained
           optimization algorithms.
*/

namespace ROL {
namespace TypeU {

template<typename Real>
class TrustRegionAlgorithm : public Algorithm<Real> {
private:
  // TRUST REGION INFORMATION
  Ptr<TrustRegion_U<Real>>      solver_; ///< Container for trust-region solver object.
  Ptr<TrustRegionModel_U<Real>> model_;  ///< Container for trust-region model.
  ETrustRegionU                 etr_;    ///< Trust-region subproblem solver type.
  Real                          delMax_; ///< Maximum trust-region radius.
  Real                          eta0_;   ///< Step acceptance threshold.
  Real                          eta1_;   ///< Radius decrease threshold.
  Real                          eta2_;   ///< Radius increase threshold.
  Real                          gamma0_; ///< Radius decrease rate (negative rho).
  Real                          gamma1_; ///< Radius decrease rate (positive rho).
  Real                          gamma2_; ///< Radius increase rate.
  Real                          TRsafe_; ///< Safeguard size for numerically evaluating ratio.
  Real                          eps_;    ///< Safeguard for numerically evaluating ratio.
  TRUtils::ETRFlag              TRflag_; ///< Trust-region exit flag.
  int                           SPflag_; ///< Subproblem solver termination flag.
  int                           SPiter_; ///< Subproblem solver iteration count.

  // NONMONOTONE INFORMATION
  bool useNM_;
  int  NMstorage_;

  // SECANT INFORMATION
  ESecant esec_; ///< Secant type.
  bool useSecantPrecond_;
  bool useSecantHessVec_;

  // INEXACT COMPUTATION PARAMETERS
  std::vector<bool> useInexact_; ///< Flags for inexact (0) objective function, (1) gradient, (2) Hessian.
  Real              scale0_;     ///< Scale for inexact gradient computation.
  Real              scale1_;     ///< Scale for inexact gradient computation.
  Real scale_, omega_, force_, forceFactor_;
  int updateIter_;
  Real gtol_;

  // VERBOSITY SETTING
  int verbosity_;    ///< Print additional information to screen if > 0.
  bool printHeader_; ///< Print header at every iteration.

  using Algorithm<Real>::state_;
  using Algorithm<Real>::status_;
  using Algorithm<Real>::initialize;

public:

  TrustRegionAlgorithm( ParameterList &parlist,
    const Ptr<Secant<Real>> &secant = nullPtr );

  using Algorithm<Real>::run;
  void run( Vector<Real>       &x,
            const Vector<Real> &g, 
            Objective<Real>    &obj,
            std::ostream       &outStream = std::cout) override;

  void writeHeader( std::ostream& os ) const override;

  void writeName( std::ostream& os ) const override;
  
  void writeOutput( std::ostream& os, const bool print_header = false ) const override;

private:

  void initialize(const Vector<Real> &x, const Vector<Real> &g, Vector<Real> &Bg,
                  Objective<Real> &obj, std::ostream &outStream = std::cout);

  Real computeValue(const Vector<Real> &x, Objective<Real> &obj, Real pRed);

  /** \brief Compute gradient to iteratively satisfy inexactness condition.

      This function attempts to ensure that the inexact gradient condition,
      \f[
         \|g_k-\nabla J(x_k)\|_{\mathcal{X}} \le \kappa_1\min\{\,\|g_k\|_{\mathcal{X}},\,\Delta_k\,\},
      \f]
      is satisfied.  This function works under the assumption that the gradient function returns 
      a gradient approximation which satisfies the error tolerance prescribed by the tol input 
      parameter.  
      @param[in]      x          is the current optimization variable.
      @param[in]      obj        is the objective function.
  */
  void computeGradient(const Vector<Real> &x, Objective<Real> &obj, bool accept);

}; // class ROL::TrustRegionAlgorithm
} // namespace TypeU
} // namespace ROL

#include "ROL_TypeU_TrustRegionAlgorithm_Def.hpp"

#endif
