// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#pragma once
#ifndef ROL2_TYPEU_TRUSTREGIONALGORITHM_DECL_H
#define ROL2_TYPEU_TRUSTREGIONALGORITHM_DECL_H

/** \class ROL2::TypeU::TrustRegionAlgorithm
    \brief Provides an interface to run trust-region methods for unconstrained
           optimization algorithms.
*/

namespace ROL2 {
namespace TypeU {

template<typename Real>
class TrustRegionAlgorithm : public ROL2::TypeU::Algorithm<Real> {
public:

  TrustRegionAlgorithm(       ParameterList&     parlist,
                        const Ptr<Secant<Real>>& secant = nullPtr );

  void run(       Vector<Real>&    x,
            const Vector<Real>&    g, 
                  Objective<Real>& obj,
                  std::ostream&    outStream = std::cout) override;

  void writeHeader( std::ostream& = std::cout ) const override;

  void writeName( std::ostream& = std::cout ) const override;
  
  void writeOutput( std::ostream& = std::cout, 
                    bool print_header = false ) const override;

  void setModel(  const Ptr<TrustRegionModel<Real>>& );
  void setSolver( const Ptr<TrustRegion<Real>>&      );
  
  using Algorithm<Real>::getState;
  using Algorithm<Real>::getStatus;
 
  const TrustRegionModel<Real>& getModel()  const;
  const Secant<Real>&           getSecant() const;
  const TrustRegion<Real>&      getSolver() const;
   
protected:
  TrustRegionModel<Real>& getModel();
  TrustRegion<Real>&      getSolver();

private:
  // TRUST REGION INFORMATION
  Ptr<TrustRegion<Real>>            solver_; ///< Container for trust-region solver object.
  Ptr<TrustRegionModel<Real>>       model_;  ///< Container for trust-region model.
  typename TrustRegion<Real>::Type  type_;   ///< Trust-region subproblem solver type.
  typename TrustRegion<Real>::Flag  flag_;   ///< Trust-region exit flag.
  Real                              delMax_; ///< Maximum trust-region radius.
  Real                              eta0_;   ///< Step acceptance threshold.
  Real                              eta1_;   ///< Radius decrease threshold.
  Real                              eta2_;   ///< Radius increase threshold.
  Real                              gamma0_; ///< Radius decrease rate (negative rho).
  Real                              gamma1_; ///< Radius decrease rate (positive rho).
  Real                              gamma2_; ///< Radius increase rate.
  Real                              TRsafe_; ///< Safeguard size for numerically evaluating ratio.
  Real                              eps_;    ///< Safeguard for numerically evaluating ratio.
  int                               SPflag_; ///< Subproblem solver termination flag.
  int                               SPiter_; ///< Subproblem solver iteration count.

  // SECANT INFORMATION
  typename Secant<Real>::Type secantType_;
  bool useSecantPrecond_;
  bool useSecantHessVec_;

  // INEXACT COMPUTATION PARAMETERS
  std::vector<bool> useInexact_; ///< Flags for inexact (0) objective function, (1) gradient, (2) Hessian.
  Real              scale0_;     ///< Scale for inexact gradient computation.
  Real              scale1_;     ///< Scale for inexact gradient computation.
  Real scale_, omega_, force_, forceFactor_;
  int updateIter_;

  // VERBOSITY SETTING
  int verbosity_;    ///< Print additional information to screen if > 0.
  bool printHeader_; ///< Print header at every iteration.



  void initialize( const Vector<Real>&    x, 
                   const Vector<Real>&    g,  
                         Vector<Real>&    Bg,
                         Objective<Real>& obj, 
                         std::ostream&    outStream = std::cout );

  Real computeValue( const Vector<Real>&    x, 
                           Objective<Real>& obj, 
                           Real             pRed );

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
  void computeGradient( const Vector<Real>&    x, 
                              Objective<Real>& obj );

}; // class ROL2::TypeU::TrustRegionAlgorithm


} // namespace TypeU
} // namespace ROL

#endif // ROL2_TYPEU_TRUSTREGIONALGORITHM_DECL_H


