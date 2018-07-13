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
#ifndef ROL_DYNAMICCONSTRAINT_HPP
#define ROL_DYNAMICCONSTRAINT_HPP

#include "ROL_DynamicFunction.hpp"
#include "ROL_TimeStamp.hpp"

/** @ingroup dynamic_group
    \class ROL::DynamicConstraint
    \brief Defines the time-dependent constraint operator interface for 
           simulation-based optimization.

    This constraint interface inherits from ROL_Constraint_SimOpt. Though 
    the interface takes two simulation space vectors from spaces
    \f$\mathcal{U_o}\times\mathcal{U_n}\f$. The space \f$\mathcal{U_o}\f$ is 
    ``old'' information that accounts for the initial condition on the time 
     interval. The space \f$\mathcal{U_n}\f$ is the ``new'' variables that can 
    be determined by satisfying constraints in the form
    \f[
      c(u_o,u_n,z,t_o,t_n) = 0 \,.
    \f]

*/


namespace ROL {

template<typename Real>
class DynamicConstraint : public DynamicFunction<Real> {
public:

  using V  = Vector<Real>;
  using PV = PartitionedVector<Real>;
  using TS = TimeStamp<Real>;
 
  virtual ~DynamicConstraint() {}

  virtual void update( const V& uo, const V& un, const V& z ) {
    update_uo( uo );
    update_un( un );
    update_z( z );
  }

  virtual void update_uo( const V& uo ) { }
  virtual void update_un( const V& un ) { }
  virtual void update_z( const V& z ) { }

  virtual void value( V& c, const V& uo, const V& un, 
                      const V& z, const TS& ts ) const = 0;

  virtual void solve( V& c, V& uo, const V& un, 
                      const V& z, const TS& ts ) const {}

  //----------------------------------------------------------------------------
  // Partial Jacobians
  virtual void applyJacobian_uo( V& jv, const V& vo, const V& uo, 
                                    const V& un, const V& z, 
                                    const TS& ts ) const {}

  virtual void applyJacobian_un( V& jv, const V& vn, const V& uo, 
                                    const V& un, const V& z, 
                                    const TS& ts ) const {}

  virtual void applyJacobian_z( V& jv, const V& vz, const V& uo, 
                                const V& un, const V& z, 
                                const TS& ts ) const {}

  //----------------------------------------------------------------------------
  // Adjoint partial Jacobians
  virtual void applyAdjointJacobian_uo( V& ajv, const V& dualv, const V& uo, 
                                           const V& un, const V& z, 
                                           const TS& ts ) const {}

  virtual void applyAdjointJacobian_un( V& ajv, const V& dualv, const V& uo, 
                                           const V& un, const V& z, 
                                           const TS& ts ) const {}

  virtual void applyAdjointJacobian_z( V& ajv, const V& dualv, const V& uo, 
                                       const V& un, const V& z, 
                                       const TS& ts ) const {}

  //----------------------------------------------------------------------------
  // Inverses
  virtual void applyInverseJacobian_un( V& ijv, const V& vn, const V& uo, 
                                           const V& un, const V& z, 
                                           const TS& ts ) const {}
    
  virtual void applyAdjointInverseJacobian_un( V& iajv, const V& vn, const V& uo, 
                                                  const V& un, const V& z, 
                                                  const TS& ts ) const {}

  //----------------------------------------------------------------------------
  // Adjoint Hessian components
  virtual void applyAdjointHessian_un_un( V& ahwv, const V& wn, const V& vn,
                                          const V& uo, const V& un, 
                                          const V& z, const TS& ts ) const {
    ahwv.zero();
  }

  virtual void applyAdjointHessian_un_uo( V& ahwv, const V& w, const V& vn,
                                          const V& uo, const V& un, 
                                          const V& z, const TS& ts ) const {
    ahwv.zero();
  }

  virtual void applyAdjointHessian_un_z( V& ahwv, const V& w, const V& vn,
                                          const V& uo, const V& un, 
                                          const V& z, const TS& ts ) const {
    ahwv.zero();
  }

  virtual void applyAdjointHessian_uo_un( V& ahwv, const V& w, const V& vo,
                                          const V& uo, const V& un, 
                                          const V& z, const TS& ts ) const {
    ahwv.zero();
  }

// This should be zero 
//  virtual void applyAdjointHessian_uo_uo( V& ahwv, const V& w, const V& v,
//                                          const V& uo, const V& un, 
//                                          const V& z, const TS& ts ) const;

  virtual void applyAdjointHessian_uo_z( V& ahwv, const V& w, const V& vo,
                                         const V& uo, const V& un, 
                                         const V& z, const TS& ts ) const {
    ahwv.zero();
  }

  virtual void applyAdjointHessian_z_un( V& ahwv, const V& w, const V& vz,
                                         const V& uo, const V& un, 
                                         const V& z, const TS& ts ) const {
    ahwv.zero();
  }

  virtual void applyAdjointHessian_z_uo( V& ahwv, const V& w, const V& vz,
                                         const V& uo, const V& un, 
                                         const V& z, const TS& ts ) const {
    ahwv.zero();
  }
  
  virtual void applyAdjointHessian_z_z( V& ahwv, const V& w, const V& vz,
                                        const V& uo, const V& un, 
                                        const V& z, const TS& ts ) const {
    ahwv.zero();
  }

}; // DynamicConstraint


} // namespace ROL


#endif // ROL_DYNAMICCONSTRAINT_HPP

