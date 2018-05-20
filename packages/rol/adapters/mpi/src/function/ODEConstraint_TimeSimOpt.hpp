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

#ifndef __ODEConstraint_TimeSimOpt_hpp__
#define __ODEConstraint_TimeSimOpt_hpp__

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_Constraint_TimeSimOpt.hpp"
#include "ROL_Vector.hpp"
#include "ROL_PartitionedVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_RandomVector.hpp"

typedef double RealT;

template <typename Real>
class ODE_Constraint : public ROL::Constraint_TimeSimOpt<Real> {
private:
  typedef ROL::StdVector<Real> VectorType;

  const std::vector<Real> & getVector( const ROL::Vector<Real> & x ) {
    return *dynamic_cast<const VectorType&>(x).getVector();
  }

  std::vector<Real> & getVector( ROL::Vector<Real>  & x ) {
    return *dynamic_cast<VectorType&>(x).getVector();
  }

  Real timestep_; 
  Real omega_; 

public:
   
  ODE_Constraint(double dt,double omega) : timestep_(dt), omega_(omega) {}

  virtual void value(ROL::Vector<Real> &c,
             const ROL::Vector<Real> &u_old,
             const ROL::Vector<Real> &u_new,
             const ROL::Vector<Real> &z,
             Real &tol) override
  {
    auto & c_data = getVector(c);
    auto & uo_data = getVector(u_old);
    auto & un_data = getVector(u_new);
    auto & z_data = getVector(z);

    // solving (u,v are states, z is control)
    // 
    //    u' = v + z, v'=-omega^2 * u, w' = sqrt(w)
    //    u(0) = 0
    //    v(0) = omega 
    //    w(0) = 1.0
    //
    //    u(t) = sin(omega*t)
    //    v(t) = omega*cos(omega*t)
    //    w(t) = 0.25*(-2.0*c1*t+c1**2+t**2), (init cond implies c1 = 2)
    //
    // using backward euler

    c_data[0] = un_data[0]-uo_data[0] - timestep_*(un_data[1] + z_data[0]);
    c_data[1] = un_data[1]-uo_data[1] + timestep_*omega_*omega_*un_data[0];
    c_data[2] = un_data[2]-uo_data[2] - timestep_*std::sqrt(un_data[2]);  // throw in some nonlinear term
  }

  virtual void solve(ROL::Vector<Real> &c,
                     const ROL::Vector<Real> &u_old, ROL::Vector<Real> &u_new,
                     const ROL::Vector<Real> &z,
                     Real &tol) override
  {
    auto & uo_data = getVector(u_old);
    auto & un_data = getVector(u_new);
    auto & z_data = getVector(z);

    Real a = omega_*omega_*timestep_;
    Real b = timestep_;
    Real gamma = (a*b+1.0);

    Real rhs_0 = uo_data[0]+timestep_*z_data[0];
    Real rhs_1 = uo_data[1];

    un_data[0] = (1.0 * rhs_0 +   b * rhs_1)/gamma;
    un_data[1] = ( -a * rhs_0 + 1.0 * rhs_1)/gamma;

    // Newton's method
    for(int i=0;i<5;i++) {
      // safety check
      TEUCHOS_ASSERT(un_data[2]>=0.0);

      double residual = un_data[2]-uo_data[2] - timestep_*std::sqrt(un_data[2]); 
      un_data[2] = un_data[2] - residual / (1.0-timestep_*(0.5/std::sqrt(un_data[2])));
    }

    value(c,u_old,u_new,z,tol);
  }

  virtual void applyJacobian_1_old(ROL::Vector<Real> &jv,
                                   const ROL::Vector<Real> &v_old,
                                   const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new,
                                   const ROL::Vector<Real> &z,
                                   Real &tol) override {
    auto & jv_data = getVector(jv);
    auto & vo_data = getVector(v_old);

    jv_data[0] = -vo_data[0];
    jv_data[1] = -vo_data[1];
    jv_data[2] = -vo_data[2];
  }

  virtual void applyJacobian_1_new(ROL::Vector<Real> &jv,
                                   const ROL::Vector<Real> &v_new,
                                   const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new,
                                   const ROL::Vector<Real> &z,
                                   Real &tol) override {
    auto & jv_data = getVector(jv);
    auto & vn_data = getVector(v_new);
    auto & un_data = getVector(u_new);

    jv_data[0] = vn_data[0] - timestep_*vn_data[1];
    jv_data[1] = vn_data[1] + timestep_*omega_*omega_*vn_data[0];
             // [      1,   -dt ]
             // [ dt*w*w,     1 ]
            
    TEUCHOS_ASSERT(un_data[2]>=0.0);
    jv_data[2] = (1.0 - timestep_*0.5/std::sqrt(un_data[2]))*vn_data[2];
  }

  virtual void applyInverseJacobian_1_new(ROL::Vector<Real> &ijv,
                                          const ROL::Vector<Real> &v_new,
                                          const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new,
                                          const ROL::Vector<Real> &z,
                                          Real &tol) override {
    auto & ijv_data = getVector(ijv);
    auto & v_data = getVector(v_new); 
    auto & un_data = getVector(u_new);
      
    Real a = omega_*omega_*timestep_;
    Real b = timestep_;
    Real gamma = 1.0/(a*b+1.0);

    ijv_data[0] = gamma*(1.0 * v_data[0] +   b * v_data[1]);
    ijv_data[1] = gamma*( -a * v_data[0] + 1.0 * v_data[1]);

    TEUCHOS_ASSERT(un_data[2]>=0.0);
    ijv_data[2] = v_data[2]/(1.0 - timestep_*0.5/std::sqrt(un_data[2]));
  }

  virtual void applyJacobian_2(ROL::Vector<Real> &jv,
                                   const ROL::Vector<Real> &v_new,
                                   const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new,
                                   const ROL::Vector<Real> &z,
                                   Real &tol) override {
    auto & jv_data = getVector(jv);
    auto & vn_data = getVector(v_new);

    jv_data[0] = -timestep_*vn_data[0];
    jv_data[1] = 0.0;
    jv_data[2] = 0.0;
  }

  virtual void applyAdjointJacobian_1_old(ROL::Vector<Real> &ajv_old,
                                      const ROL::Vector<Real> &dualv,
                                      const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new,
                                      const ROL::Vector<Real> &z,
                                      Real &tol) override {
    auto & ajv_data = getVector(ajv_old);
    auto & v_data = getVector(dualv);

    ajv_data[0] = -v_data[0];
    ajv_data[1] = -v_data[1];
    ajv_data[2] = -v_data[2];
  }

  virtual void applyAdjointJacobian_1_new(ROL::Vector<Real> &ajv_new,
                                      const ROL::Vector<Real> &dualv,
                                      const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new,
                                      const ROL::Vector<Real> &z,
                                      Real &tol) override {

    auto & ajv_data = getVector(ajv_new);
    auto & v_data = getVector(dualv);
    auto & un_data = getVector(u_new);

    ajv_data[0] = v_data[0] + timestep_*omega_*omega_*v_data[1];
    ajv_data[1] = v_data[1] - timestep_*v_data[0];

    TEUCHOS_ASSERT(un_data[2]>=0.0);
    ajv_data[2] = (1.0 - timestep_*0.5/std::sqrt(un_data[2]))*v_data[2];
  }

  virtual void applyAdjointJacobian_2_time(ROL::Vector<Real> &ajv,
                                      const ROL::Vector<Real> &dualv,
                                      const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new,
                                      const ROL::Vector<Real> &z,
                                      Real &tol) override {
    auto & ajv_data = getVector(ajv);
    auto & v_data = getVector(dualv);

    ajv_data[0] = -timestep_*v_data[0];
  }

  virtual void applyInverseAdjointJacobian_1_new(ROL::Vector<Real> &iajv,
                                                 const ROL::Vector<Real> &v_new,
                                                 const ROL::Vector<Real> &u_old, const ROL::Vector<Real> &u_new,
                                                 const ROL::Vector<Real> &z,
                                                 Real &tol) override {
    auto & iajv_data = getVector(iajv);
    auto & v_data = getVector(v_new); 
    auto & un_data = getVector(u_new);
      
    Real a = omega_*omega_*timestep_;
    Real b = timestep_;
    Real gamma = 1.0/(a*b+1.0);

    iajv_data[0] = gamma*(1.0 * v_data[0] -   a * v_data[1]);
    iajv_data[1] = gamma*(  b * v_data[0] + 1.0 * v_data[1]);

    TEUCHOS_ASSERT(un_data[2]>=0.0);
    iajv_data[2] = v_data[2]/(1.0 - timestep_*0.5/std::sqrt(un_data[2]));
  }

  virtual void applyAdjointHessian_11_old(ROL::Vector<Real> &ahwv_old,
                                          const ROL::Vector<Real> &w,
                                          const ROL::Vector<Real> &v_old,
                                          const ROL::Vector<Real> &u_old,
                                          const ROL::Vector<Real> &u_new,
                                          const ROL::Vector<Real> &z,
                                          Real &tol) override
  {
    auto & ahwv_data = getVector(ahwv_old);

    ahwv_data[0] = 0.0;
    ahwv_data[1] = 0.0;
    ahwv_data[2] = 0.0;
  }

  virtual void applyAdjointHessian_11_new(ROL::Vector<Real> &ahwv_new,
                                          const ROL::Vector<Real> &w,
                                          const ROL::Vector<Real> &v_new,
                                          const ROL::Vector<Real> &u_old,
                                          const ROL::Vector<Real> &u_new,
                                          const ROL::Vector<Real> &z,
                                          Real &tol) override
  {
    auto & ahwv_data = getVector(ahwv_new);
    auto & w_data = getVector(w); 
    auto & v_data = getVector(v_new); 
    auto & un_data = getVector(u_new);

    ahwv_data[0] = 0.0;
    ahwv_data[1] = 0.0;
    ahwv_data[2] = 0.25*timestep_*w_data[2]*v_data[2]/std::pow(un_data[2],1.5);
  }

};

#endif 
