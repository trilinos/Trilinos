// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once

#include "ROL_StdVector.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_Constraint_SimOpt.hpp"

namespace Rocket {

std::vector<double>& getVector( ROL::Vector<double>& x ) {
  return *(dynamic_cast<ROL::StdVector<double>&>(x).getVector());
}

const std::vector<double>& getVector( const ROL::Vector<double>& x ) {
  return *(dynamic_cast<const ROL::StdVector<double>&>(x).getVector());
}


class Objective : public ROL::Objective_SimOpt<double> {
private:
  using V = ROL::Vector<double>;

  int N;
  double T, dt, mt;
  double htarg, alpha;
  const ROL::Ptr<const V>  w; // Trapezoidal weights

public:
  
  Objective(  int N_, double T_, double mt_, double htarg_, double alpha_, 
              const ROL::Ptr<const V>& w_ ) :
    N(N_), T(T_), dt(T/N), mt(mt_), htarg(htarg_), alpha(alpha_), w(w_) {
  }

  double value( const V& u, const V& z, double& tol ) {
    return 0.5*std::pow(htarg-w->dot(u),2) + 0.5*alpha*dt*z.dot(z); 
  }

  void gradient_1( V& g, const V& u, const V& z, double& tol ) {
    g.set(*w);   g.scale(w->dot(u)-htarg);
  }

  void gradient_2( V& g, const V& u, const V& z, double& tol ) {
    g.set(z);    g.scale(alpha*dt);
  }

  void hessVec_11( V& hv, const V& v, const V& u, const V& z, double& tol ) {
    hv.set(*w);  hv.scale(w->dot(v));
  }

  void hessVec_12( V& hv, const V& v, const V& u, const V& z, double& tol ) {
    hv.zero();
  }

  void hessVec_21( V& hv, const V& v, const V& u, const V& z, double& tol ) {
    hv.zero(); 
  }

  void hessVec_22( V& hv, const V& v, const V& u, const V& z, double& tol ) {
    hv.set(v);  hv.scale(alpha*dt);
  }
}; // Objective




class Constraint : public ROL::Constraint_SimOpt<double> {
private:
  using V = ROL::Vector<double>;

  int N;
  double T, dt, gdt, mt, mf, ve;

  std::vector<double> mass;

public:

  Constraint( int N_, double T_, double mt_, 
    double mf_, double ve_, double g_ ) : N(N_), T(T_), 
    dt(T/N), gdt(g_*dt), mt(mt_),
    mf(mf_), ve(ve_), mass(N) {
      mass[0] = mt;
    }

  void update_2( const V& z, bool flag = true, int iter = -1 ) {
    auto& zs = getVector(z);
  
    mass[0] = mt - dt*zs[0];
    for( int i=1; i<N; ++i ) 
      mass[i] = mass[i-1] - dt*zs[i];
  }

  void solve( V& c, V& u, const V& z, double& tol ) {

    auto& us = getVector(u); 

    us[0] = -ve*std::log(mass[0]/mt) - gdt;
    for( int i=1; i<N; ++i ) 
      us[i] = us[i-1] - ve*std::log(mass[i]/mass[i-1]) - gdt;

    value(c,u,z,tol);
  }

  void value( V& c, const V& u, const V&z, double &tol ) {

    auto& cs = getVector(c); auto& us = getVector(u); 

    cs[0] = us[0] + ve*std::log(mass[0]/mt) + gdt;

    for( int i=1; i<N; ++i ) 
      cs[i] = us[i] - us[i-1] + ve*std::log(mass[i]/mass[i-1]) + gdt;
  }

  void applyJacobian_1( V& jv, const V& v, const V& u, const V& z, double& tol ) {

    auto& jvs = getVector(jv); auto& vs =  getVector(v);
 
    jvs[0] = vs[0];
    for( int i=1; i<N; ++i ) jvs[i] = vs[i] - vs[i-1];
  }

   void applyAdjointJacobian_1( V& ajv, const V& v, const V& u, const V& z, double& tol ) {

     auto& ajvs = getVector(ajv);  auto& vs = getVector(v);
 
     ajvs[N-1] = vs[N-1];
     for( int i=N-2; i>=0; --i ) ajvs[i] = vs[i] - vs[i+1];
   }

  void applyInverseJacobian_1( V& ijv, const V& v, const V& u, const V& z, double &tol ) {

     auto& ijvs = getVector(ijv);  auto& vs = getVector(v);

     ijvs[0] = vs[0];
     for( int i=1; i<N; ++i ) ijvs[i] = ijvs[i-1] + vs[i];
  }

  void applyInverseAdjointJacobian_1( V& ijv, const V& v, const V& u, const V& z, double &tol ) {

     auto& ijvs = getVector(ijv);  auto& vs = getVector(v);

     ijvs[N-1] = vs[N-1];
     for( int i=N-2; i>=0; --i ) ijvs[i] = ijvs[i+1] + vs[i];
  }
 
  void applyJacobian_2( V& jv, const V& v, const V& u, const V& z, double& tol ) {

    auto& jvs = getVector(jv);  auto& vs  = getVector(v);

    double q{-ve*dt*vs[0]};
    jvs[0] = q/mass[0];
    for( int i=1; i<N; ++i ) {
      jvs[i] = -q/mass[i-1];  q -= ve*dt*vs[i];
      jvs[i] += q/mass[i];
    }
  }

  void applyAdjointJacobian_2( V& ajv, const V& v, const V& u, const V& z, double& tol ) { 
     auto& ajvs = getVector(ajv);  auto& vs   = getVector(v);
 
     ajvs[N-1] = -ve*dt*vs[N-1]/mass[N-1];

     for( int i=N-2; i>=0; --i ) 
       ajvs[i] = ajvs[i+1]-ve*dt*(vs[i]-vs[i+1])/mass[i];
  }  
};

} // namespace Rocket
