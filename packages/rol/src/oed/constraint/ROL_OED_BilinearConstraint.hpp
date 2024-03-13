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

#ifndef ROL_OED_BILINEARCONSTRAINT_HPP
#define ROL_OED_BILINEARCONSTRAINT_HPP

#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_ProbabilityVector.hpp"
#include "ROL_OED_MomentOperator.hpp"
#include "ROL_OED_TraceSampler.hpp"
#include "ROL_OED_ProfiledClass.hpp"

namespace ROL {
namespace OED {

template<typename Real>
class BilinearConstraint : public Constraint_SimOpt<Real>, public ProfiledClass<Real,std::string> {
private:
  const Ptr<Factors<Real>>           factors_;      // Factors for regression
  const Ptr<MomentOperator<Real>>    M_;            // Specialized Moment Operator
  Ptr<TraceSampler<Real>>            traceSampler_; // Trace Sampler for A-, D-, I-optimality
  const std::string                  type_;         // Optimality type
  const Ptr<Vector<Real>>            g_;            // Storage for A-, C-, I-optimality
  Ptr<Vector<Real>>                  p_;            // Storage for A-, C-, I-optimality
  bool                               useTrace_;     // Use trace form for I-optimality
  bool                               isPinit_;

  using ProfiledClass<Real,std::string>::startTimer;
  using ProfiledClass<Real,std::string>::stopTimer;

  void computeG(Vector<Real> &g);

public:
  // C optimality
  BilinearConstraint(const Ptr<Factors<Real>>        &factors,
                     const Ptr<MomentOperator<Real>> &M,
                     const Ptr<Vector<Real>>         &c);

  // A, D, I, and R optimality
  BilinearConstraint(const Ptr<Factors<Real>>        &factors,
                     const Ptr<MomentOperator<Real>> &M,
                     const std::string               &type = "I",
                     const Ptr<TraceSampler<Real>>   &traceSampler = nullPtr);

  void update_2(const Vector<Real> &z, UpdateType type, int iter = -1 ) override;
  void value(Vector<Real> &c,const Vector<Real> &u,const Vector<Real> &z,Real &tol) override;
  void solve(Vector<Real> &c,Vector<Real> &u,const Vector<Real> &z,Real &tol) override;
  void applyJacobian_1(Vector<Real> &jv,const Vector<Real> &v,const Vector<Real> &u,
                       const Vector<Real> &z,Real &tol) override;
  void applyJacobian_2(Vector<Real> &jv,const Vector<Real> &v,const Vector<Real> &u,
                       const Vector<Real> &z,Real &tol) override;
  void applyInverseJacobian_1(Vector<Real> &ijv,const Vector<Real> &v,
                              const Vector<Real> &u,const Vector<Real> &z,Real &tol) override;
  void applyAdjointJacobian_1(Vector<Real> &ajv,const Vector<Real> &v,
                              const Vector<Real> &u,const Vector<Real> &z,Real &tol) override;
  void applyAdjointJacobian_2(Vector<Real> &ajv,const Vector<Real> &v,
                              const Vector<Real> &u,const Vector<Real> &z,Real &tol) override;
  void applyInverseAdjointJacobian_1(Vector<Real> &iajv,const Vector<Real> &v,
                                     const Vector<Real> &u,const Vector<Real> &z,Real &tol) override;
  void applyAdjointHessian_11(Vector<Real> &ahwv,const Vector<Real> &w,const Vector<Real> &v,
                              const Vector<Real> &u,const Vector<Real> &z,Real &tol) override;
  void applyAdjointHessian_12(Vector<Real> &ahwv,const Vector<Real> &w,const Vector<Real> &v,
                              const Vector<Real> &u,const Vector<Real> &z,Real &tol) override;
  void applyAdjointHessian_21(Vector<Real> &ahwv,const Vector<Real> &w,const Vector<Real> &v,
                              const Vector<Real> &u,const Vector<Real> &z,Real &tol) override;
  void applyAdjointHessian_22(Vector<Real> &ahwv,const Vector<Real> &w,const Vector<Real> &v,
                              const Vector<Real> &u,const Vector<Real> &z,Real &tol) override;

  void getFactor(Vector<Real> &F, int k) const;
  void getFactor(Vector<Real> &F, const std::vector<Real> &param) const;
  Real getNoise(int k) const;
  void getTraceSample(Vector<Real> &g, const std::vector<Real> &param) const;
  // void sumAll(Real *in, Real *out, int size) const;
  Real logDeterminant(const Vector<Real> &z);
}; // class BilinearConstraint

} // End OED Namespace
} // End ROL Namespace

#include "ROL_OED_BilinearConstraint_Def.hpp"

#endif
