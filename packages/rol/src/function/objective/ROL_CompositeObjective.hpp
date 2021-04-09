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

#ifndef ROL_COMPOSITE_OBJECTIVE_H
#define ROL_COMPOSITE_OBJECTIVE_H

#include "ROL_StdObjective.hpp"

/** @ingroup func_group
    \class ROL::CompositeObjective
    \brief Provides the interface to evaluate composite objective functions.
*/


namespace ROL {

template<typename Real>
class CompositeObjective : public Objective<Real> {
private:
  const std::vector<Ptr<Objective<Real>>> obj_vec_;
  const Ptr<StdObjective<Real>> std_obj_;

  Ptr<std::vector<Real>> obj_value_, obj_grad_, obj_gv_, obj_hess_;
  Ptr<StdVector<Real>> obj_value_vec_, obj_grad_vec_, obj_gv_vec_, obj_hess_vec_;
  std::vector<Ptr<Vector<Real>>> vec_grad_, vec_hess_;

  bool isInitialized_, isValueComputed_, isGradientComputed_;

public:
  CompositeObjective(const std::vector<Ptr<Objective<Real>>> &obj_vec,
                     const Ptr<StdObjective<Real>> &std_obj);

  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) override;
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) override;
  Real value( const Vector<Real> &x, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;

// Definitions for parametrized (stochastic) objective functions
public:
  void setParameter(const std::vector<Real> &param) override;

private:
  void initialize(const Vector<Real> &x);
  void computeValue(const Vector<Real> &x, Real &tol);
  void computeGradient(const Vector<Real> &x, Real &tol);
  void computeHessVec(const Vector<Real> &v, const Vector<Real> &x, Real &tol);
};

} // namespace ROL

#include "ROL_CompositeObjective_Def.hpp"

#endif
