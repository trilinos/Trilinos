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

#ifndef ROL_OED_QUADRATIC_OBJECTIVE_HPP
#define ROL_OED_QUADRATIC_OBJECTIVE_HPP

#include "ROL_Objective_SimOpt.hpp"
#include "ROL_SampledVector.hpp"
#include "ROL_OED_BilinearConstraint.hpp"
#include "ROL_OED_ProfiledClass.hpp"

namespace ROL {
namespace OED {

template<typename Real>
class QuadraticObjective : public Objective_SimOpt<Real>, public ProfiledClass<Real,std::string> {
private:
  const Ptr<BilinearConstraint<Real>> M_;
  Ptr<Vector<Real>>                 g_;
  bool                   isInit_, isG_;

  using ProfiledClass<Real,std::string>::startTimer;
  using ProfiledClass<Real,std::string>::stopTimer;

  void initialize(const Vector<Real> &x);
  void apply(const Vector<Real> &u, const Vector<Real> &z, Real &tol);

public:
  QuadraticObjective(const Ptr<BilinearConstraint<Real>> &M);

  Ptr<BilinearConstraint<Real>> getM() const;

  void update( const Vector<Real> &u, const Vector<Real> &z, UpdateType type, int iter = -1 ) override;
  Real value( const Vector<Real> &u, const Vector<Real> &z, Real &tol ) override;
  void gradient_1( Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) override;
  void gradient_2( Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) override;
  void hessVec_11( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) override;
  void hessVec_12( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) override;
  void hessVec_21( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) override;
  void hessVec_22( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) override;
  void setParameter( const std::vector<Real> &param ) override;

}; // class QuadraticObjective

} // End OED Namespace
} // End ROL Namespace

#include "ROL_OED_QuadraticObjective_Def.hpp"

#endif
