// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package Copyright (2014)
//               Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive license
// for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the contributors may
// be used to endorse or promote products derived from this software without
// specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
// EVENT SHALL SANDIA CORPORATION OR THE CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
// THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers: Drew Kouri   (dpkouri@sandia.gov) and
// Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER


#ifndef ROL_OED_ITRACE_HOM_OBJECTIVE_HPP
#define ROL_OED_ITRACE_HOM_OBJECTIVE_HPP

#include "ROL_OED_HomObjectiveBase.hpp"
#include "ROL_SampleGenerator.hpp"

namespace ROL {
namespace OED {
namespace Hom {

template<typename Real>
class Itrace_Objective : public ObjectiveBase<Real,std::vector<Real>> {
private:
  const Ptr<SampleGenerator<Real>> sampler_;
  const std::vector<Real> weight_;
  const Ptr<Vector<Real>> adjoint_, rhs_, sens_;
  const Ptr<VectorController<Real,int>> adjointStore_;
  std::vector<Ptr<Vector<Real>>> b_;
  Ptr<Vector<Real>> Xa_;
  const bool storage_;

  using ObjectiveBase<Real,std::vector<Real>>::setConstraint;
  using ObjectiveBase<Real,std::vector<Real>>::setObjective;
  using ObjectiveBase<Real,std::vector<Real>>::setStorage;
  using ObjectiveBase<Real,std::vector<Real>>::initialize;
  using ObjectiveBase<Real,std::vector<Real>>::getConstraint;
  using ObjectiveBase<Real,std::vector<Real>>::getObjective;
  using ObjectiveBase<Real,std::vector<Real>>::getState;
  using ObjectiveBase<Real,std::vector<Real>>::getStateSens;
  using ObjectiveBase<Real,std::vector<Real>>::solve_state_equation;
  using ObjectiveBase<Real,std::vector<Real>>::solve_state_sensitivity;

  void solveAdjointEquation(Vector<Real> &adjoint, const Vector<Real> &u,
                            const Vector<Real> &z, int i, Real &tol);

public:
  Itrace_Objective( const Ptr<BilinearConstraint<Real>> &con,
                    const Ptr<LinearObjective<Real>>  &obj,
                    const Ptr<Vector<Real>>           &state,
                    const Ptr<SampleGenerator<Real>>  &sampler,
                    const std::vector<Real>           &weight,
                    bool storage = true);

  void update( const Vector<Real> &z,
               UpdateType type,
               int iter = -1 ) override;
  Real value( const Vector<Real> &z, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) override; 
};

} // END Hom Namespace
} // END OED Namespace
} // END ROL Namespace

#include "ROL_OED_Itrace_HomObjective_Def.hpp"

#endif
