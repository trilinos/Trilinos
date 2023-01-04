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

#ifndef ROL_SERIALSTATIONARYCONTROLSOBJECTIVE_HPP
#define ROL_SERIALSTATIONARYCONTROLSOBJECTIVE_HPP

#include "ROL_Ptr.hpp"
#include "ROL_SerialObjective.hpp"
#include "ROL_PartitionedVector.hpp"

namespace ROL
{

/** @ingroup func_group
    \class ROL::SerialStationaryControlsObjective
    \brief Wrapper for SerialObjective for when the controls are stationary
           (i.e., not time-dependent).
*/
template <typename Real>
class SerialStationaryControlsObjective : public Objective_SimOpt<Real>
{
public:
  using size_type = typename std::vector<Real>::size_type;

  SerialStationaryControlsObjective(
      const Ptr<SerialObjective<Real>> &serial_obj,
      const Ptr<Vector<Real>> &z,
      const size_type Nt)
      : serial_obj_(serial_obj), Nt_(Nt)
  {
    z_dyn_ = PartitionedVector<Real>::create(*z, Nt_);
    v_dyn_ = PartitionedVector<Real>::create(*z, Nt_);
    g_dyn_ = PartitionedVector<Real>::create(z->dual(), Nt_);
  }

  virtual ~SerialStationaryControlsObjective() {}

  Real
  value(const Vector<Real> &u, const Vector<Real> &z, Real &tol)
  {
    for (size_type i = 0; i < Nt_; ++i)
      (*z_dyn_)[i].set(z);
    Real val = serial_obj_->value(u, *z_dyn_, tol);
    return val;
  }

  void
  gradient_1(Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol)
  {
    for (size_type i = 0; i < Nt_; ++i)
      (*z_dyn_)[i].set(z);
    serial_obj_->gradient_1(g, u, *z_dyn_, tol);
  }

  void
  gradient_2(Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol)
  {
    for (size_type i = 0; i < Nt_; ++i)
      (*z_dyn_)[i].set(z);
    serial_obj_->gradient_2(*g_dyn_, u, *z_dyn_, tol);
    g.zero();
    for (size_type i = 0; i < Nt_; ++i)
      g.axpy(1.0, (*g_dyn_)[i]);
  }

  // Hack in finite differences for hessVec_??.  We really should use the
  // commented out implementations later in the file, but SerialObjective()
  // doesn't implement all of the hessVec_?? functions.
  //
  // Furthermore, the finite difference perturbation size used in the
  // default implementations in Objective_SimOpt does not work for some
  // applications.
  void
  hessVec_11(Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol)
  {
    Real gtol = std::sqrt(ROL_EPSILON<Real>());
    // Compute step length
    Real h = tol;
    // if (v.norm() > std::sqrt(ROL_EPSILON<Real>())) {
    h = std::max(1.0, u.norm() / v.norm()) * tol;
    //}
    // Evaluate gradient of first component at (u+hv,z)
    ROL::Ptr<Vector<Real>> unew = u.clone();
    unew->set(u);
    unew->axpy(h, v);
    this->update(*unew, z);
    hv.zero();
    this->gradient_1(hv, *unew, z, gtol);
    // Evaluate gradient of first component at (u,z)
    ROL::Ptr<Vector<Real>> g = hv.clone();
    this->update(u, z);
    this->gradient_1(*g, u, z, gtol);
    // Compute Newton quotient
    hv.axpy(-1.0, *g);
    hv.scale(1.0 / h);
  }

  void
  hessVec_12(Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol)
  {
    // hv.zero();
    Real gtol = std::sqrt(ROL_EPSILON<Real>());
    // Compute step length
    Real h = tol;
    // if (v.norm() > std::sqrt(ROL_EPSILON<Real>())) {
    h = std::max(1.0, z.norm() / v.norm()) * tol;
    //}
    // Evaluate gradient of first component at (u,z+hv)
    ROL::Ptr<Vector<Real>> znew = z.clone();
    znew->set(z);
    znew->axpy(h, v);
    this->update(u, *znew);
    hv.zero();
    this->gradient_1(hv, u, *znew, gtol);
    // Evaluate gradient of first component at (u,z)
    ROL::Ptr<Vector<Real>> g = hv.clone();
    this->update(u, z);
    this->gradient_1(*g, u, z, gtol);
    // Compute Newton quotient
    hv.axpy(-1.0, *g);
    hv.scale(1.0 / h);
  }

  void
  hessVec_21(Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol)
  {
    // hv.zero();
    Real gtol = std::sqrt(ROL_EPSILON<Real>());
    // Compute step length
    Real h = tol;
    // if (v.norm() > std::sqrt(ROL_EPSILON<Real>())) {
    h = std::max(1.0, u.norm() / v.norm()) * tol;
    //}
    // Evaluate gradient of first component at (u+hv,z)
    ROL::Ptr<Vector<Real>> unew = u.clone();
    unew->set(u);
    unew->axpy(h, v);
    this->update(*unew, z);
    hv.zero();
    this->gradient_2(hv, *unew, z, gtol);
    // Evaluate gradient of first component at (u,z)
    ROL::Ptr<Vector<Real>> g = hv.clone();
    this->update(u, z);
    this->gradient_2(*g, u, z, gtol);
    // Compute Newton quotient
    hv.axpy(-1.0, *g);
    hv.scale(1.0 / h);
  }

  void
  hessVec_22(Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol)
  {
    // hv.zero();
    Real gtol = std::sqrt(ROL_EPSILON<Real>());
    // Compute step length
    Real h = tol;
    // if (v.norm() > std::sqrt(ROL_EPSILON<Real>())) {
    h = std::max(1.0, z.norm() / v.norm()) * tol;
    //}
    // Evaluate gradient of first component at (u,z+hv)
    ROL::Ptr<Vector<Real>> znew = z.clone();
    znew->set(z);
    znew->axpy(h, v);
    this->update(u, *znew);
    hv.zero();
    this->gradient_2(hv, u, *znew, gtol);
    // Evaluate gradient of first component at (u,z)
    ROL::Ptr<Vector<Real>> g = hv.clone();
    this->update(u, z);
    this->gradient_2(*g, u, z, gtol);
    // Compute Newton quotient
    hv.axpy(-1.0, *g);
    hv.scale(1.0 / h);
  }

  /*
  // The proper implementation of hessVec_??
  void hessVec_11(Vector<Real>& hv, const Vector<Real>& v,
                  const Vector<Real>& u, const Vector<Real>& z,
                  Real &tol)
  {
    for (size_type i=0; i<Nt_; ++i)
      (*z_dyn_)[i].set(z);
    serial_obj_->hessVec_11(hv, v, u, *z_dyn_, tol);
  }

  void hessVec_12(Vector<Real>& hv, const Vector<Real>& v,
                  const Vector<Real>& u, const Vector<Real>& z,
                  Real &tol)
  {
    for (size_type i=0; i<Nt_; ++i) {
      (*z_dyn_)[i].set(z);
      (*v_dyn_)[i].set(v);
    }
    serial_obj_->hessVec_12(hv, *v_dyn_, u, *z_dyn_, tol);
  }

  void hessVec_21(Vector<Real>& hv, const Vector<Real>& v,
                  const Vector<Real>& u, const Vector<Real>& z,
                  Real &tol)
  {
    for (size_type i=0; i<Nt_; ++i)
      (*z_dyn_)[i].set(z);
    serial_obj_->hessVec_21(*g_dyn_, v, u, *z_dyn_, tol);
    hv.zero();
    for (size_type i=0; i<Nt_; ++i)
      hv.axpy(1.0, (*g_dyn_)[i]);
  }

  void hessVec_22(Vector<Real>& hv, const Vector<Real>& v,
                  const Vector<Real>& u, const Vector<Real>& z,
                  Real &tol)
  {
    for (size_type i=0; i<Nt_; ++i) {
      (*z_dyn_)[i].set(z);
      (*v_dyn_)[i].set(v);
    }
    serial_obj_->hessVec_22(*g_dyn_, *v_dyn_, u, *z_dyn_, tol);
    hv.zero();
    for (size_type i=0; i<Nt_; ++i)
      hv.axpy(1.0, (*g_dyn_)[i]);
  }
  */

  // Bring in hidden overloads
  using Objective_SimOpt<Real>::value;

private:
  const Ptr<SerialObjective<Real>> serial_obj_;
  const size_type Nt_;
  Ptr<PartitionedVector<Real>> z_dyn_;
  Ptr<PartitionedVector<Real>> v_dyn_;
  Ptr<PartitionedVector<Real>> g_dyn_;
};
} // namespace ROL

#endif
