// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_STDOBJECTIVE_H
#define ROL_STDOBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_StdVector.hpp"

/** @ingroup func_group
    \class ROL::StdObjective
    \brief Specializes the ROL::Objective interface for objective functions
    that operate on ROL::StdVector's.
*/


namespace ROL {

template<typename Real>
class StdObjective : public virtual Objective<Real> {
public:
  virtual void update( const std::vector<Real> &x, bool flag = true, int iter = -1 ) {}

  using Objective<Real>::update;
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) override;

  virtual void update( const std::vector<Real> &x, UpdateType type, int iter = -1 ) {}

  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) override;

  virtual Real value( const std::vector<Real> &x, Real &tol ) = 0;

  using Objective<Real>::value;
  Real value( const Vector<Real> &x, Real &tol ) override;

  virtual void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol );

  using Objective<Real>::gradient;
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) override;

  virtual Real dirDeriv( const std::vector<Real> &x, const std::vector<Real> &d, Real &tol );

  using Objective<Real>::dirDeriv;
  Real dirDeriv( const Vector<Real> &x, const Vector<Real> &d, Real &tol ) override;

  virtual void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol );

  using Objective<Real>::hessVec;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;

  virtual void invHessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol );

  using Objective<Real>::invHessVec;
  void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;

  virtual void precond( std::vector<Real> &Pv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol );

  using Objective<Real>::precond;
  void precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;

private:
  Real sgn(Real x) const;

};

} // namespace ROL

#include "ROL_StdObjective_Def.hpp"

#endif
