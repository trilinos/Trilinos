// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SLACKLESSOBJECTIVE_HPP
#define ROL_SLACKLESSOBJECTIVE_HPP

#include "ROL_Objective.hpp"
#include "ROL_PartitionedVector.hpp"

/** @ingroup func_group
 *  \class ROL::SlacklessObjective
 *  \brief This class strips out the slack variables from objective evaluations
 *         to create the new objective  \f$ F(x,s) = f(x) \f$
 */

namespace ROL {

template<typename Real> 
class SlacklessObjective : public Objective<Real> {
private: 
  const Ptr<Objective<Real>> obj_;

public:
  ~SlacklessObjective() {}
  SlacklessObjective( const Ptr<Objective<Real>> &obj );

  Ptr<Objective<Real>> getObjective(void) const;
 
  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) override;
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) override;
  Real value( const Vector<Real> &x, Real &tol ) override;
  Real dirDeriv( const Vector<Real> &x, const Vector<Real> &d, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;
  void invHessVec( Vector<Real> &ihv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;
  void precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;

// Definitions for parametrized (stochastic) objective functions
public:
  void setParameter(const std::vector<Real> &param) override;

private:
  Ptr<Vector<Real>> getOpt( Vector<Real> &xs ) const;
  Ptr<const Vector<Real>> getOpt( const Vector<Real> &xs ) const;
  void zeroSlack( Vector<Real> &x ) const;
}; // class SlacklessObjective 

} // namespace ROL

#include "ROL_SlacklessObjective_Def.hpp"

#endif // ROL_SLACKLESSOBJECTIVE_HPP

