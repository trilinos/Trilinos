#pragma once
#ifndef ROL2_TYPEU_NEWTON_DEF_HPP
#define ROL2_TYPEU_NEWTON_DEF_HPP

namespace ROL2 {
namespace TypeU {

template<class Real>
void Newton<Real>::compute(       Vector<Real>&    s,
                                  Real&            snorm,
                                  Real&            sdotg,
                                  int&             iter,
                                  int&             flag,
                            const Vector<Real>&    x,
                            const Vector<Real>&    g,
                                  Objective<Real>& obj ) {

  Real tol = default_tolerance<Real>();

  // Compute unconstrained step
  obj.invHessVec(s,g,x,tol);
  sdotg = -s.apply(g);
  if (sdotg >= static_cast<Real>(0)) {
    s.set(g.dual());
    sdotg = -s.apply(g);
  }
  s.scale(static_cast<Real>(-1));
  snorm = s.norm();
  iter  = 0;
  flag  = 0; 
}

} // namespace TypeU
} // namespace ROL2








#endif //ROL2_TYPEU_NEWTON_DEF_HPP

