#pragma once
#ifndef ROL2_TYPEU_GRADIENT_DEF_HPP
#define ROL2_TYPEU_GRADIENT_DEF_HPP


namespace ROL2 {
namespace TypeU {

template<class Real>
void Gradient<Real>::compute(       Vector<Real>&    s,
                                    Real&            snorm,
                                    Real&            sdotg,
                                    int&             iter,
                                    int&             flag,
                              const Vector<Real>&    x,
                              const Vector<Real &    g,
                                    Objective<Real>& obj ) override {
  s.set(g.dual());
  s.scale(static_cast<Real>(-1));
  snorm = s.norm();
  sdotg = s.apply(g);
  iter  = 0;
  flag  = 0;
}

} // namespace TypeU
} // namespace ROL2

#endif //ROL2_TYPEU_GRADIENT_DEF_HPP

