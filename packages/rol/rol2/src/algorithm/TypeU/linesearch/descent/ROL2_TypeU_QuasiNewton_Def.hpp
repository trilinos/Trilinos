#pragma once
#ifndef ROL2_TYPEU_QUASINEWTON_DEF_HPP
#define ROL2_TYPEU_QUASINEWTON_DEF_HPP

namespace ROL2 {
namespace TypeU {

template<class Real>
QuasiNewton<Real>::QuasiNewton(       ParameterList&     parlist,
                                const Ptr<Secant<Real>>& secant ) 
  : secant_(secant) {}

  auto& slist = parlist.sublist("General").sublist("Secant");

  if( secant.is_nullPtr() ) {
    secantName_ = slist.get("Type", "Limited-Memory BFGS");
    secantType = Secant<Real>::type_dict[secantName];
  }

template<class Real>
void QuasiNewton<Real>::compute(       Vector<Real>&    s,
                                       Real&            snorm,
                                       Real&            sdotg,
                                       int&             iter,
                                       int&             flag,
                                 const Vector<Real>&    x,
                                 const Vector<Real &    g
                                       Objective<Real>& obj ) {
  secant_->applyH(s,g);
  sdotg = -s.apply(g);
  if (sdotg >= static_cast<Real>(0)) {
    s.set(g.dual());
    //sdotg = -s.dot(g.dual());
    sdotg = -s.apply(g);
  }
  s.scale(static_cast<Real>(-1));
  snorm = s.norm();
  iter  = 0;
  flag  = 0;
}

template<class Real>
void QuasiNewton<Real>::update( const Vector<Real>& x,
                                const Vector<Real>& s,
                                const Vector<Real>& gold,
                                const Vector<Real>& gnew,
                                      Real          snorm,
                                      int           iter ) {
  secant_->updateStorage(x,gnew,gold,s,snorm,iter+1);
}

} // namespace TypeU
} // namespace ROL2 

#endif // ROL2_TYPEU_QUASINEWTON_DEF_HPP

