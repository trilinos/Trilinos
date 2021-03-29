#pragma once
#ifndef ROL2_TYPEU_QUASINEWTON_DECL_HPP
#define ROL2_TYPEU_QUASINEWTON_DECL_HPP

namespace ROL2 {
namespace TypeU {

template<class Real>
class QuasiNewton : public DescentDirection<Real> {
public:

  using SecantType = typename Secant<Real>::Type;

  QuasiNewton(       ParameterList&     parlist, 
               const Ptr<Secant<Real>>& secant = nullPtr );

  void compute(       Vector<Real>&    s,
                      Real&            snorm,
                      Real&            sdotg,
                      int&             iter,
                      int&             flag,
                const Vector<Real>&    x,
                const Vector<Real>&    g,
                      Objective<Real>& obj) override;

  void update( const Vector<Real>& x,
               const Vector<Real>& s,
               const Vector<Real>& gold,
               const Vector<Real>& gnew,
                     Real          snorm,
                     int           iter ) override;
  
  void writeName( std::ostream& os ) const override {
    os << "Quasi-Newton method with " << secantName_;
  }

private:

  Ptr<Secant<Real>> secant_;
  SecantType        secantType_ = SecantType::UserDefined;
  std::string       secantName_;

}; // class QuasiNewton

} // namespace TypeU 
} // namespace ROL2 

#endif // ROL2_TYPEU_QUASINEWTON_DECL_HPP

