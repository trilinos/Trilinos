#pragma once
#ifndef ROL2_TYPEU_GRADIENT_DECL_HPP
#define ROL2_TYPEU_GRADIENT_DECL_HPP

namespace ROL2 {
namespace TypeU {

template<class Real>
class Gradient : public DescentDirection<Real> {
public:
  Gradient() = default;
  virtual ~Gradient() = default;

  void compute(       Vector<Real>&    s,
                      Real&            snorm,
                      Real&            sdotg,
                      int&             iter,
                      int&             flag,
                const Vector<Real>&    x,
                const Vector<Real &    g,
                      Objective<Real>& obj ) override; 

  void writeName( std::ostream& os ) const override { os << "Gradient Descent"; }

}; // class Gradient

} // namespace TypeU
} // namespace ROL2







#endif //ROL2_TYPEU_GRADIENT_DECL_HPP

