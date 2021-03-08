#pragma once
#ifndef ROL2_ELEMENTWISE_HPP
#define ROL2_ELEMENTWISE_HPP

namespace ROL2 {
namespace Elementwise {

template<typename Real, template<typename> class...DerivedTypes>
struct Function {
  using visitor_type = Visitor<DerivedTypes<Real>...>;
  virtual ~Function() = default;
  virtual Real apply( const Real& ) const = 0;
  virtual void accept( visitor_type& ) const {}
};

} // namespace Elementwise
} // namespace ROL2

#include "ROL2_Elementwise_UnaryFunction.hpp"
#include "ROL2_Elementwise_BinaryFunction.hpp"
#include "ROL2_Elementwise_Reduce.hpp"

#endif // ROL2_ELEMENTWISE_HPP

