#pragma once
#ifndef ROL2_ELEMENTWISE_HPP
#define ROL2_ELEMENTWISE_HPP

namespace ROL2 {
namespace Elementwise {

template<typename Real, template<typename> class...DerivedTypes>
struct Function {
  using visitor_type = Visitor<DerivedTypes<Real>...>;
  virtual ~Function() = default;
  virtual void accept( visitor_type& ) const {}
};

template<typename Real, template<typename> class...DerivedTypes>
struct UFunction : Function<Real,DerivedTypes...> {
  using visitor_type = Visitor<DerivedTypes<Real>...>;
  virtual ~UFunction() = default;
  virtual Real apply( const Real& ) const = 0;
};

template<typename Real, template<typename> class...DerivedTypes>
struct BFunction : Function<Real,DerivedTypes...> {
  using visitor_type = Visitor<DerivedTypes<Real>...>;
  virtual ~BFunction() = default;
  virtual Real apply( const Real&, const Real& ) const = 0;
};

} // namespace Elementwise
} // namespace ROL2

#include "ROL2_Elementwise_UnaryFunction.hpp"
#include "ROL2_Elementwise_BinaryFunction.hpp"
#include "ROL2_Elementwise_Reduce.hpp"

#endif // ROL2_ELEMENTWISE_HPP

