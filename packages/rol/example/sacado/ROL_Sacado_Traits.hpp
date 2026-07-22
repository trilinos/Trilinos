#ifndef ROL_SACADO_TRAITS_HPP
#define ROL_SACADO_TRAITS_HPP

#include <type_traits>
#include <utility>

#include "Sacado.hpp"
#include "ROL_Objective.hpp"

constexpr bool is_Objective( ... ) { return false; }
constexpr bool is_Constraint( ... ) { return false; }

template<class Real>
constexpr bool is_Objective( const ROL::Objective<Real>& ) { return true; }

template<class Real>
constexpr bool is_Constraint( const ROL::Constraint<Real>& ) { return true; }

constexpr bool is_DFad( ... ) { return false; }
constexpr bool is_SFad( ... ) { return false; }

template<typename T>
constexpr bool is_DFad( const Sacado::Fad::DFad<T>& ) { return true; }

template<typename T, int Num>
constexpr bool is_SFad( const Sacado::Fad::SFad<T,Num>& ) { return true; }

struct InvalidType {};

constexpr bool InvalidType Fad_value( ... );

template<typename T>
constexpr T Fad_value( const Sacado::Fad::DFad<T>& );

template<typename T, int Num>
constexpr T Fad_value( const Sacado::Fad::SFad<T,Num>& );

template<typename T>
using Fad_value_t = decltype( Fad_value( std::declval<T>() ) );

template<typename VectorType>
struct Accessor {};

template<typename Real>
struct Accessor<std::vector<Real>> {
  using container_type = std::vector<Real>;
  using value_type     = typename std::vector<Real>::value_type;
  using size_type      = typename std::vector<Real>::size_type;
  inline value_type get_value( const container_type& x, size_type i ) const {
    return x[i];
  }
  inline void set_value( container_type& x, size_type i, value_type xi ) const {
    x[i] = xi;
  }
};

// 1. Define the trait structure
template <typename T, typename = void> // SFINAE helper for more complex traits, not strictly needed here but good practice
struct is_accessor_specialization : std::false_type {};

// 2. Define the partial specialization for Accessor<VectorType>
template <typename VectorType>
struct is_accessor_specialization<Accessor<VectorType>> : std::true_type {};

// 3. Define the C++17 style _v variable template
template <typename T>
inline constexpr bool has_accessor_v = is_accessor_specialization<T>::value;

// --- Example Usage (C++17) ---

// SFINAE example for a function
template <typename T,
          typename std::enable_if_t<has_accessor_v<T>, int> = 0>
void process_if_accessor(const T& acc_specialization) {
    // This function will only be available in overload resolution if has_accessor_v<T> is true.
    // You can use acc_specialization here.
    // For demonstration:
    if constexpr (has_accessor_v<T>) { // if constexpr is C++17
        // Additional compile-time checks or type usage
    }
}


template<typename Real>
struct ObjectiveValue {
  template<typename ScalarT>
  ScalarT value  
};

#endif // ROL_SACADO_TRAITS_HPP
