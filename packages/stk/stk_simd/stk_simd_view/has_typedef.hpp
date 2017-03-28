#ifndef STK_HAS_TYPEDEF_H
#define STK_HAS_TYPEDEF_H

namespace stk {
namespace simd {

template<typename> struct void_ { typedef void type; };

template <typename T, typename=void>
struct HasUseSimd {
  typedef std::false_type type;
};

template <typename T>
struct HasUseSimd<T, typename void_<typename T::use_simd>::type> {
  typedef typename T::use_simd type;
};

}
}

#endif
