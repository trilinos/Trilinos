#ifndef STK_UTIL_UTIL_FORCEINLINE_HPP
#define STK_UTIL_UTIL_FORCEINLINE_HPP

#if !defined(FORCEINLINE) && defined( __INTEL_COMPILER )
#define FORCEINLINE  __forceinline
#endif

#if !defined(FORCEINLINE) && defined( __GNUC__ )
#define FORCEINLINE  inline __attribute__((always_inline))
#endif

#if !defined(FORCEINLINE)
#define FORCEINLINE  inline
#endif

#endif //STK_UTIL_UTIL_FORCEINLINE_HPP
