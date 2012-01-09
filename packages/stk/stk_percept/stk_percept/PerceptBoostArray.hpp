#if defined(__INTEL_COMPILER) && (__INTEL_COMPILER == 1210)

#  define PERCEPT_BOOST_DISABLE_ASSERTS 0
#  if defined(BOOST_DISABLE_ASSERTS)
#    define PERCEPT_BOOST_DISABLE_ASSERTS 1
#  endif

#  define BOOST_DISABLE_ASSERTS

#endif

#include <boost/multi_array.hpp>
#include <boost/array.hpp>

#if PERCEPT_BOOST_DISABLE_ASSERTS
#define BOOST_DISABLE_ASSERTS
#else
#undef BOOST_DISABLE_ASSERTS
#endif
