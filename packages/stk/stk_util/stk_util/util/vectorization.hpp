#ifndef STK_UTIL_UTIL_VECTORIZATION_H
#define STK_UTIL_UTIL_VECTORIZATION_H

#if defined(__INTEL_COMPILER)
#define RESTRICT_ALIAS restrict
#elif defined(__GNUC__)
#define RESTRICT_ALIAS __restrict__
#else
#define RESTRICT_ALIAS
#endif

#endif
