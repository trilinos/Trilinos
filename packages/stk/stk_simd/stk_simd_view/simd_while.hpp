#ifndef STK_SIMD_WHILE_H
#define STK_SIMD_WHILE_H

#include <stk_simd/Simd.hpp>


template <typename FuncWhile>
KOKKOS_INLINE_FUNCTION void simd_while(int numValid, FuncWhile whileFunc) {
  while ( stk::simd::are_all(whileFunc(), numValid) ) {}
}



#endif
