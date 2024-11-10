// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_util.cpp
    \brief  Performance test comparing dynrankview overhead
    \author Created by Kyungjoo Kim.
*/
#include "Kokkos_Core.hpp"
#include <Kokkos_Timer.hpp>

namespace Intrepid2 {
  
  namespace Test {
    
    template<size_t BufSize, typename ExecSpaceType>
    struct Flush {
      typedef double value_type;

      // flush a large host buffer
      Kokkos::View<value_type*,ExecSpaceType> _buf;
      Flush() : _buf("Flush::buf", BufSize) {
        Kokkos::deep_copy(_buf, 1);
      }
      
      KOKKOS_INLINE_FUNCTION
      void init(value_type &update) {
        update = 0;
      }
      
      KOKKOS_INLINE_FUNCTION
      void join(volatile value_type &update,
                const volatile value_type &input) {
        update += input;
      }
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const int i, value_type &update) const {
        update += _buf[i];
      }
      
      void run() {
        double sum = 0;

	Kokkos::RangePolicy<ExecSpaceType> policy(0, BufSize/sizeof(double));
        Kokkos::parallel_reduce(policy, *this, sum);

        FILE *fp = fopen("/dev/null", "w");
        fprintf(fp, "%f\n", sum);
        fclose(fp);
      }
    };

  } // end of namespace TEST
} // end of namespace Intrepid2
