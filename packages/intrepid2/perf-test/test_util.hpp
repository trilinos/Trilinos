// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
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
