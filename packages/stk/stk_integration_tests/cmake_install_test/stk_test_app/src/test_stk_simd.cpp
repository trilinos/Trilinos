#include "test_stk_simd.hpp"

#include <stk_util/util/ReportHandler.hpp>
#include <stk_simd/Simd.hpp>

namespace test_stk_lib {

void test_stk_simd(stk::ParallelMachine comm)
{
  stk::simd::Double simdDouble1 = 1.0;
  stk::simd::Double simdDouble2 = 2.0;

  for (int i=0; i<stk::simd::ndoubles; ++i) {
    STK_ThrowRequireMsg(simdDouble1[i] < simdDouble2[i],
                    "simdDouble1["<<i<<"] = "<<simdDouble1[i]<<" should be less than "
                    <<"simdDouble2["<<i<<"] = "<<simdDouble2[i]);
  }

  if (stk::parallel_machine_rank(comm) == 0) {
    std::cout<<"test_stk_simd, stk::simd::ndoubles = "<<stk::simd::ndoubles<<std::endl;
  }
}

}

