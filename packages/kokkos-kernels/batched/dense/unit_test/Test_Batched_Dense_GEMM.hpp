#ifndef TEST_BATCHED_DENSE_GEMM_HPP
#define TEST_BATCHED_DENSE_GEMM_HPP

// Serial kernels
#include "Test_Batched_SerialGemm.hpp"
#include "Test_Batched_SerialGemm_Real.hpp"
#include "Test_Batched_SerialGemm_Complex.hpp"
#include "Test_Batched_BatchedGemm.hpp"
#include "Test_Batched_BatchedGemm_Real.hpp"
#include "Test_Batched_BatchedGemm_Complex.hpp"

// Team Kernels
#include "Test_Batched_TeamGemm.hpp"
#include "Test_Batched_TeamGemm_Real.hpp"
#include "Test_Batched_TeamGemm_Complex.hpp"

// TeamVector Kernels
#include "Test_Batched_TeamVectorGemm.hpp"
#include "Test_Batched_TeamVectorGemm_Real.hpp"
#include "Test_Batched_TeamVectorGemm_Complex.hpp"

#endif  // TEST_BATCHED_DENSE_GEMM_HPP
