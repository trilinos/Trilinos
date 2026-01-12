// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <gtest/gtest.h>
#include <stk_util/stk_config.h>
#include <stk_mesh/base/Ngp.hpp>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/NgpField.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/PerformanceTester.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_unit_test_utils/timer.hpp>
#include <stk_performance_tests/stk_mesh/calculate_centroid.hpp>
#include <stk_performance_tests/stk_mesh/multi_block.hpp>

namespace {

class FieldBLAS : public stk::unit_test_util::MeshFixture
{
public:
  FieldBLAS()
    : batchTimer(get_comm())
  { }

protected:
  stk::mesh::Field<double>& declare_element_vector_field(const std::string& fieldName)
  {
    stk::mesh::Field<double>& newField = get_meta().declare_field<double>(stk::topology::ELEM_RANK, fieldName);
    stk::mesh::put_field_on_mesh(newField, get_meta().universal_part(), 3, nullptr);
    return newField;
  }

  void verify_result(const stk::mesh::Field<double>& field, double expectedResult)
  {
    auto fieldData = field.data();
    stk::mesh::BucketVector buckets = get_bulk().buckets(stk::topology::ELEM_RANK);

    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::ScalarIdx scalar : bucketValues.scalars()) {
          EXPECT_DOUBLE_EQ(bucketValues(entity, scalar), expectedResult);
        }
      }
    }
  }

  void verify_result(const stk::mesh::Field<double>& field, const std::array<double, 3>& expectedResult)
  {
    auto fieldData = field.data();
    stk::mesh::BucketVector buckets = get_bulk().buckets(stk::topology::ELEM_RANK);

    for (stk::mesh::Bucket* bucket : buckets) {
      auto bucketValues = fieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx entity : bucket->entities()) {
        for (stk::mesh::ScalarIdx scalar : bucketValues.scalars()) {
          EXPECT_DOUBLE_EQ(bucketValues(entity, scalar), expectedResult[scalar]);
        }
      }
    }
  }

  stk::unit_test_util::BatchTimer batchTimer;
  stk::mesh::Field<double> *m_centroidField;
  stk::mesh::Field<double, stk::mesh::Layout::Left> *m_centroidFieldLeft;
  stk::mesh::Field<double, stk::mesh::Layout::Right> *m_centroidFieldRight;
};

constexpr int NUM_RUNS = 5;

//------------------------------------------------------------------------------
// field_axpy: y[i] = a*x[i] + y[i]
//
TEST_F(FieldBLAS, field_axpy)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int NUM_ITERS = 500;

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<double>& fieldX = declare_element_vector_field("fieldX");
  stk::mesh::Field<double>& fieldY = declare_element_vector_field("fieldY");
  stk::io::fill_mesh("generated:100x100x100", get_bulk());

  const double a = 2.0;
  const double initX = 1.0;
  const double initY = 0.0;

  for (int run = 0; run < NUM_RUNS; ++run) {
    stk::mesh::field_fill(initX, fieldX);
    stk::mesh::field_fill(initY, fieldY);
    double expected = initY;

    batchTimer.start_batch_timer();
    for (int iter = 0; iter < NUM_ITERS; ++iter) {
      stk::mesh::field_axpy(a, fieldX, fieldY);
      expected = a*initX + expected;
    }
    batchTimer.stop_batch_timer();

    verify_result(fieldY, expected);
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

//------------------------------------------------------------------------------
// field_axpby: y[i] = a*x[i] + b*y[i]
//
TEST_F(FieldBLAS, field_axpby)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int NUM_ITERS = 500;

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<double>& fieldX = declare_element_vector_field("fieldX");
  stk::mesh::Field<double>& fieldY = declare_element_vector_field("fieldY");
  stk::io::fill_mesh("generated:100x100x100", get_bulk());

  const double a = 2.0;
  const double b = 1.001;
  const double initX = 1.0;
  const double initY = 0.0;

  for (int run = 0; run < NUM_RUNS; ++run) {
    stk::mesh::field_fill(initX, fieldX);
    stk::mesh::field_fill(initY, fieldY);
    double expected = initY;

    batchTimer.start_batch_timer();
    for (int iter = 0; iter < NUM_ITERS; ++iter) {
      stk::mesh::field_axpby(a, fieldX, b, fieldY);
      expected = a*initX + b*expected;
    }
    batchTimer.stop_batch_timer();

    verify_result(fieldY, expected);
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

//------------------------------------------------------------------------------
// field_product: z[i] = x[i] * y[i]
//
TEST_F(FieldBLAS, field_product)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int NUM_ITERS = 500;

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<double>& fieldX = declare_element_vector_field("fieldX");
  stk::mesh::Field<double>& fieldY = declare_element_vector_field("fieldY");
  stk::mesh::Field<double>& fieldZ = declare_element_vector_field("fieldZ");
  stk::io::fill_mesh("generated:100x100x100", get_bulk());

  const double initX = 2.0;
  const double initY = 3.0;

  for (int run = 0; run < NUM_RUNS; ++run) {
    stk::mesh::field_fill(initX, fieldX);
    stk::mesh::field_fill(initY, fieldY);

    batchTimer.start_batch_timer();
    for (int iter = 0; iter < NUM_ITERS; ++iter) {
      stk::mesh::field_product(fieldX, fieldY, fieldZ);
    }
    batchTimer.stop_batch_timer();

    double expected = initX * initY;
    verify_result(fieldZ, expected);
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

//------------------------------------------------------------------------------
// field_copy: y[i] = x[i]
//
TEST_F(FieldBLAS, field_copy)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int NUM_ITERS = 500;

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<double>& fieldX = declare_element_vector_field("fieldX");
  stk::mesh::Field<double>& fieldY = declare_element_vector_field("fieldY");
  stk::io::fill_mesh("generated:100x100x100", get_bulk());

  const double initX = 10.0;
  const double initY = 0.0;

  for (int run = 0; run < NUM_RUNS; ++run) {
    stk::mesh::field_fill(initX, fieldX);
    stk::mesh::field_fill(initY, fieldY);

    batchTimer.start_batch_timer();
    for (int iter = 0; iter < NUM_ITERS; ++iter) {
      stk::mesh::field_copy(fieldX, fieldY);
    }
    batchTimer.stop_batch_timer();

    double expected = initX;
    verify_result(fieldY, expected);
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

//------------------------------------------------------------------------------
// field_dot: global_sum( sum_i( x[i]*y[i] ) )
//
TEST_F(FieldBLAS, field_dot)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int NUM_ITERS = 500;

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<double>& fieldX = declare_element_vector_field("fieldX");
  stk::mesh::Field<double>& fieldY = declare_element_vector_field("fieldY");
  stk::io::fill_mesh("generated:100x100x100", get_bulk());
  const size_t numElems = 100*100*100;

  const double initX = 1.01;
  const double initY = 1.02;

  for (int run = 0; run < NUM_RUNS; ++run) {
    stk::mesh::field_fill(initX, fieldX);
    stk::mesh::field_fill(initY, fieldY);
    double expected = 0.0;

    batchTimer.start_batch_timer();
    double result  = 0.0;
    for (int iter = 0; iter < NUM_ITERS; ++iter) {
      result = stk::mesh::field_dot(fieldX, fieldY);
    }
    batchTimer.stop_batch_timer();

    expected = numElems*(initX*initY + initX*initY + initX*initY);
    EXPECT_NEAR(result, expected, 1.e-6);
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

//------------------------------------------------------------------------------
// field_nrm2: sqrt( global_sum( sum_i( x[i]*x[i] )))
//
TEST_F(FieldBLAS, field_nrm2)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int NUM_ITERS = 500;

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<double>& fieldX = declare_element_vector_field("fieldX");
  stk::io::fill_mesh("generated:100x100x100", get_bulk());
  const size_t numElems = 100*100*100;

  const double initX = 1.01;

  for (int run = 0; run < NUM_RUNS; ++run) {
    stk::mesh::field_fill(initX, fieldX);
    double expected = 0.0;

    batchTimer.start_batch_timer();
    double result  = 0.0;
    for (int iter = 0; iter < NUM_ITERS; ++iter) {
      result = stk::mesh::field_nrm2(fieldX);
    }
    batchTimer.stop_batch_timer();

    expected = std::sqrt(numElems*(initX*initX + initX*initX + initX*initX));
    EXPECT_NEAR(result, expected, 1.e-6);
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

//------------------------------------------------------------------------------
// field_scale: x[i] = alpha * x[i]
//
TEST_F(FieldBLAS, field_scale)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int NUM_ITERS = 500;

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<double>& fieldX = declare_element_vector_field("fieldX");
  stk::io::fill_mesh("generated:100x100x100", get_bulk());

  const double alpha = 1.01;
  const double initX = 1.0;

  for (int run = 0; run < NUM_RUNS; ++run) {
    stk::mesh::field_fill(initX, fieldX);
    double expected = initX;

    batchTimer.start_batch_timer();
    for (int iter = 0; iter < NUM_ITERS; ++iter) {
      stk::mesh::field_scale(alpha, fieldX);
      expected *= alpha;
    }
    batchTimer.stop_batch_timer();

    verify_result(fieldX, expected);
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

//------------------------------------------------------------------------------
// field_fill: x[i] = alpha
//
TEST_F(FieldBLAS, field_fill)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int NUM_ITERS = 500;

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<double>& fieldX = declare_element_vector_field("fieldX");
  stk::io::fill_mesh("generated:100x100x100", get_bulk());

  const double alpha = 5.0;
  const double initX = 0.0;

  for (int run = 0; run < NUM_RUNS; ++run) {
    stk::mesh::field_fill(initX, fieldX);

    batchTimer.start_batch_timer();
    for (int iter = 0; iter < NUM_ITERS; ++iter) {
      stk::mesh::field_fill(alpha, fieldX);
    }
    batchTimer.stop_batch_timer();

    double expected = alpha;
    verify_result(fieldX, expected);
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

//------------------------------------------------------------------------------
// field_fill_component: x[i] = alpha[i]
//
TEST_F(FieldBLAS, field_fill_component)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int NUM_ITERS = 500;

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<double>& fieldX = declare_element_vector_field("fieldX");
  stk::io::fill_mesh("generated:100x100x100", get_bulk());

  const std::array<double, 3> alpha {1.0, 2.0, 3.0};
  const double initX = 0.0;

  for (int run = 0; run < NUM_RUNS; ++run) {
    stk::mesh::field_fill(initX, fieldX);

    batchTimer.start_batch_timer();
    for (int iter = 0; iter < NUM_ITERS; ++iter) {
      stk::mesh::field_fill_component(alpha.data(), fieldX);
    }
    batchTimer.stop_batch_timer();

    verify_result(fieldX, alpha);
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

//------------------------------------------------------------------------------
// field_swap: x[i] = y[i]; y[i] = x[i]
//
TEST_F(FieldBLAS, field_swap)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int NUM_ITERS = 501;  // Odd number so that the final is different from initial

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<double>& fieldX = declare_element_vector_field("fieldX");
  stk::mesh::Field<double>& fieldY = declare_element_vector_field("fieldY");
  stk::io::fill_mesh("generated:100x100x100", get_bulk());

  const double initX = 1.0;
  const double initY = 2.0;

  for (int run = 0; run < NUM_RUNS; ++run) {
    stk::mesh::field_fill(initX, fieldX);
    stk::mesh::field_fill(initY, fieldY);

    batchTimer.start_batch_timer();
    for (int iter = 0; iter < NUM_ITERS; ++iter) {
      stk::mesh::field_swap(fieldX, fieldY);
    }
    batchTimer.stop_batch_timer();

    verify_result(fieldX, initY);
    verify_result(fieldY, initX);
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

//------------------------------------------------------------------------------
// field_asum: global_sum( sum( abs(x[i]) ) )
//
TEST_F(FieldBLAS, field_asum)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int NUM_ITERS = 500;

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<double>& fieldX = declare_element_vector_field("fieldX");
  stk::io::fill_mesh("generated:100x100x100", get_bulk());
  const size_t numElems = 100*100*100;

  const double initX = -1.0;

  for (int run = 0; run < NUM_RUNS; ++run) {
    stk::mesh::field_fill(initX, fieldX);

    batchTimer.start_batch_timer();
    double result  = 0.0;
    for (int iter = 0; iter < NUM_ITERS; ++iter) {
      result = stk::mesh::field_asum(fieldX);
    }
    batchTimer.stop_batch_timer();

    const double expected = numElems*(std::abs(initX) + std::abs(initX) + std::abs(initX));
    EXPECT_NEAR(result, expected, 1.e-6);
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

//------------------------------------------------------------------------------
// field_amax: global_max( max_i( abs(x[i]) ) )
//
TEST_F(FieldBLAS, field_amax)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int NUM_ITERS = 500;

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<double>& fieldX = declare_element_vector_field("fieldX");
  stk::io::fill_mesh("generated:100x100x100", get_bulk());

  const std::array<double, 3> initX { -1.0, -2.0, -3.0 };

  for (int run = 0; run < NUM_RUNS; ++run) {
    stk::mesh::field_fill_component(initX.data(), fieldX);

    batchTimer.start_batch_timer();
    double result  = 0.0;
    for (int iter = 0; iter < NUM_ITERS; ++iter) {
      result = stk::mesh::field_amax(fieldX);
    }
    batchTimer.stop_batch_timer();

    const double expected = std::max({std::abs(initX[0]), std::abs(initX[1]), std::abs(initX[2])});
    EXPECT_NEAR(result, expected, 1.e-6);
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

//------------------------------------------------------------------------------
// field_amin: global_min( min( abs(x[i]) ) )
//
TEST_F(FieldBLAS, field_amin)
{
  if (get_parallel_size() != 1) GTEST_SKIP();
  const int NUM_ITERS = 500;

  setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::Field<double>& fieldX = declare_element_vector_field("fieldX");
  stk::io::fill_mesh("generated:100x100x100", get_bulk());

  const std::array<double, 3> initX { -1.0, -2.0, -3.0 };

  for (int run = 0; run < NUM_RUNS; ++run) {
    stk::mesh::field_fill_component(initX.data(), fieldX);

    batchTimer.start_batch_timer();
    double result  = 0.0;
    for (int iter = 0; iter < NUM_ITERS; ++iter) {
      result = stk::mesh::field_amin(fieldX);
    }
    batchTimer.stop_batch_timer();

    const double expected = std::min({std::abs(initX[0]), std::abs(initX[1]), std::abs(initX[2])});
    EXPECT_NEAR(result, expected, 1.e-6);
  }

  batchTimer.print_batch_timing(NUM_ITERS);
}

}
