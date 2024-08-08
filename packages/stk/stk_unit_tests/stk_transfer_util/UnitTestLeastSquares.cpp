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

#include <math.h>                                     // for sqrt
#include <cstddef>                                    // for size_t
#include <cstdint>                                    // for int64_t, uint64_t
#include <limits>                                     // for numeric_limits
#include <stdexcept>                                  // for logic_error
#include <string>                                     // for string, basic_s...
#include <typeinfo>                                   // for type_info
#include <utility>                                    // for move, pair
#include <vector>                                     // for vector, swap
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include "mpi.h"                                      // for MPI_COMM_WORLD

#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/Field.hpp>
#include "stk_mesh/base/Entity.hpp"                   // for Entity
#include "stk_mesh/base/FieldBase.hpp"                // for field_data, Fie...
#include "stk_mesh/base/FieldState.hpp"               // for FieldState
#include "stk_mesh/base/Part.hpp"                     // for Part
#include "stk_mesh/base/Selector.hpp"                 // for Selector, opera...
#include "stk_mesh/base/SkinBoundary.hpp"             // for create_all_sides
#include "stk_mesh/base/Types.hpp"                    // for EntityRank, Ent...
#include "stk_topology/topology.hpp"                  // for topology, topol...
#include "stk_unit_test_utils/FieldEvaluator.hpp"
#include "stk_unit_test_utils/MeshFixture.hpp"        // for MeshFixtureNoTest
#include "stk_unit_test_utils/TextMesh.hpp"           // for setup_text_mesh
#include "stk_unit_test_utils/getOption.h"            // for get_command_lin...
#include "stk_util/diag/String.hpp"                   // for String
#include "stk_util/parallel/Parallel.hpp"             // for parallel_machin...
#include "stk_util/util/ReportHandler.hpp"            // for ThrowRequireMsg
#include <stk_transfer_util/EntityCentroidRecoverField.hpp>
#include <stk_transfer_util/LeastSquares.hpp>
#include "stk_transfer_util/LeastSquaresInterpolation.hpp"
#include <stk_transfer_util/Patch.hpp>
#include <stk_transfer_util/RecoverField.hpp>

namespace {

class LeastSquaresTester : public stk::unit_test_util::MeshFixtureNoTest, public ::testing::Test {
 public:
  LeastSquaresTester()
    : stk::unit_test_util::MeshFixtureNoTest(3)
  {
    m_scaleFactorX = stk::unit_test_util::get_command_line_option("-sx", 1.0);
    m_scaleFactorY = stk::unit_test_util::get_command_line_option("-sy", 1.0);
    m_stretchFactorZ = stk::unit_test_util::get_command_line_option("-sz", 1.0);

    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  }

  ~LeastSquaresTester()
  {
    delete m_patch;
    delete m_patchFilter;
    delete m_selector;
  }

 protected:
  double m_scaleFactorX = 1.0;
  double m_scaleFactorY = 1.0;
  double m_stretchFactorZ = 1.0;

  stk::mesh::FieldBase* m_field = nullptr;

  stk::mesh::Entity m_patchSeed;
  stk::transfer::EntityPatchFilter* m_patchFilter = nullptr;
  stk::transfer::Patch<stk::transfer::EntityPatchFilter>* m_patch = nullptr;

  stk::mesh::Part* m_part = nullptr;
  stk::mesh::Selector *m_selector = nullptr;

  std::vector<double> m_evalPoint;

  unsigned m_nComponents = 0;
  unsigned m_nSamples = 0;
  unsigned m_basisSize = 0;


  void stretch_mesh(stk::mesh::BulkData& bulk, const double scaleFactorX, const double scaleFactorY,
                    const double stretchFactorZ)
  {
    double delta = stretchFactorZ - 1.0;

    stk::mesh::MetaData& meta = bulk.mesh_meta_data();
    const stk::mesh::FieldBase* coordinates = meta.coordinate_field();

    stk::mesh::EntityVector nodes;
    stk::mesh::get_entities(bulk, stk::topology::NODE_RANK, nodes);

    for(stk::mesh::Entity node : nodes) {
      double* data = static_cast<double*>(stk::mesh::field_data(*coordinates, node));

      data[0] *= scaleFactorX;
      data[1] *= scaleFactorY;

      if(data[2] > 0.0) {
        data[2] += delta;
      }
    }
  }

  void generate_stretched_mesh(stk::mesh::BulkData& bulk, unsigned numElemPerDim)
  {
    double scaleFactorX = stk::unit_test_util::get_command_line_option("-sx", 1.0);
    double scaleFactorY = stk::unit_test_util::get_command_line_option("-sy", 1.0);
    double stretchFactorZ = stk::unit_test_util::get_command_line_option("-sz", 1.0);

    ASSERT_TRUE(scaleFactorX >= 0.0);
    ASSERT_TRUE(scaleFactorY >= 0.0);
    ASSERT_TRUE(stretchFactorZ >= 1.0);

    std::string meshSpec("generated:" + std::to_string(numElemPerDim) + "x" + std::to_string(numElemPerDim) + "x" +
                         std::to_string(numElemPerDim));
    stk::io::fill_mesh(meshSpec, bulk);
    stretch_mesh(bulk, scaleFactorX, scaleFactorY, stretchFactorZ);
  }

  void create_mesh(unsigned numElemsPerAxis = 3)
  {
    if(m_field == nullptr) {
      m_field = &get_meta().declare_field<double>(stk::topology::ELEM_RANK, "transfer_field", 1);
      stk::mesh::put_field_on_mesh(*m_field, get_meta().universal_part(), 1, nullptr);

      generate_stretched_mesh(get_bulk(), numElemsPerAxis);
    }
  }

  void create_patch_filter()
  {
    delete m_patchFilter;
    delete m_patch;
    delete m_selector;

    m_part = get_meta().get_part("block_1");
    m_selector = new stk::mesh::Selector(get_meta().universal_part());

    m_patchFilter = new stk::transfer::EntityPatchFilter(m_part, m_selector);
  }

  void create_patch_seed(stk::mesh::EntityId seedId)
  {
    m_patchSeed = get_bulk().get_entity(stk::topology::ELEM_RANK, seedId);
    ASSERT_TRUE(get_bulk().is_valid(m_patchSeed));
  }

  void create_patch(stk::mesh::EntityId seedId)
  {
    create_patch_filter();
    create_patch_seed(seedId);
    m_patch = new stk::transfer::LinearPatch<stk::transfer::EntityPatchFilter>(get_bulk(), m_patchSeed, *m_patchFilter, *m_selector);
  }

  void create_cubic_patch(stk::mesh::EntityId seedId)
  {
    create_patch_filter();
    create_patch_seed(seedId);
    m_patch = new stk::transfer::CubicPatch<stk::transfer::EntityPatchFilter>(get_bulk(), m_patchSeed, *m_patchFilter, *m_selector);
  }

  std::vector<double> get_default_evaluation_point() const
  {
    double x = (3 + 1.0 / 30) * m_scaleFactorX;
    double y = 1.5 * m_scaleFactorY;
    double z = m_stretchFactorZ;

    return std::vector<double>{ x, y, z };
  }

  void finish_extrapolation_setup(stk::unit_test_util::FieldEvaluator& eval,
      stk::transfer::RecoverField::RecoveryType recoveryType)
  {
    stk::unit_test_util::set_element_field(get_bulk(), *m_field, eval);

    double delta = 0.05;
    m_evalPoint = get_default_evaluation_point();
    m_evalPoint[2] += delta;

    m_nComponents = 1;
    m_nSamples = m_patch->size();
    m_basisSize = (unsigned)recoveryType;
  }

  void setup_extrapolation(stk::mesh::EntityId seedId, stk::unit_test_util::FieldEvaluator& eval,
      stk::transfer::RecoverField::RecoveryType recoveryType)
  {
    create_patch(seedId);
    finish_extrapolation_setup(eval, recoveryType);
  }

  void setup_cubic_extrapolation(stk::mesh::EntityId seedId, stk::unit_test_util::FieldEvaluator& eval,
      stk::transfer::RecoverField::RecoveryType recoveryType)
  {
    create_cubic_patch(seedId);
    finish_extrapolation_setup(eval, recoveryType);
  }

  std::vector<std::vector<double> > get_trilinear_coeffs() const
  {
    std::vector<double> xCoeffs{ 1.0, 1.0 };
    std::vector<double> yCoeffs{ 1.0, 1.0 };
    std::vector<double> zCoeffs{ 1.0, 1.0 };
    return std::vector<std::vector<double> >{ xCoeffs, yCoeffs, zCoeffs };
  }

  std::vector<std::vector<double> > get_triquadratic_coeffs() const
  {
    std::vector<double> xCoeffs{ 1.0, 1.0, 1.0 };
    std::vector<double> yCoeffs{ 1.0, 1.0, 1.0 };
    std::vector<double> zCoeffs{ 1.0, 1.0, 1.0 };
    return std::vector<std::vector<double> >{ xCoeffs, yCoeffs, zCoeffs };
  }

  std::vector<std::vector<double> > get_tricubic_coeffs() const
  {
    std::vector<double> xCoeffs{ 1.0, 1.0, 1.0, 1.0 };
    std::vector<double> yCoeffs{ 1.0, 1.0, 1.0, 1.0 };
    std::vector<double> zCoeffs{ 1.0, 1.0, 1.0, 1.0 };
    return std::vector<std::vector<double> >{ xCoeffs, yCoeffs, zCoeffs };
  }

  double get_linear_evaluation(stk::transfer::LeastSquares& leastSquaresCalculator, const std::vector<double>& evalPoint)
  {
    double interpolatedValue;
    stk::mesh::FieldBase const* coordField = get_meta().coordinate_field();

    const stk::transfer::EntityInterpolationData interpData(get_bulk(), m_field, 0, 0, coordField,
                                                             m_patch->get_patch_seed(), m_part, m_selector);
    bool solvable = stk::transfer::least_squares_linear_interpolation(*m_patch, leastSquaresCalculator, evalPoint,
                                                                      interpData, 1u, &interpolatedValue);

    EXPECT_TRUE(solvable);
    return interpolatedValue;
  }

  void test_linear_evaluation(stk::transfer::LeastSquares& leastSquaresCalculator,
                              stk::unit_test_util::FieldEvaluator& evaluator,
                              const std::vector<double>& evalPoint)
  {
    double interpolatedValue = get_linear_evaluation(leastSquaresCalculator, evalPoint);
    double expectedValue = evaluator(stk::mesh::Entity(), evalPoint[0], evalPoint[1], evalPoint[2]);
    EXPECT_NEAR(expectedValue, interpolatedValue, 1.0e-8);
  }

  double
  get_quadratic_evaluation(stk::transfer::LeastSquares& leastSquaresCalculator, const std::vector<double>& evalPoint)
  {
    double interpolatedValue;
    stk::mesh::FieldBase const* coordField = get_meta().coordinate_field();

    const stk::transfer::EntityInterpolationData interpData(get_bulk(), m_field, 0, 0, coordField,
                                                             m_patch->get_patch_seed(), m_part, m_selector);
    bool solvable = stk::transfer::least_squares_quadratic_interpolation(*m_patch, leastSquaresCalculator, evalPoint,
                                                                         interpData, 1u, &interpolatedValue);

    EXPECT_TRUE(solvable);
    return interpolatedValue;
  }

  void test_quadratic_evaluation(stk::transfer::LeastSquares& leastSquaresCalculator,
                                 stk::unit_test_util::FieldEvaluator& evaluator,
                                 const std::vector<double>& evalPoint)
  {
    double interpolatedValue = get_quadratic_evaluation(leastSquaresCalculator, evalPoint);
    double expectedValue = evaluator(stk::mesh::Entity(), evalPoint[0], evalPoint[1], evalPoint[2]);
    EXPECT_NEAR(expectedValue, interpolatedValue, 1.0e-6);
  }

  double
  get_cubic_evaluation(stk::transfer::LeastSquares& leastSquaresCalculator, const std::vector<double>& evalPoint)
  {
    double interpolatedValue;
    stk::mesh::FieldBase const* coordField = get_meta().coordinate_field();

    const stk::transfer::EntityInterpolationData interpData(get_bulk(), m_field, 0, 0, coordField,
                                                             m_patch->get_patch_seed(), m_part, m_selector);
    bool solvable = stk::transfer::least_squares_cubic_interpolation(*m_patch, leastSquaresCalculator, evalPoint,
                                                                     interpData, 1u, &interpolatedValue);

    EXPECT_TRUE(solvable);
    return interpolatedValue;
  }

  void test_cubic_evaluation(stk::transfer::LeastSquares& leastSquaresCalculator,
                             stk::unit_test_util::FieldEvaluator& evaluator,
                             const std::vector<double>& evalPoint)
  {
    double interpolatedValue = get_cubic_evaluation(leastSquaresCalculator, evalPoint);
    double expectedValue = evaluator(stk::mesh::Entity(), evalPoint[0], evalPoint[1], evalPoint[2]);
    EXPECT_NEAR(expectedValue, interpolatedValue, 5.0e-4);
  }

  void test_mls_weights(const std::vector<double>& evalPoint, const unsigned expectedMaxIndex)
  {
    const stk::mesh::BulkData& bulk = get_bulk();
    const stk::mesh::MetaData& meta = get_meta();

    stk::mesh::EntityVector elements;
    stk::mesh::get_selected_entities(meta.universal_part(), bulk.buckets(stk::topology::ELEM_RANK), elements);

    const stk::mesh::FieldBase* coordField = meta.coordinate_field();
    unsigned numComponents = 1;

    stk::transfer::GeometricMovingLeastSquares leastSquaresCalculator(numComponents, elements.size(),
                                                                      stk::transfer::RecoverField::TRILINEAR, get_bulk(),
                                                                      elements, coordField, evalPoint);

    std::vector<double> weights = leastSquaresCalculator.compute_weights();

    EXPECT_EQ(weights.size(), elements.size());

    unsigned maxIndex = 0;
    double maxWeight = -1.0;

    for(unsigned i = 0; i < weights.size(); ++i) {
      EXPECT_TRUE(weights[i] >= 0.0);

      if(weights[i] > maxWeight) {
        maxIndex = i + 1;
        maxWeight = weights[i];
      }
    }

    EXPECT_EQ(expectedMaxIndex, maxIndex);
  }
};

TEST_F(LeastSquaresTester, computeWeights)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) != 1) {
    return;
  }

  create_mesh();

  double delta = 0.05;
  std::vector<double> evalPoint = get_default_evaluation_point();
  evalPoint[2] += delta;
  unsigned expectedMaxIndex = 15u;
  test_mls_weights(evalPoint, expectedMaxIndex);

  evalPoint = get_default_evaluation_point();
  evalPoint[2] -= delta;
  expectedMaxIndex = 6u;
  test_mls_weights(evalPoint, expectedMaxIndex);
}

TEST_F(LeastSquaresTester, extrapolateLinear_UsingLeastSquares)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) != 1) {
    return;
  }

  create_mesh();

  std::vector<std::vector<double> > coeffs = get_trilinear_coeffs();
  stk::unit_test_util::LinearFieldEvaluatorWithCoefficients eval(coeffs[0], coeffs[1], coeffs[2]);
  setup_extrapolation(stk::mesh::EntityId(14u), eval, stk::transfer::RecoverField::TRILINEAR);

  stk::transfer::LeastSquares leastSquaresCalculator(m_nComponents, m_nSamples, m_basisSize);

  test_linear_evaluation(leastSquaresCalculator, eval, m_evalPoint);
}

TEST_F(LeastSquaresTester, extrapolateLinear_UsingMovingLeastSquares)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) != 1) {
    return;
  }

  create_mesh();

  std::vector<std::vector<double> > coeffs = get_trilinear_coeffs();
  stk::unit_test_util::LinearFieldEvaluatorWithCoefficients eval(coeffs[0], coeffs[1], coeffs[2]);
  setup_extrapolation(stk::mesh::EntityId(14u), eval, stk::transfer::RecoverField::TRILINEAR);

  stk::transfer::GeometricMovingLeastSquares leastSquaresCalculator(m_nComponents, m_nSamples, m_basisSize, get_bulk(),
                                                                   m_patch->get_patch_entities(),
                                                                   get_meta().coordinate_field(), m_evalPoint);

  test_linear_evaluation(leastSquaresCalculator, eval, m_evalPoint);
}

TEST_F(LeastSquaresTester, extrapolateQuadratic_AndCompareLinearLeastSquaresVsLinearMovingLeastSquares)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) != 1) {
    return;
  }

  create_mesh();

  std::vector<std::vector<double> > coeffs = get_triquadratic_coeffs();
  stk::unit_test_util::QuadraticFieldEvaluatorWithCoefficients evaluator(coeffs[0], coeffs[1], coeffs[2]);
  setup_extrapolation(stk::mesh::EntityId(14u), evaluator, stk::transfer::RecoverField::TRILINEAR);

  stk::transfer::LeastSquares lsCalculator(m_nComponents, m_nSamples, m_basisSize);

  stk::transfer::GeometricMovingLeastSquares mlsCalculator(m_nComponents, m_nSamples, m_basisSize, get_bulk(),
                                                          m_patch->get_patch_entities(), get_meta().coordinate_field(),
                                                          m_evalPoint);

  double lsEvaluation = get_linear_evaluation(lsCalculator, m_evalPoint);
  double mlsEvaluation = get_linear_evaluation(mlsCalculator, m_evalPoint);

  double expectedValue = evaluator(stk::mesh::Entity(), m_evalPoint[0], m_evalPoint[1], m_evalPoint[2]);

  double lsError = std::abs(lsEvaluation - expectedValue);
  double mlsError = std::abs(mlsEvaluation - expectedValue);

  EXPECT_TRUE(mlsError < lsError);
}

TEST_F(LeastSquaresTester, extrapolateQuadratic_UsingLeastSquares)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) != 1) {
    return;
  }

  create_mesh();

  std::vector<std::vector<double> > coeffs = get_triquadratic_coeffs();
  stk::unit_test_util::QuadraticFieldEvaluatorWithCoefficients evaluator(coeffs[0], coeffs[1], coeffs[2]);
  setup_extrapolation(stk::mesh::EntityId(14u), evaluator, stk::transfer::RecoverField::TRIQUADRATIC);

  stk::transfer::LeastSquares leastSquaresCalculator(m_nComponents, m_nSamples, m_basisSize);

  test_quadratic_evaluation(leastSquaresCalculator, evaluator, m_evalPoint);
}

TEST_F(LeastSquaresTester, extrapolateQuadratic_UsingMovingLeastSquares)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) != 1) {
    return;
  }

  create_mesh();

  std::vector<std::vector<double> > coeffs = get_triquadratic_coeffs();
  stk::unit_test_util::QuadraticFieldEvaluatorWithCoefficients evaluator(coeffs[0], coeffs[1], coeffs[2]);
  setup_extrapolation(stk::mesh::EntityId(14u), evaluator, stk::transfer::RecoverField::TRIQUADRATIC);

  stk::transfer::GeometricMovingLeastSquares leastSquaresCalculator(m_nComponents, m_nSamples, m_basisSize, get_bulk(),
                                                                   m_patch->get_patch_entities(),
                                                                   get_meta().coordinate_field(), m_evalPoint);

  test_quadratic_evaluation(leastSquaresCalculator, evaluator, m_evalPoint);
}

TEST_F(LeastSquaresTester, interpolateQuadratic_UsingLeastSquares_AtCentroid)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) != 1) {
    return;
  }

  create_mesh();

  unsigned numElements = stk::mesh::count_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::ELEM_RANK));

  std::vector<std::vector<double> > coeffs = get_triquadratic_coeffs();
  stk::unit_test_util::QuadraticFieldEvaluatorWithCoefficients evaluator(coeffs[0], coeffs[1], coeffs[2]);

  const unsigned spatialDimension = get_meta().spatial_dimension();
  stk::mesh::FieldBase const* coord = get_meta().coordinate_field();

  for(unsigned patchSeed = 1; patchSeed <= numElements; ++patchSeed) {
    setup_extrapolation(stk::mesh::EntityId(patchSeed), evaluator, stk::transfer::RecoverField::TRIQUADRATIC);

    stk::transfer::LeastSquares leastSquaresCalculator(m_nComponents, m_nSamples, m_basisSize);

    stk::mesh::EntityVector elements = m_patch->get_patch_entities();

    for(stk::mesh::Entity e : elements) {
      std::vector<double> centroid;
      stk::unit_test_util::determine_centroid(spatialDimension, e, *coord, centroid);

      leastSquaresCalculator.resize_data(m_nComponents, m_nSamples, m_basisSize);
      test_quadratic_evaluation(leastSquaresCalculator, evaluator, centroid);
    }
  }
}




TEST_F(LeastSquaresTester, extrapolateCubic_UsingLeastSquares)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) != 1) {
    return;
  }

  create_mesh(5);

  std::vector<std::vector<double> > coeffs = get_tricubic_coeffs();
  stk::unit_test_util::CubicFieldEvaluatorWithCoefficients evaluator(coeffs[0], coeffs[1], coeffs[2]);
  setup_cubic_extrapolation(stk::mesh::EntityId(63u), evaluator, stk::transfer::RecoverField::TRICUBIC);

  stk::transfer::LeastSquares leastSquaresCalculator(m_nComponents, m_nSamples, m_basisSize);

  test_cubic_evaluation(leastSquaresCalculator, evaluator, m_evalPoint);
}

TEST_F(LeastSquaresTester, extrapolateCubic_UsingMovingLeastSquares)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) != 1) {
    return;
  }

  create_mesh(5);

  std::vector<std::vector<double> > coeffs = get_tricubic_coeffs();
  stk::unit_test_util::CubicFieldEvaluatorWithCoefficients evaluator(coeffs[0], coeffs[1], coeffs[2]);
  setup_cubic_extrapolation(stk::mesh::EntityId(63u), evaluator, stk::transfer::RecoverField::TRICUBIC);

  stk::transfer::GeometricMovingLeastSquares leastSquaresCalculator(m_nComponents, m_nSamples, m_basisSize, get_bulk(),
                                                                    m_patch->get_patch_entities(),
                                                                    get_meta().coordinate_field(), m_evalPoint);

  test_cubic_evaluation(leastSquaresCalculator, evaluator, m_evalPoint);
}

}


