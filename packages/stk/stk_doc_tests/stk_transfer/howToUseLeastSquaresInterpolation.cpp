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
#include <stddef.h>
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

#include <stk_io/FillMesh.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/Field.hpp>
#include "stk_mesh/base/Entity.hpp"                   // for Entity
#include "stk_mesh/base/FieldBase.hpp"                // for field_data, Fie...
#include "stk_mesh/base/FieldState.hpp"               // for FieldState
#include "stk_mesh/base/ForEachEntity.hpp"
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
#include <stk_transfer_util/LeastSquares.hpp>
#include "stk_transfer_util/LeastSquaresInterpolation.hpp"
#include <stk_transfer_util/Patch.hpp>
#include <stk_transfer_util/RecoverField.hpp>

namespace {
double evaluate_polynomial(double x, const std::vector<double>& coeffs)
{
  double value = 0.0;
  unsigned dim = coeffs.size();

  for(unsigned i = 0; i < dim; ++i) {
    value += coeffs[i] * std::pow(x, i);
  }

  return value;
}

double evaluate_tripolynomial_function(const std::vector<std::vector<double> >& coeffs,
                                       double x, double y, double z)
{
  EXPECT_EQ(3u, coeffs.size());
  EXPECT_EQ(coeffs[0].size(), coeffs[1].size());
  EXPECT_EQ(coeffs[0].size(), coeffs[2].size());

  const std::vector<double>& xcoeffs = coeffs[0];
  const std::vector<double>& ycoeffs = coeffs[1];
  const std::vector<double>& zcoeffs = coeffs[2];

  double xValue = evaluate_polynomial(x, xcoeffs);
  double yValue = evaluate_polynomial(y, ycoeffs);
  double zValue = evaluate_polynomial(z, zcoeffs);

  return xValue * yValue * zValue;
}

std::vector<double> compute_centroid(const stk::mesh::BulkData& mesh, stk::mesh::Entity elem, const stk::mesh::FieldBase& nodalCoordField)
{
  const unsigned spatialDimension = mesh.mesh_meta_data().spatial_dimension();

  std::vector<double> centroid(spatialDimension, 0.0);

  const stk::mesh::Entity* const nodes = mesh.begin_nodes(elem);
  const unsigned numNodes = mesh.num_nodes(elem);

  for(unsigned i = 0; i < spatialDimension; ++i) {
    for(unsigned iNode = 0; iNode < numNodes; ++iNode) {
      stk::mesh::Entity node = nodes[iNode];
      double* coor = static_cast<double*>(stk::mesh::field_data(nodalCoordField, node));
      STK_ThrowRequire(coor);
      centroid[i] += coor[i];
    }
  }
  for(unsigned i = 0; i < spatialDimension; ++i) {
    centroid[i] /= numNodes;
  }

  return centroid;
}

void set_element_field_vals_from_polynomial_coefficients(const stk::mesh::BulkData& mesh,
                                                         const stk::mesh::Field<double>& field,
                                                         const std::vector<std::vector<double> >& coeffs)
{
  STK_ThrowRequireMsg(field.entity_rank() == stk::topology::ELEM_RANK, "Input field: " << field.name() << " must be ELEMENT_RANK");
  stk::mesh::FieldBase const* nodalCoordField = mesh.mesh_meta_data().coordinate_field();

  stk::mesh::for_each_entity_run(mesh, stk::topology::ELEM_RANK,
                                 mesh.mesh_meta_data().locally_owned_part(),
                                 [&field, &nodalCoordField, &coeffs](const stk::mesh::BulkData& mesh, const stk::mesh::Entity& elem)
  {
    std::vector<double> centroid = compute_centroid(mesh, elem, *nodalCoordField);
    double x = centroid[0];
    double y = centroid[1];
    double z = (mesh.mesh_meta_data().spatial_dimension() == 2 ? 0.0 : centroid[2]);

    double * scalar = stk::mesh::field_data(field, elem);
    *scalar = evaluate_tripolynomial_function(coeffs, x, y, z);
  });
}

void set_node_field_vals_from_polynomial_coefficients(const stk::mesh::BulkData& mesh,
                                                      const stk::mesh::Field<double>& field,
                                                      const std::vector<std::vector<double> >& coeffs)
{
  STK_ThrowRequireMsg(field.entity_rank() == stk::topology::NODE_RANK, "Input field: " << field.name() << " must be NODE_RANK");
  stk::mesh::FieldBase const* nodalCoordField = mesh.mesh_meta_data().coordinate_field();

  stk::mesh::for_each_entity_run(mesh, stk::topology::NODE_RANK,
                                 mesh.mesh_meta_data().locally_owned_part(),
                                 [&field, &nodalCoordField, &coeffs](const stk::mesh::BulkData& mesh, const stk::mesh::Entity& node)
  {
    double * coords = (double*)stk::mesh::field_data(*nodalCoordField, node);
    double x = coords[0];
    double y = coords[1];
    double z = (mesh.mesh_meta_data().spatial_dimension() == 2 ? 0.0 : coords[2]);

    double * scalar = stk::mesh::field_data(field, node);
    *scalar = evaluate_tripolynomial_function(coeffs, x, y, z);
  });
}

//BEGIN
TEST(StkTransferHowTo, useNodeLinearLeastSquaresInterpolation)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 1) { GTEST_SKIP(); }

  const std::string meshSpec("generated:2x2x2");
  double initVals = std::numeric_limits<double>::max();
  const unsigned spatialDim = 3;

  stk::mesh::MeshBuilder builder(communicator);
  builder.set_spatial_dimension(spatialDim);
  std::shared_ptr<stk::mesh::BulkData> mesh = builder.create();
  stk::mesh::MetaData& meta = mesh->mesh_meta_data();
  stk::mesh::Field<double> &transferField = meta.declare_field<double>(stk::topology::NODE_RANK, "transfer_field", 1);
  stk::mesh::put_field_on_mesh(transferField, meta.universal_part(), &initVals);
  stk::io::fill_mesh(meshSpec, *mesh);

  // Populate node field with tri-linear spatial function based on node coordinates
  const std::vector<std::vector<double> > triLinearCoeffs = {{1, 1}, {1, 1}, {1, 1}};
  set_node_field_vals_from_polynomial_coefficients(*mesh, transferField, triLinearCoeffs);

  stk::mesh::EntityId patchSeedId = 14u; // Middle node
  stk::mesh::Entity patchSeed = mesh->get_entity(stk::topology::NODE_RANK, patchSeedId);
  stk::mesh::Part* part = meta.get_part("block_1");
  stk::mesh::Selector selector(meta.universal_part());
  stk::transfer::EntityPatchFilter patchFilter(part, &selector);
  stk::transfer::LinearPatch<stk::transfer::EntityPatchFilter> patch(*mesh, patchSeed, patchFilter, selector);

  unsigned nComponents = 1; // scalar field
  unsigned nSamples = patch.get_patch_entities().size(); // All nodes
  unsigned basisSize = (unsigned)stk::transfer::RecoverField::TRILINEAR;  // 8 basis function evaluations for linear

  EXPECT_EQ(27u, nSamples); // 27 nodes

  stk::transfer::LeastSquares leastSquaresCalculator(nComponents, nSamples, basisSize);

  double x = (2 + 1.0 / 20);
  double y = 1;
  double z = 1.05;
  std::vector<double> evalPoint{x, y, z}; // Slightly outside domain so it's an extrapolation

  double interpolatedValue;
  const stk::transfer::EntityInterpolationData interpData(*mesh, &transferField, 0, 0, meta.coordinate_field(),
                                                           patchSeed, part, &selector);
  bool solvable = stk::transfer::least_squares_linear_interpolation(patch, leastSquaresCalculator, evalPoint,
                                                                    interpData, 1u, &interpolatedValue);

  EXPECT_TRUE(solvable);
  double expectedValue = evaluate_tripolynomial_function(triLinearCoeffs, x, y, z);
  EXPECT_NEAR(expectedValue, interpolatedValue, 1.0e-6);
}
//END

//BEGIN
TEST(StkTransferHowTo, useElementCentroidLinearLeastSquaresInterpolation)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 1) { GTEST_SKIP(); }

  const std::string meshSpec("generated:3x3x3");
  double initVals = std::numeric_limits<double>::max();
  const unsigned spatialDim = 3;

  stk::mesh::MeshBuilder builder(communicator);
  builder.set_spatial_dimension(spatialDim);
  std::shared_ptr<stk::mesh::BulkData> mesh = builder.create();
  stk::mesh::MetaData& meta = mesh->mesh_meta_data();
  stk::mesh::Field<double> &transferField = meta.declare_field<double>(stk::topology::ELEM_RANK, "transfer_field", 1);
  stk::mesh::put_field_on_mesh(transferField, meta.universal_part(), &initVals);
  stk::io::fill_mesh(meshSpec, *mesh);

  // Populate element field with tri-linear spatial function based on element centroid
  const std::vector<std::vector<double> > triLinearCoeffs = {{1, 1}, {1, 1}, {1, 1}};
  set_element_field_vals_from_polynomial_coefficients(*mesh, transferField, triLinearCoeffs);

  stk::mesh::EntityId patchSeedId = 14u; // Middle element
  stk::mesh::Entity patchSeed = mesh->get_entity(stk::topology::ELEM_RANK, patchSeedId);
  stk::mesh::Part* part = meta.get_part("block_1");
  stk::mesh::Selector selector(meta.universal_part());
  stk::transfer::EntityPatchFilter patchFilter(part, &selector);
  stk::transfer::LinearPatch<stk::transfer::EntityPatchFilter> patch(*mesh, patchSeed, patchFilter, selector);

  unsigned nComponents = 1; // scalar field
  unsigned nSamples = patch.size();
  unsigned basisSize = (unsigned)stk::transfer::RecoverField::TRILINEAR;  // 8 basis function evaluations for linear

  stk::transfer::LeastSquares leastSquaresCalculator(nComponents, nSamples, basisSize);

  double x = (3 + 1.0 / 30);
  double y = 1.5;
  double z = 1.05;
  std::vector<double> evalPoint{x, y, z}; // Slightly outside domain so it's an extrapolation

  double interpolatedValue;
  const stk::transfer::EntityInterpolationData interpData(*mesh, &transferField, 0, 0, meta.coordinate_field(),
                                                           patchSeed, part, &selector);
  bool solvable = stk::transfer::least_squares_linear_interpolation(patch, leastSquaresCalculator, evalPoint,
                                                                    interpData, 1u, &interpolatedValue);

  EXPECT_TRUE(solvable);
  double expectedValue = evaluate_tripolynomial_function(triLinearCoeffs, x, y, z);
  EXPECT_NEAR(expectedValue, interpolatedValue, 1.0e-6);
}
//END


//BEGIN
TEST(StkTransferHowTo, useElementCentroidLinearMovingLeastSquaresInterpolation)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 1) { GTEST_SKIP(); }

  const std::string meshSpec("generated:3x3x3");
  double initVals = std::numeric_limits<double>::max();
  const unsigned spatialDim = 3;

  stk::mesh::MeshBuilder builder(communicator);
  builder.set_spatial_dimension(spatialDim);
  std::shared_ptr<stk::mesh::BulkData> mesh = builder.create();
  stk::mesh::MetaData& meta = mesh->mesh_meta_data();
  stk::mesh::Field<double> &transferField = meta.declare_field<double>(stk::topology::ELEM_RANK, "transfer_field", 1);
  stk::mesh::put_field_on_mesh(transferField, meta.universal_part(), &initVals);
  stk::io::fill_mesh(meshSpec, *mesh);

  // Populate element field with tri-linear spatial function based on element centroid
  const std::vector<std::vector<double> > triLinearCoeffs = {{1, 1}, {1, 1}, {1, 1}};
  set_element_field_vals_from_polynomial_coefficients(*mesh, transferField, triLinearCoeffs);

  stk::mesh::EntityId patchSeedId = 14u; // Middle element
  stk::mesh::Entity patchSeed = mesh->get_entity(stk::topology::ELEM_RANK, patchSeedId);
  stk::mesh::Part* part = meta.get_part("block_1");
  stk::mesh::Selector selector(meta.universal_part());
  stk::transfer::EntityPatchFilter patchFilter(part, &selector);
  stk::transfer::LinearPatch<stk::transfer::EntityPatchFilter> patch(*mesh, patchSeed, patchFilter, selector);

  unsigned nComponents = 1; // scalar field
  unsigned nSamples = patch.size();
  unsigned basisSize = (unsigned)stk::transfer::RecoverField::TRILINEAR;  // 8 basis function evaluations for linear

  double x = (3 + 1.0 / 30);
  double y = 1.5;
  double z = 1.05;
  std::vector<double> evalPoint{x, y, z}; // Slightly outside domain so it's an extrapolation

  stk::transfer::GeometricMovingLeastSquares leastSquaresCalculator(nComponents, nSamples, basisSize,
                                                                    *mesh, patch.get_patch_entities(),
                                                                    meta.coordinate_field(), evalPoint);


  double interpolatedValue;
  const stk::transfer::EntityInterpolationData interpData(*mesh, &transferField, 0, 0, meta.coordinate_field(),
                                                           patchSeed, part, &selector);
  bool solvable = stk::transfer::least_squares_linear_interpolation(patch, leastSquaresCalculator, evalPoint,
                                                                    interpData, 1u, &interpolatedValue);

  EXPECT_TRUE(solvable);
  double expectedValue = evaluate_tripolynomial_function(triLinearCoeffs, x, y, z);
  EXPECT_NEAR(expectedValue, interpolatedValue, 1.0e-6);
}
//END

//BEGIN
TEST(StkTransferHowTo, useElementCentroidQuadraticLeastSquaresInterpolation)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 1) { GTEST_SKIP(); }

  const std::string meshSpec("generated:3x3x3");
  double initVals = std::numeric_limits<double>::max();
  const unsigned spatialDim = 3;

  stk::mesh::MeshBuilder builder(communicator);
  builder.set_spatial_dimension(spatialDim);
  std::shared_ptr<stk::mesh::BulkData> mesh = builder.create();
  stk::mesh::MetaData& meta = mesh->mesh_meta_data();
  stk::mesh::Field<double> &transferField = meta.declare_field<double>(stk::topology::ELEM_RANK, "transfer_field", 1);
  stk::mesh::put_field_on_mesh(transferField, meta.universal_part(), &initVals);
  stk::io::fill_mesh(meshSpec, *mesh);

  // Populate element field with tri-quadratic spatial function based on element centroid
  const std::vector<std::vector<double> > triQuadraticCoeffs = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
  set_element_field_vals_from_polynomial_coefficients(*mesh, transferField, triQuadraticCoeffs);

  stk::mesh::EntityId patchSeedId = 14u; // Middle element
  stk::mesh::Entity patchSeed = mesh->get_entity(stk::topology::ELEM_RANK, patchSeedId);
  stk::mesh::Part* part = meta.get_part("block_1");
  stk::mesh::Selector selector(meta.universal_part());
  stk::transfer::EntityPatchFilter patchFilter(part, &selector);
  stk::transfer::QuadraticPatch<stk::transfer::EntityPatchFilter> patch(*mesh, patchSeed, patchFilter, selector);

  unsigned nComponents = 1; // scalar field
  unsigned nSamples = patch.size();
  unsigned basisSize = (unsigned)stk::transfer::RecoverField::TRIQUADRATIC;  // 27 basis function evaluations for quadratic

  stk::transfer::LeastSquares leastSquaresCalculator(nComponents, nSamples, basisSize);

  double x = (3 + 1.0 / 30);
  double y = 1.5;
  double z = 1.05;
  std::vector<double> evalPoint{x, y, z}; // Slightly outside domain so it's an extrapolation

  double interpolatedValue;
  const stk::transfer::EntityInterpolationData interpData(*mesh, &transferField, 0, 0, meta.coordinate_field(),
                                                           patchSeed, part, &selector);
  bool solvable = stk::transfer::least_squares_quadratic_interpolation(patch, leastSquaresCalculator, evalPoint,
                                                                       interpData, 1u, &interpolatedValue);

  EXPECT_TRUE(solvable);
  double expectedValue = evaluate_tripolynomial_function(triQuadraticCoeffs, x, y, z);
  EXPECT_NEAR(expectedValue, interpolatedValue, 1.0e-6);
}
//END

//BEGIN
TEST(StkTransferHowTo, useElementCentroidCubicLeastSquaresInterpolation)
{
  MPI_Comm communicator = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(communicator) != 1) { GTEST_SKIP(); }

  const std::string meshSpec("generated:5x5x5");
  double initVals = std::numeric_limits<double>::max();
  const unsigned spatialDim = 3;

  stk::mesh::MeshBuilder builder(communicator);
  builder.set_spatial_dimension(spatialDim);
  std::shared_ptr<stk::mesh::BulkData> mesh = builder.create();
  stk::mesh::MetaData& meta = mesh->mesh_meta_data();
  stk::mesh::Field<double> &transferField = meta.declare_field<double>(stk::topology::ELEM_RANK, "transfer_field", 1);
  stk::mesh::put_field_on_mesh(transferField, meta.universal_part(), &initVals);
  stk::io::fill_mesh(meshSpec, *mesh);

  // Populate element field with tri-cubic spatial function based on element centroid
  const std::vector<std::vector<double> > triCubicCoeffs = {{1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}};
  set_element_field_vals_from_polynomial_coefficients(*mesh, transferField, triCubicCoeffs);

  stk::mesh::EntityId patchSeedId = 63u; // Middle element
  stk::mesh::Entity patchSeed = mesh->get_entity(stk::topology::ELEM_RANK, patchSeedId);
  stk::mesh::Part* part = meta.get_part("block_1");
  stk::mesh::Selector selector(meta.universal_part());
  stk::transfer::EntityPatchFilter patchFilter(part, &selector);
  stk::transfer::CubicPatch<stk::transfer::EntityPatchFilter> patch(*mesh, patchSeed, patchFilter, selector);

  unsigned nComponents = 1; // scalar field
  unsigned nSamples = patch.size();
  unsigned basisSize = (unsigned)stk::transfer::RecoverField::TRICUBIC;  // 64 basis function evaluations for quadratic

  stk::transfer::LeastSquares leastSquaresCalculator(nComponents, nSamples, basisSize);

  double x = (5 + 1.0 / 50);
  double y = 1.5;
  double z = 1.05;
  std::vector<double> evalPoint{x, y, z}; // Slightly outside domain so it's an extrapolation

  double interpolatedValue;
  const stk::transfer::EntityInterpolationData interpData(*mesh, &transferField, 0, 0, meta.coordinate_field(),
                                                           patchSeed, part, &selector);
  bool solvable = stk::transfer::least_squares_cubic_interpolation(patch, leastSquaresCalculator, evalPoint,
                                                                   interpData, 1u, &interpolatedValue);

  EXPECT_TRUE(solvable);
  double expectedValue = evaluate_tripolynomial_function(triCubicCoeffs, x, y, z);
  EXPECT_NEAR(expectedValue, interpolatedValue, 5.0e-4);
}
//END


} // namespace


