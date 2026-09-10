// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Seth-Hill hyperelastic mesh smoothing prototype.
// Ports the mesh smoothing method from Norma.jl
// (https://github.com/sandialabs/Norma.jl) into the krino unit test
// framework, using MiniTensor for tensor algebra.
//
// For each tet4 element, an ideal regular tetrahedron reference configuration
// is constructed from the initial mesh (using Norma's "max(equal_volume,
// average_edge_length)" rule), and the Seth-Hill hyperelastic energy
// measures how far the current element is from that ideal shape.  The
// gradient of the total energy with respect to nodal coordinates gives
// internal forces that drive nodes toward better element quality.
//
// Configuration is read from a Norma-format YAML input file.  Pass the
// input path via "--input=<path>" on the krino_unit command line; default
// is "awful-cube-norma-original.yaml" in the current working directory.
//
// Available solvers (selected by "solver.type" in the YAML):
//   "steepest descent" - Matches Norma.jl exactly (direct port).
//   "rol"              - ROL Quasi-Newton L-BFGS with backtracking line
//                        search.  Uses initial energy for scaling.
//
// YAML structure parsed (Norma.jl format):
//   input mesh file:  <path>
//   output mesh file: <path>
//   model:
//     smooth reference: max | equal volume | average edge length
//     material:
//       blocks: { <block-name>: <material-name> }
//       <material-name>:
//         m: <int>
//         n: <int>
//         bulk modulus:  <double>
//         shear modulus: <double>
//   time integrator:
//     initial time, final time, time step  (produces num time steps)
//   solver:
//     type: steepest descent | rol
//     minimum iterations, maximum iterations
//     absolute tolerance, relative tolerance
//     step length
//     line search backtrack factor, decrease factor, maximum iterations
//
// Boundary conditions are hardcoded in the test body for the awful-cube
// node sets because Teuchos' YAML parser does not support YAML sequences
// of maps (the "- node set: ..." syntax Norma uses).

#include <Akri_MeshFromFileFixture.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_MeshQuality.hpp>
#include <Akri_QualityMetric.hpp>
#include <Akri_SmoothMesh.hpp>
#include <Akri_LogRedirecter.hpp>
#include <Akri_OutputUtils.hpp>
#include <Akri_AllReduce.hpp>
#include <stk_unit_test_utils/CommandLineArgs.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <vector>

#include <Akri_Objective_SethHill.hpp>
#include <Akri_Objective_SizeShape.hpp>
#include <Akri_Smoothing_Utils.hpp>
#include <Akri_SethHillConfig.hpp>
#include <Akri_NormaOptimize.hpp>
#include <Akri_Objective_MeanRatioNorm.hpp>
#include <Akri_Optimize.hpp>
#include <Akri_ROLOptimize.hpp>
#include <stk_util/parallel/ParallelVectorConcat.hpp>

namespace krino {

namespace {

// ==========================================================================
// Boundary constraint: zero gradient in constrained directions only
// Matches Norma.jl: each Dirichlet BC constrains one component (x=0, y=1, z=2)
// on a node set.  A node on a face has 1 component fixed, on an edge 2, on a
// corner 3.
// ==========================================================================

std::function<void(const stk::mesh::Entity, const stk::math::Vector3d &, stk::math::Vector3d &)>
build_boundary_node_filter(const stk::mesh::BulkData & mesh,
    const std::vector<DirichletBC> & bcs)
{
  // For each node, store a bitmask of constrained components (bit 0=x, 1=y, 2=z)
  std::map<stk::mesh::Entity, unsigned> constrainedComponents;
  const stk::mesh::MetaData & meta = mesh.mesh_meta_data();
  for (const auto & bc : bcs)
  {
    const stk::mesh::Part * part = meta.get_part(bc.nodeSet);
    if (part != nullptr)
    {
      std::vector<stk::mesh::Entity> nodes;
      stk::mesh::get_entities(
          mesh, stk::topology::NODE_RANK, stk::mesh::Selector(*part), nodes);
      for (const auto & node : nodes)
        constrainedComponents[node] |= (1u << bc.component);
    }
  }

  return [constrainedComponents](const stk::mesh::Entity node, const stk::math::Vector3d & x, stk::math::Vector3d & dir)
  {
    auto it = constrainedComponents.find(node);
    if (it != constrainedComponents.end())
    {
      const unsigned mask = it->second;
      for (int c = 0; c < 3; ++c)
        if (mask & (1u << c))
          dir[c] = 0.0;
    }
  };
}

// ==========================================================================
// Solver dispatch table — extensible registry of available solvers
// ==========================================================================

using SolverFn = std::function<void(
    DistributedVector &,
    const ObjectiveInterface &,
    const SolverConfig &,
    const bool doOutput)>;

void seth_hill_rol_solve(
    DistributedVector & x,
    const ObjectiveInterface & objFn,
    const SolverConfig & cfg,
    const bool doOutput)
{
  const double xTol = 1.e-14;
  const int totalIterations = cfg.numTimeSteps * cfg.maxIterations;
  const double initialEnergy = objFn.compute_value(x);
  rol_optimize(objFn, x, xTol, cfg.absoluteTolerance, initialEnergy, totalIterations, doOutput);
}

void lbfgs_solve(
    DistributedVector & x,
    const ObjectiveInterface & objFn,
    const SolverConfig & cfg,
    const bool doOutput)
{
  const double xTol = 1.e-4;
  const int totalIterations = cfg.numTimeSteps * cfg.maxIterations;
  const unsigned numLevels = 10;
  const double initialEnergy = objFn.compute_value(x);
  lbfgs(objFn, x, xTol, cfg.absoluteTolerance, initialEnergy, totalIterations, numLevels, doOutput);
}

void steepest_descent_solve(
    DistributedVector & x,
    const ObjectiveInterface & objFn,
    const SolverConfig & cfg,
    const bool doOutput)
{
  const double xTol = 1.e-4;
  const int totalIterations = cfg.numTimeSteps * cfg.maxIterations;
  const double initialEnergy = objFn.compute_value(x);
  steepest_descent(objFn, x, xTol, cfg.absoluteTolerance, initialEnergy, totalIterations, doOutput);
}

const std::map<OptimizerType, SolverFn> & solver_registry()
{
  static const std::map<OptimizerType, SolverFn> registry = {
    { OptimizerType::NormaSteepestDescent,
      [](DistributedVector & x, const ObjectiveInterface & objFn,
         const SolverConfig & cfg, const bool doOutput)
      {
        norma_steepest_descent_solve(x, objFn, cfg, doOutput);
      }
    },
    { OptimizerType::ROL, seth_hill_rol_solve },
    { OptimizerType::SteepestDescent, steepest_descent_solve },
    { OptimizerType::LBFGS, lbfgs_solve },
  };
  return registry;
}

// Look for a "--input=path" argument in the global CLI args.
// Returns the path if found, otherwise returns the provided default.
// Note: GlobalCommandLineArguments captures argc before gtest's InitGoogleTest
// shifts the array, so we walk until argv[i] is null rather than trusting argc.
std::string input_path_from_cli(const std::string & defaultPath)
{
  auto & cliArgs = stk::unit_test_util::GlobalCommandLineArguments::self();
  const int argc = cliArgs.get_argc();
  char ** argv = cliArgs.get_argv();
  const std::string prefix = "--input=";
  for (int i = 1; i < argc && argv[i] != nullptr; ++i)
  {
    const std::string arg = argv[i];
    if (arg.compare(0, prefix.size(), prefix) == 0)
      return arg.substr(prefix.size());
  }
  return defaultPath;
}

} // anonymous namespace

// ==========================================================================
// Test: Smooth awful-cube.g using Seth-Hill energy, reproducing Norma.jl
// ==========================================================================

using SethHillSmoother = MeshFromFileFixture;

TEST_F(SethHillSmoother, smoothAwfulCubeUsingSethHill)
{
  // Load configuration from a Norma.jl YAML input file.  Default is
  // awful-cube-norma-original.yaml, override via "--input=path" on the
  // krino_unit command line.
  const std::string inputFile = input_path_from_cli("awful-cube-norma-original.yaml");
  if (!does_file_exist(inputFile))
  {
    std::cout << "Skipping test: " << inputFile << " not found.\n";
    return;
  }
  SethHillConfig cfg = parse_norma_yaml(inputFile);
  if (0 == parallel_rank())
    std::cout << "Loaded config from " << inputFile << std::endl;

  // Read input mesh (name from YAML)
  if (cfg.inputMeshFile.empty())
  {
    std::cout << "Skipping test: input mesh file not specified in YAML.\n";
    return;
  }
  if (!read_mesh_if_present_and_supported(cfg.inputMeshFile))
  {
    std::cout << "Skipping test: " << cfg.inputMeshFile << " not found.\n";
    return;
  }

  LogRedirecter log;

  // Boundary conditions (Teuchos cannot parse YAML sequences of maps)
  cfg.boundaryConditions = {
    {"nsx-", 0}, {"nsx+", 0},
    {"nsy-", 1}, {"nsy+", 1},
    {"nsz-", 2}, {"nsz+", 2}
  };

  const stk::mesh::BulkData & mesh = get_mesh();
  const CoordinatesFieldRef coordsField = get_coordinates_field();
  const stk::mesh::Selector elemSelector(get_aux_meta().active_part());

  // Report configuration
  if (0 == parallel_rank())
  {
    std::cout << "Seth-Hill parameters: mu=" << cfg.params.mu
              << " kappa=" << cfg.params.kappa
              << " m=" << cfg.params.m
              << " n=" << cfg.params.n << std::endl;
    std::cout << "Solver: " << solver_name(cfg.solverConfig.solverType)
              << ", step_length=" << cfg.solverConfig.stepLength
              << " max_iter=" << cfg.solverConfig.maxIterations
              << " time_steps=" << cfg.solverConfig.numTimeSteps << std::endl;
    std::cout << "Output mesh: " << cfg.outputMeshFile << std::endl;
  }

  SethHillSmoothingObjective elemSmoothingObj(cfg.params);
  KinematicMeshObjective meshSmoothingObj(mesh, coordsField, elemSelector, elemSmoothingObj);

  // Constrain boundary nodes (normal direction only, per Norma.jl convention)
  const auto boundaryFilter = build_boundary_node_filter(mesh,
      cfg.boundaryConditions);
  meshSmoothingObj.set_gradient_filter(boundaryFilter);
  meshSmoothingObj.set_element_ref_size(coordsField, cfg.smoothReference);

  // Dispatch to selected solver via the solver registry
  const auto & registry = solver_registry();
  const auto solverIt = registry.find(cfg.solverConfig.solverType);
  ASSERT_NE(solverIt, registry.end()) << "Solver type not registered";
  if (0 == parallel_rank())
    std::cout << "Using " << solver_name(cfg.solverConfig.solverType) << " solver" << std::endl;

  const bool doSolverOutput = (0 == parallel_rank());
  const auto optimizer = [solverIt, &cfg, doSolverOutput](const ObjectiveInterface &objFn, DistributedVector &x)
    {
      solverIt->second(x, objFn, cfg.solverConfig, doSolverOutput);
    };

  const bool doPrint = true;
  smooth_mesh_by_optimization(get_mesh(),
          get_coordinates_field(),
          optimizer,
          meshSmoothingObj,
          doPrint);

  const std::string outputFile =
      cfg.outputMeshFile.empty() ? "awful-cube-smoothed.e" : cfg.outputMeshFile;
  output_composed_mesh_with_fields(
      get_mesh(), mesh.mesh_meta_data().universal_part(),
      outputFile, 1, 0.);

  if (0 == parallel_rank())
    std::cout << log.get_log() << std::endl;
}

} // namespace krino
