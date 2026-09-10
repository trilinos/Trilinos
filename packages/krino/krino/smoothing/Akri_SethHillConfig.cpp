// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_SethHillConfig.hpp>

#include <Teuchos_YamlParameterListCoreHelpers.hpp>
#include <Teuchos_ParameterList.hpp>
#include <fstream>

namespace krino {

static std::map<OptimizerType, std::string> validSolvers = {
    {OptimizerType::NormaSteepestDescent, "NORMA STEEPEST DESCENT"},
    {OptimizerType::ROL, "ROL L-BFGS"},
    {OptimizerType::SteepestDescent, "KRINO STEEPEST DESCENT"},
    {OptimizerType::LBFGS, "KRINO L-BFGS"}
};

std::string solver_name(OptimizerType type)
{
  return validSolvers.at(type);
}

OptimizerType solver_type(const std::string & name)
{
  std::string upperName = name;
  std::transform(upperName.begin(), upperName.end(), upperName.begin(), ::toupper);
  for (const auto & entry : validSolvers)
    if (entry.second == upperName)
      return entry.first;
  std::cout << "Unrecognized solver \"" << upperName << "\"\n";
  std::cout << "Supported solvers:\n";
  for (const auto & entry : validSolvers)
    std::cout << " " << entry.second << "\n";
  return OptimizerType::NormaSteepestDescent;
}

// ==========================================================================
// YAML configuration reader — parses Norma.jl format
// ==========================================================================

// Teuchos infers scalar types from YAML values: "1.0e+03" becomes int 1000,
// "0.0" becomes int 0.  These helpers read the correct type with a default.
double get_double(const Teuchos::ParameterList & pl, const std::string & name,
    double def_value = 0.0)
{
  if (!pl.isParameter(name))
    return def_value;
  if (pl.isType<double>(name))
    return pl.get<double>(name);
  return static_cast<double>(pl.get<int>(name));
}

int get_int(const Teuchos::ParameterList & pl, const std::string & name,
    int def_value = 0)
{
  if (!pl.isParameter(name))
    return def_value;
  return pl.get<int>(name);
}

SethHillConfig parse_norma_yaml(const std::string & filename)
{
  // Teuchos requires a single top-level ParameterList and cannot parse
  // YAML sequences of maps (the "- key: val" syntax used for BCs).
  // Wrap under a root key and skip the "boundary conditions" block.
  std::ifstream infile(filename);
  std::string wrapped = "norma:\n";
  std::string line;
  bool skipping = false;
  while (std::getline(infile, line))
  {
    // Skip the boundary conditions block (contains YAML sequences)
    if (line.find("boundary conditions:") != std::string::npos)
    {
      skipping = true;
      continue;
    }
    // Stop skipping when we hit the next top-level key
    if (skipping && !line.empty() && line[0] != ' ' && line[0] != '#')
      skipping = false;
    if (!skipping)
      wrapped += "  " + line + "\n";
  }

  Teuchos::ParameterList params;
  Teuchos::updateParametersFromYamlString(wrapped, Teuchos::inoutArg(params), true);
  SethHillConfig cfg;

  if (params.isParameter("input mesh file"))
    cfg.inputMeshFile = params.get<std::string>("input mesh file");
  if (params.isParameter("output mesh file"))
    cfg.outputMeshFile = params.get<std::string>("output mesh file");

  if (params.isSublist("model"))
  {
    const auto & model = params.sublist("model");
    if (model.isParameter("smooth reference"))
    {
      const std::string sr = model.get<std::string>("smooth reference");
      if (sr == "max") cfg.smoothReference = SmoothReference::Max;
      else if (sr == "average edge length") cfg.smoothReference = SmoothReference::AvgEdgeLength;
      else cfg.smoothReference = SmoothReference::EqualVolume;
    }
    if (model.isSublist("material"))
    {
      const auto & material = model.sublist("material");
      // Find elastic material name from blocks map
      std::string elasticName = "elastic";
      if (material.isSublist("blocks"))
      {
        const auto & blocks = material.sublist("blocks");
        for (auto it = blocks.begin(); it != blocks.end(); ++it)
          elasticName = blocks.get<std::string>(it->first);
      }
      if (material.isSublist(elasticName))
      {
        const auto & elastic = material.sublist(elasticName);
        cfg.params.m     = get_int(elastic, "m");
        cfg.params.n     = get_int(elastic, "n");
        cfg.params.kappa = get_double(elastic, "bulk modulus");
        cfg.params.mu    = get_double(elastic, "shear modulus");
      }
    }
  }

  if (params.isSublist("time integrator"))
  {
    const auto & ti = params.sublist("time integrator");
    const double t0 = get_double(ti, "initial time");
    const double tf = get_double(ti, "final time");
    const double dt = get_double(ti, "time step", 1.0);
    cfg.solverConfig.numTimeSteps = static_cast<int>((tf - t0) / dt);
  }

  // Note: boundary conditions use YAML sequences of maps which Teuchos
  // does not support.  BCs are geometry-specific and set in the test body.

  if (params.isSublist("solver"))
  {
    const auto & solver = params.sublist("solver");
    if (solver.isParameter("type"))
    {
      const std::string solverType = solver.get<std::string>("type");
      cfg.solverConfig.solverType = solver_type(solverType);
    }
    cfg.solverConfig.minIterations             = get_int(solver, "minimum iterations");
    cfg.solverConfig.maxIterations             = get_int(solver, "maximum iterations");
    cfg.solverConfig.relativeTolerance         = get_double(solver, "relative tolerance");
    cfg.solverConfig.absoluteTolerance         = get_double(solver, "absolute tolerance");
    cfg.solverConfig.stepLength                = get_double(solver, "step length");
    cfg.solverConfig.lineSearchBacktrackFactor = get_double(solver, "line search backtrack factor");
    cfg.solverConfig.lineSearchDecreaseFactor  = get_double(solver, "line search decrease factor");
    cfg.solverConfig.lineSearchMaxIterations   = get_int(solver, "line search maximum iterations");
  }

  return cfg;
}

}
