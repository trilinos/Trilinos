// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_SMOOTHING_AKRI_SETHHILLCONFIG_HPP_
#define KRINO_KRINO_SMOOTHING_AKRI_SETHHILLCONFIG_HPP_
#include <string>
#include <vector>

#include <Akri_OptimizeTypes.hpp>

namespace krino {

struct SethHillParams
{
  double mu;      // Shear modulus
  double kappa;   // Bulk modulus
  int m;          // Bulk exponent
  int n;          // Shear exponent
};

enum class SmoothReference { EqualVolume, AvgEdgeLength, Max };

struct DirichletBC
{
  std::string nodeSet;
  int component;  // 0=x, 1=y, 2=z
};


struct SethHillConfig
{
  // I/O
  std::string inputMeshFile;
  std::string outputMeshFile;

  // Material
  SethHillParams params{0.0, 0.0, 0, 0};
  SmoothReference smoothReference = SmoothReference::EqualVolume;

  SolverConfig solverConfig;

  // Boundary conditions (node set + constrained component)
  // Set in the test body since Teuchos cannot parse YAML sequences of maps
  std::vector<DirichletBC> boundaryConditions;
};

std::string solver_name(OptimizerType type);
SethHillConfig parse_norma_yaml(const std::string & filename);

}

#endif /* KRINO_KRINO_SMOOTHING_AKRI_SETHHILLCONFIG_HPP_ */
