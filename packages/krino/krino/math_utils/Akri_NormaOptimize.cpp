// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_NormaOptimize.hpp>

#include <Akri_DistributedVector.hpp>
#include <Akri_OptimizeTypes.hpp>

namespace krino {

// ==========================================================================
// Steepest descent with backtracking line search
// Matches Norma.jl src/solver.jl lines 478-668
// ==========================================================================

// Backtracking line search — Norma.jl solver.jl:478-519
// Returns the accepted step vector.
// On entry, grad contains the gradient at the current position.
// On exit, grad contains the gradient at the accepted trial position,
// and nodeCoords is restored to its original value.
std::vector<double> backtrack_line_search(
    DistributedVector & nodeCoords,
    DistributedVector & grad,
    const std::vector<double> & direction,
    const ObjectiveInterface & objFn,
    const SolverConfig & cfg)
{
  const size_t n = nodeCoords.size();
  const double meritInit = 0.5 * Dot(grad, grad);

  double stepLength = cfg.stepLength;
  DistributedVector xSaved = nodeCoords;

  for (int lsIter = 0; lsIter < cfg.lineSearchMaxIterations; ++lsIter)
  {
    // Trial: x = x_saved + stepLength * direction
    for (size_t i = 0; i < n; ++i)
      nodeCoords[i] = xSaved[i] + stepLength * direction[i];

    // Evaluate gradient at trial point
    objFn.fill_gradient(nodeCoords, grad);

    // Merit = 0.5 * ||g_new||^2
    const double merit = 0.5 * Dot(grad, grad);

    // dot(g_new, direction)
    double dotGradDir = 0.0;
    for (size_t i = 0; i < grad.local_size(); ++i)
      dotGradDir += grad[i] * direction[i];

    // Armijo-like acceptance: merit ≤ merit_init + c * alpha * dot(g_new, d)
    if (merit <= meritInit + cfg.lineSearchDecreaseFactor * stepLength * dotGradDir)
    {
      // Accepted — build step, restore position
      std::vector<double> step(n);
      for (size_t i = 0; i < n; ++i)
        step[i] = stepLength * direction[i];
      nodeCoords = xSaved;
      return step;
    }

    stepLength *= cfg.lineSearchBacktrackFactor;
  }

  // Line search exhausted — return last tried step
  std::vector<double> step(n);
  for (size_t i = 0; i < n; ++i)
    step[i] = stepLength * direction[i];
  nodeCoords = xSaved;
  return step;
}

// Main solver — Norma.jl solver.jl:608-668
void norma_steepest_descent_solve(
    DistributedVector & nodeCoords,
    const ObjectiveInterface & objFn,
    const SolverConfig & cfg,
    const bool doOutput)
{
  for (int timeStep = 0; timeStep < cfg.numTimeSteps; ++timeStep)
  {
    // Evaluate initial gradient (Norma: predict + evaluate)
    DistributedVector grad;
    objFn.fill_gradient(nodeCoords, grad);

    double normResidual = std::sqrt(Dot(grad, grad));
    const double initialNorm = normResidual;

    if (doOutput)
      printf("Time step %d Iteration [%d] |R| = %.2e |r| = %.2e [WAIT]\n",
          timeStep, 0, normResidual, 1.0);

    int iterationNumber = 1;
    while (true)
    {
      // Compute descent direction = normalize(-gradient)
      const size_t n = grad.size();
      const double gradNorm = std::sqrt(Dot(grad, grad));
      std::vector<double> direction(n, 0.0);
      if (gradNorm > 0.0)
        for (size_t i = 0; i < n; ++i)
          direction[i] = -grad[i] / gradNorm;

      // Backtracking line search
      auto step = backtrack_line_search(
          nodeCoords, grad, direction, objFn, cfg);

      // Apply step (Norma: solver.solution[free] += step)
      for (size_t i = 0; i < nodeCoords.size(); ++i)
        nodeCoords[i] += step[i];

      // Re-evaluate gradient (Norma: correct + evaluate)
      objFn.fill_gradient(nodeCoords, grad);

      normResidual = std::sqrt(Dot(grad, grad));
      const double relativeError =
          (initialNorm > 0.0) ? normResidual / initialNorm : normResidual;
      const bool converged =
          (normResidual <= cfg.absoluteTolerance) ||
          (relativeError <= cfg.relativeTolerance);

      if (doOutput)
      {
        const char * status = converged ? "[DONE]" : "[WAIT]";
        printf("Time step %d Iteration [%d] |R| = %.2e |r| = %.2e %s\n",
            timeStep, iterationNumber, normResidual, relativeError, status);
      }

      iterationNumber++;

      // Stop conditions (Norma: stop_solve)
      if (iterationNumber > cfg.minIterations)
      {
        if (normResidual == 0.0 || iterationNumber > cfg.maxIterations || converged)
          break;
      }
    }
  }
}

std::function<void(const ObjectiveInterface &objFn, DistributedVector &x)>
build_norma_steepest_descent_optimizer_with_parameters(const SolverConfig & cfg,
    const bool doOutput)
{
  auto fn = [=](const ObjectiveInterface &objFn, DistributedVector &x) { norma_steepest_descent_solve(x, objFn, cfg, doOutput); };
  return fn;
}

}


