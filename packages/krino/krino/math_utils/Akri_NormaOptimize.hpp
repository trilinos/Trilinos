// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_MATH_UTILS_AKRI_NORMAOPTIMIZE_HPP_
#define KRINO_KRINO_MATH_UTILS_AKRI_NORMAOPTIMIZE_HPP_

#include <functional>
#include <Akri_OptimizeTypes.hpp>
#include <Akri_ObjectiveInterface.hpp>

namespace krino {

void norma_steepest_descent_solve(
    DistributedVector & nodeCoords,
    const ObjectiveInterface & objFn,
    const SolverConfig & cfg,
    const bool doOutput = false);

std::function<void(const ObjectiveInterface &objFn, DistributedVector &x)>
build_norma_steepest_descent_optimizer_with_parameters(const SolverConfig & cfg,
    const bool doOutput = false);

}

#endif /* KRINO_KRINO_MATH_UTILS_AKRI_NORMAOPTIMIZE_HPP_ */
