// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_SMOOTHING_AKRI_KINEMATICOBJECTIVES_HPP_
#define KRINO_KRINO_SMOOTHING_AKRI_KINEMATICOBJECTIVES_HPP_

#include <memory>

#include <Akri_Objective_KinematicMesh.hpp>
#include <Akri_SethHillConfig.hpp>

namespace krino
{

// The kinematic element objectives of the material models: the element energy and its sensitivity
// to the node coordinates, for tri3 elements in 2D and tet4 elements in 3D.
std::shared_ptr<KinematicElementObjective> make_size_shape_objective();
std::shared_ptr<KinematicElementObjective> make_modified_neo_hookean_objective();
std::shared_ptr<KinematicElementObjective> make_seth_hill_objective(const SethHillParams & params);
std::shared_ptr<KinematicElementObjective> make_simplified_seth_hill_objective();

} // namespace krino

#endif /* KRINO_KRINO_SMOOTHING_AKRI_KINEMATICOBJECTIVES_HPP_ */
