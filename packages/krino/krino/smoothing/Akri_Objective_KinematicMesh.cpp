// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#include <Akri_Objective_KinematicMesh.hpp>

#include <Akri_KinematicUtils.hpp>
#include <Akri_Smoothing_Utils.hpp>

namespace krino {

void KinematicMeshObjective::set_element_ref_size(const CoordinatesFieldRef coordsField, const SmoothReference & reference)
{
  // Precompute ideal reference configurations from initial mesh (computed once)
  DistributedVector xInitial = gather_node_coordinates(get_mesh(), coordsField, get_nodes());

  std::vector<const double *> elemNodeCoords;
  const size_t numElems = get_owned_elements().size();
  myElemRefSizes.resize(numElems, 0.);
  for (size_t e = 0; e < numElems; ++e)
  {
    fill_elem_node_locations(get_spatial_dim(), xInitial, get_element_node_indices(), e, elemNodeCoords);
    const double * const * elemNodeCoordData = elemNodeCoords.data();
    myElemRefSizes[e] = kinematicUtils::compute_ideal_tet4_h(elemNodeCoordData, reference);
  }
}

void KinematicMeshObjective::set_constant_element_ref_size(const double refSize)
{
  const size_t numElems = get_owned_elements().size();
  myElemRefSizes.assign(numElems, refSize);
}

void KinematicMeshObjective::set_element_ref_sizes(const std::vector<double> & elemRefSizes)
{
  STK_ThrowRequireMsg(elemRefSizes.size() == get_owned_elements().size(), "Element reference size vector has size that does not match the number of locally owned elements.");
  myElemRefSizes = elemRefSizes;
}

double KinematicMeshObjective::compute_element_objective(
    const unsigned spatialDim,
    const unsigned npe,
    const size_t elemIndex,
    const double * const * elemNodeCoords) const
{
  return myElemObj.compute_element_objective(get_spatial_dim(), npe, myElemRefSizes[elemIndex], elemNodeCoords);
}

void KinematicMeshObjective::fill_element_sensitivity(
    const unsigned spatialDim,
    const unsigned npe,
    const size_t elemIndex,
    const double * const * elemNodeCoords,
    std::vector<double> & elemGradContrib) const
{
  myElemObj.fill_element_sensitivity(get_spatial_dim(), npe, myElemRefSizes[elemIndex], elemNodeCoords, elemGradContrib);
}

}
