// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#include <Akri_Objective_Mesh.hpp>

#include <Akri_AllReduce.hpp>
#include <Akri_Smoothing_Utils.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace krino {

NodeCoordinateObjective::NodeCoordinateObjective(const NodeCoordinateObjectiveInterface & nodeCoordObj,
      const stk::mesh::BulkData & mesh,
      const CoordinatesFieldRef coordsField,
      const stk::mesh::Selector & elemSelector,
      const GradientFilter & gradientFilter,
      const stk::mesh::Entity node)
: myNodeCoordObj(nodeCoordObj),
  myMesh(mesh),
  myCoordsField(coordsField),
  myElemSelector(elemSelector),
  myGradientFilter(gradientFilter),
  myNode(node)
{
}

double NodeCoordinateObjective::compute_value(const stk::math::Vector3d & x) const
{
  return myNodeCoordObj.compute_value(myMesh, myCoordsField, myElemSelector, myNode, x);
}

void NodeCoordinateObjective::fill_gradient(const stk::math::Vector3d & x,
    stk::math::Vector3d & g) const
{
  myNodeCoordObj.fill_gradient(myMesh, myCoordsField, myElemSelector, myGradientFilter, myNode, x, g);
}

MeshObjectiveBase::MeshObjectiveBase(const stk::mesh::BulkData &mesh,
      const CoordinatesFieldRef coordsField,
      const stk::mesh::Selector & elemSelector)
: myMesh(mesh),
  myCoordsField(coordsField),
  myElemSelector(elemSelector),
  mySpatialDim(mesh.mesh_meta_data().spatial_dimension())
{
  setup_nodes_and_elements(myElemSelector);
}

void MeshObjectiveBase::project_solution(DistributedVector & x) const
{
  apply_projection(mySpatialDim, mySolutionProjector, myNodes, x);
}

void MeshObjectiveBase::setup_nodes_and_elements(const stk::mesh::Selector & elementSelector)
{
  myNodes = get_sorted_owned_nodes_and_unowned_shared_nodes(myMesh, elementSelector);

  const stk::mesh::Selector ownedElementSelector = elementSelector & myMesh.mesh_meta_data().locally_owned_part();
  stk::mesh::get_entities(myMesh, stk::topology::ELEMENT_RANK, ownedElementSelector, myOwnedElems);

  myMeshElemNodeIndices.build_for_mesh(myMesh, myOwnedElems, myNodes);
}

double MeshObjectiveBase::compute_value(const DistributedVector & x) const
{
  std::vector<const double *> elemNodeCoords;

  double objective = 0.0;
  for (size_t e = 0; e < myMeshElemNodeIndices.num_elements(); ++e)
  {
    fill_elem_node_locations(mySpatialDim, x, myMeshElemNodeIndices, e, elemNodeCoords);
    objective += float(compute_element_objective(mySpatialDim, elemNodeCoords.size(), e, elemNodeCoords.data())); // reduced order addend to produce parallel consistency
  }
  all_reduce_sum(myMesh.parallel(), objective);

  return objective;
}

void MeshObjectiveBase::fill_gradient(
    const DistributedVector & x,
    DistributedVector & g) const
{
  g.assign(x.comm(), x.size(), x.local_size(), 0.0);

  std::vector<const double *> elemNodeCoords;
  std::vector<double> elemGradContrib;

  for (size_t e = 0; e < myMeshElemNodeIndices.num_elements(); ++e)
  {
    const ElemNodeIndices elemNodeIndices = myMeshElemNodeIndices.element_node_indices(e);
    fill_elem_node_locations(mySpatialDim, x, elemNodeIndices, elemNodeCoords);
    fill_element_sensitivity(mySpatialDim, elemNodeCoords.size(), e, elemNodeCoords.data(), elemGradContrib);
    unsigned i=0;
    for (auto elemNodeIndex : elemNodeIndices)
      for (unsigned j = 0; j < mySpatialDim; ++j)
        g[mySpatialDim * elemNodeIndex + j] += float(elemGradContrib[i++]); // reduced order addend to produce parallel consistency
  }

  communicate_gradient(myMesh, myNodes, g);
  apply_gradient_filter(mySpatialDim, myGradientFilter, myNodes, x, g);
}

double MeshObjective::compute_element_objective(
    const unsigned spatialDim,
    const unsigned npe,
    const size_t elemIndex,
    const double * const * elemNodeCoords) const
{
  return myElemObj.compute_element_objective(get_spatial_dim(), npe, elemNodeCoords);
}

void MeshObjective::fill_element_sensitivity(
    const unsigned spatialDim,
    const unsigned npe,
    const size_t elemIndex,
    const double * const * elemNodeCoords,
    std::vector<double> & elemGradContrib) const
{
  myElemObj.fill_element_sensitivity(get_spatial_dim(), npe, elemNodeCoords, elemGradContrib);
}

}


