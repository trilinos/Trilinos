// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#ifndef KRINO_KRINO_SMOOTHING_AKRI_OBJECTIVE_MESH_HPP_
#define KRINO_KRINO_SMOOTHING_AKRI_OBJECTIVE_MESH_HPP_

#include <stk_math/StkVector.hpp>
#include <Akri_ObjectiveInterface.hpp>
#include <Akri_Smoothing_Utils.hpp>

namespace krino {

class NodeCoordinateObjectiveInterface
{
public:
  virtual ~NodeCoordinateObjectiveInterface() {}
  virtual double compute_value(const stk::mesh::BulkData & mesh,
      const CoordinatesFieldRef coordsField,
      const stk::mesh::Selector & elemSelector,
      const stk::mesh::Entity node,
      const stk::math::Vector3d & x) const = 0;

  virtual void fill_gradient(const stk::mesh::BulkData & mesh,
      const CoordinatesFieldRef coordsField,
      const stk::mesh::Selector & elemSelector,
      const GradientFilter & gradientFilter,
      const stk::mesh::Entity node,
      const stk::math::Vector3d & x,
      stk::math::Vector3d & g) const = 0;
};

class NodeCoordinateObjective : public Objective3DInterface
{
public:
  NodeCoordinateObjective(const NodeCoordinateObjectiveInterface & nodeCoordObj,
    const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const GradientFilter & gradientFilter,
    const stk::mesh::Entity node);

  virtual double compute_value(const stk::math::Vector3d & x) const override;
  virtual void fill_gradient(const stk::math::Vector3d & x,
      stk::math::Vector3d & g) const override;

private:
  const NodeCoordinateObjectiveInterface & myNodeCoordObj;
  const stk::mesh::BulkData & myMesh;
  const CoordinatesFieldRef myCoordsField;
  const stk::mesh::Selector & myElemSelector;
  GradientFilter myGradientFilter;
  stk::mesh::Entity myNode;
};

class ElementObjective
{
public:
  virtual ~ElementObjective() {}
  virtual double compute_element_objective(
      const unsigned spatialDim,
      const unsigned npe,
      const double * const * elemNodeCoords) const = 0;

  virtual void fill_element_sensitivity(
      const unsigned spatialDim,
      const unsigned npe,
      const double * const * elemNodeCoords,
      std::vector<double> & elemGradContrib) const = 0;
};

class MeshObjectiveBase : public ObjectiveInterface
{
public:
  MeshObjectiveBase(const stk::mesh::BulkData &mesh,
      const CoordinatesFieldRef coordsField,
      const stk::mesh::Selector & elemSelector);

  virtual void project_solution(DistributedVector & x) const override;
  virtual double compute_value(const DistributedVector & x) const override;
  virtual void fill_gradient(const DistributedVector & x, DistributedVector & g) const override;

  virtual double compute_element_objective(
      const unsigned spatialDim,
      const unsigned npe,
      const size_t elemIndex,
      const double * const * elemNodeCoords) const = 0;
  virtual void fill_element_sensitivity(
      const unsigned spatialDim,
      const unsigned npe,
      const size_t elemIndex,
      const double * const * elemNodeCoords,
      std::vector<double> & elemGradContrib) const = 0;

  void set_gradient_filter(const GradientFilter & gradientFilter)
    { myGradientFilter = gradientFilter; }
  void set_solution_projector(const SolutionProjector & solutionProjector)
      { mySolutionProjector = solutionProjector; }

  unsigned get_spatial_dim() const { return mySpatialDim; }
  const DistributedEntityVector & get_nodes() const { return myNodes; }
  const std::vector<stk::mesh::Entity> & get_owned_elements() const { return myOwnedElems; }

protected:
  const stk::mesh::BulkData & get_mesh() const { return myMesh; }
  const MeshElemNodeIndices & get_element_node_indices() const { return myMeshElemNodeIndices; }

private:
  void setup_nodes_and_elements(const stk::mesh::Selector & elementSelector);

  const stk::mesh::BulkData & myMesh;
  CoordinatesFieldRef myCoordsField;
  const stk::mesh::Selector & myElemSelector;
  unsigned mySpatialDim;
  MeshElemNodeIndices myMeshElemNodeIndices;
  DistributedEntityVector myNodes;
  std::vector<stk::mesh::Entity> myOwnedElems;
  GradientFilter myGradientFilter;
  SolutionProjector mySolutionProjector;
};

class MeshObjective : public MeshObjectiveBase
{
public:
  MeshObjective(const stk::mesh::BulkData &mesh,
     const CoordinatesFieldRef coordsField,
     const stk::mesh::Selector & elemSelector,
     const ElementObjective & elemSmoothingObj)
  : MeshObjectiveBase(mesh, coordsField, elemSelector), myElemObj(elemSmoothingObj) {}

  virtual double compute_element_objective(
          const unsigned spatialDim,
          const unsigned npe,
          const size_t elemIndex,
          const double * const * elemNodeCoords) const override;
  virtual void fill_element_sensitivity(
      const unsigned spatialDim,
      const unsigned npe,
      const size_t elemIndex,
      const double * const * elemNodeCoords,
      std::vector<double> & elemGradContrib) const override;

private:
  const ElementObjective & myElemObj;
};

}



#endif /* KRINO_KRINO_SMOOTHING_AKRI_OBJECTIVE_MESH_HPP_ */
