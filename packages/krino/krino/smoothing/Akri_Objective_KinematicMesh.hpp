// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#ifndef KRINO_KRINO_SMOOTHING_AKRI_OBJECTIVE_KINEMATICMESH_HPP_
#define KRINO_KRINO_SMOOTHING_AKRI_OBJECTIVE_KINEMATICMESH_HPP_

#include <Akri_Objective_Mesh.hpp>

#include <vector>

namespace krino {

class KinematicElementObjective
{
public:
  virtual ~KinematicElementObjective() {}

  virtual double compute_element_objective(
      const unsigned spatialDim,
      const unsigned npe,
      const double refSize,
      const double * const * elemNodeCoords) const = 0;

  virtual void fill_element_sensitivity(
      const unsigned spatialDim,
      const unsigned npe,
      const double refSize,
      const double * const * elemNodeCoords,
      std::vector<double> & elemGradContrib) const = 0;

  // Anisotropic forms, with the element's ideal shape and size given by a metric M = F_M^T F_M
  // through its metric mapping F_M (metric::compute_metric_mapping).  The energy is
  // integrated over the element's ideal reference configuration X = F_M^-1 Y.  Tri3 in 2D and tet4 in 3D.
  virtual double compute_element_objective(
      const unsigned spatialDim,
      const unsigned npe,
      const double * metricMapping,
      const double * const * elemNodeCoords) const = 0;

  virtual void fill_element_sensitivity(
      const unsigned spatialDim,
      const unsigned npe,
      const double * metricMapping,
      const double * const * elemNodeCoords,
      std::vector<double> & elemGradContrib) const = 0;
};

class KinematicMeshObjective : public MeshObjectiveBase
{
public:
  KinematicMeshObjective(
      const stk::mesh::BulkData & mesh,
      const CoordinatesFieldRef coordsField,
      const stk::mesh::Selector & elemSelector,
      const KinematicElementObjective & elemSmoothingObj)
    : MeshObjectiveBase(mesh, coordsField, elemSelector),
      myElemObj(elemSmoothingObj)
  {}

  virtual ~KinematicMeshObjective() {}

protected:
  const KinematicElementObjective & get_element_objective() const { return myElemObj; }

private:
  const KinematicElementObjective & myElemObj;
};

class KinematicIsotropicMeshObjective final : public KinematicMeshObjective
{
public:
  KinematicIsotropicMeshObjective(
      const stk::mesh::BulkData & mesh,
      const CoordinatesFieldRef coordsField,
      const stk::mesh::Selector & elemSelector,
      const KinematicElementObjective & elemSmoothingObj);

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

  // Ideal size of each locally owned element, in the order of get_owned_elements().  See
  // compute_element_ref_sizes_from_reference.
  void set_element_ref_sizes(const std::vector<double> & elemRefSizes);

private:
  std::vector<double> myElemRefSizes;
};

class KinematicConstantIsotropicMeshObjective final : public KinematicMeshObjective
{
public:
  KinematicConstantIsotropicMeshObjective(
      const stk::mesh::BulkData & mesh,
      const CoordinatesFieldRef coordsField,
      const stk::mesh::Selector & elemSelector,
      const KinematicElementObjective & elemSmoothingObj,
      const double refSize);

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
  double myRefSize;
};

class KinematicAnisotropicMeshObjective final : public KinematicMeshObjective
{
public:
  KinematicAnisotropicMeshObjective(
      const stk::mesh::BulkData & mesh,
      const CoordinatesFieldRef coordsField,
      const stk::mesh::Selector & elemSelector,
      const KinematicElementObjective & elemSmoothingObj);

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

  // Metric mapping F_M of each locally owned element, in the order of get_owned_elements(), with the
  // spatialDim*spatialDim components of each F_M in row-major order.  See
  // compute_element_metric_mappings_from_nodal_metric_field.
  void set_element_metric_map_data(const std::vector<double> & elemMetricMapData);

  // As above, from the metric M at each of the nodes of the objective (get_nodes()), in that order, with
  // the spatialDim*spatialDim components of each M in row-major order.  See
  // compute_element_metric_mappings_from_nodal_metrics.
  void set_element_metric_map_data_from_nodal_metrics(const std::vector<double> & nodalMetrics);

private:
  std::vector<double> myElemMetricMapData;
};


class KinematicConstantAnisotropicMeshObjective final : public KinematicMeshObjective
{
public:
  KinematicConstantAnisotropicMeshObjective(
      const stk::mesh::BulkData & mesh,
      const CoordinatesFieldRef coordsField,
      const stk::mesh::Selector & elemSelector,
      const KinematicElementObjective & elemSmoothingObj,
      const double * metricMapping);

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
  std::vector<double> myMetricMapData;
};

}

#endif /* KRINO_KRINO_SMOOTHING_AKRI_OBJECTIVE_KINEMATICMESH_HPP_ */
