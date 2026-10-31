// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#include <Akri_Objective_KinematicMesh.hpp>

#include <Akri_KinematicMeshUtils.hpp>
#include <Akri_Smoothing_Utils.hpp>

namespace krino {
namespace {

void require_anisotropic_tri3_or_tet4(const unsigned spatialDim, const unsigned npe)
{
  STK_ThrowRequireMsg((spatialDim == 2 && npe == 3) || (spatialDim == 3 && npe == 4),
      "Anisotropic KinematicElementObjective metric mapping form is only supported for tri3 elements in 2D and tet4 elements in 3D.");
}

}

KinematicConstantIsotropicMeshObjective::KinematicConstantIsotropicMeshObjective(
    const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const KinematicElementObjective & elemSmoothingObj,
    const double refSize)
  : KinematicMeshObjective(mesh, coordsField, elemSelector, elemSmoothingObj),
    myRefSize(refSize)
{
}

double KinematicConstantIsotropicMeshObjective::compute_element_objective(
    const unsigned spatialDim,
    const unsigned npe,
    const size_t /*elemIndex*/,
    const double * const * elemNodeCoords) const
{
  return get_element_objective().compute_element_objective(
      spatialDim,
      npe,
      myRefSize,
      elemNodeCoords);
}

void KinematicConstantIsotropicMeshObjective::fill_element_sensitivity(
    const unsigned spatialDim,
    const unsigned npe,
    const size_t /*elemIndex*/,
    const double * const * elemNodeCoords,
    std::vector<double> & elemGradContrib) const
{
  get_element_objective().fill_element_sensitivity(
      spatialDim,
      npe,
      myRefSize,
      elemNodeCoords,
      elemGradContrib);
}

KinematicIsotropicMeshObjective::KinematicIsotropicMeshObjective(
    const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const KinematicElementObjective & elemSmoothingObj)
  : KinematicMeshObjective(mesh, coordsField, elemSelector, elemSmoothingObj)
{
}

void KinematicIsotropicMeshObjective::set_element_ref_sizes(const std::vector<double> & elemRefSizes)
{
  STK_ThrowRequireMsg(elemRefSizes.size() == get_owned_elements().size(), "Element reference size vector has size that does not match the number of locally owned elements.");
  myElemRefSizes = elemRefSizes;
}

double KinematicIsotropicMeshObjective::compute_element_objective(
    const unsigned spatialDim,
    const unsigned npe,
    const size_t elemIndex,
    const double * const * elemNodeCoords) const
{
  STK_ThrowRequireMsg(!myElemRefSizes.empty(),
      "Element sizes have not been set.");
  STK_ThrowRequireMsg(elemIndex < myElemRefSizes.size(),
      "Element index is out of range for element reference sizes.");

  return get_element_objective().compute_element_objective(
      spatialDim,
      npe,
      myElemRefSizes[elemIndex],
      elemNodeCoords);
}

void KinematicIsotropicMeshObjective::fill_element_sensitivity(
    const unsigned spatialDim,
    const unsigned npe,
    const size_t elemIndex,
    const double * const * elemNodeCoords,
    std::vector<double> & elemGradContrib) const
{
  STK_ThrowRequireMsg(elemIndex < myElemRefSizes.size(),
      "Element index is out of range for element reference sizes.");

  get_element_objective().fill_element_sensitivity(
      spatialDim,
      npe,
      myElemRefSizes[elemIndex],
      elemNodeCoords,
      elemGradContrib);
}

KinematicConstantAnisotropicMeshObjective::KinematicConstantAnisotropicMeshObjective(
    const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const KinematicElementObjective & elemSmoothingObj,
    const double * metricMapData)
  : KinematicMeshObjective(mesh, coordsField, elemSelector, elemSmoothingObj)
{
  myMetricMapData.assign(metricMapData, metricMapData+coordsField.dim()*coordsField.dim());
}

double KinematicConstantAnisotropicMeshObjective::compute_element_objective(
    const unsigned spatialDim,
    const unsigned npe,
    const size_t /*elemIndex*/,
    const double * const * elemNodeCoords) const
{
  require_anisotropic_tri3_or_tet4(spatialDim, npe);

  return get_element_objective().compute_element_objective(
      spatialDim,
      npe,
      myMetricMapData.data(),
      elemNodeCoords);
}

void KinematicConstantAnisotropicMeshObjective::fill_element_sensitivity(
    const unsigned spatialDim,
    const unsigned npe,
    const size_t /*elemIndex*/,
    const double * const * elemNodeCoords,
    std::vector<double> & elemGradContrib) const
{
  require_anisotropic_tri3_or_tet4(spatialDim, npe);

  get_element_objective().fill_element_sensitivity(
      spatialDim,
      npe,
      myMetricMapData.data(),
      elemNodeCoords,
      elemGradContrib);
}

KinematicAnisotropicMeshObjective::KinematicAnisotropicMeshObjective(
    const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const stk::mesh::Selector & elemSelector,
    const KinematicElementObjective & elemSmoothingObj)
  : KinematicMeshObjective(mesh, coordsField, elemSelector, elemSmoothingObj)
{
}

void KinematicAnisotropicMeshObjective::set_element_metric_map_data(const std::vector<double> & elemMetricMapData)
{
  const unsigned spatialDim = get_spatial_dim();
  STK_ThrowRequireMsg(elemMetricMapData.size() == get_owned_elements().size()*spatialDim*spatialDim,
      "Element metric mapping vector has size that does not match the number of locally owned elements.");
  myElemMetricMapData = elemMetricMapData;
}

void KinematicAnisotropicMeshObjective::set_element_metric_map_data_from_nodal_metrics(const std::vector<double> & nodalMetrics)
{
  set_element_metric_map_data(compute_element_metric_mappings_from_nodal_metrics(get_mesh(), get_nodes(), get_element_node_indices(), nodalMetrics));
}

double KinematicAnisotropicMeshObjective::compute_element_objective(
    const unsigned spatialDim,
    const unsigned npe,
    const size_t elemIndex,
    const double * const * elemNodeCoords) const
{
  STK_ThrowRequireMsg(elemIndex*spatialDim*spatialDim < myElemMetricMapData.size(),
      "Element index is out of range for element metric mappings.");

  require_anisotropic_tri3_or_tet4(spatialDim, npe);

  return get_element_objective().compute_element_objective(
      spatialDim,
      npe,
      &myElemMetricMapData[elemIndex*spatialDim*spatialDim],
      elemNodeCoords);
}

void KinematicAnisotropicMeshObjective::fill_element_sensitivity(
    const unsigned spatialDim,
    const unsigned npe,
    const size_t elemIndex,
    const double * const * elemNodeCoords,
    std::vector<double> & elemGradContrib) const
{
  STK_ThrowRequireMsg(elemIndex*spatialDim*spatialDim < myElemMetricMapData.size(),
      "Element index is out of range for element metric mappings.");

  require_anisotropic_tri3_or_tet4(spatialDim, npe);

  get_element_objective().fill_element_sensitivity(
      spatialDim,
      npe,
      &myElemMetricMapData[elemIndex*spatialDim*spatialDim],
      elemNodeCoords,
      elemGradContrib);
}

}
