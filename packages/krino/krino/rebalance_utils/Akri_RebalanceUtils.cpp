// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_RebalanceUtils.hpp>
#include <Akri_RebalanceUtils_Impl.hpp>
#include <Akri_CDMesh.hpp>
#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_mesh/base/Field.hpp>
#include <string>

namespace krino {
namespace rebalance_utils {

namespace {

class MultipleCriteriaSettings : public stk::balance::GraphCreationSettings
{
public:
  MultipleCriteriaSettings(stk::mesh::BulkData & stkMeshBulkData,
      const std::vector<stk::mesh::Field<double> *> critFields,
      const double default_weight = 0.0)
      : m_stkMeshBulkData(stkMeshBulkData),
        m_critFields(critFields),
        m_defaultWeight(default_weight)
  {
    method = "rcb";
    setUseNodeBalancer(true);
    setNodeBalancerTargetLoadBalance(getImbalanceTolerance());
    setNodeBalancerMaxIterations(10);
  }
  virtual ~MultipleCriteriaSettings() = default;

  virtual double
  getGraphEdgeWeight(stk::topology element1Topology, stk::topology element2Topology) const override
  {
    return 1.0;
  }
  virtual bool areVertexWeightsProvidedViaFields() const override { return true; }
  virtual bool includeSearchResultsInGraph() const override { return false; }
  virtual int getGraphVertexWeight(stk::topology type) const override { return 1; }
  virtual double getImbalanceTolerance() const override { return 1.05; }
  virtual void setDecompMethod(const std::string & input_method) override { method = input_method; }
  virtual std::string getDecompMethod() const override { return method; }
  virtual int getNumCriteria() const override { return m_critFields.size(); }
  virtual bool isMultiCriteriaRebalance() const override { return true; }
  virtual bool shouldFixMechanisms() const override { return false; }

  virtual double getGraphVertexWeight(stk::mesh::Entity entity, int criteria_index) const override
  {
    ThrowRequireWithSierraHelpMsg(
        criteria_index >= 0 && static_cast<size_t>(criteria_index) < m_critFields.size());
    const double * weight = stk::mesh::field_data(*m_critFields[criteria_index], entity);
    if (weight != nullptr)
    {
      ThrowRequireWithSierraHelpMsg(*weight >= 0);
      return *weight;
    }
    else
    {
      return m_defaultWeight;
    }
  }

protected:
  MultipleCriteriaSettings() = delete;
  MultipleCriteriaSettings(const MultipleCriteriaSettings &) = delete;
  MultipleCriteriaSettings & operator=(const MultipleCriteriaSettings &) = delete;

  const stk::mesh::BulkData & m_stkMeshBulkData;
  const std::vector<stk::mesh::Field<double> *> m_critFields;
  const double m_defaultWeight;
};

class CDFEMRebalance final : public MultipleCriteriaSettings
{
public:
  CDFEMRebalance(stk::mesh::BulkData & bulk_data,
      CDMesh * cdmesh,
      const std::string & coordinates_field_name,
      const std::vector<stk::mesh::Field<double> *> & weights_fields,
      const double imbalance_threshold,
      const double default_weight = 0.)
      : MultipleCriteriaSettings(bulk_data, weights_fields, default_weight),
        my_cdmesh(cdmesh),
        my_bulk_data(bulk_data),
        my_coordinates_field_name(coordinates_field_name),
        my_imbalance_threshold(imbalance_threshold)
  {}

  ~CDFEMRebalance() = default;

  double getImbalanceTolerance() const override { return my_imbalance_threshold; }
  std::string getCoordinateFieldName() const override { return my_coordinates_field_name; }

  void modifyDecomposition(stk::balance::DecompositionChangeList & decomp_changes) const override;
  bool shouldPrintMetrics() const override { return true; }
  bool isIncrementalRebalance() const override { return true; }

  virtual double getGraphVertexWeight(stk::mesh::Entity entity, int criteria_index) const override
  {
    double scaleVertexWeightForTestingDueToSmallMesh = 12;
    return scaleVertexWeightForTestingDueToSmallMesh*MultipleCriteriaSettings::getGraphVertexWeight(entity, criteria_index);
  }

private:
  CDMesh * my_cdmesh;
  stk::mesh::BulkData & my_bulk_data;
  std::string my_coordinates_field_name;
  double my_imbalance_threshold;
};

void CDFEMRebalance::modifyDecomposition(stk::balance::DecompositionChangeList & decomp_changes) const
{
  /* Algorithm:
   * 1) Go through change list and update destinations of all adaptivity children to match their
   *    root parent. Also move family trees with the parent+child entities.
   * 2) Go through change list and remove any CDFEM child elements whose parents aren't being moved,
   *    as well as update the destination of all CDFEM children with parents that are being moved to
   *    the parent destination.
   */

  impl::update_rebalance_for_adaptivity(decomp_changes, my_bulk_data);

  if(my_cdmesh)
  {
    impl::update_rebalance_for_cdfem(decomp_changes, my_bulk_data, *my_cdmesh);
  }
}

} // anonymous namespace

void
update_parent_child_rebalance_weights(const stk::mesh::BulkData & bulk_data,
    stk::mesh::Field<double> & element_weights_field,
    const CDMesh * cdmesh)
{
  // First sum CDFEM child weights to their CDFEM parents, then adaptivity children to
  // their parents. Adaptivity intermediate parents should not have their weights
  // summed, they should be 0.
  if(cdmesh)
  {
    impl::accumulate_cdfem_child_weights_to_parents(bulk_data, element_weights_field, *cdmesh);
  }

  impl::accumulate_adaptivity_child_weights_to_parents(bulk_data, element_weights_field);
}

bool rebalance_mesh(stk::mesh::BulkData & bulk_data,
    CDMesh * cdmesh,
    const std::string & element_weights_field_name,
    const std::string & coordinates_field_name,
    const std::vector<stk::mesh::Selector> & selections_to_rebalance_separately,
    const std::string & decomp_method,
    const double imbalance_threshold)
{
  const auto & meta = bulk_data.mesh_meta_data();
  auto weights_base = meta.get_field(stk::topology::ELEMENT_RANK, element_weights_field_name);
  ThrowRequireMsg(weights_base,
      "Failed to find element rank field " << element_weights_field_name
      << " to use for rebalance weights.");
  const auto element_weights_field = static_cast<stk::mesh::Field<double> *>(weights_base);

  update_parent_child_rebalance_weights(bulk_data, *element_weights_field, cdmesh);

  ThrowAssert(impl::check_family_tree_element_and_side_ownership(bulk_data));

  CDFEMRebalance balancer(
      bulk_data, cdmesh, coordinates_field_name, {element_weights_field}, imbalance_threshold);
  balancer.setDecompMethod(decomp_method);
  const bool rebalanced =
          stk::balance::balanceStkMesh(balancer, bulk_data, selections_to_rebalance_separately);

  ThrowAssert(impl::check_family_tree_element_and_side_ownership(bulk_data));

  if(cdmesh)
  {
    cdmesh->rebuild_after_rebalance();
  }

  return rebalanced;
}

bool rebalance_mesh(stk::mesh::BulkData & bulk_data,
    CDMesh * cdmesh,
    const std::vector<std::string> & element_weights_field_names,
    const std::string & coordinates_field_name,
    const std::string & decomp_method,
    const double imbalance_threshold)
{
  const auto & meta = bulk_data.mesh_meta_data();

  std::vector<stk::mesh::Field<double> *> weights_fields;
  for (auto && field_name : element_weights_field_names)
  {
    auto weights_base = meta.get_field(stk::topology::ELEMENT_RANK, field_name);
    ThrowRequireMsg(weights_base,
        "Failed to find element rank field " << field_name << " to use for rebalance weights.");
    const auto element_weights_field = static_cast<stk::mesh::Field<double> *>(weights_base);

    update_parent_child_rebalance_weights(bulk_data, *element_weights_field, cdmesh);
    weights_fields.push_back(element_weights_field);
  }

  ThrowAssert(impl::check_family_tree_element_and_side_ownership(bulk_data));

  CDFEMRebalance balancer(
      bulk_data, cdmesh, coordinates_field_name, weights_fields, imbalance_threshold);
  balancer.setDecompMethod(decomp_method);
  const bool rebalanced = stk::balance::balanceStkMesh(balancer, bulk_data);
  if (rebalanced)
    stk::balance::balanceStkMeshNodes(balancer, bulk_data);

  ThrowAssert(impl::check_family_tree_element_and_side_ownership(bulk_data));

  if (cdmesh)
  {
    cdmesh->rebuild_after_rebalance();
  }

  return rebalanced;
}
}
}

