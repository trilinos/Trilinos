#include <Akri_ClosestPointRedistance.hpp>
#include <Akri_LevelSetInterfaceGeometry.hpp>
#include <Akri_LevelSetSurfaceInterfaceGeometry.hpp>
#include <Akri_MeshHelpers.hpp>

#include <stk_mesh/base/GetEntities.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include "Akri_AnalyticSurfaceInterfaceGeometry.hpp"
#include "Akri_CDFEM_Parent_Edges.hpp"
#include "Akri_ContourElement.hpp"
#include "Akri_FieldRef.hpp"
#include "Akri_LevelSet.hpp"

namespace krino {

LevelSetSurfaceInterfaceGeometry::LevelSetSurfaceInterfaceGeometry(const int dim,
    const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport,
    const std::vector<LS_Field> & LSFields)
: AnalyticSurfaceInterfaceGeometry(activePart, cdfemSupport, phaseSupport),
  myLSFields(LSFields)
{
  for (auto && lsField : myLSFields)
  {
    mySurfaceIdentifiers.push_back(lsField.identifier);
    myLSSurfaces.emplace_back(FacetedSurfaceBase::build(dim));
  }

  const std::vector<stk::mesh::Selector> levelsetElementSelectors = get_levelset_element_selectors();
  for (size_t i=0; i<mySurfaceIdentifiers.size(); ++i)
    add_surface(mySurfaceIdentifiers[i], *myLSSurfaces[i], levelsetElementSelectors[i]);
}

std::vector<stk::mesh::Selector> LevelSetSurfaceInterfaceGeometry::get_levelset_element_selectors() const
{
  const bool isCdfemUseCase = get_phase_support().is_cdfem_use_case();

  std::vector<stk::mesh::Selector> levelsetElementSelector;
  for (auto && lsField : myLSFields)
  {
    if (isCdfemUseCase)
      levelsetElementSelector.push_back(get_phase_support().get_levelset_decomposed_blocks_selector(lsField.identifier));
    else
      levelsetElementSelector.emplace_back(stk::mesh::selectField(lsField.isovar));
  }

  return levelsetElementSelector;
}

static bool is_node_snapped_to_levelset(const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    const stk::mesh::Entity node,
    const int lsIndex)
{
  const auto iter = nodesToCapturedDomains.find(node);
  if (iter == nodesToCapturedDomains.end())
    return false;
  const std::vector<int> & nodeSortedDomains = iter->second;
  return std::binary_search(nodeSortedDomains.begin(), nodeSortedDomains.end(), lsIndex);
}

static int determine_node_sign_from_nodal_levelset(const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    const stk::mesh::Entity node,
    const std::vector<LS_Field> & lsFields,
    const int lsIndex)
{
  // In general maybe it makes sense to use ls_node_value here
  if (is_node_snapped_to_levelset(nodesToCapturedDomains, node, lsIndex))
    return 0;
  const double * ls = field_data<double>(lsFields[lsIndex].isovar, node);
  if (ls == nullptr)
    return 0;
  if (*ls < lsFields[lsIndex].isoval)
    return -1;
  return 1;
}

static void fill_node_signs_from_nodal_levelsets(const stk::mesh::BulkData & mesh,
  const stk::mesh::Entity node,
  const std::vector<stk::mesh::Selector> & perSurfaceElementSelector,
  const NodeToCapturedDomainsMap & nodesToCapturedDomains,
  const std::vector<LS_Field> & lsFields,
  std::vector<int8_t> & nodeSigns)
{
  const size_t numLS = perSurfaceElementSelector.size();
  nodeSigns.assign(numLS, -2);

  for (size_t lsIndex=0; lsIndex<numLS; ++lsIndex)
    if (is_entity_selected(mesh, perSurfaceElementSelector[lsIndex], node))
      nodeSigns[lsIndex] = determine_node_sign_from_nodal_levelset(nodesToCapturedDomains, node, lsFields, lsIndex);
}

void LevelSetSurfaceInterfaceGeometry::set_node_signs_from_nodal_levelsets(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const
{
  NodeToSignsMap nodesToSigns;

  const stk::mesh::Selector elementSelector = get_mesh_parent_element_selector();
  const std::vector<stk::mesh::Selector> levelsetElementSelectors = get_levelset_element_selectors();

  for(const auto & bucketPtr : mesh.get_buckets(stk::topology::NODE_RANK, elementSelector))
    for (const auto & node : *bucketPtr)
      fill_node_signs_from_nodal_levelsets(mesh, node, levelsetElementSelectors, nodesToCapturedDomains, myLSFields, nodesToSigns[node]);

  set_node_signs(nodesToSigns);
}

void LevelSetSurfaceInterfaceGeometry::prepare_to_decompose_elements(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const
{
  build_levelset_facets_if_needed(mesh);

  set_elements_to_intersect_and_prepare_to_compute_with_surfaces(mesh, get_mesh_parent_elements(mesh));

  const bool doDetermineNodeSignFromNodalLevelSets = false;
  if (doDetermineNodeSignFromNodalLevelSets)
    set_node_signs_from_nodal_levelsets(mesh, nodesToCapturedDomains);
  else
    set_node_signs_from_surfaces(mesh, get_levelset_element_selectors(), nodesToCapturedDomains);
}


void LevelSetSurfaceInterfaceGeometry::prepare_to_intersect_elements(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const
{
  build_levelset_facets_if_needed(mesh);
  AnalyticSurfaceInterfaceGeometry::prepare_to_intersect_elements(mesh, nodesToCapturedDomains);
}

void LevelSetSurfaceInterfaceGeometry::prepare_to_intersect_elements(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & elementsToIntersect,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const
{
  build_levelset_facets_if_needed(mesh);
  AnalyticSurfaceInterfaceGeometry::prepare_to_intersect_elements(mesh, elementsToIntersect, nodesToCapturedDomains);
}

void LevelSetSurfaceInterfaceGeometry::build_levelset_facets_if_needed(const stk::mesh::BulkData & mesh) const
{
  if (myDoUpdateFacetsWhenMeshChanges && mesh.synchronized_count() != myLastMeshSyncCount)
  {
    build_levelset_facets(mesh);
    myLastMeshSyncCount = mesh.synchronized_count();
  }
}

void LevelSetSurfaceInterfaceGeometry::build_levelset_facets(const stk::mesh::BulkData & mesh) const
{
  if (myLSFields.empty())
    return;

  const stk::mesh::Selector elementSelector = (get_phase_support().is_cdfem_use_case()) ?
      stk::mesh::Selector(get_active_part() & get_phase_support().get_all_decomposed_blocks_selector()) :
      stk::mesh::Selector(get_active_part());
  const stk::mesh::Selector ownedElementSelector = elementSelector & mesh.mesh_meta_data().locally_owned_part();

  std::vector<stk::mesh::Entity> elementsToIntersect;
  const FieldRef coordsField = get_coordinates_field(mesh);
  stk::mesh::get_selected_entities( ownedElementSelector, mesh.buckets(stk::topology::ELEMENT_RANK), elementsToIntersect, false );
  const double avgEdgeLength = compute_global_average_edge_length_for_elements(mesh, coordsField, elementsToIntersect);

  const std::vector<stk::mesh::Selector> levelsetElementSelectors = get_levelset_element_selectors();
  for (size_t i=0; i<myLSFields.size(); ++i)
  {
    const stk::mesh::Selector elemFieldSelector = ownedElementSelector & levelsetElementSelectors[i];
    stk::mesh::get_selected_entities( elemFieldSelector, mesh.buckets(stk::topology::ELEMENT_RANK), elementsToIntersect, false );

    if (myLSFields[i].ptr && myLSFields[i].ptr->is_transient())
    {
      if (2 == mesh.mesh_meta_data().spatial_dimension())
        myLSSurfaces[i]->get_facets_2d() = myLSFields[i].ptr->get_facets().get_facets_2d();
      else
        myLSSurfaces[i]->get_facets_3d() = myLSFields[i].ptr->get_facets().get_facets_3d();
    }
    else
    {
      ClosestPointRedistance::build_facets_for_elements(mesh, coordsField, myLSFields[i].isovar, elementsToIntersect, avgEdgeLength, *myLSSurfaces[i]);
    }
  }
}

std::vector<stk::mesh::Entity> LevelSetSurfaceInterfaceGeometry::get_possibly_cut_elements(const stk::mesh::BulkData & mesh) const
{
  // NOTE: Uses levelset field directly, not facetted interface, because it needs the analog, not discrete version
  return LevelSetInterfaceGeometry::get_active_elements_that_may_be_cut_by_levelsets(mesh, get_active_part(), myLSFields);
}


}
