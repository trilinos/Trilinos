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

  for (size_t i=0; i<mySurfaceIdentifiers.size(); ++i)
    add_surface(mySurfaceIdentifiers[i], *myLSSurfaces[i]);
}

std::vector<stk::mesh::Selector> LevelSetSurfaceInterfaceGeometry::get_levelset_element_selectors() const
{
  const bool isCdfemUseCase = is_cdfem_use_case(get_phase_support());

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

void LevelSetSurfaceInterfaceGeometry::prepare_to_decompose_elements(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const
{
  build_levelset_facets_if_needed(mesh);

  set_elements_to_intersect_and_prepare_to_compute_with_surfaces(mesh, get_mesh_parent_elements(mesh));

  const std::vector<stk::mesh::Selector> levelsetElementSelector = get_levelset_element_selectors();
  set_element_signs(mesh, levelsetElementSelector);
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

  const stk::mesh::Selector elementSelector = (is_cdfem_use_case(get_phase_support())) ?
      stk::mesh::Selector(get_active_part() & get_phase_support().get_all_decomposed_blocks_selector()) :
      stk::mesh::Selector(get_active_part());

  std::vector<stk::mesh::Entity> elementsToIntersect;
  const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());
  stk::mesh::get_selected_entities( elementSelector, mesh.buckets(stk::topology::ELEMENT_RANK), elementsToIntersect, false );
  const double avgEdgeLength = compute_global_average_edge_length_for_elements(mesh, coordsField, elementsToIntersect);

  const std::vector<stk::mesh::Selector> levelsetElementSelectors = get_levelset_element_selectors();
  for (size_t i=0; i<myLSFields.size(); ++i)
  {
    const stk::mesh::Selector elemFieldSelector = elementSelector & levelsetElementSelectors[i];
    stk::mesh::get_selected_entities( elemFieldSelector, mesh.buckets(stk::topology::ELEMENT_RANK), elementsToIntersect, false );

    LevelSet::build_facets_for_elements(mesh, coordsField, myLSFields[i].isovar, elementsToIntersect, avgEdgeLength, *myLSSurfaces[i]);
  }
}

std::vector<stk::mesh::Entity> LevelSetSurfaceInterfaceGeometry::get_possibly_cut_elements(const stk::mesh::BulkData & mesh) const
{
  // NOTE: Uses levelset field directly, not facetted interface, because it needs the analog, not discrete version
  return LevelSetInterfaceGeometry::get_active_elements_that_may_be_cut_by_levelsets(mesh, get_active_part(), myLSFields);
}


}
