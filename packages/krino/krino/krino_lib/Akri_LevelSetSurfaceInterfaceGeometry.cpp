#include <Akri_LevelSetSurfaceInterfaceGeometry.hpp>

#include <stk_mesh/base/GetEntities.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include "Akri_AnalyticSurfaceInterfaceGeometry.hpp"
#include "Akri_CDFEM_Parent_Edges.hpp"
#include "Akri_ContourElement.hpp"
#include "Akri_FieldRef.hpp"
#include "Akri_LevelSet.hpp"

namespace krino {

LevelSetSurfaceInterfaceGeometry::LevelSetSurfaceInterfaceGeometry(const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport,
    const std::vector<LS_Field> & LSFields)
: AnalyticSurfaceInterfaceGeometry(activePart, cdfemSupport, phaseSupport),
  myLSFields(LSFields)
{
  for (auto && lsField : myLSFields)
  {
    mySurfaceIdentifiers.push_back(lsField.identifier);
    myLSSurfaces.emplace_back("My Facets");
  }

  for (size_t i=0; i<mySurfaceIdentifiers.size(); ++i)
    add_surface(mySurfaceIdentifiers[i], myLSSurfaces[i]);
}

void LevelSetSurfaceInterfaceGeometry::prepare_to_process_elements(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const
{
  build_levelset_facets_if_needed(mesh);
  AnalyticSurfaceInterfaceGeometry::prepare_to_process_elements(mesh, nodesToCapturedDomains);
}

void LevelSetSurfaceInterfaceGeometry::prepare_to_process_elements(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & elementsToIntersect,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const
{
  build_levelset_facets_if_needed(mesh);
  AnalyticSurfaceInterfaceGeometry::prepare_to_process_elements(mesh, elementsToIntersect, nodesToCapturedDomains);
}

void LevelSetSurfaceInterfaceGeometry::build_levelset_facets_if_needed(const stk::mesh::BulkData & mesh) const
{
  if (mesh.synchronized_count() != myLastMeshSyncCount)
  {
    build_levelset_facets(mesh);
    myLastMeshSyncCount = mesh.synchronized_count();
  }
}

void LevelSetSurfaceInterfaceGeometry::build_levelset_facets(const stk::mesh::BulkData & mesh) const
{
  stk::mesh::Selector conformingElementSelector = get_active_part() & get_phase_support().get_all_decomposed_blocks_selector();

  std::vector<stk::mesh::Entity> elementsToIntersect;
  stk::mesh::get_selected_entities( conformingElementSelector, mesh.buckets(stk::topology::ELEMENT_RANK), elementsToIntersect, false );

  const FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());
  const FieldRef firstIsoField = myLSFields[0].isovar;

  const double avgEdgeLength = LevelSet::compute_global_average_edge_length_for_elements(mesh, coordsField, firstIsoField, elementsToIntersect);
  for (size_t i=0; i<myLSFields.size(); ++i)
    LevelSet::build_facets_for_elements(mesh, coordsField, myLSFields[i].isovar, elementsToIntersect, avgEdgeLength, myLSSurfaces[i]);
}

}
