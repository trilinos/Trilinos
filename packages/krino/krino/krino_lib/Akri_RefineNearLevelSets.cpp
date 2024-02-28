#include <Akri_RefineNearLevelSets.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <Akri_AdaptivityHelpers.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_CDFEM_Support.hpp>
#include <Akri_CDMesh_Refinement.hpp>
#include <Akri_Faceted_Surface.hpp>
#include <Akri_LevelSetSurfaceInterfaceGeometry.hpp>
#include <Akri_Phase_Support.hpp>
#include <Akri_RefinementInterface.hpp>
#include <Akri_RefinementSupport.hpp>
#include <stk_util/diag/Timer.hpp>

namespace krino {

std::vector<LS_Field> build_levelset_fields(const std::vector<stk::mesh::Field<double>*> & stkLevelSetFields)
{
  unsigned lsId = 0;
  std::vector<LS_Field> lsFields;
  for (auto stkLevelSetField : stkLevelSetFields)
    lsFields.emplace_back(stkLevelSetField->name(),Surface_Identifier(lsId++),FieldRef(stkLevelSetField),0.,nullptr);
  return lsFields;
}

void refine_elements_that_intersect_distance_interval_from_levelsets(stk::mesh::BulkData & mesh,
    const stk::mesh::Part & activePart,
    const std::vector<stk::mesh::Field<double>*> & stkLevelSetFields,
    const std::function<void()> & initialize_levelsets,
    const std::array<double,2> & refinementDistanceInterval,
    const unsigned numRefinementLevels)
{
  RefinementSupport & refinementSupport = RefinementSupport::get(mesh.mesh_meta_data());
  refinementSupport.activate_nonconformal_adaptivity(numRefinementLevels);
  refinementSupport.set_refinement_interval(refinementDistanceInterval);
  RefinementInterface & refinement = krino::KrinoRefinement::create(mesh.mesh_meta_data());
  refinementSupport.set_non_interface_conforming_refinement(refinement);

  const CDFEM_Support & cdfemSupport = CDFEM_Support::get(mesh.mesh_meta_data());
  const Phase_Support & phaseSupport = Phase_Support::get(mesh.mesh_meta_data());
  const auto lsFields = build_levelset_fields(stkLevelSetFields);
  const LevelSetSurfaceInterfaceGeometry interfaceGeom(mesh.mesh_meta_data().spatial_dimension(), activePart, cdfemSupport, phaseSupport, lsFields);

  const int numRefinementSteps = 2*refinementSupport.get_interface_maximum_refinement_level(); // Make sure refinement completes so that elements touching interval are fully refined

  std::function<void(int)> initialize_levelsets_and_compute_distance_and_mark_elements_that_intersect_interval =
    [&mesh, &interfaceGeom, &refinementSupport, &initialize_levelsets, numRefinementSteps](int refinementIterCount)
    {
      initialize_levelsets();
      if (refinementIterCount < numRefinementSteps)
      {
        krino::mark_elements_that_intersect_interval(mesh,
            refinementSupport.get_non_interface_conforming_refinement(),
            interfaceGeom,
            refinementSupport,
            refinementIterCount);
      }
    };

  perform_multilevel_adaptivity(refinementSupport.get_non_interface_conforming_refinement(),
      mesh,
      initialize_levelsets_and_compute_distance_and_mark_elements_that_intersect_interval,
      refinementSupport.get_do_not_refine_or_unrefine_selector());
}

}
