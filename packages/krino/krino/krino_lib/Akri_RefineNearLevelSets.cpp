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

std::vector<LS_Field> build_levelset_fields(const stk::mesh::Field<double> & stkLevelSetField)
{
  unsigned lsId = 0;
  std::vector<LS_Field> lsFields;
  lsFields.emplace_back(stkLevelSetField.name(),Surface_Identifier(lsId++),FieldRef(stkLevelSetField),0.,nullptr);
  return lsFields;
}

void refine_elements_that_intersect_distance_interval_from_levelset(stk::mesh::BulkData & mesh,
    const stk::mesh::Part & activePart,
    const stk::mesh::Field<double> & stkLevelSetField,
    const std::function<void()> & initialize_levelset,
    const std::array<double,2> & refinementDistanceInterval,
    const unsigned numRefinementLevels)
{
  RefinementInterface & refinement = krino::KrinoRefinement::get_or_create(mesh.mesh_meta_data());

  const CDFEM_Support & cdfemSupport = CDFEM_Support::get(mesh.mesh_meta_data());
  const Phase_Support & phaseSupport = Phase_Support::get(mesh.mesh_meta_data());
  const auto lsFields = build_levelset_fields(stkLevelSetField);
  const LevelSetSurfaceInterfaceGeometry interfaceGeom(mesh.mesh_meta_data().spatial_dimension(), activePart, cdfemSupport, phaseSupport, lsFields);

  const int numRefinementSteps = 2*numRefinementLevels; // Make sure refinement completes so that elements touching interval are fully refined
  constexpr bool isDefaultCoarsen = false; // Refinement here will be "cumulative"

  std::function<void(int)> initialize_levelset_and_compute_distance_and_mark_elements_that_intersect_interval =
    [&mesh, &initialize_levelset, &refinement, &interfaceGeom, &refinementDistanceInterval, numRefinementSteps, numRefinementLevels](int refinementIterCount)
    {
      initialize_levelset();
      if (refinementIterCount < numRefinementSteps)
      {
        krino::mark_elements_that_intersect_interval(mesh,
            refinement,
            interfaceGeom,
            refinementDistanceInterval,
            numRefinementLevels,
            isDefaultCoarsen);
      }
    };

  stk::mesh::Selector emptySelector;
  perform_multilevel_adaptivity(refinement,
      mesh,
      initialize_levelset_and_compute_distance_and_mark_elements_that_intersect_interval,
      emptySelector);
}

}
