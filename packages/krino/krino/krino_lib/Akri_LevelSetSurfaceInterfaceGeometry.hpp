#ifndef KRINO_KRINO_KRINO_LIB_AKRI_LEVELSETSURFACEINTERFACEGEOMETRY_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_LEVELSETSURFACEINTERFACEGEOMETRY_HPP_
#include <stk_mesh/base/Part.hpp>
#include <Akri_AnalyticSurfaceInterfaceGeometry.hpp>
#include <Akri_Surface.hpp>
#include <Akri_Faceted_Surface.hpp>
#include "Akri_Phase_Support.hpp"

namespace krino {

class CDFEM_Support;
class Phase_Support;

class LevelSetSurfaceInterfaceGeometry : public AnalyticSurfaceInterfaceGeometry {

public:
  LevelSetSurfaceInterfaceGeometry(const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport,
    const std::vector<LS_Field> & LSFields);

  virtual void prepare_to_process_elements(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const override;
  virtual void prepare_to_process_elements(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & elementsToIntersect,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const override;

private:
  void build_levelset_facets_if_needed(const stk::mesh::BulkData & mesh) const;
  void build_levelset_facets(const stk::mesh::BulkData & mesh) const;

  std::vector<LS_Field> myLSFields;
  std::vector<Surface_Identifier> mySurfaceIdentifiers;
  mutable std::vector<Faceted_Surface> myLSSurfaces;
  mutable size_t myLastMeshSyncCount{0};
};

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_LEVELSETSURFACEINTERFACEGEOMETRY_HPP_ */
