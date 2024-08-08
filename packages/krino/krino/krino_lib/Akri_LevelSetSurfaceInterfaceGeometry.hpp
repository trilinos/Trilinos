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
  LevelSetSurfaceInterfaceGeometry(const int dim,
    const stk::mesh::Part & activePart,
    const CDFEM_Support & cdfemSupport,
    const Phase_Support & phaseSupport,
    const std::vector<LS_Field> & LSFields);

  using AnalyticSurfaceInterfaceGeometry::prepare_to_intersect_elements;

  virtual bool might_have_interior_or_face_intersections() const override { return mySurfaceIdentifiers.size() > 1; }
  virtual void prepare_to_decompose_elements(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const override;
  virtual void prepare_to_intersect_elements(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const override;
  virtual void prepare_to_intersect_elements(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & elementsToIntersect,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const override;
  virtual void set_do_update_geometry_when_mesh_changes(const bool flag) const override { myDoUpdateFacetsWhenMeshChanges = flag; }

  // Methods that use levelset directly and not facetted levelset surface, because they need the analog, not discrete version
  virtual std::vector<stk::mesh::Entity> get_possibly_cut_elements(const stk::mesh::BulkData & mesh) const override;

private:
  std::vector<stk::mesh::Selector> get_levelset_element_selectors() const;
  void build_levelset_facets_if_needed(const stk::mesh::BulkData & mesh) const;
  void build_levelset_facets(const stk::mesh::BulkData & mesh) const;

  std::vector<LS_Field> myLSFields;
  std::vector<Surface_Identifier> mySurfaceIdentifiers;
  mutable std::vector<std::unique_ptr<FacetedSurfaceBase>> myLSSurfaces;
  mutable size_t myLastMeshSyncCount{0};
  mutable bool myDoUpdateFacetsWhenMeshChanges{true};
};

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_LEVELSETSURFACEINTERFACEGEOMETRY_HPP_ */
