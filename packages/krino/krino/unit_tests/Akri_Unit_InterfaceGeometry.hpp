/*
 * Akri_Unit_InterfaceGeometry.hpp
 *
 *  Created on: May 4, 2023
 *      Author: drnoble
 */

#ifndef KRINO_KRINO_UNIT_TESTS_AKRI_UNIT_INTERFACEGEOMETRY_HPP_
#define KRINO_KRINO_UNIT_TESTS_AKRI_UNIT_INTERFACEGEOMETRY_HPP_
#include <Akri_InterfaceGeometry.hpp>
#include <Akri_PhaseTag.hpp>

namespace krino {

class IntersectionPointFromNodalLevelsetInterfaceGeometry : public InterfaceGeometry
{
public:
  virtual ~IntersectionPointFromNodalLevelsetInterfaceGeometry() {}
  virtual void prepare_to_decompose_elements(const stk::mesh::BulkData & mesh, const NodeToCapturedDomainsMap & nodesToCapturedDomains) const override {}
  virtual void prepare_to_intersect_elements(const stk::mesh::BulkData & mesh) const override {}
  virtual void prepare_to_intersect_elements(const stk::mesh::BulkData & mesh, const NodeToCapturedDomainsMap & nodesToCapturedDomains) const override {}
  virtual void prepare_to_intersect_elements(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & elementsToIntersect,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains) const override {}

  virtual std::vector<stk::mesh::Entity> get_possibly_cut_elements(const stk::mesh::BulkData & mesh) const override
    { STK_ThrowRequireMsg(false, "Unimplemented"); std::vector<stk::mesh::Entity> empty; return empty; }
  virtual void fill_elements_that_intersect_distance_interval(const stk::mesh::BulkData & mesh, const Surface_Identifier surfaceIdentifier, const std::array<double,2> loAndHi, std::vector<stk::mesh::Entity> & elementsThaIntersectInterval) const override
    { STK_ThrowRequireMsg(false, "Unimplemented"); }

  virtual bool might_have_interior_or_face_intersections() const override { STK_ThrowRequireMsg(false, "Unimplemented"); return false; }
  virtual bool snapped_elements_may_have_new_intersections() const override { return false; }

  virtual void append_element_intersection_points(const stk::mesh::BulkData & mesh,
    const NodeToCapturedDomainsMap & nodesToCapturedDomains,
    const std::vector<stk::mesh::Entity> & elementsToIntersect,
    const IntersectionPointFilter & intersectionPointFilter,
    std::vector<IntersectionPoint> & intersectionPoints) const override;

  virtual const std::vector<Surface_Identifier> & get_surface_identifiers() const override
    { STK_ThrowRequireMsg(false, "Unimplemented"); static std::vector<Surface_Identifier> empty; return empty; }

  virtual void store_phase_for_uncut_elements(const stk::mesh::BulkData & mesh) const override
    { STK_ThrowRequireMsg(false, "Unimplemented"); }
  virtual void store_phase_for_elements_that_will_be_uncut_after_snapping(const stk::mesh::BulkData & mesh,
      const std::vector<IntersectionPoint> & intersectionPoints,
      const std::vector<SnapInfo> & independentSnapInfos,
      const NodeToCapturedDomainsMap & nodesToCapturedDomains) const override
    { STK_ThrowRequireMsg(false, "Unimplemented"); }

  virtual const ElementToDomainMap & get_phase_for_uncut_elements() const override
    { STK_ThrowRequireMsg(false, "Unimplemented"); static ElementToDomainMap empty; return empty; }

  virtual std::unique_ptr<ElementCutter> build_element_cutter(const stk::mesh::BulkData & mesh,
    stk::mesh::Entity element,
    const std::function<bool(const std::array<unsigned,4> &)> & intersectingPlanesDiagonalPicker) const override
    { STK_ThrowRequireMsg(false, "Unimplemented"); return nullptr; }

  virtual PhaseTag get_starting_phase(const ElementCutter * cutter) const override
    { STK_ThrowRequireMsg(false, "Unimplemented"); PhaseTag empty; return empty; }

  virtual std::vector<IntersectionPoint> get_edge_intersection_points(const stk::mesh::BulkData & mesh,
      const NodeToCapturedDomainsMap & nodesToCapturedDomains) const override
    { STK_ThrowRequireMsg(false, "Unimplemented"); static std::vector<IntersectionPoint> empty; return empty; }

  void set_nodal_levelset(const stk::mesh::BulkData & mesh, const std::vector<stk::mesh::EntityId> & nodeIds, const std::vector<double> & nodeLs);

private:
  std::map<stk::mesh::Entity, double> nodeLSValues;
};

}



#endif /* KRINO_KRINO_UNIT_TESTS_AKRI_UNIT_INTERFACEGEOMETRY_HPP_ */
