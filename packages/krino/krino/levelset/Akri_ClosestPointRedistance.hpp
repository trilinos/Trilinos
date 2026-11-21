#ifndef KRINO_KRINO_LEVELSET_AKRI_CLOSESTPOINTREDISTANCE_HPP_
#define KRINO_KRINO_LEVELSET_AKRI_CLOSESTPOINTREDISTANCE_HPP_
#include <stk_mesh/base/Types.hpp>
#include <stk_util/diag/Timer.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_FacetAppender.hpp>

namespace krino {

class MeshFacetedSurface;

class ClosestPointRedistance {
public:
ClosestPointRedistance(const stk::mesh::BulkData & mesh,
    const FieldRef& coordinates,
    const FieldRef& distance);
ClosestPointRedistance(const stk::mesh::BulkData & mesh,
    const FieldRef& coordinates,
    const FieldRef& distance,
    stk::diag::Timer & parentTimer);

static void build_facets_for_elements(const stk::mesh::BulkData & mesh,
    const FieldRef coordinates,
    const FieldRef distance,
    const std::vector<stk::mesh::Entity> & elementsToIntersect,
    const double elementLengthScale,
    FacetAppender & facets);
static void build_facets_for_selected_elements(const stk::mesh::BulkData & mesh,
    const FieldRef coordinates,
    const FieldRef distance,
    const stk::mesh::Selector & elementSelector,
    const double elementLengthScale,
    FacetAppender & facets);

// Methods for building facets and redistancing using them
void build_isosurface_facets(const stk::mesh::Selector & activeElementSelector,
    const double elementLengthScale,
    FacetAppender & facets) const;
void redistance_using_facets(const stk::mesh::Selector & activeElementSelector,
    MeshFacetedSurface & facets,
    const double narrowBandDistance = 0.,
    const bool doEnforceSignAwayFromInterface = false,
    const double elementLengthScale = 0.) const;
void redistance_given_nodes_using_facets(MeshFacetedSurface & facets,
    const std::vector<stk::mesh::Entity> & nodesToRedistance,
    const double narrowBandDistance = 0.,
    const bool doEnforceSignAwayFromInterface = false,
    const double elementLengthScale = 0.) const;

private:
static stk::mesh::Selector elements_to_contour_and_redistance(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & activeElementSelector, const FieldRef distance);
void redistance_selected_element_nodes_using_prepared_facets(
    const stk::mesh::Selector & activeElementSelector,
    const MeshFacetedSurface & facets,
    const double narrowBandDistance,
    const bool doEnforceSignAwayFromInterface,
    const double elementLengthScale) const;

const stk::mesh::BulkData & myMesh;
FieldRef myCoordinates;
FieldRef myDistance;
mutable stk::diag::Timer myTimer;

};

}



#endif /* KRINO_KRINO_LEVELSET_AKRI_CLOSESTPOINTREDISTANCE_HPP_ */
