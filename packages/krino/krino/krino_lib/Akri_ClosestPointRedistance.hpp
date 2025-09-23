#ifndef KRINO_KRINO_KRINO_LIB_AKRI_CLOSESTPOINTREDISTANCE_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_CLOSESTPOINTREDISTANCE_HPP_
#include <stk_mesh/base/Types.hpp>
#include <stk_util/diag/Timer.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_Faceted_Surface.hpp>

namespace krino {

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
    FacetedSurfaceBase & facets);
static void build_facets_for_selected_elements(const stk::mesh::BulkData & mesh,
    const FieldRef coordinates,
    const FieldRef distance,
    const stk::mesh::Selector & elementSelector,
    const double elementLengthScale,
    FacetedSurfaceBase & facets);

// Complete method that contours the elements and then uses resulting facets to redistance
void redistance(const stk::mesh::Selector & activeElementSelector,
    const double elementLengthScale,
    const double narrowBandDistance=0.,
    const bool doEnforceSignAwayFromInterface=false) const;

// Separate methods for building facets and redistancing
void build_isosurface_facets(const stk::mesh::Selector & activeElementSelector, const double elementLengthScale, FacetedSurfaceBase & facets) const;
void redistance_using_facets(const stk::mesh::Selector & activeElementSelector,
    FacetedSurfaceBase & facets,
    const double narrowBandDistance=0.,
    const bool doEnforceSignAwayFromInterface=false,
    const double elementLengthScale=0.) const;
void redistance_given_nodes_using_facets(FacetedSurfaceBase & facets,
    const std::vector<stk::mesh::Entity> & nodesToRedistance,
    const double narrowBandDistance=0.,
    const bool doEnforceSignAwayFromInterface=false,
    const double elementLengthScale=0.) const;

private:
static stk::mesh::Selector elements_to_contour_and_redistance(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & activeElementSelector, const FieldRef distance);
void redistance_selected_element_nodes_using_prepared_facets(const stk::mesh::Selector & activeElementSelector,
    const FacetedSurfaceBase & facets,
    const double narrowBandDistance,
    const bool doEnforceSignAwayFromInterface,
    const double elementLengthScale) const;

const stk::mesh::BulkData & myMesh;
FieldRef myCoordinates;
FieldRef myDistance;
mutable stk::diag::Timer myTimer;

};

}



#endif /* KRINO_KRINO_KRINO_LIB_AKRI_CLOSESTPOINTREDISTANCE_HPP_ */
