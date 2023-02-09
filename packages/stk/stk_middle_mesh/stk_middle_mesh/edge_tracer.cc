#include "edge_tracer.h"
#include "field.h"
#include "mesh_entity.h"
#include "predicates/edge_intersection_primitive.h"
#include <limits>

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

using predicates::impl::PointClassification;

void EdgeTracer::trace_edge(mesh::MeshEntityPtr edge, const predicates::impl::PointRecord& recordStart,
                            const predicates::impl::PointRecord& recordEnd,
                            std::vector<EdgeIntersection>& intersections)
{
  assert(edge->get_type() == stk::middle_mesh::mesh::MeshEntityType::Edge);
  assert(is_entity_on_mesh(edge));

  auto& normalField        = *m_normalField;
  utils::Point normalStart = normalField(edge->get_down(0), 0, 0);
  utils::Point normalEnd   = normalField(edge->get_down(1), 0, 0);

  trace_edge(edge->get_down(0)->get_point_orig(0), normalStart, recordStart, edge->get_down(1)->get_point_orig(0),
             normalEnd, recordEnd, intersections);
}

void EdgeTracer::trace_edge(const utils::Point& ptStart, const utils::Point& normalStart,
                            const predicates::impl::PointRecord& recordStart, const utils::Point& ptEnd,
                            const utils::Point& normalEnd, const predicates::impl::PointRecord& recordEnd,
                            std::vector<EdgeIntersection>& intersections)
{
  assert(!is_entity_on_mesh(recordStart.el));
  assert(!is_entity_on_mesh(recordEnd.el));

  intersections.clear();
  bool atEndOfLine = false;
  // std::cout << "\ncomputing first intersection" << std::endl;

  // std::cout << "tracing edge from " << pt_start << " to " << pt_end << std::endl;
  // std::cout << "start_record = " << record_start << ", el = " << record_start.el << std::endl;
  // std::cout << "end_record = " << record_end << ", el = " << record_end.el << std::endl;
  EdgeIntersection currentIntersection = compute_next_intersection(
      predicates::impl::PointRecord(), recordStart, recordEnd, ptStart, ptEnd, normalStart, normalEnd, atEndOfLine);
  predicates::impl::PointRecord prevRecord    = recordStart;
  predicates::impl::PointRecord currentRecord = currentIntersection.record;
  while (!atEndOfLine)
  {
    intersections.push_back(currentIntersection);

    currentIntersection = compute_next_intersection(prevRecord, currentRecord, recordEnd, ptStart, ptEnd, normalStart,
                                                    normalEnd, atEndOfLine);

    prevRecord    = currentRecord;
    currentRecord = currentIntersection.record;
  }

  redistribute_points_near_endpoints(intersections, m_tolerances.endpointRedistributionDistance);
}

EdgeIntersection EdgeTracer::compute_next_intersection(const predicates::impl::PointRecord& prevIntersection,
                                                       const predicates::impl::PointRecord& currentIntersection,
                                                       const predicates::impl::PointRecord& endVert,
                                                       const utils::Point& ptStart, const utils::Point& ptEnd,
                                                       const utils::Point& normalStart, const utils::Point& normalEnd,
                                                       bool& atEndOfLine)
{
  assert(endVert.type != PointClassification::Exterior);
  assert(currentIntersection.type != PointClassification::Exterior);
  if (prevIntersection.type != PointClassification::Exterior)
  {
    assert(currentIntersection.type == PointClassification::Vert ||
           currentIntersection.type == PointClassification::Edge);
  }

  // std::cout << "current intersection = " << current_intersection << ", el = " << current_intersection.el <<
  // std::endl;

  mesh::MeshEntityPtr entity = get_entity(currentIntersection);
  // std::cout << "current intersection is on " << entity << std::endl;
  // TODO: in most cases there is only 1 excluded element, the current_intersection.el.  For
  //       performance, it might be worthwhile to special case it
  std::vector<mesh::MeshEntityPtr> includedElements =
      get_included_elements(entity, prevIntersection, currentIntersection);

  // std::cout << "all elements = " << std::endl;
  // for (auto& el : all_elements)
  //   std::cout << el << std::endl;
  //
  // std::cout << "excluded elements = " << std::endl;
  // for (auto& el : excluded_elements)
  //  std::cout << el << std::endl;

  if (at_end_of_line(includedElements, endVert))
  {
    atEndOfLine = true;
    return {endVert, 1, -1};
  }

  atEndOfLine = false;

  std::vector<mesh::MeshEntityPtr> includedEdges = get_included_edges(includedElements, currentIntersection);
  // std::cout << "initially, included_edges.size = " << included_edges.size() << std::endl;
  // std::cout << "\nincluded edges: " << std::endl;
  // for (auto& included_edge : included_edges)
  //   std::cout << included_edge << std::endl;

  /*
    // get included entities for other end of the edge
    // TODO: when traversing an edge, the included_edges_endpoint gets recomputed several times.  Consider caching it
    std::vector<mesh::MeshEntityPtr> included_elements_endpoint = getIncludedElements(get_entity(end_vert),
    predicates::impl::PointRecord(), end_vert); std::vector<mesh::MeshEntityPtr> included_edges_endpoint =
    getIncludedEdges(included_elements_endpoint, end_vert); std::vector<mesh::MeshEntityPtr> common_edges =
    get_commonEntities(included_edges, included_edges_endpoint);

    std::cout << "\nendpoint included edges: " << std::endl;
    for (auto& included_edge : included_edges_endpoint)
      std::cout << included_edge << std::endl;

    std::cout << "\ncommon edges: " << std::endl;
    for (auto& common_edge : common_edges)
      std::cout << common_edge << std::endl;

    std::cout << "common_edges.size = " << common_edges.size() << ", included_edges.size = " << included_edges.size() <<
    std::endl; std::vector<mesh::MeshEntityPtr>& edges_to_intersect = common_edges.size() > 0 ? common_edges :
    included_edges;
  */
  int edgeIdx;
  double alpha, beta;
  edgeIdx = compute_best_intersection(ptStart, ptEnd, normalStart, normalEnd, includedEdges, alpha, beta);
  // std::cout << "found intersection on edge " << edge_idx << ", alpha = " << alpha << ", beta = " << beta <<
  // std::endl;
  return create_edge_intersection(includedEdges[edgeIdx], includedElements, alpha, beta);
}

std::vector<mesh::MeshEntityPtr>
EdgeTracer::get_included_elements(mesh::MeshEntityPtr entity, const predicates::impl::PointRecord& prevIntersection,
                                  const predicates::impl::PointRecord& currentIntersection)
{
  std::vector<mesh::MeshEntityPtr> allElements      = get_elements(entity);
  std::vector<mesh::MeshEntityPtr> excludedElements = get_excluded_elements(prevIntersection, currentIntersection);
  return get_included_entities(allElements, excludedElements);
}

std::vector<mesh::MeshEntityPtr>
EdgeTracer::get_included_edges(const std::vector<mesh::MeshEntityPtr>& includedElements,
                               const predicates::impl::PointRecord& currentIntersection)
{
  std::vector<mesh::MeshEntityPtr> allEdges      = get_edges(includedElements);
  std::vector<mesh::MeshEntityPtr> excludedEdges = get_excluded_edges(currentIntersection);
  return get_included_entities(allEdges, excludedEdges);
}

std::vector<mesh::MeshEntityPtr> EdgeTracer::get_elements(mesh::MeshEntityPtr entity)
{
  if (entity->get_type() == stk::middle_mesh::mesh::MeshEntityType::Vertex)
  {
    std::vector<mesh::MeshEntityPtr> els;
    get_upward(entity, 2, els);
    return els;
  } else if (entity->get_type() == stk::middle_mesh::mesh::MeshEntityType::Edge)
  {
    std::vector<mesh::MeshEntityPtr> els(entity->count_up());
    for (int i = 0; i < entity->count_up(); ++i)
      els[i] = entity->get_up(i);

    return els;
  } else
    return {entity};
}

std::vector<mesh::MeshEntityPtr>
EdgeTracer::get_excluded_elements(const predicates::impl::PointRecord& prevIntersection,
                                  const predicates::impl::PointRecord& currentIntersection)
{
  assert(currentIntersection.type != PointClassification::Exterior);

  if (currentIntersection.type == PointClassification::Interior ||
      prevIntersection.type == PointClassification::Exterior)
  {
    assert(prevIntersection.type == PointClassification::Exterior);
    return {};
  }

  if (prevIntersection.type == PointClassification::Interior)
    return {prevIntersection.el};

  return get_common_elements(get_entity(prevIntersection), get_entity(currentIntersection));
}

std::vector<mesh::MeshEntityPtr> EdgeTracer::get_common_elements(mesh::MeshEntityPtr entity1,
                                                                 mesh::MeshEntityPtr entity2)
{
  assert(entity1->get_type() == stk::middle_mesh::mesh::MeshEntityType::Vertex ||
         entity1->get_type() == stk::middle_mesh::mesh::MeshEntityType::Edge);
  assert(entity2->get_type() == stk::middle_mesh::mesh::MeshEntityType::Vertex ||
         entity2->get_type() == stk::middle_mesh::mesh::MeshEntityType::Edge);

  std::vector<mesh::MeshEntityPtr> els1, els2, elsCommon;
  get_upward(entity1, 2, els1);
  get_upward(entity2, 2, els2);

  return get_common_entities(els1, els2);
  /*
  std::sort(els1.begin(), els1.end(), mesh::is_less);
  std::sort(els2.begin(), els2.end(), mesh::is_less);
  std::set_intersection(els1.begin(), els1.end(), els2.begin(), els2.end(), std::back_inserter(els_common),
  mesh::is_less);
  */
  return elsCommon;
}

std::vector<mesh::MeshEntityPtr> EdgeTracer::get_common_entities(std::vector<mesh::MeshEntityPtr>& entities1,
                                                                 std::vector<mesh::MeshEntityPtr>& entities2)
{
  std::vector<mesh::MeshEntityPtr> commonEntities;
  std::sort(entities1.begin(), entities1.end(), mesh::is_less);
  std::sort(entities2.begin(), entities2.end(), mesh::is_less);
  std::set_intersection(entities1.begin(), entities1.end(), entities2.begin(), entities2.end(),
                        std::back_inserter(commonEntities), mesh::is_less);

  return commonEntities;
}

std::vector<mesh::MeshEntityPtr>
EdgeTracer::get_included_entities(std::vector<mesh::MeshEntityPtr>& allEntities,
                                  const std::vector<mesh::MeshEntityPtr>& excludedEntities)
{
  assert(std::is_sorted(excludedEntities.begin(), excludedEntities.end(), mesh::is_less));
  std::sort(allEntities.begin(), allEntities.end(), mesh::is_less);

  std::vector<mesh::MeshEntityPtr> includedEntities;
  std::set_difference(allEntities.begin(), allEntities.end(), excludedEntities.begin(), excludedEntities.end(),
                      std::back_inserter(includedEntities), mesh::is_less);

  return includedEntities;
}

bool EdgeTracer::at_end_of_line(const std::vector<mesh::MeshEntityPtr>& includedElements,
                                const predicates::impl::PointRecord& endVert)
{
  mesh::MeshEntityPtr entity = get_entity(endVert);
  return contains_entity(includedElements, entity);
}

bool EdgeTracer::contains_entity(const std::vector<mesh::MeshEntityPtr>& elements, mesh::MeshEntityPtr entity)
{
  if (get_type_dimension(entity->get_type()) == 2)
  {
    return std::find(elements.begin(), elements.end(), entity) != elements.end();
  } else
  {
    int entityDim = get_type_dimension(entity->get_type());
    std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> downEntities;
    for (auto& el : elements)
    {
      int ndown = get_downward(el, entityDim, downEntities.data());
      for (int i = 0; i < ndown; ++i)
        if (downEntities[i] == entity)
          return true;
    }

    return false;
  }

  throw std::runtime_error("unreachable");
}

std::vector<mesh::MeshEntityPtr> EdgeTracer::get_edges(const std::vector<mesh::MeshEntityPtr>& elements)
{
  std::vector<mesh::MeshEntityPtr> edges;
  for (auto& el : elements)
    for (int i = 0; i < el->count_down(); ++i)
      edges.push_back(el->get_down(i));

  std::sort(edges.begin(), edges.end(), mesh::is_less);
  auto it = std::unique(edges.begin(), edges.end());
  edges.erase(it, edges.end());

  return edges;
}

std::vector<mesh::MeshEntityPtr> EdgeTracer::get_excluded_edges(const predicates::impl::PointRecord& record)
{
  assert(record.type != PointClassification::Exterior);

  if (record.type == PointClassification::Interior)
    return {};
  else if (record.type == PointClassification::Edge)
    return {get_entity(record)};
  else
  {
    mesh::MeshEntityPtr vert = get_entity(record);
    std::vector<mesh::MeshEntityPtr> edges(vert->count_up());
    for (int i = 0; i < vert->count_up(); ++i)
      edges[i] = vert->get_up(i);

    std::sort(edges.begin(), edges.end(), mesh::is_less);
    return edges;
  }
}

int EdgeTracer::compute_best_intersection(const utils::Point& startPt, const utils::Point& endPt,
                                          const utils::Point& startNormal, const utils::Point& endNormal,
                                          const std::vector<mesh::MeshEntityPtr>& edges, double& alpha, double& beta)
{
  // std::cout << "\nEntered computeBestIntersection" << std::endl;
  // std::cout << "input edge = " << start_pt << " to " << end_pt << std::endl;
  // std::cout << "normal vectors = " << start_normal << " and " << end_normal << std::endl;
  // std::cout << "second point on normal vectors = " << start_pt + start_normal << ", " << end_pt + end_normal <<
  // std::endl;
  assert(edges.size() > 0);
  // alphas are for the line (start_pt, end_pt), betas are for the edges
  std::vector<double> alphas, betas;
  std::vector<int> intersectionIdxToEdgeIdx;
  for (size_t i = 0; i < edges.size(); ++i)
  {
    // std::cout << "\ncomputing intersection with " << edges[i]->get_down(0)->get_point_orig(0) << " to " <<
    // edges[i]->get_down(1)->get_point_orig(0) << std::endl;
    predicates::impl::EdgeIntersectionResult result = predicates::impl::compute_edge_intersection(
        edges[i]->get_down(0)->get_point_orig(0), edges[i]->get_down(1)->get_point_orig(0), startPt, endPt, startNormal,
        endNormal, m_tolerances.edgeIntersectionPrimitiveZeroTol, m_tolerances.edgeIntersectionPrimitiveErrorTol);

    if (result.intersection1_found())
    {
      // std::cout << std::boolalpha << "found first intersection: " << result.intersection1_found() << std::endl;
      //  the edge intersection primitive defines alpha and beta backwards compared to
      //  what this class expects
      alphas.push_back(result.get_beta1());
      betas.push_back(result.get_alpha1());
      // std::cout << "alpha = " << alphas.back() << ", beta = " << betas.back() << std::endl;
      intersectionIdxToEdgeIdx.push_back(i);
    }

    if (result.intersection2_found())
    {
      // std::cout << std::boolalpha << "found second intersection: " << result.intersection2_found() << std::endl;
      //  the edge intersection primitive defines alpha and beta backwards compared to
      //  what this class expects
      alphas.push_back(result.get_beta2());
      betas.push_back(result.get_alpha2());
      // std::cout << "alpha = " << alphas.back() << ", beta = " << betas.back() << std::endl;
      intersectionIdxToEdgeIdx.push_back(i);
    }
  }

  // for (size_t i=0; i < alphas.size(); ++i)
  //   std::cout << "possible intersection " << i << ", alpha = " << alphas[i] << ", beta = " << betas[i] << std::endl;

  // std::cout << "alphas.size = " << alphas.size() << std::endl;
  assert(alphas.size() > 0);

  int intersection = choose_best_intersection(alphas, betas, edges, intersectionIdxToEdgeIdx, startPt, endPt);
  // std::cout << "chose intersection " << intersection << " on edge " << intersection_idx_to_edge_idx[intersection] <<
  // std::endl;
  alpha = alphas[intersection];
  beta  = betas[intersection];
  // std::cout << "chosen alpha, beta = " << alpha << ", " << beta << std::endl;
  return intersectionIdxToEdgeIdx[intersection];
}

int EdgeTracer::choose_best_intersection(const std::vector<double>& alphas, const std::vector<double>& betas,
                                         const std::vector<mesh::MeshEntityPtr>& edges,
                                         const std::vector<int>& intersectionIdxToEdgeIdx,
                                         const utils::Point& edgeStartPt, const utils::Point& edgeEndPt)
{
  assert(alphas.size() == betas.size());
  std::vector<int> intersectionsInRange;

  for (size_t i = 0; i < alphas.size(); ++i)
    if (predicates::impl::in_range(alphas[i], 0, 1, m_tolerances.inRangeTol) &&
        predicates::impl::in_range(betas[i], 0, 1, m_tolerances.inRangeTol))
      intersectionsInRange.push_back(i);

  if (intersectionsInRange.size() == 1)
  {
    return intersectionsInRange[0];
  } else if (intersectionsInRange.size() == 0)
  {
    for (size_t i = 0; i < alphas.size(); ++i)
      if (predicates::impl::in_range(alphas[i], 0, 1, m_tolerances.rangeNarrowingTol) &&
          predicates::impl::in_range(betas[i], 0, 1, m_tolerances.rangeNarrowingTol))
        intersectionsInRange.push_back(i);

    if (intersectionsInRange.size() == 0)
      throw std::runtime_error("could not find an edge intersection anywhere close to within range");

    return choose_intersection_with_minimum_angle(alphas, betas, intersectionsInRange, edges, intersectionIdxToEdgeIdx,
                                                  edgeStartPt, edgeEndPt);
  } else
  {
    return choose_intersection_with_minimum_angle(alphas, betas, intersectionsInRange, edges, intersectionIdxToEdgeIdx,
                                                  edgeStartPt, edgeEndPt);
  }
}

// returns a value describing how far outside the range [val_min, val_max] the
// given value is.  The value is always positive
double EdgeTracer::compute_deviation(double val, double valMin, double valMax)
{
  if (val >= valMin && val <= valMax)
    return 0;
  else
  {
    double devMin = std::abs(val - valMin);
    double devMax = std::abs(val - valMax);
    return std::min(devMin, devMax);
  }
}

// TODO: this takes a lot of arguments: refactor somehow?
int EdgeTracer::choose_intersection_with_minimum_angle(const std::vector<double>& alphas,
                                                       const std::vector<double>& betas,
                                                       const std::vector<int>& intersectionsInRange,
                                                       const std::vector<mesh::MeshEntityPtr>& edges,
                                                       const std::vector<int>& intersectionIdxToEdgeIdx,
                                                       const utils::Point& edgeStartPt, const utils::Point& edgeEndPt)
{
  // std::cout << "found multiple solutions in range" << std::endl;
  assert(intersectionsInRange.size() > 0);
  utils::Point edgeDirection = edgeEndPt - edgeStartPt;
  edgeDirection              = edgeDirection / std::sqrt(dot(edgeDirection, edgeDirection));

  double maxCosTheta  = std::numeric_limits<double>::min();
  int maxIntersection = -1;
  for (int idx : intersectionsInRange)
  {
    // compute vector from start_pt to intersection point (alpha[i], beta[i])
    mesh::MeshEntityPtr edge1 = edges[intersectionIdxToEdgeIdx[idx]];
    utils::Point ptStart      = edge1->get_down(0)->get_point_orig(0);
    utils::Point ptEnd        = edge1->get_down(1)->get_point_orig(0);

    double beta                     = betas[idx];
    utils::Point intersectionPt     = ptStart + beta * (ptEnd - ptStart); // TODO: numerically accuracy - subtraction?
    utils::Point intersectionVector = intersectionPt - edgeStartPt;
    intersectionVector              = intersectionVector / std::sqrt(dot(intersectionVector, intersectionVector));

    // compute angle between vector and edge_direction
    double cosTheta = dot(edgeDirection, intersectionVector);

    // choose minimum angle
    if (cosTheta > maxCosTheta)
    {
      maxCosTheta     = cosTheta;
      maxIntersection = idx;
    }
  }

  assert(maxIntersection != -1);

  return maxIntersection;
}

EdgeIntersection EdgeTracer::create_edge_intersection(mesh::MeshEntityPtr edge,
                                                      const std::vector<mesh::MeshEntityPtr>& includedElements,
                                                      double alpha, double beta)
{
  mesh::MeshEntityPtr el = get_parent_element(edge, includedElements);

  alpha                      = clamp(alpha, 0, 1);
  beta                       = clamp(beta, 0, 1);
  mesh::MeshEntityPtr entity = get_intersected_entity(edge, beta);

  assert(entity->get_type() == stk::middle_mesh::mesh::MeshEntityType::Vertex ||
         entity->get_type() == stk::middle_mesh::mesh::MeshEntityType::Edge);
  PointClassification type = entity->get_type() == stk::middle_mesh::mesh::MeshEntityType::Vertex ?
                                 PointClassification::Vert :
                                 PointClassification::Edge;
  beta                     = type == PointClassification::Edge ? beta : -1;

  int id = predicates::impl::get_entity_id(el, entity);

  return {predicates::impl::PointRecord(type, id, el), alpha, beta};
}

mesh::MeshEntityPtr EdgeTracer::get_intersected_entity(mesh::MeshEntityPtr edge, double beta)
{
  // figure out where beta is on the edge
  // TODO: make this branchless.
  mesh::MeshEntityPtr entity;
  if (std::abs(beta) < m_tolerances.vertexClassificationTol)
  {
    entity = edge->get_down(0);
  } else if (std::abs(beta - 1) < m_tolerances.vertexClassificationTol)
  {
    entity = edge->get_down(1);
  } else
  {
    entity = edge;
  }

  return entity;
}

double EdgeTracer::clamp(double val, double valMin, double valMax)
{
  double val1 = std::max(val, valMin);
  return std::min(val1, valMax);
}

mesh::MeshEntityPtr EdgeTracer::get_parent_element(mesh::MeshEntityPtr edge,
                                                   const std::vector<mesh::MeshEntityPtr>& includedElements)
{
#ifndef NDEBUG
  int nfound = 0;
  for (auto& el : includedElements)
  {
    for (int i = 0; i < el->count_down(); ++i)
      if (el->get_down(i) == edge)
        nfound++;
  }

  assert(nfound == 1);
#endif
  for (auto& el : includedElements)
  {
    for (int i = 0; i < el->count_down(); ++i)
      if (el->get_down(i) == edge)
        return el;
  }

  throw std::runtime_error("could not find parent element");
}

bool EdgeTracer::is_entity_on_mesh(mesh::MeshEntityPtr entity)
{
  int id          = entity->get_id();
  auto& entityVec = m_mesh->get_mesh_entities(get_type_dimension(entity->get_type()));

  return id < int(entityVec.size()) && entityVec[id] == entity;
}

// finds edge intersections within distance eps of the endpoints and uniformly distributes them, while
// preserving ordering.  This avoids the creation of zero length edges when, for example, an intersection
// is found with alpha < 0 or alpha > 1, which is then clamped into the range [0, 1]
void EdgeTracer::redistribute_points_near_endpoints(std::vector<EdgeIntersection>& intersections, double eps)
{
  redistribute_points_near_endpoint(intersections, true, eps);
  redistribute_points_near_endpoint(intersections, false, eps);
}

void EdgeTracer::redistribute_points_near_endpoint(std::vector<EdgeIntersection>& intersections, bool isLeftEndpoint,
                                                   double eps)
{
  double alphaStart = isLeftEndpoint ? 0 : 1 - eps;
  double alphaEnd   = isLeftEndpoint ? eps : 1;
  std::vector<int> indicesToRedistribute;
  for (size_t i = 0; i < intersections.size(); ++i)
    if (std::abs(intersections[i].alpha - alphaStart) < eps)
      indicesToRedistribute.push_back(i);

  if (indicesToRedistribute.size() == 0)
    return;

  // we don't want there to be an intersection on the endpoint, so pretend
  // there is one more point than there actually is, and don't assign
  // the endpoint
  int npts          = indicesToRedistribute.size();
  double deltaAlpha = (alphaEnd - alphaStart) / (npts + 1);
  for (int i = 0; i < npts; ++i)
  {
    double alphaNew                               = alphaStart + deltaAlpha * (i + 1);
    intersections[indicesToRedistribute[i]].alpha = alphaNew;
  }
}

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
