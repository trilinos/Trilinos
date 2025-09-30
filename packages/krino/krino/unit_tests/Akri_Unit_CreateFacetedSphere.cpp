#include <Akri_Unit_CreateFacetedSphere.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <Akri_OutputUtils.hpp>

namespace krino {

static void scale_vertex_onto_sphere(const double radius, stk::math::Vector3d & vertex)
{
  vertex *= radius / vertex.length();
}

static stk::math::Vector3d midpoint_on_sphere(const double radius,
    const stk::math::Vector3d & vertexA,
    const stk::math::Vector3d & vertexB)
{
  stk::math::Vector3d mid = 0.5*(vertexA+vertexB);
  scale_vertex_onto_sphere(radius, mid);
  return mid;
}

static void fill_icosahedron(const double radius,
  std::vector<stk::math::Vector3d> & vertices,
  std::vector<std::array<unsigned,3>> & facetConnectivity)
{
  const double phi = (1 + std::sqrt(5)) / 2; // Golden ratio
  vertices =
  {
      {-1,  phi, 0},
      { 1,  phi, 0},
      {-1, -phi, 0},
      { 1, -phi, 0},
      {0, -1,  phi},
      {0,  1,  phi},
      {0, -1, -phi},
      {0,  1, -phi},
      { phi, 0, -1},
      { phi, 0,  1},
      {-phi, 0, -1},
      {-phi, 0,  1},
  };

  for (auto& vertex : vertices)
    scale_vertex_onto_sphere(radius, vertex);

  facetConnectivity =
  {
    {{0, 11, 5}},
    {{0, 5, 1}},
    {{0, 1, 7}},
    {{0, 7, 10}},
    {{0, 10, 11}},
    {{1, 5, 9}},
    {{5, 11, 4}},
    {{11, 10, 2}},
    {{10, 7, 6}},
    {{7, 1, 8}},
    {{3, 9, 4}},
    {{3, 4, 2}},
    {{3, 2, 6}},
    {{3, 6, 8}},
    {{3, 8, 9}},
    {{4, 9, 5}},
    {{2, 4, 11}},
    {{6, 2, 10}},
    {{8, 6, 7}},
    {{9, 8, 1}},
  };
}

static unsigned add_or_retrieve_child_vertex(const double radius,
  const unsigned indexA,
  const unsigned indexB,
  std::vector<stk::math::Vector3d>& vertices,
  std::map<std::pair<unsigned,unsigned>,unsigned>& childMap)
{
  const std::pair<unsigned,unsigned> key = (indexA<indexB) ?
    std::make_pair(indexA,indexB) :
    std::make_pair(indexB,indexA);
  auto iter = childMap.find(key);

  if (iter != childMap.end())
    return iter->second;

  const unsigned childIndex = vertices.size();
  vertices.push_back(midpoint_on_sphere(radius, vertices[indexA], vertices[indexB]));
  childMap.emplace_hint(iter, key, childIndex);
  return childIndex;
}

static void subdivide_facets_on_sphere(const double radius,
  std::vector<stk::math::Vector3d>& vertices,
  std::vector<std::array<unsigned,3>>& facetConnectivity)
{
  std::map<std::pair<unsigned,unsigned>,unsigned> childMap;
  std::vector<std::array<unsigned,3>> newFacetConn;
  newFacetConn.reserve(4*facetConnectivity.size());

  for (const auto& conn : facetConnectivity)
  {
    const std::array<unsigned,6> refinedFacet =
    {
      conn[0], conn[1], conn[2],
      add_or_retrieve_child_vertex(radius, conn[0], conn[1], vertices, childMap),
      add_or_retrieve_child_vertex(radius, conn[1], conn[2], vertices, childMap),
      add_or_retrieve_child_vertex(radius, conn[2], conn[0], vertices, childMap)
    };

    newFacetConn.emplace_back(std::array<unsigned,3>{refinedFacet[0],refinedFacet[3],refinedFacet[5]});
    newFacetConn.emplace_back(std::array<unsigned,3>{refinedFacet[1],refinedFacet[4],refinedFacet[3]});
    newFacetConn.emplace_back(std::array<unsigned,3>{refinedFacet[2],refinedFacet[5],refinedFacet[4]});
    newFacetConn.emplace_back(std::array<unsigned,3>{refinedFacet[3],refinedFacet[4],refinedFacet[5]});
  }

  facetConnectivity = newFacetConn;
}

void fill_sphere_vertices_and_connectivity(const double radius,
  const double meshSize,
  std::vector<stk::math::Vector3d> & vertices,
  std::vector<std::array<unsigned,3>> & facetConnectivity)
{
  fill_icosahedron(radius, vertices, facetConnectivity);

  double facetSize = (vertices[1]-vertices[0]).length();
  const unsigned maxRefineLevels = 9;
  STK_ThrowRequireMsg(meshSize > facetSize*std::pow(0.5,maxRefineLevels), "Too small of mesh size requested.");

  while(meshSize < 0.75*facetSize)
  {
    subdivide_facets_on_sphere(radius, vertices, facetConnectivity);
    facetSize *= 0.5;
  }
}


void write_stl_for_sphere(const std::string &filename,
  const double radius,
  const double meshSize,
  const stk::ParallelMachine comm)
{
  std::vector<stk::math::Vector3d> vertices;
  std::vector<std::array<unsigned,3>> facetConnectivity;

  if (0 == stk::parallel_machine_rank(comm))
  {
    fill_sphere_vertices_and_connectivity(radius, meshSize, vertices, facetConnectivity);
    write_stl(filename, vertices, facetConnectivity);
  }
}

}
