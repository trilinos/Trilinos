#ifndef __IBMCPP__
#include <sierra/mesh/fixture/csr_mesh_factory.hpp>
#include <sierra/mesh/fixture/hex_fixture.hpp>
#include <sierra/mesh/csr/csr_mesh.hpp>
#include <sierra/mesh/fixture/csr_mesh_factory.hpp>

#include <sierra/mesh/mesh_traits.hpp>
#include <sierra/mesh/csr_mesh_api.hpp>

#include <stk_performance_test_includes/generic_gather.hpp>
#include <stk_performance_test_includes/bucket_traverse.hpp>

#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <iostream>

#include <boost/range.hpp>

using namespace sierra::mesh;

namespace stk {
namespace performance_tests {

STKUNIT_UNIT_TEST( csr_mesh, generic_hex_gather )
{
  double start_time = stk::cpu_time();
  unsigned ex=100, ey=100, ez=100;
  unsigned num_elems = ex*ey*ez;
  sierra::mesh::fixture::hex_fixture fixture(ex,ey,ez);
  fixture.generate_mesh();
  double end_time = stk::cpu_time() - start_time;

  std::cout << "hex_fixture: "<<std::endl;
  std::cout << "\tNum Nodes: " << fixture.m_num_nodes<<std::endl;
  std::cout << "\tNum Elements: " << fixture.m_num_elements << std::endl;
  std::cout << "Time to create hex mesh: " << end_time << std::endl;

  start_time = stk::cpu_time();
  boost::shared_ptr<csr_mesh> mesh = csr_mesh_factory::create_from_modifiable(fixture.m_mesh);
  end_time = stk::cpu_time() - start_time;

  std::cout << "Time to convert mesh: " << end_time << std::endl;

  std::vector<double> avg_centroid(3,0);
  const double tolerance = 1.e-6;
  const double expected = ((double)ex)/2;

  start_time = stk::cpu_time();

  const int num_iters = 100;
  for(int t=0; t<num_iters; ++t) {
    mesh_traits<csr_mesh>::selector select(fixture.m_hex_part&!fixture.m_node_part);
    gather_centroid_algorithm(select, fixture.m_coordinates, *mesh, avg_centroid);

    for(size_t d=0; d<3u; ++d) {
      EXPECT_LT(std::abs(avg_centroid[d]/num_elems - expected), tolerance);
    }

    avg_centroid[0] = 0;
    avg_centroid[1] = 0;
    avg_centroid[2] = 0;
  }

  end_time = stk::cpu_time();

  double test_time = end_time - start_time;
  std::cout << "Num centroid iterations: " << num_iters << std::endl;
  std::cout << "Time to compute centroids: " << test_time << std::endl;
}

template <class Mesh, class Field, int spatial_dim>
class CentroidVisitor1: public BucketVisitor<Mesh>
{
public:
  typedef typename BucketVisitor<Mesh>::entity_descriptor       entity_descriptor;
  typedef typename BucketVisitor<Mesh>::relation_range          relation_range;
  typedef typename BucketVisitor<Mesh>::relation_descriptor     relation_descriptor;

  CentroidVisitor1(const Field &coordinatefield)
    : m_coordinateField(coordinatefield),
      m_numElements(0),
      m_numNodes(0)
  {}

  void initialize(Mesh & /*mesh */) {
    m_numElements = 0;

    for (int i = 0; i < spatial_dim; ++i)
      m_averageCentroid[i] = 0.0;
  }

  void discover_entity(entity_descriptor entity, Mesh & /* mesh */) {
    ++m_numElements;
    m_numNodes = 0;
    for (int i = 0; i < spatial_dim; ++i)
      m_coordinate[i] = 0.0;
  }

  void examine_relation(relation_descriptor rd, Mesh &mesh) {
    entity_descriptor node = target_entity(rd, mesh);
    double * node_coords = sierra::mesh::get_field_data(m_coordinateField, node, mesh);

    ++m_numNodes;
    for (int i = 0; i < spatial_dim; ++i)
      m_coordinate[i] += node_coords[i];
  }

  void finish_entity(entity_descriptor /* entity */, Mesh & /* mesh */) {
    for (int i = 0; i < spatial_dim; ++i)
      m_averageCentroid[i] += m_coordinate[i]/m_numNodes;
  }

private:
  const Field &         m_coordinateField;
  size_t                m_numElements;
  size_t                m_numNodes;
  double                m_coordinate[spatial_dim];

public:
  double                m_averageCentroid[spatial_dim];
};


STKUNIT_UNIT_TEST(csr_mesh, generic_hex_gather_bucket_traverse )
{
  double start_time = stk::cpu_time();
  unsigned ex=100, ey=100, ez=100;
  unsigned num_elems = ex*ey*ez;
  sierra::mesh::fixture::hex_fixture fixture(ex,ey,ez);
  fixture.generate_mesh();
  double end_time = stk::cpu_time() - start_time;

  std::cout << "hex_fixture: "<<std::endl;
  std::cout << "\tNum Nodes: " << fixture.m_num_nodes<<std::endl;
  std::cout << "\tNum Elements: " << fixture.m_num_elements << std::endl;
  std::cout << "Time to create hex mesh: " << end_time << std::endl;

  start_time = stk::cpu_time();
  boost::shared_ptr<csr_mesh> mesh = csr_mesh_factory::create_from_modifiable(fixture.m_mesh);
  end_time = stk::cpu_time() - start_time;

  std::cout << "Time to convert mesh: " << end_time << std::endl;

  const double tolerance = 1.e-6;
  const double expected = ((double)ex)/2;

  start_time = stk::cpu_time();

  CentroidVisitor1<csr_mesh, sierra::mesh::fixture::hex_fixture::CoordinateField, 3> centroid_visitor(fixture.m_coordinates);

  const int num_iters = 100;
  for(int t=0; t<num_iters; ++t) {
    mesh_traits<csr_mesh>::selector select(fixture.m_hex_part&!fixture.m_node_part);

    const mesh_traits<csr_mesh>::entity_rank node_rank(0);

    bucket_traverse(select, node_rank, centroid_visitor, *mesh);

    for(size_t d=0; d<3u; ++d) {
      EXPECT_LT(std::abs(centroid_visitor.m_averageCentroid[d]/num_elems - expected), tolerance);
    }
  }

  end_time = stk::cpu_time();

  double test_time = end_time - start_time;
  std::cout << "Num centroid iterations: " << num_iters << std::endl;
  std::cout << "Time to compute centroids: " << test_time << std::endl;
}


template <class Mesh, class Field, int spatial_dim>
struct CentroidVisitor2: public BucketVisitor<Mesh>
{
  typedef typename BucketVisitor<Mesh>::entity_descriptor       entity_descriptor;
  typedef typename BucketVisitor<Mesh>::relation_range          relation_range;
  typedef typename BucketVisitor<Mesh>::relation_descriptor     relation_descriptor;

  CentroidVisitor2(const Field &coordinatefield)
    : m_coordinateField(coordinatefield),
      m_numElements(0),
      m_numNodes(0),
      m_elementCentroid(spatial_dim)
  {}

  void initialize(Mesh & /*mesh */) {
    m_numElements = 0;

    for (int i = 0; i < spatial_dim; ++i)
      m_averageCentroid[i] = 0.0;
  }

  void discover_entity(entity_descriptor entity, Mesh & /* mesh */) {
    ++m_numElements;
  }

  void examine_relation_range(const relation_range &elem_node_relation_range, Mesh &mesh) {
    m_numNodes = boost::size(elem_node_relation_range);
    m_offset = 0;

    m_elementNodeCoordinate.resize(spatial_dim*m_numNodes);
  }

  void examine_relation(relation_descriptor rd, Mesh &mesh) {
    entity_descriptor node = target_entity(rd, mesh);

    const double * node_coords = sierra::mesh::get_field_data(m_coordinateField, node, mesh);

    for (int i = 0; i < spatial_dim; ++i)
      m_elementNodeCoordinate[m_offset++] = node_coords[i];
  }

  void finish_entity(entity_descriptor /* entity */, Mesh & /* mesh */) {
    stk::performance_tests::calculate_centroid_3d(m_numNodes, &m_elementNodeCoordinate[0], &m_elementCentroid[0]);

    for (int i = 0; i < spatial_dim; ++i) {
      m_averageCentroid[i] += m_elementCentroid[i];
      m_elementCentroid[i] = 0.0;
    }
  }

private:
  const Field &         m_coordinateField;
  size_t                m_numElements;
  size_t                m_numNodes;
  size_t                m_offset;
  std::vector<double>   m_elementNodeCoordinate;
  std::vector<double>   m_elementCentroid;


public:
  double                m_averageCentroid[spatial_dim];
};


STKUNIT_UNIT_TEST(csr_mesh, generic_hex_gather_bucket_traverse_2)
{
  double start_time = stk::cpu_time();
  unsigned ex=100, ey=100, ez=100;
  unsigned num_elems = ex*ey*ez;
  sierra::mesh::fixture::hex_fixture fixture(ex,ey,ez);
  fixture.generate_mesh();
  double end_time = stk::cpu_time() - start_time;

  std::cout << "hex_fixture: "<<std::endl;
  std::cout << "\tNum Nodes: " << fixture.m_num_nodes<<std::endl;
  std::cout << "\tNum Elements: " << fixture.m_num_elements << std::endl;
  std::cout << "Time to create hex mesh: " << end_time << std::endl;

  start_time = stk::cpu_time();
  boost::shared_ptr<csr_mesh> mesh = csr_mesh_factory::create_from_modifiable(fixture.m_mesh);
  end_time = stk::cpu_time() - start_time;

  std::cout << "Time to convert mesh: " << end_time << std::endl;

  const double tolerance = 1.e-6;
  const double expected = ((double)ex)/2;

  start_time = stk::cpu_time();

  CentroidVisitor2<csr_mesh, sierra::mesh::fixture::hex_fixture::CoordinateField, 3> centroid_visitor(fixture.m_coordinates);

  const int num_iters = 100;
  for(int t=0; t<num_iters; ++t) {
    mesh_traits<csr_mesh>::selector select(fixture.m_hex_part&!fixture.m_node_part);

    const mesh_traits<csr_mesh>::entity_rank node_rank(0);

    bucket_traverse(select, node_rank, centroid_visitor, *mesh);

    for(size_t d=0; d<3u; ++d) {
      EXPECT_LT(std::abs(centroid_visitor.m_averageCentroid[d]/num_elems - expected), tolerance);
    }
  }

  end_time = stk::cpu_time();

  double test_time = end_time - start_time;
  std::cout << "Num centroid iterations: " << num_iters << std::endl;
  std::cout << "Time to compute centroids: " << test_time << std::endl;
}

}
}
#endif
