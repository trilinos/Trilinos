#ifndef __IBMCPP__

/*--------------------------------------------------------------------*/
/*    Copyright 2009, 2011 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#include <sstream>
#include <gtest/gtest.h>

#include <stk_util/environment/CPUTime.hpp>
#include <sierra/io/modifiable_mesh_reader.hpp>
#include <sierra/io/io_fixture.hpp>
#include <sierra/mesh/csr_mesh_api.hpp>
#include <sierra/mesh/fixture/csr_mesh_factory.hpp>
#include <sierra/mesh/details/cell_topology.hpp>
#include <sierra/mesh/details/selected_buckets.hpp>

#include <Geom_Vec3d.h>

extern int* STKUNIT_ARGC;
extern char** STKUNIT_ARGV;

namespace Fperf_generic
{

void get_mesh_file_name(int argc, char** argv, const std::string& flag, std::string& mesh_file_name)
{
  mesh_file_name.clear();

  for(int i=0; i<argc; ++i) {
    if (argv[i] == NULL) continue;
    if (flag != argv[i]) continue;
    if (i < argc-1) {
      mesh_file_name = argv[i+1];
    }
  }
}

sierra::mesh::csr_mesh::part_key get_rank_part(const sierra::mesh::details::entity_rank& rank, const sierra::mesh::csr_mesh& csr)
{
  using namespace sierra::mesh;
  using namespace sierra::mesh::details;

  csr_mesh::part_range parts = csr.get_parts();

  for(csr_mesh::part_iterator part_iter=boost::begin(parts); part_iter!=boost::end(parts); ++part_iter) {
    csr_mesh::part_key part = *part_iter;
    const csr_mesh::part_property& part_prop = csr[part];

    if (!part_prop.has_property<sierra::mesh::entity_rank>()) continue;
    //we want a part with an entity-rank property but not a cell_topology property:
    if (part_prop.has_property<cell_topology>()) continue;

    if (part_prop.get_property<sierra::mesh::entity_rank>() == rank) {
      return part;
    }
  }

  csr_mesh::part_key invalid_part;
  return invalid_part;
}

void get_tet_parts(sierra::mesh::csr_mesh& csr,
                   std::vector<sierra::mesh::csr_mesh::part_key>& tet_parts)
{
  using namespace sierra::mesh;
  using namespace sierra::mesh::details;

  tet_parts.clear();

  csr_mesh::part_range parts = csr.get_parts();

  const cell_topology tet_topology = shards::getCellTopologyData<shards::Tetrahedron<4> >();
  for(csr_mesh::part_iterator part_iter=boost::begin(parts); part_iter!=boost::end(parts); ++part_iter) {
    csr_mesh::part_key part = *part_iter;
    const csr_mesh::part_property& part_prop = csr[part];
    std::cout << "Part: " << part_prop.name() << std::endl;

    if (!part_prop.has_property<cell_topology>()) continue;

    const cell_topology& cell_top = part_prop.get_property<cell_topology>();
    if (cell_top == tet_topology) {
      tet_parts.push_back(part);
      std::cout << "     tet part!" << std::endl;
    }
  }
}

template < class Mesh, class Field, class Selector >
void ExtractVertexGradientOperator(
    const sierra::mesh::entity_descriptor & node,
    const Field & gradop_field,
    const Selector & tet_selector,
    std::vector<geometry::Vec3d> & gradient_vector,
    const Mesh & mesh
    )
{
//This test mimics mesh-traversal that is done by the
//function ExtractVertexGradientOperator in strumento/src/element/ug3dt4n.C.
//This test only does the stencilSize==1 case, which is the case that is exercised by
//the heavy-test a021_snow_grosh_wheel.
  using namespace sierra::mesh;
  using namespace geometry;

  gradient_vector.clear();
  std::vector<entity_descriptor> node_descriptors;
  node_descriptors.reserve(32);

  const unsigned nodes_per_tet = 4;
  entity_rank node_rank(0);//\TODO fix this hardwired magic number!
  entity_rank element_rank(3);//\TODO fix this hardwired magic number!

  //get element-relations for node:
  relation_range elem_range = get_relations(node, element_rank, mesh);
  for(relation_iterator elem_rel=boost::begin(elem_range), elem_end=boost::end(elem_range); elem_rel!=elem_end; ++elem_rel) {
    entity_descriptor elem = target_entity(*elem_rel,mesh);

    //skip elements that are not tet elements:
    if (!is_selected(get_bucket(elem,mesh),tet_selector,mesh)) { continue; }

    const double curScale = 0.25;
    double * gradop = get_field_data(gradop_field, elem, mesh);

    //get node-relations for elem:
    relation_range node_range = get_relations(elem, node_rank, mesh);
    relation_iterator node_rel=boost::const_begin(node_range);

    for(unsigned inode=0; inode < nodes_per_tet; ++inode) {
      entity_descriptor node_local = target_entity(*node_rel,mesh);
      ++node_rel;
      bool need_to_push_back = true;

      for(unsigned n=0, sz=node_descriptors.size(); n<sz; ++n) {
        if (node_descriptors[n] == node_local) {
          gradient_vector[n] += Vec3d(curScale*gradop[inode*3+0], curScale*gradop[inode*3+1], curScale*gradop[inode*3+2]);
          need_to_push_back = false;
          break;
        }
      }

      if (need_to_push_back) {
        node_descriptors.push_back(node_local);
        gradient_vector.push_back(Vec3d(curScale*gradop[inode*3+0], curScale*gradop[inode*3+1], curScale*gradop[inode*3+2]));
      }
    }
  }
}

//This test will visit all nodes that are connected to tet elements,
//and for each node call the ExtractVertexGradientOperator function,
//and use the resulting gradient_vector to update a nodal field.
template<class ElemField, class NodeField>
void do_nodal_based_tet_test(sierra::mesh::details::selector& tet_selector,
                             const ElemField& gradop_field,
                             const NodeField& nodal_aspect_field,
                             const sierra::mesh::csr_mesh& csr)
{
  using namespace sierra::mesh;
  using namespace geometry;

  std::vector<Vec3d> gradient_vector;

  selector tet_node_select = tet_selector & get_rank_part(csr.node_rank(), csr);
  selected_bucket_range buckets = sierra::mesh::get_buckets(tet_node_select, csr);
  for(selected_bucket_iterator b_iter=boost::begin(buckets); b_iter!=boost::end(buckets); ++b_iter) {
    const bucket_key b = *b_iter;

    double* node_aspect = sierra::mesh::get_field_data(nodal_aspect_field, b, csr);

    csr_mesh::bucket_entity_range nodes = sierra::mesh::get_entities(b, csr);
    unsigned i = 0;
    for(csr_mesh::bucket_entity_iterator n_iter=boost::begin(nodes); n_iter!=boost::end(nodes); ++n_iter, ++i) {
      entity_descriptor node = *n_iter;
      ExtractVertexGradientOperator(node, gradop_field, tet_selector, gradient_vector, csr);

      double grad_inner_product = 0;
      for(std::vector<Vec3d>::const_iterator gv=gradient_vector.begin(), gv_end=gradient_vector.end(); gv!=gv_end; ++gv) {
        const Vec3d& grad_val = *gv;
        grad_inner_product += Dot(4.0*grad_val, 4.0*grad_val);
      }

      node_aspect[i] = grad_inner_product;
    }
  }
}

//Here's this test-program's "main":
TEST(CSR, Nodal_Based_Tet_meshAPI_generic)
{
  using namespace sierra::mesh;
  using details::constant_size_field;

  double start_time = stk::cpu_time();

  std::string mesh_file_name;
  get_mesh_file_name(*STKUNIT_ARGC, STKUNIT_ARGV, "-mesh", mesh_file_name);
  if (mesh_file_name.empty()) {
    std::cout << "no mesh file specified." << std::endl;
    return;
  }

  std::cout << "Attempting to read mesh file: " << mesh_file_name << std::endl;

  sierra::mesh::io::io_fixture fixture(MPI_COMM_WORLD, "exodus", mesh_file_name);
  boost::shared_ptr<csr_mesh> csr = csr_mesh_factory::create_from_modifiable(fixture.mmesh);

  std::vector<csr_mesh::part_key> tet_parts;
  get_tet_parts(*csr, tet_parts);
  std::cout << "num tet parts: " << tet_parts.size()<<std::endl;

  selector tet_selector;
  for(size_t i=0; i<tet_parts.size(); ++i) tet_selector |= tet_parts[i];

  const unsigned nodes_per_tet = 4;
  const unsigned spatial_dim = 3;

  selector tet_elem_selector = tet_selector & get_rank_part(csr->element_rank(), *csr);
  constant_size_field<double,nodes_per_tet*spatial_dim> gradop_field;
  gradop_field.update_from_mesh(tet_elem_selector, *csr);

  selector tet_node_selector = tet_selector & get_rank_part(csr->node_rank(), *csr);
  constant_size_field<double,1> nodal_aspect_field;
  nodal_aspect_field.update_from_mesh(tet_node_selector, *csr);

  double mesh_create_time = stk::cpu_time() - start_time;


  //initialize the contents of the gradient_operator field to an arbitrary value.
  std::vector<double>& gradop_contents = gradop_field.field_data_flat_array();
  std::fill(gradop_contents.begin(), gradop_contents.end(), 0.5);

  start_time = stk::cpu_time();

  const int num_iters = 100;
  for(int t=0; t<num_iters; ++t) {

    do_nodal_based_tet_test(tet_selector, gradop_field, nodal_aspect_field, *csr);

  }

  double test_time = stk::cpu_time() - start_time;

  std::cout << "Time to read/create mesh: " << mesh_create_time << std::endl;
  std::cout << "Time to do Nodal-Based-Tet test: " << test_time << std::endl;

  std::cout << "PASSED" << std::endl;
}

} //namespace Fperf_generic

#endif
