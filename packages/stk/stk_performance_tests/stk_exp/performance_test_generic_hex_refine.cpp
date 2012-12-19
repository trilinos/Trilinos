#ifndef __IBMCPP__
#include <boost/unordered_map.hpp>

#include <sierra/mesh/modifiable/modifiable_mesh.hpp>
#include <sierra/mesh/details/cell_topology.hpp>
#include <sierra/mesh/fixture/hex_fixture.hpp>
#include <sierra/mesh/mesh_traits.hpp>

#include <stk_performance_test_includes/hex_refine_info.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/environment/CPUTime.hpp>

#include <iostream>

#include <boost/range.hpp>

using namespace sierra::mesh;
using namespace sierra::mesh::details;

namespace stk {
namespace performance_tests {

namespace {

void create_entities( modifiable_mesh & mesh,
                      details::part_key& node_part,
                      details::part_key& hex_part,
                      // FIXME Part& active_elements_part,
                      HexRefineInfo& refine_info)

{

  HexRefineInfo refine_info_half(refine_info.m_level-1, refine_info.m_nx, refine_info.m_ny, refine_info.m_nz, 0);

  unsigned eid_start = 1 + refine_info.elem_id_offset();
  unsigned eid_end = eid_start + refine_info.num_elems();

  unsigned nid_start = 1 + refine_info.node_id_offset();
  unsigned nid_end = nid_start + refine_info.num_nodes();

  std::cout << "nid_start = " << nid_start << " nid_end= " << nid_end << " diff= " << nid_end - nid_start << std::endl;
  std::cout << "eid_start = " << eid_start << " eid_end= " << eid_end << " diff= " << eid_end - eid_start << std::endl;

  entity_rank node_rank(0);
  entity_rank elem_rank(3);

  // FIXME - need modifiable_mesh::add_n_entities();
  for(unsigned nid=nid_start; nid<nid_end; ++nid) {
    entity_key key = mesh.add_entity(entity_property(node_rank, entity_id(nid)));
    (void) key; // avoid an unused but set variable warning from gcc 4.6.3
  }

  for(unsigned eid=eid_start; eid<eid_end; ++eid) {
    entity_key elem_key = mesh.add_entity(entity_property(elem_rank, entity_id(eid)));
    mesh.change_entity_parts(elem_key, &hex_part, &hex_part+1);

    unsigned ix = 0, iy = 0, iz = 0;
    refine_info.elem_x_y_z(eid, ix, iy, iz);
    unsigned ie_check = refine_info.elem_id(ix, iy, iz);
    EXPECT_EQ(ie_check, eid);

    unsigned elem_node[8] ;

    elem_node[0] = refine_info.node_id( ix   , iy   , iz   );
    elem_node[1] = refine_info.node_id( ix+1 , iy   , iz   );
    elem_node[2] = refine_info.node_id( ix+1 , iy   , iz+1 );
    elem_node[3] = refine_info.node_id( ix   , iy   , iz+1 );
    elem_node[4] = refine_info.node_id( ix   , iy+1 , iz   );
    elem_node[5] = refine_info.node_id( ix+1 , iy+1 , iz   );
    elem_node[6] = refine_info.node_id( ix+1 , iy+1 , iz+1 );
    elem_node[7] = refine_info.node_id( ix   , iy+1 , iz+1 );

    // check if a parent node
    if (1)
    {
      for (unsigned i = 0; i<8; ++i) {
        unsigned ixn = 0, iyn = 0, izn = 0;
        refine_info.node_x_y_z(elem_node[i], ixn, iyn, izn);
        unsigned in_check = refine_info.node_id(ixn, iyn, izn);
        EXPECT_EQ(in_check, elem_node[i]);

        if (
          ((ixn - 1) % 2 == 0) &&
          ((iyn - 1) % 2 == 0) &&
          ((izn - 1) % 2 == 0))
        {
          elem_node[i] = refine_info_half.node_id(ixn/2, iyn/2, izn/2);
        }
      }
    }


    for(size_t j=0; j<8u; ++j) {

      // FIXME - add_relation expects an entity_key for the "to" arg, so how do we get a key from an id?
      entity_key node = static_cast<entity_key>(elem_node[j]);
      mesh.add_relation(elem_key, node, relation_position(0,j));
      //mesh.change_entity_parts(node_map[elem_node[j]], &node_part, &node_part+1);
      //mesh.change_entity_parts(node_map[elem_node[j]], &hex_part, &hex_part+1);
    }

    // FIXME: what's the equivalent here?  how do we check if a node is in the database given its id?
    // in the above code, I just used a local map, but in general we will need this lookup
#if 0

    for (unsigned i = 0; i<8; ++i) {
      stk::mesh::Entity const node = bulk.get_entity( fem::FEMMetaData::NODE_RANK , elem_node[i] );
      bulk.change_entity_parts(*node, add_parts);

      ThrowRequireMsg( node != NULL, "found null node in create_entities");

      // Compute and assign coordinates to the node
      unsigned nx = 0, ny = 0, nz = 0;
      refine_info.node_x_y_z(elem_node[i], nx, ny, nz);

      Scalar * data = stk::mesh::field_data( m_coord_field , *node );

      data[0] = nx ;
      data[1] = ny ;
      data[2] = -(Scalar)nz ;
    }
#endif
  }
}

} // empty namespace

STKUNIT_UNIT_TEST(modifiable_mesh, hex_refine_generic)
{
  double start_time = stk::cpu_time();
  //unsigned ex=100, ey=100, ez=100;
  const int max_levels = 2;
  unsigned nn = 50/(1 << (max_levels-1));   // use 4 here for a quick running callgrind case
  unsigned ex=nn, ey=nn, ez=nn;

  //unsigned num_elems = ex*ey*ez;
  sierra::mesh::fixture::hex_fixture fixture(ex,ey,ez);
  fixture.generate_mesh();
  double end_time = stk::cpu_time() - start_time;

  std::cout << "Time to create hex mesh: " << end_time << std::endl;

  //const double tolerance = 1.e-6;
  //const double expected = ((double)ex)/2;

  start_time = stk::cpu_time();

  for (int level=1; level <= max_levels; ++level) {

    HexRefineInfo refine_info(level, ex, ey, ez, 0);

    //Selector hex_elem_selector(fixture.m_hex_part & !fixture.m_node_part);

    unsigned num_elems_new = refine_info.num_elems();
    std::cout << "num_elems_new for level = " << level << " = " << num_elems_new << std::endl;
    //modification_begin(fixture.m_mesh);
    create_entities(fixture.m_mesh, fixture.m_node_part, fixture.m_hex_part, refine_info);
    //modification_end(fixture.m_mesh);
  }

  end_time = stk::cpu_time();

  double test_time = end_time - start_time;
  std::cout << "Num refines: " << max_levels << std::endl;
  std::cout << "Time to refine: " << test_time << std::endl;
}

}
}

#endif
