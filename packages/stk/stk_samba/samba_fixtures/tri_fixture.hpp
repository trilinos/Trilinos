#ifndef SAMBA_FIXTURES_TRI_FIXTURE_HPP
#define SAMBA_FIXTURES_TRI_FIXTURE_HPP

#include <samba/mesh.hpp>
#include <samba/field.hpp>
#include <samba/rank_field.hpp>

namespace samba { namespace fixtures {


#include <iostream>
#include <iterator>
#include <algorithm>
#include <sstream>

#include <boost/range/algorithm/copy.hpp>
#include <boost/range/algorithm/equal.hpp>


class tri_fixture {

  static const uint32_t NX_default = 10;
  static const uint32_t NY_default = 10;
  //static const uint32_t NZ = 0;

  uint32_t NX;
  uint32_t NY;

  inline size_t node_offset(uint32_t x, uint32_t y)
  { return (x + ( NX + 1 ) * ( y  ));  }

  inline size_t quad_offset(uint32_t x, uint32_t y)
  { return (x + NX * ( y ));  }


public:

  tri_fixture(uint32_t nx = NX_default, uint32_t ny = NY_default) : NX(nx), NY(ny) {}

  void samba_mesh_create(samba::mesh mesh)
  {
    const uint32_t num_nodes = (NX+1)*(NY+1); //*(NZ+1);
    const uint32_t num_quads = NX*NY; //*NZ;
    const uint32_t num_tris = 2*num_quads;

    mesh.begin_modification();

    samba::entity_key_interval nodes = mesh.add_entities( samba::entity_topology::node(),  num_nodes);
    samba::entity_key_interval tris = mesh.add_entities( samba::entity_topology::triangle_3(), num_tris);

    for (uint32_t x=0; x<NX; ++x) {
      for (uint32_t y=0; y<NY; ++y) {


        samba::entity_key tri0 = tris[2*quad_offset(x,y)+0];
        samba::entity_key tri1 = tris[2*quad_offset(x,y)+1];
        samba::entity_key node;

        // ---
        samba::connectivity_ordinal ordinal0 = {0};
        node = nodes[node_offset( x   , y     )];
        mesh.add_connectivity(tri0,node,ordinal0);
        mesh.add_connectivity(node,tri0,ordinal0);
        ++ordinal0;

        node = nodes[node_offset( x+1 , y     )];
        mesh.add_connectivity(tri0,node,ordinal0);
        mesh.add_connectivity(node,tri0,ordinal0);
        ++ordinal0;

        node = nodes[node_offset( x+1 , y+1   )];
        mesh.add_connectivity(tri0,node,ordinal0);
        mesh.add_connectivity(node,tri0,ordinal0);
        ++ordinal0;

        samba::connectivity_ordinal ordinal1 = {0};

        // ---
        node = nodes[node_offset( x   , y     )];
        mesh.add_connectivity(tri1,node,ordinal1);
        mesh.add_connectivity(node,tri1,ordinal1);
        ++ordinal1;

        node = nodes[node_offset( x+1 , y+1   )];
        mesh.add_connectivity(tri1,node,ordinal1);
        mesh.add_connectivity(node,tri1,ordinal1);
        ++ordinal1;

        node = nodes[node_offset( x   , y+1   )];
        mesh.add_connectivity(tri1,node,ordinal1);
        mesh.add_connectivity(node,tri1,ordinal1);
        ++ordinal1;


      }
    }

    mesh.end_modification();
  }

};


}} //namespace samba::fixtures

#endif //SAMBA_FIXTURES_TRI_FIXTURE_HPP
