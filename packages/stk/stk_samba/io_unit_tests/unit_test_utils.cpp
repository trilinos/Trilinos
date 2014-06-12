#include "unit_test_utils.hpp"

Ioss::Region *make_region(const std::string mesh_type, const std::string file_name)
{
  Ioss::Init::Initializer init_db;
  Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(mesh_type, file_name, Ioss::READ_MODEL, 0 );
  return new Ioss::Region(dbi, "input_model");
}

bool is_nodally_equal(samba::entity_key entity1, samba::mesh mesh1, samba::io::ioss_mapper &mapper1,
                      samba::entity_key entity2, samba::mesh mesh2, samba::io::ioss_mapper &mapper2,
                      double epsilon)
{
  samba::entity_topology topology = entity1.topology();
  if (entity2.topology() != topology)
  {
    return false;
  }

  samba::io::entity_key_to_ioss_id_field to_global1 = mapper1.m_to_ioss_global_id;
  samba::io::entity_key_to_ioss_id_field to_global2 = mapper2.m_to_ioss_global_id;
  coords_field_type coords_field1 = mapper1.m_node_coordinates;
  coords_field_type coords_field2 = mapper2.m_node_coordinates;

  if (samba::entity_topology::node() == topology)
  {
    if (to_global1[entity1] != to_global2[entity2])
    {
      return false;
    }

    double coords1[3], coords2[3];
    get_coords(entity1, coords_field1, coords1);
    get_coords(entity2, coords_field2, coords2);
    return is_equal_3d(coords1, coords2, epsilon);
  }

  samba::entity_proxy e_proxy1 = mesh1[entity1];
  samba::entity_proxy e_proxy2 = mesh2[entity2];

  size_t num_nodes1 = e_proxy1.num_nodes();
  if (num_nodes1 != e_proxy2.num_nodes())
  {
    return false;
  }

  samba::entity_key_iterator nodes1_begin, nodes2_begin;
  nodes1_begin = e_proxy1.begin_nodes<samba::entity_key>();
  nodes2_begin = e_proxy2.begin_nodes<samba::entity_key>();
  for (size_t i = 0; i < num_nodes1; ++i)
  {
    samba::entity_key node1, node2;
    node1 = *(nodes1_begin + i);
    node2 = *(nodes2_begin + i);

    if (to_global1[node1] != to_global2[node2])
    {
      return false;
    }

    double coords1[3], coords2[3];
    get_coords(node1, coords_field1, coords1);
    get_coords(node2, coords_field2, coords2);
    if (!is_equal_3d(coords1, coords2, epsilon))
    {
      return false;
    }
  }

  return true;
}
