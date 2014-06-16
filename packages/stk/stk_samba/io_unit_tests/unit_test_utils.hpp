#include <gtest/gtest.h>

#include <Ioss_SubSystem.h>

#include <samba_io/mesh_reader.hpp>
#include <samba_io/mesh_reader_fixture.hpp>

#include <init/Ionit_Initializer.h>
#include <stk_util/use_cases/UseCaseEnvironment.hpp>

typedef samba::io::coordinates_field_type coords_field_type;

Ioss::Region *make_region(const std::string mesh_type, const std::string file_name);

template <typename Index>
size_t get_nodes(samba::entity_proxy eproxy, std::vector<Index> &nodes_out)
{
  nodes_out.clear();

  for (samba::entity_key_iterator nodes_iter = eproxy.begin_nodes<Index>(), nodes_end = eproxy.end_nodes<Index>();
       nodes_iter != nodes_end;
       ++nodes_iter)
  {
    nodes_out.push_back(*nodes_iter);
  }
  return nodes_out.size();
}

bool is_nodally_equal(samba::entity_key entity1, samba::mesh mesh1, samba::io::ioss_mapper &mapper1,
                      samba::entity_key entity2, samba::mesh mesh2, samba::io::ioss_mapper &mapper2,
                      double epsilon = 0.000001);

template <typename Index>
inline
void get_coords(Index node, coords_field_type node_coords_field,
                double coords[3])
{
  coords_field_type::const_reference node_coords_ref  = node_coords_field[node];
  for (size_t i = 0; i < 3; ++i)
  {
    coords[i] = node_coords_ref[i];
  }
}


inline bool is_equal_3d(const double vec1[3], const double vec2[3], double epsilon = 0.000001)
{
  double diff[3];
  for (size_t i = 0; i < 3; ++i)
  {
    diff[i] = vec1[i] - vec2[i];
  }
  double distsq = (diff[0] * diff[0]) + (diff[1] * diff[1]) + (diff[2] * diff[2]);
  return (distsq <= epsilon * epsilon);
}
