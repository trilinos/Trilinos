#include <gtest/gtest.h>

#include <Ioss_SubSystem.h>

#include <samba_io/mesh_reader.hpp>
#include <samba_io/mesh_reader_fixture.hpp>
#include <io_unit_tests/unit_test_utils.hpp>

#include <init/Ionit_Initializer.h>
#include <stk_util/use_cases/UseCaseEnvironment.hpp>

typedef samba::io::coordinates_field_type coords_field_type;

Ioss::Region *make_region(const std::string mesh_type, const std::string file_name)
{
  Ioss::Init::Initializer init_db;
  Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(mesh_type, file_name, Ioss::READ_MODEL, 0 );

  return new Ioss::Region(dbi, "input_model");
}


bool check_aa_coplanar(const std::vector<samba::entity_key> &nodes,
                       const coords_field_type node_coords_field)
{
  if (nodes.empty())
  {
    return false;
  }
  const samba::entity_key anchor_node = nodes[0];
  if (anchor_node.topology() != samba::entity_topology::node())
  {
    return false;
  }
  coords_field_type::const_reference anchor_ncr = node_coords_field[anchor_node];
  double anchor_coords[3] = {0, 0, 0};
  anchor_coords[0] = *anchor_ncr;
  anchor_coords[1] = *(anchor_ncr + 1);
  anchor_coords[2] = *(anchor_ncr + 2);

  bool same[3] = {true, true, true};
  size_t num_nodes = nodes.size();
  for (size_t i = 1; i < num_nodes; ++i)
  {
    const samba::entity_key node = nodes[i];
    if (node.topology() != samba::entity_topology::node())
    {
      return false;
    }
    coords_field_type::const_reference node_coords_ref = node_coords_field[node];

    for (size_t j = 0; j < 3; ++j)
    {
      if (std::abs(anchor_coords[j] - node_coords_ref[j]) > 0.000001)
      {
        same[j] = false;
      }
    }
  }

  size_t same_dims = 0;
  size_t not_same_dims = 0;
  for (size_t j = 0; j < 3; ++j)
  {
    if (same[j])
      ++same_dims;
    else
      ++not_same_dims;
  }

  // The nodes lie on a plane perpendicular to an axis, and they
  // do not lie on a line parallel to an axis.
  return ((same_dims == 1) && (not_same_dims == 2) );
}

bool update_cobound_count(samba::entity_key face,
                          samba::mesh &mesh_arg,
                          std::map<samba::entity_key, size_t> &cobound_counter)
{
  samba::entity_proxy face_proxy = mesh_arg[face];

  samba::entity_key_iterator connectivity_iter, connectivity_end;
  connectivity_iter = face_proxy.begin_elements<samba::entity_key>();
  connectivity_end  = face_proxy.end_elements<samba::entity_key>();
  if (connectivity_iter == connectivity_end)
  {
    return false;
  }

  for ( ; connectivity_iter != connectivity_end; ++connectivity_iter)
  {
    samba::entity_key curr_elt = *connectivity_iter;
    samba::entity_proxy elt_proxy = mesh_arg[curr_elt];

    bool found_face = false;
    for (samba::entity_key_iterator faces_i = elt_proxy.begin_faces<samba::entity_key>(), faces_end = elt_proxy.end_faces<samba::entity_key>();
         faces_i != faces_end;
         ++faces_i)
    {
      if (*faces_i == face)
      {
        ++cobound_counter[curr_elt];
        found_face = true;
      }
    }
    if (!found_face)
    {
      return false;
    }
  }
  return true;
}


bool check_cobound_count(std::map<samba::entity_key, size_t> &cobound_counter,
                         const size_t single_side_elt_count_arg,
                         const size_t edge_elt_count_arg,
                         const size_t corner_elt_count_arg)
{
  size_t single_side_elt_count =0, edge_elt_count =0, corner_elt_count =0;
  for (std::map<samba::entity_key, size_t>::iterator cbc_i = cobound_counter.begin();
       cbc_i != cobound_counter.end();
       ++cbc_i)
  {
    switch (cbc_i->second)
    {
    case 1:
      ++single_side_elt_count;
      break;
    case 2:
      ++edge_elt_count;
      break;
    case 3:
      ++ corner_elt_count;
      break;
    default:
      return false;
      break;
    }
  }

  return ((single_side_elt_count_arg == single_side_elt_count)
          && (edge_elt_count_arg == edge_elt_count)
          && (corner_elt_count_arg == corner_elt_count));
}


TEST(samba_io, mesh_reader_exodus_basic)
{
  Ioss::Region *region = make_region("exodus", "cube.par");

  samba::mesh mesh;
  samba::io::mesh_reader reader(mesh, region);

  // Read in the entities and connectivity.
  reader.process_mesh();

  EXPECT_EQ(mesh.num_elements(), 1000u);
  EXPECT_EQ(mesh.num_nodes(), 1331u);

  delete region;
}

#if 0
TEST(samba_io, mesh_reader_exodus_coords_fields)
{
  Ioss::Region *region = make_region("exodus", "cube.par");

  samba::mesh mesh(samba::spatial_dimension::create(3));
  samba::io::mesh_reader reader(mesh, region);

  // Right now, samba seems buggy in handling fields constructed after a mesh is populated,
  // we constuct these before calling reader.process_mesh(.).
  desc_idx_coords_field_type  node_coords_by_eDesc(mesh, samba::entity_rank::node());
  key_idx_coords_field_type node_coords_by_eKey(mesh, samba::entity_topology::node());

  //  Read in the entities and connectivity.
  reader.process_mesh();

  Ioss::NodeBlockContainer node_blocks = region->get_node_blocks();

  // Read in the coordinates field to descriptor-indexed version.
  reader.read_nodal_desc_field(node_coords_by_eDesc, "mesh_model_coordinates");

  // Read in the coordinates field to key-indexed version.
  reader.read_nodal_key_field(node_coords_by_eKey, "mesh_model_coordinates");

  // Go through the nodes and do a indexing-invariant fields content check.
  samba::mesh::entity_key_vector node_keys;
  mesh.get_entities(samba::set_expression(samba::entity_rank::node()), node_keys);
  size_t num_nodes = node_keys.size();
  EXPECT_EQ(num_nodes, 1331u);
  for (size_t i = 0; i < num_nodes; ++i)
  {
    samba::entity_key nk_i = node_keys[i];
    key_idx_coords_field_type::reference coords_by_k = node_coords_by_eKey[nk_i];
    desc_idx_coords_field_type::reference coords_by_d = node_coords_by_eDesc[mesh.convert(nk_i)];
    EXPECT_EQ(*coords_by_k,  *coords_by_d);
    EXPECT_EQ(*(coords_by_k + 1), *(coords_by_d + 1));
    EXPECT_EQ(*(coords_by_k + 2), *(coords_by_d + 2));
  }

  delete region;
}
#endif

TEST(samba_io, mesh_reader_exodus_sidesets)
{
  Ioss::Region *region = make_region("exodus", "cube.par");

  // Need face-element back-connectivity for certain sideset capabilities.
  samba::connectivity_map dflt_relmap_plus_face_element = samba::connectivity_map::default_map();
  dflt_relmap_plus_face_element(samba::entity_rank::face(), samba::entity_rank::element())
    = samba::connectivity_kind::dynamic();

  samba::mesh mesh(dflt_relmap_plus_face_element);
  samba::io::mesh_reader reader(mesh, region);

  EXPECT_TRUE(samba::io::mesh_reader::supports_direct_get_elements(mesh));

  coords_field_type node_coords(mesh);

  reader.process_mesh();

  reader.read_nodal_field(node_coords, "mesh_model_coordinates");

  // Read in sidesets.
  EXPECT_TRUE(samba::io::mesh_reader::supports_sidesets(mesh));
  std::vector<samba::io::SidesetPtr> sidesets = reader.process_sidesets();

  //
  // Now do some error checking on the sidesets.
  //
  bool entered = false;
  std::map<samba::entity_key, size_t> cobound_counter;

  size_t num_sidesets = sidesets.size();

  EXPECT_EQ(num_sidesets, 6u);
  for (size_t i = 0; i < num_sidesets; ++i)
  {
    samba::io::SidesetPtr sideset = sidesets[i];
    // EXPECT_EQ(expected_id, sideset->m_id);
    EXPECT_EQ(samba::entity_rank::face()(), sideset->m_rank());

    samba::entity_block_key sideset_key = sideset->m_entity_block;

    samba::mesh::entity_key_vector side_entities;
    mesh.get_entities(sideset_key & samba::entity_rank::face(), side_entities);

    // Go through all the entities on the sideset.
    size_t num_side_entities = side_entities.size();
    EXPECT_GT(num_side_entities, 0u);
    for (size_t k = 0; k < num_side_entities; ++k)
    {
      samba::entity_key side_face = side_entities[k];
      samba::entity_proxy side_face_proxy = mesh[side_face];

      //  std::cout << "face " << side_face << std::endl;

      std::vector<samba::entity_key> face_nodes;
      size_t num_face_nodes = get_nodes(side_face_proxy, face_nodes);
      EXPECT_EQ(num_face_nodes, 4u);
      EXPECT_TRUE(check_aa_coplanar(face_nodes, node_coords));

      bool cobound_elt_found = update_cobound_count(side_face, mesh, cobound_counter);
      EXPECT_TRUE(cobound_elt_found);

      entered = true;
    }
  }

  EXPECT_TRUE(entered);
  EXPECT_TRUE(check_cobound_count(cobound_counter, 384, 96, 8) );

  delete region;
}

////
////  NEED A REGRESSION TEST TO TEST nodesets.
////
