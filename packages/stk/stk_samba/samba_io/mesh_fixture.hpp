#ifndef STK_SAMBA_MESH_FIXTURE_HPP
#define STK_SAMBA_MESH_FIXTURE_HPP

#include <samba_io/mesh_reader.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <init/Ionit_Initializer.h>

namespace samba {
namespace io {


class samba_mesh_fixture {
 public:

  samba_mesh_fixture(stk::ParallelMachine comm, const std::string& mesh_type, const std::string& file_name)
   : m_mesh(), m_coords(m_mesh, samba::entity_topology::node())
  {
    Ioss::Init::Initializer init_db;
    Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(mesh_type, file_name, Ioss::READ_MODEL, comm);
    Ioss::Region *region = new Ioss::Region(dbi, "input_model");

    mesh_reader reader(m_mesh);
    reader.process_mesh(region);

    Ioss::NodeBlockContainer node_blocks = region->get_node_blocks();

    reader.process_field(m_coords, node_blocks.begin(), node_blocks.end(),"mesh_model_coordinates");

    reader.process_ent_field(m_coords, node_blocks.begin(), node_blocks.end(),"mesh_model_coordinates");

    sidesets_manager sidesets_mgr(m_mesh);
    sidesets_mgr.process_region(region, reader);

    delete region;
  }

  ~samba_mesh_fixture(){}

  mesh m_mesh;

  coords_field_type m_coords;
};
}}
#endif
