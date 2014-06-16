#ifndef SAMBA_SAMBA_IO_MESH_READER_FIXTURE_HPP
#define SAMBA_SAMBA_IO_MESH_READER_FIXTURE_HPP

#include <samba_io/mesh_reader.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <init/Ionit_Initializer.h>

namespace samba {
namespace io {

class mesh_reader_fixture
{
public:

  static samba::connectivity_map make_connectivity_map_supporting_sidesets()
  {
    // Need face-element back-connectivity for certain sideset capabilities.
    samba::connectivity_map retval = samba::connectivity_map::default_map();
    retval (samba::entity_rank::face(), samba::entity_rank::element()) = samba::connectivity_kind::dynamic();

    return retval;
  }

  /** Constructor uses the mesh specified by mesh_type and file_name
  *   to populate the modifiable_mesh 'mesh'.
  */
  mesh_reader_fixture(stk::ParallelMachine comm,
                      const std::string& mesh_type, const std::string& file_name)
    : m_mesh(make_connectivity_map_supporting_sidesets())
    , m_region(0)
    , m_node_coords(m_mesh)
  {
    Ioss::Init::Initializer init_db;
    Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(mesh_type, file_name, Ioss::READ_MODEL, comm);

    m_region = new Ioss::Region(dbi, "input_model");

    m_reader = mesh_reader::Ptr(new mesh_reader(m_mesh, m_region));

    m_reader->process_mesh();

    m_node_coords = m_reader->get_ioss_mapper()->m_node_coordinates;

    m_reader->process_sidesets();

    m_ioss_mapper_ptr = m_reader->get_ioss_mapper();
  };

  ~mesh_reader_fixture()
  {
    delete m_region;
  }

  samba::mesh m_mesh;

  samba::io::mesh_reader::Ptr m_reader;

  ioss_mapper::Ptr m_ioss_mapper_ptr;

  Ioss::Region *m_region;

  samba::io::coordinates_field_type  m_node_coords;
};

} // namespace io
} // namespace samba

#endif
