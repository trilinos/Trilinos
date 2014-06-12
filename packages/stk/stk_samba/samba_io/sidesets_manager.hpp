#ifndef SAMBA_SAMBA_IO_SIDESETS_MANAGER_HPP
#define SAMBA_SAMBA_IO_SIDESETS_MANAGER_HPP
#include <Ioss_SubSystem.h>

#include <samba/mesh.hpp>
#include <samba/field.hpp>
#include <samba_io/mesh_reader.hpp>


namespace samba {
namespace io {

class sidesets_manager
{
public:

  static std::vector<SidesetPtr> process_sidesets(mesh_reader &reader);

  sidesets_manager( mesh &mesh_arg)
    : m_mesh(mesh_arg)
    , m_new_sideset_id(0)
  { }

  std::vector<SidesetPtr> process_region(Ioss::Region *io_region, const mesh_reader &reader);

private:

  SidesetPtr add_sideset(const std::string &name);

  int get_new_sideset_id() { return m_new_sideset_id++; }

  SidesetPtr process_sideset(const Ioss::SideSet* sset, const entity_key_vector &elem_keys);

  mesh m_mesh;

  int  m_new_sideset_id;
};

}  // namespace io
}  // namespace samba

#endif
