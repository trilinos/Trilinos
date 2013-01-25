#include <stk_io/util/IO_Fixture.hpp>
#include <stk_io/MeshReadWriteUtils.hpp>

#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>

namespace stk {
namespace io {
namespace util {

IO_Fixture::IO_Fixture(stk::ParallelMachine comm)
  : m_comm(comm),
    m_mesh_data()
{}

IO_Fixture::~IO_Fixture()
{}

void IO_Fixture::output_ioss_region(Teuchos::RCP<Ioss::Region> ioss_output_region) {
  m_mesh_data.set_output_io_region(ioss_output_region);
}

void IO_Fixture::create_output_mesh(
                                    const std::string & base_exodus_filename,
                                    bool  add_transient,
                                    bool  add_all_fields
                                   )
{
  // TODO: Check that the meta-data has a coordinates field and IO parts defined.
  m_mesh_data.create_output_mesh(base_exodus_filename);

  if (add_transient) {
    m_mesh_data.define_output_fields(add_all_fields);
  }
}

void IO_Fixture::add_timestep_to_output_mesh( double time )
{
  m_mesh_data.process_output_request(time);
}

void IO_Fixture::initialize_meta_data( const std::string & base_filename, const std::string & mesh_type)
{
  m_mesh_data.create_input_mesh(mesh_type, base_filename, m_comm);
}

void IO_Fixture::initialize_bulk_data()
{
  m_mesh_data.populate_bulk_data();
}

void IO_Fixture::set_input_ioss_region( Teuchos::RCP<Ioss::Region> input_region )
{
  m_mesh_data.set_input_io_region(input_region);
}

} // namespace util
} // namespace io
} // namespace stk
