#include <Ioss_SubSystem.h>
#include <init/Ionit_Initializer.h>

#include <samba_io/io_fixture.hpp>

namespace samba {
namespace io {

void io_fixture::set_input_ioss_region(Ioss::Region *input_region)
{
  // WRITE ME
}


int io_fixture::set_input_step(int step)
{
  // WRITE ME
  return -1;
}


int io_fixture::set_input_step_from_time(double time_arg)
{
  // WRITE ME
  return -1;
}


std::vector<Sideset::Ptr> io_fixture::populate_mesh()
{
  // WRITE ME
  std::vector<Sideset::Ptr> retval;

  return retval;
}


void io_fixture::set_output_ioss_region(Ioss::Region *output_region)
{
  // WRITE ME
}


int io_fixture::set_output_step(int step)
{
  // WRITE ME
  return -1;
}


int io_fixture::set_output_step_from_time(double time_arg)
{
  // WRITE ME
  return -1;
}

void io_fixture::set_io_entity_block(entity_block_key ebk, bool val)
{
  // WRITE ME
}


bool io_fixture::is_io_entity_block(entity_block_key ebk)
{
  // WRITE ME
  return false;
}

void io_fixture::set_output_filter(set_expression anded_set_sexpression)
{
  // WRITE ME
}

void io_fixture::output_mesh()
{
  // WRITE ME
}


}
}
