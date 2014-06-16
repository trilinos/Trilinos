#include <samba_io/read_write_utils.hpp>

namespace samba {
namespace io {


void show_mesh_help()
{

}

void create_input_mesh(const std::string &type
                       , const std::string &filename
                       , MPI_Comm comm
                       , mesh mesh_arg
                       , Ioss::Region *input_region)
{


}


void populate_bulk_data(mesh mesh_arg, Ioss::Region *input_region)
{


}

void define_input_fields(mesh mesh_arg, Ioss::Region *input_region)
{

}


void process_input_request(mesh mesh_arg, Ioss::Region *input_region, int step)
{

}

void process_input_request(mesh mesh_arg, Ioss::Region *input_region, double time)
{

}


void create_output_mesh(const std::string &filename
                        , MPI_Comm comm
                        , mesh mesh_arg
                        , Ioss::Region *input_region
                        , Ioss::Region *output_region )
{

}

void define_output_fields(Ioss::Region *output_region, mesh mesh_arg, bool add_all_fields)
{

}

int process_output_request(Ioss::Region *output_region, mesh mesh_arg, double time)
{

  return -1;
}

}
}
