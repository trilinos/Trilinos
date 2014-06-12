#ifndef SAMBA_SAMBA_IO_READ_WRITE_UTILS_HPP
#define SAMBA_SAMBA_IO_READ_WRITE_UTILS_HPP

#include <Ioss_SubSystem.h>

#include <samba/mesh.hpp>
#include <samba/field.hpp>


namespace samba {
namespace io {

void show_mesh_help();


Ioss::Region *
create_input_mesh_region(const std::string &type
                         , const std::string &filename
                         , MPI_Comm comm);


// This is analagous to process_mesh.
void populate_bulk_data(mesh mesh_arg, Ioss::Region *input_region);


void define_input_fields(mesh mesh_arg, Ioss::Region *input_region);


// The guts of these two is basically process_fields.
void process_input_request(mesh mesh_arg, Ioss::Region *input_region, int step);
void process_input_request(mesh mesh_arg, Ioss::Region *input_region, double time);



void create_output_mesh(const std::string &filename
                        , MPI_Comm comm
                        , mesh mesh_arg
                        , Ioss::Region *input_region
                        , Ioss::Region *output_region );

void define_output_fields(Ioss::Region *output_region, mesh mesh_arg, bool add_all_fields = false);

int process_output_request(Ioss::Region *output_region, mesh mesh_arg, double time);

template <typename INT>
void get_element_block_sizes(Ioss::Region *input_region, std::vector<INT>& el_blocks)
{


}


} // namespace io
} // namespace samba


#endif
