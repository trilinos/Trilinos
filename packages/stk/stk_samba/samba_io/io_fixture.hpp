#ifndef SAMBA_SAMBA_IO_IO_FIXTURE_HPP
#define SAMBA_SAMBA_IO_IO_FIXTURE_HPP

#include <Ioss_SubSystem.h>

#include <samba/mesh.hpp>
#include <samba/field.hpp>
#include <samba_io/ioss_mapper.hpp>

namespace samba {
namespace io {


/**
*  Simple class for reading/writing and exodus file; ported/adapted from stk_io
*/
class io_fixture
{
public:

  typedef mesh::entity_key_vector                                     entity_key_vector;

  typedef boost::unordered_map<entity_block_key, Ioss::EntityBlock *> entity_block_to_iossblock_map;
  typedef boost::unordered_map<entity_block_key,                      entity_rank> entity_block_to_rank_map;

  io_fixture(mesh mesh_arg)
    : m_mesh(mesh_arg)
    , m_ioss_mapper(mesh_arg)
    , m_input_region(NULL)
    , m_output_region(NULL)
    , m_ioss_input_step(0)
    , m_ioss_output_step(0)
  {  }

  void set_input_ioss_region(Ioss::Region *input_region);

  // Set the database step for which field data will be read.
  int set_input_step(int step);

  // Set the database step for which field data will be read,
  // computing it for the given time.
  int set_input_step_from_time(double time_arg);

  // Populate mesh from the current m_input_region.
  std::vector<Sideset::Ptr> populate_mesh();

  // Populate the given entity descriptor indexed field with data from the
  // database at the database step that is currently set.
  template <class Field_T, typename IossPartIterator>
  void populate_descriptor_indexed_field(Field_T field_arg,
                                         IossPartIterator part_begin,
                                         IossPartIterator part_end,
                                         const std::string &io_field_name);

  // Populate the given  entity key indexed field with data from the database
  // at the database step that is currently set.
  template <class Field_T, typename IossPartIterator>
  void populate_field(Field_T field_arg,
                       IossPartIterator part_begin,
                       IossPartIterator part_end,
                       const std::string &io_field_name);




  void set_output_ioss_region(Ioss::Region *output_region);

  // Set the database step to which field data will be written.
  int set_output_step(int step);

  // Set the database step to which field data will be written,
  // computing it for the given time.
  int set_output_step_from_time(double time_arg);

  // Set whether mesh (field) data for the given block will written out.
  void set_io_entity_block(entity_block_key ebk, bool val);

  // Query whether mesh (field) data for the given block will be written out.
  bool is_io_entity_block(entity_block_key ebk);

  // Additional requirements ("anded on") for a block to have its
  // data included.
  void set_output_filter(set_expression anded_set_sexpression);

  // Output the mesh to the current m_output_region.
  void output_mesh();

  // Populate the given field data to the database at the database step 
  // that is currently set.
  template <class Field_T, typename IossPartIterator>
  void output_field(Field_T field_arg, const std::string &io_field_name);

private:

  mesh m_mesh;

  ioss_mapper m_ioss_mapper;

  Ioss::Region *m_input_region;
  Ioss::Region *m_output_region;

  int m_ioss_input_step;
  int m_ioss_output_step;

  std::set<entity_block_key> m_io_entity_blocks;

  set_expression m_anded_set_expression;

};


} // namespace io
} // namespace samba


#endif
