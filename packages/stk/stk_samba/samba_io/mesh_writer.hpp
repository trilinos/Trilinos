#ifndef SAMBA_SAMBA_IO_MESH_WRITER_HPP
#define SAMBA_SAMBA_IO_MESH_WRITER_HPP

#include <Ioss_SubSystem.h>
#include <samba_io/mesh_reader.hpp>
#include <samba_io/ioss_field_helper.hpp>
#include <stk_util/environment/ReportHandler.hpp>

namespace samba {
namespace io {
namespace detail {

template<typename DimFunctorT>
struct copy_data_from_field
{
  template<typename FieldT, typename IndexT, typename OutputIteratorT>
  inline
  void operator() (FieldT field_arg, IndexT idx, OutputIteratorT first_out)
  {
    std::copy(field_arg[idx], field_arg[idx] +  field_arg.dimension(idx), first_out);
  }
};


template <>
struct copy_data_from_field<scalar_functor>
{
  template<typename FieldT, typename IndexT, typename OutputIteratorT>
  inline
  void operator() (FieldT field_arg, IndexT idx, OutputIteratorT first_out)
  {
    *first_out = field_arg[idx];
  }
};

} // namespace detail inside samba::io

class mesh_writer
{
public:

  // Simple constructor.  Will add filtering samba::set_expression argument later.
  mesh_writer(MPI_Comm comm,
              const std::string& file_name,
              mesh mesh_arg,
              ioss_mapper::Ptr ioss_mapper_arg = ioss_mapper::Ptr() );

  ~mesh_writer();

  // Needs to be tested.
  template <typename Field_T>
  void write_nodal_field(Field_T data_values, const std::string &name, double time);

  void write_data_vector(const std::vector<double>& data, size_t dim, const std::string &name, double time);

private:

  void define_output_db(Ioss::Region& io_region);
  void write_output_db(Ioss::Region& io_region);

  void define_node_block(entity_block_key node_block_key, Ioss::Region& io_region);
  void define_element_block(entity_block_key elt_block_key, Ioss::Region &io_region);
  void define_sideset(entity_block_key sideset_key, Ioss::Region &io_region);
  void define_nodeset(entity_block_key nodeset_key, Ioss::Region &io_region);

  void write_node_block(entity_block_key node_block_key, Ioss::NodeBlock& nb);
  void write_element_block(entity_block_key elt_block_key, Ioss::ElementBlock &eb);
  void write_sideset(Ioss::SideSet &io_sideset);

  // This one still needs to be written.
  void write_nodeset(entity_block_key nodeset_key, Ioss::Region &io_region);

  void write_node_coordinates( Ioss::NodeBlock& nb);

  Ioss::Region* m_io_region;
  Ioss::PropertyManager* m_property_manager;

  mesh m_mesh;

  ioss_mapper::Ptr m_ioss_mapper_ptr;

  int m_next_element_local_id;

  std::set<std::string> m_transient_fields;
};

inline void mesh_writer::write_data_vector(const std::vector<double>& data, size_t dim, const std::string &name, double time)
{
  Ioss::NodeBlock* nb = m_io_region->get_node_blocks()[0];

  std::set<std::string>::iterator iter = m_transient_fields.find(name);

  if (iter == m_transient_fields.end())
  {
    Ioss::Field::BasicType basic_type = Ioss::Field::REAL;

    std::string storage_type;

    if ( dim == 2 )
    {
        storage_type = "vector_2d";
    }
    else if ( dim == 3 )
    {
        storage_type = "vector_3d";
    }
    else if (dim == 1)
    {
        storage_type = "REAL[1]";
    }
    else
    {
        ThrowRequireMsg(false, "Unsupported dim " << dim);
    }

    size_t num_entities = data.size()/dim;

    m_transient_fields.insert(name);
    m_io_region->begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
    Ioss::Field field(name, basic_type, storage_type,  Ioss::Field::TRANSIENT, num_entities);
    nb->field_add(field);
    m_io_region->end_mode(Ioss::STATE_DEFINE_TRANSIENT);
  }

  m_io_region->begin_mode(Ioss::STATE_TRANSIENT);
  int step = m_io_region->add_state(time);
  m_io_region->begin_state(step);
  nb->put_field_data(name, const_cast<std::vector<double> &>(data));
  m_io_region->end_state(step);
  m_io_region->end_mode(Ioss::STATE_TRANSIENT);
}

template <class Field_T>
void mesh_writer::write_nodal_field(Field_T field_arg, const std::string &name, double time)
{
  detail::copy_data_from_field<typename Field_T::dimension_functor> copy_data_out;

  const size_t num_nodes =  m_mesh.num_nodes();
  if (num_nodes == 0)
  {
    std::cout << "samba::io::mesh_writer::write_nodal_field  --- field has no data; skipping"
              << std::endl;
    return;
  }

  Ioss::NodeBlock* nb = m_io_region->get_node_blocks()[0];
  std::set<std::string>::iterator iter = m_transient_fields.find(name);

  detail::IossFieldHelper<typename Field_T::data_type, typename Field_T::dimension_functor> configure_field;
  Ioss::Field::BasicType basic_type = configure_field.basicType();
  std::string storage_type = configure_field.storage(m_mesh);

  // std::cout << "field basic_type is " << basic_type << "  storage_type is " << storage_type << std::endl;

  if (storage_type == "UNKNOWN")
  {
    throw std::runtime_error("Unsupported field.");
  }

  if (iter == m_transient_fields.end())
  {
    m_transient_fields.insert(name);
    m_io_region->begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
    Ioss::Field field(name, basic_type, storage_type, Ioss::Field::TRANSIENT, num_nodes);
    nb->field_add(field);
    m_io_region->end_mode(Ioss::STATE_DEFINE_TRANSIENT);
  }

  entity_key_vector nodes;
  m_mesh.get_entities(entity_rank::node(), nodes);

  if (num_nodes != nodes.size())
  {
    throw std::runtime_error("number of nodes from m_mesh.get_entities(..) does not match m_mesh.num_nodes()");
  }

  ioss_mapper &mapper = *m_ioss_mapper_ptr;
  entity_key_to_ioss_id_field to_local_id = mapper.m_to_ioss_local_id;
  size_t dimension = field_arg.dimension(nodes[0]);

  std::vector<typename Field_T::data_type> field_data(num_nodes * dimension);
  for (size_t i = 0; i < num_nodes; ++i)
  {
    entity_key node_key = nodes[i];
    int local_id = to_local_id[node_key];
    typename Field_T::data_type *entry_begin = &field_data[dimension * local_id];
    copy_data_out(field_arg, node_key, entry_begin);
  }

  m_io_region->begin_mode(Ioss::STATE_TRANSIENT);
  int step = m_io_region->add_state(time);
  m_io_region->begin_state(step);
  nb->put_field_data(name, const_cast<std::vector<double> &>(field_data));
  m_io_region->end_state(step);
  m_io_region->end_mode(Ioss::STATE_TRANSIENT);
}

} // namespace io
} // namespace samba

#endif
