#ifndef SAMBA_SAMBA_IO_MESH_READER_HPP
#define SAMBA_SAMBA_IO_MESH_READER_HPP

#include <Ioss_SubSystem.h>

#include <samba/mesh.hpp>
#include <samba/field.hpp>

#include <samba_io/ioss_mapper.hpp>
#include <samba_io/detail/sidesets_info.hpp>


#define MESH_READER_BABBLES

// This enables sidesets to be handled faster.  We can work on an efficient implementation
// that supports each sideblock (in a sideset) having more than one side topology in it
// if/when it matters.
#define SAMBA_IO_ASSUME_SIDEBLOCK_HAS_SINGLE_TOPOLOGY


namespace samba {
namespace io {

namespace detail {

template<typename DimFunctorT>
struct copy_data_into_field
{
  template<typename InputIteratorT, typename FieldT, typename IndexT>
  void operator()(InputIteratorT first, InputIteratorT last, FieldT field_arg, IndexT idx)
  {
    std::copy(first, last, field_arg[idx]);
  }
};

template <>
struct copy_data_into_field<scalar_functor>
{
  template<typename InputIteratorT, typename FieldT, typename IndexT>
  inline
  void operator()(InputIteratorT first, InputIteratorT last, FieldT field_arg, IndexT idx)
  {
    BOOST_ASSERT(last - first == 1);
    field_arg[idx] = *first;
  }
};

} // namespace detail inside samba::io

class mesh_reader
{
public:

  // Can get rid of this one when mesh::get_rank(entity_block_key) works.
  typedef boost::unordered_map<entity_block_key, entity_rank> entity_block_to_rank_map;

  typedef boost::shared_ptr<mesh_reader> Ptr;

  static bool supports_direct_get_elements(const mesh &mesh_arg);
  static bool supports_sidesets(const mesh &mesh_arg);

  mesh_reader( mesh mesh_arg, Ioss::Region *io_region,
               ioss_mapper::Ptr ioss_mapper_arg = ioss_mapper::Ptr() )
    : m_mesh(mesh_arg)
    , m_io_region(io_region)
    , m_ioss_mapper_ptr(ioss_mapper_arg)
    , m_ioss_mapper_rawptr(0)
    , m_meshRankProcessed(-1)
  {
    if (!m_ioss_mapper_ptr)
    {
      m_ioss_mapper_ptr = ioss_mapper::Ptr(new ioss_mapper(m_mesh));
    }
    m_ioss_mapper_rawptr = m_ioss_mapper_ptr.get();
  }

  // We can drop this if we get the okay to assume that a mesh reader will
  // only be attached to (constructed from) a virgin mesh.
  const entity_key_vector &get_entity_keys(const entity_rank rank) const
  {
    std::map<entity_rank, entity_key_vector>::const_iterator oek_probe = m_ordered_entity_keys.find(rank);
    BOOST_ASSERT(m_ordered_entity_keys.find(rank) != m_ordered_entity_keys.end());
    return oek_probe->second;
  }

  /// Load the entity and connectivity data for the region.
  ///
  void process_mesh();

  std::vector<SidesetPtr> process_sidesets();

  /// Populate the given nodal field with data from the named Ioss field.
  template <class Field_T>
  void read_nodal_field(Field_T field_arg,
                        const std::string &io_field_name);

  template <class Field_T>
  void read_nodal_field(Field_T field_arg,
                        const std::string &io_field_name,
                        double time);

  /// Populate the given field with data from the named Ioss field for the given Ioss parts.
  template <typename Field_T, typename IossPartIterator>
  void process_field_by_parts(Field_T field_arg,
                              IossPartIterator part_begin,
                              IossPartIterator part_end,
                              const std::string &io_field_name);

  /// Needed to use a mesh_writer.
  ioss_mapper::Ptr get_ioss_mapper() { return m_ioss_mapper_ptr;}

private:

  template <typename T>
  void process_blocks(const std::vector<T*> &blocks, const entity_rank rank);

  void process_block(Ioss::EntityBlock *block, const entity_rank rank);

  /// \todo Need a unit test or regression test to hit this.
  void process_nodesets();

  /// \todo Need a unit test or regression test to hit this.
  void process_nodeset(const Ioss::NodeSet *nset);

  /// The mesh into which this is reading.
  mesh m_mesh;

  Ioss::Region *m_io_region;

  /// Shared pointer to structure that augments the mesh to facilitate using, reading, and writing
  /// a model.
  ioss_mapper::Ptr m_ioss_mapper_ptr;

  /// Raw pointer cached from m_ioss_mapper_ptr.
  ioss_mapper *m_ioss_mapper_rawptr;

  /// For each rank, map from local ids in the input file to samba entity keys.
  std::map<entity_rank, entity_key_vector> m_ordered_entity_keys;

  /// Map from Samba entity blocks (sets of entities) to Ioss::EntityBlocks (*s) that have been read
  /// in.
  entity_block_to_iossblock_map  m_entityBlocks;

  /// Used internally, somewhat an artifact from ModifiableMesh heritage.
  int m_meshRankProcessed;
};


template <class Field_T>
void mesh_reader::read_nodal_field(Field_T field_arg, const std::string &io_field_name)
{
  Ioss::NodeBlockContainer node_blocks = m_io_region->get_node_blocks();
  process_field_by_parts(field_arg, node_blocks.begin(), node_blocks.end(),
                         io_field_name);
}


template <class Field_T>
void mesh_reader::read_nodal_field(Field_T field_arg, const std::string &io_field_name,
                                        double time)
{
  // Get state count and all states...
  int step_count = m_io_region->get_property("state_count").get_int();

  int step = -1;
  for (int istep = 1; istep <= step_count; istep++)
  {
    double state_time = m_io_region->get_state_time(istep);
    if (std::abs(state_time - time) < 1.e-12)
    {
      step = istep;
    }
  }

  if (step > -1)
  {
    m_io_region->begin_state(step);

    Ioss::NodeBlockContainer node_blocks = m_io_region->get_node_blocks();
    process_field_by_parts(field_arg, node_blocks.begin(), node_blocks.end(),
                                    io_field_name);
    m_io_region->end_state(step);
  }
  else
  {
    std::cout<<"samba::io::read_nodal_key_field WARNING, requested time ("
             << time << ") not found in input-mesh-file, skipping this read."<<std::endl;
  }
}

template <typename Field_T, typename IossPartIterator>
void mesh_reader::process_field_by_parts(Field_T field_arg,
                                         IossPartIterator part_begin,
                                         IossPartIterator part_end,
                                         const std::string &io_field_name)
{
  detail::copy_data_into_field<typename Field_T::dimension_functor> copy_data_in;

  typedef typename Field_T::index_type IndexType;

  std::vector<entity_block_key> io_parts;
  for( ; part_begin != part_end; ++part_begin)
  {
    entity_block_key part_ebk = m_mesh.find_entity_block((*part_begin)->name());
    BOOST_ASSERT(is_valid(part_ebk));

    io_parts.push_back(part_ebk);
  }

  size_t num_io_parts = io_parts.size();
  for (size_t i = 0; i < num_io_parts; ++i)
  {
    entity_block_key part_ebk = io_parts[i];
    entity_block_to_iossblock_map::iterator esb_probe = m_entityBlocks.find(part_ebk);
    BOOST_ASSERT(esb_probe != m_entityBlocks.end());

    Ioss::EntityBlock const *block = esb_probe->second;
    if (!block->field_exists(io_field_name))
    {
      continue;
    }

    entity_rank rank = m_mesh.get_rank(part_ebk);
    const entity_key_vector &entities = m_ordered_entity_keys[rank];

    std::vector<double> field_data;
    block->get_field_data(io_field_name, field_data);

    const int entity_count_int = block->get_property("entity_count").get_int();
    const int offset = block->get_offset();
    BOOST_ASSERT((entity_count_int >= 0) && (offset >= 0));

    // In general, we have to write the field data one entity at at time because induced
    // part membership can cause blocks to be split over partitions (buckets).
    const size_t entity_count = static_cast<size_t>(entity_count_int);
    for (size_t j = 0; j < entity_count; ++j)
    {
      IndexType ed_j = m_mesh.convert<IndexType>(entities[offset+j]);
      size_t dimension = field_arg.dimension(ed_j);

      int dbegin = (j + offset) * dimension;
      int dend   = dbegin + dimension;
      copy_data_in(&field_data[dbegin], &field_data[dend], field_arg, ed_j);
    }
  }
}

template <typename T>
void mesh_reader::process_blocks(const std::vector<T*> &blocks, const entity_rank rank)
{
  BOOST_ASSERT(m_meshRankProcessed < static_cast<int>(rank()));

  if (m_meshRankProcessed < rank())
  {
    m_meshRankProcessed = rank();
  }

  for(size_t i=0; i < blocks.size(); i++) {
    T* block = blocks[i];
    process_block(block, rank);
  }
}



}  // namespace io
}  // namespace samba

#endif
