#include <sstream>
#include <samba_io/sidesets_manager.hpp>
#include <samba_io/ioss_topology_map.hpp>

using namespace samba;
using namespace samba::io;


SidesetPtr sidesets_manager::add_sideset(const std::string &name)
{
  SidesetPtr sideset_ptr(new Sideset);

  sideset_ptr->m_rank = (m_mesh.spatial_dimension() == spatial_dimension::create(3)
                         ? entity_rank::face()
                         : entity_rank::edge());

  sideset_ptr->m_name = name;

  const bool inducible = true;
  sideset_ptr->m_entity_block = m_mesh.add_entity_block(name, sideset_ptr->m_rank, inducible);

  if (sideset_ptr->m_entity_block == entity_block_key::invalid())
  {
    return SidesetPtr();
  }

  return sideset_ptr;
}


std::vector<SidesetPtr>  sidesets_manager::process_region(Ioss::Region *io_region, const mesh_reader &reader)
{
  std::vector<SidesetPtr> sidesets_out;

  m_mesh.begin_modification();

  const Ioss::SideSetContainer& side_sets = io_region->get_sidesets();

  const entity_key_vector &element_keys = reader.get_entity_keys(entity_rank::element());

  for(Ioss::SideSetContainer::const_iterator sideset_iterator=side_sets.begin(),
        sideset_end=side_sets.end()
        ; sideset_iterator != sideset_end
        ; ++sideset_iterator)
  {
    SidesetPtr sideset = process_sideset(*sideset_iterator, element_keys);
    sidesets_out.push_back(sideset);
  }

  m_mesh.end_modification();

  return sidesets_out;
}


SidesetPtr 
sidesets_manager::process_sideset(const Ioss::SideSet* sset, const entity_key_vector &element_keys)
{
  if (!mesh_reader::supports_sidesets(m_mesh))
  {
    return SidesetPtr();
  }

  const bool use_face_element_connectivity = mesh_reader::supports_direct_get_elements(m_mesh);

  SidesetPtr sideset_ptr = add_sideset(sset->name());
  size_t block_count = sset->block_count();

  // std::cout << "sideset \"" << sset->name() << " has " << block_count << " blocks" << std::endl;

  for(size_t i=0; i<block_count; ++i)
  {
    Ioss::SideBlock* block = sset->get_block(i);

#ifdef SAMBA_IO_ASSUME_SIDEBLOCK_HAS_SINGLE_TOPOLOGY
    entity_topology sideblock_topology = map_topology_ioss_to_cell(block->topology());
    // std::cout << "sideblock_topology = " << sideblock_topology << std::endl;
#endif

    detail::SideBlock ss_block;
    ss_block.m_block_name = block->name();

    // std::cout << "Sideset block named \"" << ss_block.m_block_name << std::endl;

    std::vector<int> element_side_data;
    block->get_field_data("element_side_raw", element_side_data);
    size_t element_side_data_sz = element_side_data.size();

#ifdef SAMBA_IO_ASSUME_SIDEBLOCK_HAS_SINGLE_TOPOLOGY
    size_t num_sides = element_side_data_sz / 2;
    const entity_key_interval side_keys =
      m_mesh.add_entities(sideblock_topology, num_sides, &sideset_ptr->m_entity_block,
                          &sideset_ptr->m_entity_block + 1);
#endif

    size_t ns = 0;
    size_t nn = 0;

    size_t ii = 0;
    while (ii < element_side_data_sz )
    {
      int elem_num = element_side_data[ii++] - 1;    // Because Exodus starts numbering at 1...
      int local_side = element_side_data[ii++] - 1;  // Because Exodus starts numbering at 1..

      entity_key elem_key      = element_keys[elem_num];
      entity_topology elem_topology = elem_key.topology();
      entity_proxy elem_proxy = m_mesh[elem_key];

      entity_topology side_topology = samba::side_topology(elem_topology, local_side);

#ifdef SAMBA_IO_ASSUME_SIDEBLOCK_HAS_SINGLE_TOPOLOGY
      if (sideblock_topology != side_topology)
      {
        throw std::runtime_error("ERROR:  samba::io currently expects a SideBlock to have a single topology.");
      }
      entity_key side_key = side_keys[ns];
#else
      const entity_key_interval side_keys =
        m_mesh.add_entities(side_topology, 1, &sideset_ptr->m_entity_block,
                            &sideset_ptr->m_entity_block + 1);
      entity_key side_key = side_keys[0];
#endif

      connectivity_ordinal side_ordinal = {local_side};
      m_mesh.add_connectivity(elem_key, side_key, side_ordinal);
      if (use_face_element_connectivity)
      {
        m_mesh.add_connectivity(side_key, elem_key, side_ordinal);
      }

      entity_key_iterator elem_nodes_begin = elem_proxy.begin_nodes<entity_key>();
      connectivity_ordinal const* const side_nodes = samba::side_nodes(elem_topology, local_side);

      const size_t nodes_per_side = num_nodes(side_topology);
      for (size_t nj = 0; nj < nodes_per_side; ++nj)
      {
        const entity_key node_to = *(elem_nodes_begin + side_nodes[nj]());
        const connectivity_ordinal ro_j = {nj};
        m_mesh.add_connectivity(side_key, node_to, ro_j);

        ++nn;
      }
      ++ns;
    }
    // std::cout << "  processed " << ns << " sides facets " << " with " << nn << " nodes."
    //           << std::endl;
    sideset_ptr->m_blocks.push_back(ss_block);
  }

  return sideset_ptr;
}

