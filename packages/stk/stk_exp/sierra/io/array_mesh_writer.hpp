#ifndef ARRAY_MESH_WRITER_HPP_
#define ARRAY_MESH_WRITER_HPP_

#include <sierra/mesh/array_mesh/array_mesh.hpp>
#include <sierra/io/array_mesh_ioss_topology.hpp>
#include <Ioss_SubSystem.h>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <set>

namespace sierra
{
namespace mesh
{
namespace io
{

class array_mesh_writer
{
public:
    array_mesh_writer(MPI_Comm comm, const std::string& file_name, array_mesh& mesh, const std::vector<double>& coordinates);

    ~array_mesh_writer()
    {
        if(m_created_io_region)
            delete m_io_region;
        if(m_property_manager)
            delete m_property_manager;
    }

    void write_nodal_field(const std::string& name, const std::vector<double>& field_data, double time);

    Ioss::Region* get_ioss_region()
    {
        return m_io_region;
    }

private:
    void define_output_db(Ioss::Region& io_region);
    void write_output_db(Ioss::Region& io_region, const std::vector<double>& coordinates);

    Ioss::Region* m_io_region;
    Ioss::PropertyManager* m_property_manager;
    bool m_created_io_region;
    array_mesh& m_mesh;
    std::set<std::string> m_transient_fields;
};
//class ArrayMeshWriter

inline array_mesh_writer::array_mesh_writer(MPI_Comm comm, const std::string& file_name, array_mesh& mesh, const std::vector<double>& coordinates) :
        m_io_region(NULL), m_created_io_region(true), m_mesh(mesh), m_transient_fields()
{
    m_property_manager = new Ioss::PropertyManager;
    Ioss::DatabaseIO *dbo = Ioss::IOFactory::create("exodusII", file_name, Ioss::WRITE_RESULTS, comm, *m_property_manager);
    if(dbo == NULL || !dbo->ok())
    {
        std::ostringstream oss;
        oss << "ERROR: Could not open results exodusII database '" << file_name;
        throw std::runtime_error(oss.str());
    }

    // NOTE: 'out_region' owns 'dbo' pointer at this time...
    m_io_region = new Ioss::Region(dbo, "results_output");

    define_output_db(*m_io_region);
    write_output_db(*m_io_region, coordinates);
}

inline
void array_mesh_writer::define_output_db(Ioss::Region& io_region)
{
    io_region.begin_mode(Ioss::STATE_DEFINE_MODEL);

    const int spatial_dimension = 3; //TODO: get this from the mesh
    m_io_region->property_add(Ioss::Property("spatial_dimension", spatial_dimension));

    array_mesh::BlockRange blocks = m_mesh.get_blocks();
    for(array_mesh::BlockIterator block_iterator = blocks.first, block_end = blocks.second; block_iterator != block_end; ++block_iterator)
    {
        int rank = m_mesh.get_rank(*block_iterator);
        int topology = m_mesh.get_topology(*block_iterator);
        size_t block_size = m_mesh.get_num_elems(*block_iterator);
        std::string topology_name = map_array_mesh_topology_to_ioss(topology);
        const std::string& name = m_mesh.get_name(*block_iterator);

        if(rank == array_mesh::Node)
        {
            Ioss::NodeBlock* nb = new Ioss::NodeBlock(io_region.get_database(), name, block_size, spatial_dimension);
            io_region.add(nb);
        }

        if(rank == array_mesh::Element)
        {
            Ioss::ElementBlock *eb = new Ioss::ElementBlock(io_region.get_database(), name, topology_name, block_size);
            io_region.add(eb);
        }
    }

    array_mesh::NodesetRange nodesets = m_mesh.get_nodesets();
    for(array_mesh::NodesetIterator nodeset_iterator = nodesets.first, nodeset_end = nodesets.second; nodeset_iterator != nodeset_end; ++nodeset_iterator)
    {
        Ioss::NodeSet* ns = new Ioss::NodeSet(io_region.get_database(), m_mesh.get_name(*nodeset_iterator), m_mesh.get_size(*nodeset_iterator));
        io_region.add(ns);
    }

    array_mesh::SidesetRange sidesets = m_mesh.get_sidesets();
    for(array_mesh::SidesetIterator sideset_iterator = sidesets.first, sideset_end = sidesets.second; sideset_iterator != sideset_end; ++sideset_iterator)
    {
        Ioss::SideSet* ss = new Ioss::SideSet(io_region.get_database(), m_mesh.get_name(*sideset_iterator));
        io_region.add(ss);

        size_t num_sides = m_mesh.get_sideset_elem_nums(*sideset_iterator).size();
        Ioss::SideBlock* sb = new Ioss::SideBlock(ss->get_database(), m_mesh.get_name(*sideset_iterator), "unknown", "unknown", num_sides);
        ss->add(sb);
    }

    io_region.end_mode(Ioss::STATE_DEFINE_MODEL);
}

inline
void array_mesh_writer::write_output_db(Ioss::Region& io_region, const std::vector<double>& coordinates)
{
    io_region.begin_mode(Ioss::STATE_MODEL);

    Ioss::NodeBlock& nb = *io_region.get_node_blocks()[0]; //TODO hard-coded that only 1 node-block exists?
    if(io_region.get_node_blocks().size() > 1)
    {
        throw std::runtime_error("ERROR in ArrayMeshWriter, assumed only 1 node-block, found more...");
    }
    std::vector<int> node_ids(m_mesh.get_node_ids().first, m_mesh.get_node_ids().second);
    size_t num_ids_written = nb.put_field_data("ids", node_ids);
    if(node_ids.size() != num_ids_written)
    {
        throw std::runtime_error("FAILED in Ioss::NodeBlock::put_field_data");
    }

    const Ioss::NodeSetContainer& node_sets = io_region.get_nodesets();
    for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin(); it != node_sets.end(); ++it)
    {
        std::string name = (*it)->name();
        array_mesh::NodesetIndex nodeset_index = m_mesh.get_nodeset(name);
        array_mesh::ConstIntRange node_range = m_mesh.get_nodes(nodeset_index);
        std::vector<int> nodes(node_range.first, node_range.second);
        (*it)->put_field_data("ids_raw", nodes);
    }

    const Ioss::SideSetContainer& side_sets = io_region.get_sidesets();
    for(Ioss::SideSetContainer::const_iterator it = side_sets.begin(); it != side_sets.end(); ++it)
    {
        Ioss::SideSet* ss = *it;
        const std::string& sideset_name = ss->name();
        array_mesh::SidesetIndex sideset_index = m_mesh.get_sideset(sideset_name);
        const std::vector<int>& elem_nums = m_mesh.get_sideset_elem_nums(sideset_index);
        const std::vector<int>& local_sides = m_mesh.get_sideset_elem_local_sides(sideset_index);
        std::vector<int> raw_side_data;
        raw_side_data.reserve(elem_nums.size() + local_sides.size());
        for(size_t i = 0; i < elem_nums.size(); ++i)
        {
            //add 1 because exodus/ioss is 1-based
            raw_side_data.push_back(elem_nums[i] + 1);
            raw_side_data.push_back(local_sides[i] + 1);
        }

        //we know there is only 1 side-block (we set up the sideset in the above define_output_db method)
        Ioss::SideBlock* sb = ss->get_side_blocks()[0];

        const size_t num_sides_written = sb->put_field_data("element_side_raw", raw_side_data);
        if(num_sides_written != raw_side_data.size() / 2)
        {
            throw std::runtime_error("ArrayMeshWriter ERROR, failed in sideset output.");
        }
    }

    //TODO: why does put_field_data take a non-const std::vector?
    size_t num_coords_written = nb.put_field_data("mesh_model_coordinates", const_cast<std::vector<double>&>(coordinates));
    if(coordinates.size() / 3 != num_coords_written)
    {
        throw std::runtime_error("FAILED in Ioss::NodeBlock::put_field_data for coordinates.");
    }

    std::vector<int> elem_ids(m_mesh.get_element_ids().first, m_mesh.get_element_ids().second);

    const Ioss::ElementBlockContainer& elem_blocks = io_region.get_element_blocks();
    for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin(), end = elem_blocks.end(); it != end; ++it)
    {
        array_mesh::BlockIndex block_index = m_mesh.get_block((*it)->name());
        const std::vector<int>& elem_nums = m_mesh.get_block_elem_nums(block_index);
        std::vector<int> block_elem_ids(elem_nums.size());
        for(size_t i=0; i<elem_nums.size(); i++)
        {
            block_elem_ids[i] = elem_ids[elem_nums[i]];
        }

        (*it)->put_field_data("ids", block_elem_ids);

        std::vector<int> connectivity = m_mesh.get_block_connectivity(block_index);

        for(size_t i = 0; i < connectivity.size(); ++i) {
            connectivity[i] += 1;
        }

        (*it)->put_field_data("connectivity_raw", connectivity);
    }

    io_region.end_mode(Ioss::STATE_MODEL);
}

inline
void array_mesh_writer::write_nodal_field(const std::string& name, const std::vector<double>& field_data, double time)
{
    Ioss::NodeBlock* nb = m_io_region->get_node_blocks()[0];
    const size_t num_nodes = m_mesh.get_num_nodes();
    if(num_nodes > 0 && field_data.size() == 0)
    {
        std::cout << "ArrayMeshWriter::write_nodal_data WARNING, no field-data provided, skipping this write." << std::endl;
        return;
    }

    std::set<std::string>::iterator iter = m_transient_fields.find(name);
    if(iter == m_transient_fields.end())
    {
        m_transient_fields.insert(name);
        m_io_region->begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
        if(field_data.size() == num_nodes * 3)
        {
            Ioss::Field field(name, Ioss::Field::REAL, "vector_3d", Ioss::Field::TRANSIENT, num_nodes);
            nb->field_add(field);
        }
        else if(field_data.size() == num_nodes)
        {
            Ioss::Field field(name, Ioss::Field::REAL, "REAL[1]", Ioss::Field::TRANSIENT, num_nodes);
            nb->field_add(field);
        }
        else
        {
            //TODO: use ThrowRequire or BOOST_ASSERT
            throw std::runtime_error("Unsupported field: not 3D-vector or scalar!!");
        }
        m_io_region->end_mode(Ioss::STATE_DEFINE_TRANSIENT);
    }

    m_io_region->begin_mode(Ioss::STATE_TRANSIENT);
    int step = m_io_region->add_state(time);
    m_io_region->begin_state(step);
    nb->put_field_data(name, const_cast<std::vector<double>&>(field_data));
    m_io_region->end_state(step);
    m_io_region->end_mode(Ioss::STATE_TRANSIENT);
}

} //namespace io
} //namespace mesh
} //namespace sierra

#endif /* ARRAY_MESH_WRITER_HPP_ */
