#ifndef Stk_Mesh_Use_Cases_UseCase_4_hpp
#define Stk_Mesh_Use_Cases_UseCase_4_hpp

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/fem/FieldTraits.hpp>
#include <stk_mesh/fem/EntityTypes.hpp>
#include <stk_mesh/fem/Stencils.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

/** stk_mesh Use Case 4
 */

namespace stk {
namespace mesh {
namespace use_cases {



typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType ;
typedef stk::mesh::Field<double>                      ScalarFieldType ;
typedef stk::mesh::Field<double*,stk::mesh::ElementNode> ElementNodePointerFieldType ;

// Two part MetaData with four entity types:  
// Node, Edge, Face, Element 
// and two parts (partLeft and partRight)
// and three fields (coordinates, temperature, and volume)
class UseCase_4_MetaData
{
  public:

    const stk::mesh::MetaData & metaData() const
    { return m_metaData; }

    Part & block_hex20() const 
      { return m_block_hex20; }

    Part & block_wedge15() const 
      { return m_block_wedge15; }

    Part & part_vertex_nodes() const 
      { return m_part_vertex_nodes; }

    Part & side_part() const 
      { return m_side_part; }

    VectorFieldType & coordinates_field()
    { return m_coordinates_field; }

    const VectorFieldType & const_coordinates_field() const
    { return m_coordinates_field; }

    VectorFieldType & velocity_field() 
    { return m_velocity_field; }

    const VectorFieldType & const_velocity_field() const
    { return m_velocity_field; }

    VectorFieldType & centroid_field() 
    { return m_centroid_field; }

    const VectorFieldType & const_centroid_field() const
    { return m_centroid_field; }

    ScalarFieldType & temperature_field()
    { return m_temperature_field; }

    const ScalarFieldType & const_temperature_field() const
    { return m_temperature_field; }

    ScalarFieldType & pressure_field()
    { return m_pressure_field; }

    const ScalarFieldType & const_pressure_field() const
    { return m_pressure_field; }

    VectorFieldType & boundary_field() 
    { return m_boundary_field; }

    const VectorFieldType & const_boundary_field() const
    { return m_boundary_field; }


    ElementNodePointerFieldType & element_node_coordinates_field()
    { return m_element_node_coordinates_field; }

    const ElementNodePointerFieldType & const_element_node_coordinates_field() const
    { return m_element_node_coordinates_field; }

  protected:
    UseCase_4_MetaData(const std::vector<std::string> & entity_type_names);

    stk::mesh::MetaData m_metaData;

    Part & m_block_hex20;
    Part & m_block_wedge15;
    Part & m_part_vertex_nodes;
    Part & m_side_part;

    VectorFieldType & m_coordinates_field;
    VectorFieldType & m_velocity_field;
    VectorFieldType & m_centroid_field;
    ScalarFieldType & m_temperature_field;
    ScalarFieldType & m_pressure_field;
    VectorFieldType & m_boundary_field;
    ElementNodePointerFieldType & m_element_node_coordinates_field;

};


// Mesh from two part Meta Data above plus 
// field_data_chunk_size = 1000
class UseCase_4_Mesh : public UseCase_4_MetaData 
{
  public:
    enum {
      field_data_chunk_size = 1000
    };

    UseCase_4_Mesh( stk::ParallelMachine comm );

    ~UseCase_4_Mesh();

    const stk::mesh::BulkData & bulkData() const
    { return m_bulkData; }

    stk::mesh::BulkData       & modifiableBulkData()
    { return m_bulkData; }

  private:

    stk::mesh::BulkData m_bulkData;

};

void populate( UseCase_4_Mesh & mesh );
void runAlgorithms( UseCase_4_Mesh & mesh );
bool verifyMesh( const UseCase_4_Mesh & mesh );

} //namespace use_cases 
} //namespace mesh 
} //namespace stk 

#endif // Stk_Mesh_Use_Cases_UseCase_4_hpp
