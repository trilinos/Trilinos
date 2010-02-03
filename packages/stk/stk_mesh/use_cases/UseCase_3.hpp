#ifndef Stk_Mesh_Use_Cases_UseCase_3_hpp
#define Stk_Mesh_Use_Cases_UseCase_3_hpp

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

/** stk_mesh Use Case 3
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
class UseCase_3_MetaData
{
  public:

    const stk::mesh::MetaData & metaData() const
    { return m_metaData; }

    stk::mesh::Part & block_hex() const
    { return m_block_hex; }
    
    stk::mesh::Part & block_wedge() const
    { return m_block_wedge; }
    
    stk::mesh::Part & block_tet() const
    { return m_block_tet; }
    
    stk::mesh::Part & block_pyramid() const
    { return m_block_pyramid; }
    
    stk::mesh::Part & block_quad_shell() const
    { return m_block_quad_shell; }
    
    stk::mesh::Part & block_tri_shell() const
    { return m_block_tri_shell; }
    
    VectorFieldType & coordinates_field()
    { return m_coordinates_field; }

    const VectorFieldType & const_coordinates_field() const
    { return m_coordinates_field; }

    VectorFieldType & centroid_field() 
    { return m_centroid_field; }

    const VectorFieldType & const_centroid_field() const
    { return m_centroid_field; }

    ScalarFieldType & temperature_field()
    { return m_temperature_field; }

    const ScalarFieldType & const_temperature_field() const
    { return m_temperature_field; }

    ScalarFieldType & volume_field()
    { return m_volume_field; }

    const ScalarFieldType & const_volume_field() const
    { return m_volume_field; }

    ElementNodePointerFieldType & element_node_coordinates_field()
    { return m_element_node_coordinates_field; }

    const ElementNodePointerFieldType & const_element_node_coordinates_field() const
    { return m_element_node_coordinates_field; }

  protected:
    UseCase_3_MetaData(const std::vector<std::string> & entity_type_names);

    stk::mesh::MetaData m_metaData;

    stk::mesh::Part & m_block_hex;
    stk::mesh::Part & m_block_wedge;
    stk::mesh::Part & m_block_tet;
    stk::mesh::Part & m_block_pyramid;
    stk::mesh::Part & m_block_quad_shell;
    stk::mesh::Part & m_block_tri_shell;

    VectorFieldType & m_coordinates_field;
    VectorFieldType & m_centroid_field;
    ScalarFieldType & m_temperature_field;
    ScalarFieldType & m_volume_field;
    ElementNodePointerFieldType & m_element_node_coordinates_field;

};


// Mesh from two part Meta Data above plus 
// field_data_chunk_size = 1000
class UseCase_3_Mesh : public UseCase_3_MetaData 
{
  public:
    enum {
      field_data_chunk_size = 1000
    };

    UseCase_3_Mesh( stk::ParallelMachine comm );

    ~UseCase_3_Mesh();

    const stk::mesh::BulkData & bulkData() const
    { return m_bulkData; }

    stk::mesh::BulkData       & modifiableBulkData()
    { return m_bulkData; }

  private:

    stk::mesh::BulkData m_bulkData;

};

void populate( UseCase_3_Mesh & mesh );
bool verifyMesh( const UseCase_3_Mesh & mesh );

} //namespace use_cases 
} //namespace mesh 
} //namespace stk 

#endif // Stk_Mesh_Use_Cases_UseCase_3_hpp
