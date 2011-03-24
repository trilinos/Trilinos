/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_percept_HeterogeneousFixture_hpp
#define stk_percept_HeterogeneousFixture_hpp

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/Stencils.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/DefaultFEM.hpp>

#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/fem/CoordinateSystems.hpp>

#include <stk_percept/fixtures/HeterogeneousFixture.hpp>

/** stk_mesh Use Case 3 - copied and modified here */

#define HET_FIX_INCLUDE_EXTRA_ELEM_TYPES 0

namespace stk {
  namespace percept {

    using namespace stk::mesh;

    typedef Field<double,Cartesian>    VectorFieldType ;
    typedef Field<double>              ScalarFieldType ;
    typedef Field<double*,ElementNode> ElementNodePointerFieldType ;

    /** Use case with mixed element topologies and
     *  field relations to provide fast access to node field data
     *  from an element.
     *
     *  copied from stk_mesh and modified
     */

    class HeterogeneousFixture {
    public:


      ~HeterogeneousFixture();

      HeterogeneousFixture( stk::ParallelMachine comm, bool doCommit = true);

      void populate();

      const int m_spatial_dimension;
      mesh::fem::FEMMetaData m_metaData;
      BulkData m_bulkData;

      Part & m_block_hex;
      Part & m_block_wedge;
      Part & m_block_tet;
#if HET_FIX_INCLUDE_EXTRA_ELEM_TYPES
      Part & m_block_pyramid;
      Part & m_block_quad_shell;
      Part & m_block_tri_shell;
#endif
      const EntityRank m_elem_rank;

      VectorFieldType & m_coordinates_field;
      VectorFieldType & m_centroid_field;
      ScalarFieldType & m_temperature_field;
      ScalarFieldType & m_volume_field;
      ElementNodePointerFieldType & m_element_node_coordinates_field;
    };

    bool verifyMesh( const HeterogeneousFixture & mesh );

  } //namespace percept
} //namespace stk

#endif // Stk_Mesh_Use_Cases_UseCase_3_hpp
