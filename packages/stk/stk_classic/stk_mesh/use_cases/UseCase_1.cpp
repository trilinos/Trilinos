/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/// doxygen tutorial start #includes
#include <use_cases/UseCase_1.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>
/// end code snippet

//----------------------------------------------------------------------

namespace stk {
namespace mesh {
namespace use_cases {

/// doxygen tutorial start source
enum { field_data_chunk_size = 10 };

UseCase_1_Mesh::UseCase_1_Mesh( stk::ParallelMachine comm )
  : m_metaData( stk::mesh::fem::entity_rank_names(SpatialDim) )
  , m_bulkData(  m_metaData , comm , field_data_chunk_size )
{
  /// Done populating the mesh meta data.
  /// Commit the meta data: this locks out changes,
  /// verifies consistency of potentially complex meta data relationships,
  /// and allows the internal data structures to be optimized
  /// for subsquent use by mesh bulk data.
  m_metaData.commit();
}

UseCase_1_Mesh::~UseCase_1_Mesh()
{ }
/// end code snippet

} //namespace use_cases
} //namespace mesh
} //namespace stk


