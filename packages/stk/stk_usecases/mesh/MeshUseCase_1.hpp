/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef Stk_Mesh_Use_Cases_UseCase_1_hpp
#define Stk_Mesh_Use_Cases_UseCase_1_hpp

/// doxygen tutorial start #includes
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
/// end code snippet

/** stk_mesh Use Case 1
 * This is basically hello world for STK Mesh
 *
 */

namespace stk {
namespace mesh {
namespace use_cases {

class UseCase_1_Mesh
{
public:
  ~UseCase_1_Mesh();

  UseCase_1_Mesh( stk::ParallelMachine comm );

  enum { SpatialDim = 1};

  stk::mesh::MetaData m_metaData;
  stk::mesh::BulkData m_bulkData;
};


} //namespace use_cases
} //namespace mesh
} //namespace stk

#endif // Stk_Mesh_Use_Cases_UseCase_1_hpp

