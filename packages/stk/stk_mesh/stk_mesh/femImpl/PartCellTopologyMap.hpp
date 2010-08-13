/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_femImpl_PartCellTopologyMap_hpp
#define stk_femImpl_PartCellTopologyMap_hpp

#include <Shards_CellTopologyData.h>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace stk {
namespace mesh {
namespace impl {

class PartCellTopologyMap {
public:

  /** \brief  This method is static to allow free-functions with
   *          just a part, bucket, or entity to query cell topologies,
   *          without having to carry around a reference to the
   *          FiniteElementMesh object.  The part and then the part's
   *          supersets are queried for an associated topology.
   */
  static const CellTopologyData * get_cell_topology( const Part & , const char * const );
  static const CellTopologyData * get_cell_topology( const Bucket & , const char * const );
  
  /** \brief  Query the part uniquely associated with the topology. */
  Part * get_part( const CellTopologyData & , const char * const ) const ;
  
  /** \brief  Query or declare a part to be uniquely associated
   *          with the topology.  Throw an exception if the cell topology
   *          is not appropriate for the spatial dimension.
   */
  Part & declare_part( const CellTopologyData & , unsigned entity_rank );
  
  /** \brief  Declare parts for each predefined, standard cell topology
   *          that is defined within the Shards library and is
   *          appropriate for the spatial dimension.
   */
  PartCellTopologyMap( MetaData & , unsigned spatial_dimension );

  ~PartCellTopologyMap();

private:
  PartVector     m_parts ;
  MetaData     & m_metaData ;
  const unsigned m_spatial_dimension ;

  PartCellTopologyMap();
  PartCellTopologyMap( const PartCellTopologyMap & );
  PartCellTopologyMap & operator = ( const PartCellTopologyMap & );
};

}
}
}

#endif

