// @HEADER
// *****************************************************************************
//                Shards : Shared Discretization Tools
//
// Copyright 2008-2011 NTESS and the Shards contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//#define HAVE_SHARDS_DEBUG

#include <stdexcept>
#include <sstream>
#include <Shards_CellTopology.hpp>
#include <Shards_BasicTopologies.hpp>
#include <iostream>

namespace shards {

typedef CellTopologyData_Subcell Subcell ;

class CellTopologyManagedData : public CellTopologyData
{
public:
  /** \brief  Empty topology */
  CellTopologyManagedData( const std::string & name );

  /** \brief  1D topology */
  CellTopologyManagedData( const std::string & name , const unsigned nodeCount );

  /** \brief  2D topology */
  CellTopologyManagedData(
    const std::string                             & name,
    const unsigned                                  vertexCount,
    const unsigned                                  nodeCount,
    const std::vector< const CellTopologyData * > & edges ,
    const std::vector< unsigned >                 & edge_node_map ,
    const CellTopologyData                        * base );

  /** \brief  3D topology */
  CellTopologyManagedData(
    const std::string                             & name,
    const unsigned                                  vertexCount,
    const unsigned                                  nodeCount,
    const std::vector< const CellTopologyData * > & edges ,
    const std::vector< unsigned >                 & edge_node_map ,
    const std::vector< const CellTopologyData * > & faces ,
    const std::vector< unsigned >                 & face_node_map ,
    const CellTopologyData                        * base );

private:
  CellTopologyManagedData();
  CellTopologyManagedData( const CellTopologyManagedData & );
  CellTopologyManagedData & operator = ( const CellTopologyManagedData & );

private:
  const std::string             m_name ;
  std::vector< Subcell >        m_subcell ;
  std::vector< unsigned >       m_node_map ;
};

CellTopologyManagedData *
createCellTopology(
  const std::string & name );

CellTopologyManagedData *
createCellTopology(
  const std::string & name,
  const unsigned      node_count );

CellTopologyManagedData *
createCellTopology(
  const std::string                             & name,
  const unsigned                                  vertex_count ,
  const unsigned                                  node_count,
  const std::vector< const CellTopologyData * > & edges ,
  const std::vector< unsigned >                 & edge_node_map ,
  const CellTopologyData                        * base );

CellTopologyManagedData *
createCellTopology( const std::string                             & name,
                    const unsigned                                  vertex_count,
                    const unsigned                                  node_count,
                    const std::vector< const CellTopologyData * > & edges ,
                    const std::vector< unsigned >                 & edge_node_map ,
                    const std::vector< const CellTopologyData * > & faces ,
                    const std::vector< unsigned >                 & face_node_map ,
                    const CellTopologyData                        * base);


} // namespace shards
