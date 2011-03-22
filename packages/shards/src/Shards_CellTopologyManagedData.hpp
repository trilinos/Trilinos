/*
//@HEADER
// ************************************************************************
//
//                Shards : Shared Discretization Tools
//                 Copyright 2008, 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Carter Edwards (hcedwar@sandia.gov),
//                    Pavel Bochev (pbboche@sandia.gov), or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
//@HEADER
*/

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
