/*------------------------------------------------------------------------*/
/*                  shards : Shared Discretization Tools                  */
/*                Copyright (2008) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/

#ifndef Shards_CellTopology_hpp
#define Shards_CellTopology_hpp

#include <Shards_CellTopologyData.h>

namespace shards {

/** \addtogroup  shards_package_cell_topology
 *  \{
 */

/** \brief  Generate integer key from topological dimensions
 *  \param dimension  maximum value = 7
 *  \param face_count maximum value = 63
 *  \param edge_count maximum value = 63
 *  \param vertex_count maximum value = 63
 *  \param node_count maximum value = 1023
 *  \return -1 if any count exceeds its limit
 */
inline
int cell_topology_key( const int dimension ,
                       const int face_count ,
                       const int edge_count ,
                       const int vertex_count ,
                       const int node_count )
{
  const bool bad = ( dimension    >> 3 ) ||
                   ( face_count   >> 6 ) ||
                   ( edge_count   >> 6 ) ||
                   ( vertex_count >> 6 ) ||
                   ( node_count   >> 10 );

  const int key = bad ? -1 : ( dimension    << 28  ) |
                             ( face_count   << 22  ) |
                             ( edge_count   << 16  ) |
                             ( vertex_count << 10  ) |
                             ( node_count          ) ;

  return key ;
}

/** \} */
} // namespace shards

#endif // Shards_CellTopology_hpp

