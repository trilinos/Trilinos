// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.



#include <limits>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include <stdexcept>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
//#include <stk_mesh/base/FieldData.hpp>

#include <stk_util/util/ReportHandler.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <percept/stk_rebalance/GeomDecomp.hpp>

static const stk::mesh::EntityRank NODE_RANK = stk::topology::NODE_RANK;

namespace stk {
namespace rebalance {

std::vector< mesh::Entity > GeomDecomp::entity_coordinates(stk::mesh::BulkData& bulk_data, const mesh::Entity                 & entity,
                                                                  const stk::mesh::FieldBase            & nodal_coor,
                                                                  std::vector<std::vector<double> >  & coordinates)
{
  coordinates.clear();
  std::vector< mesh::Entity > mesh_nodes;

  const mesh::EntityRank enttype   = bulk_data.entity_rank(entity);
  if ( enttype == NODE_RANK )
  {
    throw std::runtime_error("GeomDecomp::entity_coordinates Error: Can not be called for nodal entities.");
  } else {

    // Loop over node relations in mesh entities
    const percept::MyPairIterRelation nr   (bulk_data, entity , NODE_RANK);

    for (unsigned inr=0; inr < nr.size(); ++inr)
    {
      const percept::MyPairIterRelation::MyRelation  &rel = nr[inr];
      //if (rel.entity_rank() ==  NODE_RANK) { // %fixme: need to check for USES relation
      if (bulk_data.entity_rank(rel.entity()) ==  NODE_RANK) { // %fixme: need to check for USES relation
        const mesh::Entity nent = rel.entity();
        const unsigned ndim = bulk_data.mesh_meta_data().spatial_dimension();
        double * coor = static_cast<double*>(mesh::field_data(nodal_coor, nent));
        if (!coor) {
          throw std::runtime_error("GeomDecomp::entity_coordinates Error: The coordinate field does not exist.");
        }
        std::vector<double> temp(ndim);
        for ( unsigned i = 0; i < ndim; ++i ) { temp[i] = coor[i]; }
        coordinates.push_back(temp);
        mesh_nodes.push_back(nent);
      }
    }
  }
  return mesh_nodes;
}

std::vector<std::vector<double> > GeomDecomp::compute_entity_centroid(stk::mesh::BulkData& bulk_data, const mesh::Entity & entity,
                                                                   const stk::mesh::FieldBase & nodal_coor_ref,
                                                                   std::vector<double>   & centroid)
{
  std::vector<std::vector<double> > coordinates;
  stk::mesh::EntityRank entity_rank = bulk_data.entity_rank(entity);
  if (entity_rank == stk::topology::ELEMENT_RANK + 1)
    {
      for (stk::mesh::EntityRank irank=stk::topology::NODE_RANK; irank <= stk::topology::ELEMENT_RANK; ++irank)
        {
          const percept::MyPairIterRelation nr(bulk_data, entity , irank);
          for (unsigned ii=0; ii < nr.size(); ++ii)
            {
              stk::mesh::Entity elem = nr[ii].entity();
              std::vector<std::vector<double> > coordinates_1;
              entity_coordinates(bulk_data, elem, nodal_coor_ref, coordinates_1);
              coordinates.insert(coordinates.end(), coordinates_1.begin(), coordinates_1.end());
            }
        }
    }
  else
    {
      entity_coordinates(bulk_data, entity, nodal_coor_ref, coordinates);
    }

  const int num_nodes = coordinates.size();
  const int ndim      = coordinates.front().size();

  centroid.resize(ndim);
  for (int i=0; i<ndim; ++i) { centroid[i] = 0; }
  for ( int j = 0; j < num_nodes; ++j ) {
    for ( int i = 0; i < ndim; ++i ) { centroid[i] += coordinates[j][i]; }
  }
  if (1 != num_nodes) {
    for (int i=0; i<ndim; ++i) { centroid[i] /= num_nodes; }
  }
  return coordinates;
}

namespace {
void apply_rotation (std::vector<double> &coor)
{
  // Apply slight transformation to "disalign" RCB coordinates
  // from the model coordinates.  This causes the RCB axis cuts
  // to be "disaligned" from straight lines of the model.

  static const double tS = 0.0001 ; /* sin( angle / 2 ), angle = 0.012 deg */
  static const double tC = sqrt( (double)( 1.0 - tS * tS ) );
  static const double tQ = tS / sqrt( (double) 3.0 );
  static const double t1 = tC * tC - tQ * tQ ;
  static const double t2 =  2.0 * tQ * ( tC + tQ );
  static const double t3 = -2.0 * tQ * ( tC - tQ );


  std::vector<double> temp(coor);
  const size_t nd = temp.size();

  // Apply minute transformation to the coordinate
  // to rotate the RCB axis slightly away from the model axis.

  if ( nd == 3 ) {
    coor[0] = t1 * temp[0] + t3 * temp[1] + t2 * temp[2] ;
    coor[1] = t2 * temp[0] + t1 * temp[1] + t3 * temp[2] ;
    coor[2] = t3 * temp[0] + t2 * temp[1] + t1 * temp[2] ;
  }
  else if ( nd == 2 ) {
    coor[0] = tC * temp[0] - tS * temp[1] ;
    coor[1] = tS * temp[0] + tC * temp[1] ;
  }
  else if ( nd == 1 ) {
    coor[0] = temp[0] ;
  }
  else {
    STK_ThrowRequireMsg(false, "Spatial Dimention not 1, 2, or 3, can not apply rotation."); // Should never make it here
  }
  return;
}
}

//: Convert a mesh entity to a single point
//: in cartesian coordinates (x,y,z)
  void GeomDecomp::entity_to_point (stk::mesh::BulkData& bulk_data, const mesh::Entity            & entity,
                               const stk::mesh::FieldBase & nodeCoord,
                               std::vector<double>           & coor)
{
  compute_entity_centroid(bulk_data, entity, nodeCoord, coor);
  apply_rotation (coor);
}
} // namespace rebalance
} // namespace sierra
