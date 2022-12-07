// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <cmath>
#include <cstdlib>

#include <percept/xfer/TransferHelper.hpp>

#include <stk_mesh/base/MetaData.hpp>

#include "Shards_CellTopology.hpp"

namespace percept {

double 
parametricDistanceToEntity(const double*                 p,
			   const shards::CellTopology &  cellTopo)
{
  double dist = -1.0;
  double dist2 = -1.0;
  double X,Y,Z;
  const double oneThird = 1./3.;
  const double oneFourth = 1./4.;

  unsigned key = cellTopo.getBaseCellTopologyData() -> key ;

  switch( key ) {

    case shards::Hexahedron<>::key :
      // parametric coords are -/+1
      dist = std::max(std::fabs(p[0]), 
		      std::max(std::fabs(p[1]),
			       std::fabs(p[2])));
      break;
    case shards::Quadrilateral<>::key :
      // parametric coords are -/+1
      dist = std::max(std::fabs(p[0]), 
		      std::fabs(p[1]));
      break;
    case shards::Triangle<>::key :
      X=p[0] - oneThird;
      Y=p[1] - oneThird;
      dist = std::max(std::max(-3*X,-3*Y),3*(X+Y));
      break;
    case shards::Tetrahedron<>::key :
      X=p[0] - oneFourth;
      Y=p[1] - oneFourth;
      Z=p[2] - oneFourth;
      dist  = std::max(std::max(-4*X,-4*Y),std::max(-4*Z,4*(X+Y+Z)));
      break;
    case shards::Wedge<>::key :
      X=p[0] - oneThird;
      Y=p[1] - oneThird;
      Z=p[2];
      dist2 = std::max(std::max(-3*X,-3*Y),3*(X+Y));
      dist  = std::max(dist2,std::fabs(Z));
      break;
    case shards::Pyramid<>::key :
      X=p[0];
      Y=p[1];
      Z=p[2] - oneThird;
      dist2 = (3./2.)*(Z + std::max(std::fabs(X),std::fabs(Y)));
      dist  = std::max(dist2, -3*Z);
      break;
    default:
      std::exit(EXIT_FAILURE);    
  }

  return dist;
}

void
compute_element_centroid(const stk::mesh::FieldBase & coordinates,
                         const stk::mesh::Entity & entity, 
                         double * coords) {

  const unsigned nDim = coordinates.get_mesh().mesh_meta_data().spatial_dimension();

  for (unsigned i=0; i<nDim; i++) {
    coords[i] = 0.0;
  }

  unsigned num_nodes = coordinates.get_mesh().num_connectivity(entity, stk::topology::NODE_RANK);
  const stk::mesh::Entity *elem_nodes = coordinates.get_mesh().begin(entity, stk::topology::NODE_RANK);

  for (unsigned inode=0; inode < num_nodes; inode++) {
    stk::mesh::Entity node = elem_nodes[inode];

    double * const node_data = static_cast<double*>(stk::mesh::field_data(coordinates, node));
    for (unsigned i=0; i<nDim; i++) {
      coords[i] += node_data[i] / (double) num_nodes;
    }
  }
}

void
compute_nodal_coords(const stk::mesh::FieldBase & coordinates,
                     const stk::mesh::Entity & entity, 
                     double * coords) {
  stk::mesh::MetaData& meta = coordinates.get_mesh().mesh_meta_data();
  stk::mesh::FieldBase *ucf = meta.get_field(stk::topology::NODE_RANK, "unprojected_coordinates");
  double * coord = static_cast<double*>(stk::mesh::field_data(coordinates, entity));
  if (ucf) coord = (double *)stk::mesh::field_data(*ucf, entity);
  const unsigned nDim = coordinates.get_mesh().mesh_meta_data().spatial_dimension();
  for (unsigned i=0; i<nDim; i++) {
    coords[i] = coord[i];
  }
}

}
