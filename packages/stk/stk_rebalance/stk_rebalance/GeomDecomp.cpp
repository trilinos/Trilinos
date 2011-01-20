/*--------------------------------------------------------------------*/
/*    Copyright 2002, 2010 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

// Copyright 2001, 2002 Sandia Corporation, Albuquerque, NM.

#include <limits>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include <stdexcept>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>

#include <stk_util/environment/ReportHandler.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_rebalance/GeomDecomp.hpp>

using stk::mesh::fem::NODE_RANK;

namespace stk {
namespace rebalance {

std::vector<const mesh::Entity *> GeomDecomp::entity_coordinates(const mesh::Entity                 & obj,
                                                                  const VectorField            & nodal_coor,
                                                                  std::vector<std::vector<double> >  & coordinates)
{
  coordinates.clear();
  std::vector<const mesh::Entity *> mesh_nodes;

  const mesh::EntityRank objtype   = obj.entity_rank();
  if ( objtype == NODE_RANK )
  {
    throw std::runtime_error("GeomDecomp::entity_coordinates Error: Can not be called for nodal objects.");
  } else {

    // Loop over node relations in mesh object
    mesh::PairIterRelation nr   = obj.relations( NODE_RANK );

    for ( ; nr.first != nr.second; ++nr.first )
    {
      const mesh::Relation &rel = *nr.first;
      if (rel.entity_rank() ==  NODE_RANK) { // %fixme: need to check for USES relation
        const mesh::Entity *nobj = rel.entity();
        const unsigned ndim(field_data_size(nodal_coor, *nobj)/sizeof(double)); // TODO - is there a better way to get this info?
        double * coor = mesh::field_data(nodal_coor, *nobj);
        if (!coor) {
          throw std::runtime_error("GeomDecomp::entity_coordinates Error: The coordinate field does not exist.");
        }
        std::vector<double> temp(ndim);
        for ( unsigned i = 0; i < ndim; ++i ) { temp[i] = coor[i]; }
        coordinates.push_back(temp);
        mesh_nodes.push_back(nobj);
      }
    }
  }
  return mesh_nodes;
}

std::vector<std::vector<double> > GeomDecomp::compute_obj_centroid(const mesh::Entity & obj,
                                                                   const VectorField & nodal_coor_ref,
                                                                   std::vector<double>   & centroid)
{
  std::vector<std::vector<double> > coordinates;
  entity_coordinates(obj, nodal_coor_ref, coordinates);

  const int ndim      = coordinates.front().size();
  const int num_nodes = coordinates.size();

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
    assert(0);
  }
  return;
}
}

//: Convert a mesh entity to a single point
//: in cartesian coordinates (x,y,z)
void GeomDecomp::obj_to_point (const mesh::Entity            & obj,
                               const VectorField & nodeCoord,
                               std::vector<double>           & coor)
{
  compute_obj_centroid(obj, nodeCoord, coor);
  apply_rotation (coor);
}
} // namespace rebalance
} // namespace sierra
