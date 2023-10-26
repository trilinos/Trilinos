/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_SEARCH_UTIL_OBJECT_COORDINATES_HPP
#define STK_SEARCH_UTIL_OBJECT_COORDINATES_HPP

#include <map>
#include <ostream>                                   // for operator<<, etc
#include <stk_mesh/base/Types.hpp>                   // for EntityRank
#include <stk_util/diag/String.hpp>                  // for operator<<, etc
#include <string>                                    // for allocator, etc
#include <utility>                                   // for pair, make_pair
#include <vector>                                    // for vector

#include "stk_mesh/base/MetaData.hpp"                // for MetaData
#include "stk_mesh/base/BulkData.hpp"                // for BulkData
#include "stk_mesh/base/Entity.hpp"                  // for Entity
#include "stk_mesh/base/FieldBase.hpp"
#include "stk_topology/topology.hpp"                 // for topology, etc
#include "stk_util/util/ReportHandler.hpp"    // for ThrowRequireMsg

namespace stk {
  namespace search {

  /** Convert a single mesh object to a single point.
   * This is used in the case where a mesh object is
   * an element with many nodes.  Then something like
   * the element centroid can be used to define
   * a single coordinate point for it.
   * coor is the output coordinate and is assumed to
   * be large enough to hold a single coordinate.
   * Note the maximum needed in all cases is length 3.
   *
   * entity_coordinates is used to return the nodal coordinates
   * that compute_entity_centroid averages.  The return value is
   * the mesh objects from which the coordinates were obtained.
   * compute_entity_centroid returns a vector of vectors containing
   * the coordinates of the nodes that were used to compute the
   * centroid.
   */
  inline std::vector<stk::mesh::Entity>
  entity_coordinates(stk::mesh::Entity entity, const stk::mesh::FieldBase* stkField,
                     std::vector<std::vector<double> >& coordinates)
  {
    coordinates.clear();
    std::vector<stk::mesh::Entity> meshNodes;
    const unsigned ndim = stkField->mesh_meta_data().spatial_dimension();

    const stk::mesh::BulkData& meshBulk = stkField->get_mesh();

    const stk::mesh::EntityRank fieldType = stkField->entity_rank();
    const stk::mesh::EntityRank objType = meshBulk.entity_rank(entity);
    const bool fieldExistsOnEntity = (fieldType == objType);

    if(objType == stk::topology::NODE_RANK || fieldExistsOnEntity) {
      const double* const coor = (double*)stk::mesh::field_data(*stkField, entity);
      if(ndim) {
        STK_ThrowRequireMsg(coor, " Error: The field does not exist, NULL pointer found.\n"
                                  << " Field Name:" << stkField->name() << "\n"
                                  << " Mesh Object number:" << meshBulk.identifier(entity) << "\n");
        coordinates.emplace_back(coor, coor + ndim);
        meshNodes.push_back(entity);
      }
    }
    else {
      // Loop over node relations in mesh object
      stk::mesh::Entity const* nodes = meshBulk.begin_nodes(entity);
      for(int i = 0, ie = meshBulk.num_nodes(entity); i < ie; ++i) {
        stk::mesh::Entity node = nodes[i];
        double* coor = (double*)stk::mesh::field_data(*stkField, node);

        STK_ThrowRequireMsg(coor, " Error: The field does not exist, NULL pointer found.\n"
                                  << " Field Name:" << stkField->name() << "\n"
                                  << " Mesh Object number:" << meshBulk.identifier(entity) << "\n");

        coordinates.emplace_back(coor, coor + ndim);
        meshNodes.push_back(node);
      }
    }
    return meshNodes;
}

inline void
average_coordinates_for_centroid(const std::vector<std::vector<double> >& coordinates, std::vector<double>& centroid)
{
  const int ndim = coordinates.front().size();
  const int num_nodes = coordinates.size();

  centroid.assign(ndim, 0.0);

  for(int j = 0; j < num_nodes; ++j) {
    for(int i = 0; i < ndim; ++i) {
      centroid[i] += coordinates[j][i];
    }
  }
  if(1 != num_nodes) {
    for(int i = 0; i < ndim; ++i) {
      centroid[i] /= num_nodes;
    }
  }
}

inline std::vector<std::vector<double> >
compute_entity_centroid(stk::mesh::Entity entity, const stk::mesh::FieldBase& nodalCoorRef, std::vector<double>& centroid)
{
  std::vector<std::vector<double> > coordinates;
  entity_coordinates(entity, &nodalCoorRef, coordinates);

  average_coordinates_for_centroid(coordinates, centroid);
  return coordinates;
}

} // end namespace search
} // end namespace stk

#endif

