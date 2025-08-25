/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "SecondShortestEdgeFaceSearchTolerance.hpp"
#include "stk_mesh/base/Field.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include <cmath>
#include <iterator>

namespace stk {
namespace balance {

template <typename EntityValuesType1, typename EntityValuesType2>
double distanceBetweenNodes(EntityValuesType1 firstNodeCoords, EntityValuesType2 secondNodeCoords, const int dimension)
{
  double sum = 0.0;
  for (stk::mesh::ComponentIdx component=0_comp; component < dimension; ++component) {
    sum += std::pow(firstNodeCoords(component) - secondNodeCoords(component), 2);
  }
  return std::sqrt(sum);
}

double SecondShortestEdgeFaceSearchTolerance::compute(const stk::mesh::BulkData& mesh,
                                                      const stk::mesh::FieldBase& coordField,
                                                      const stk::mesh::Entity* faceNodes,
                                                      const unsigned numFaceNodes) const
{
  std::vector<double> edgeLengthVector(numFaceNodes);
  const int dimension = mesh.mesh_meta_data().spatial_dimension();
  stk::mesh::Entity lastNode = faceNodes[numFaceNodes-1];

  stk::mesh::field_data_execute<double, stk::mesh::ReadOnly>(coordField,
    [&](auto& coordFieldData) {
      auto oldNodePosition = coordFieldData.entity_values(lastNode);

      for (size_t inode = 0; inode < numFaceNodes; ++inode) {
        auto nodePosition = coordFieldData.entity_values(faceNodes[inode]);
        edgeLengthVector[inode] = distanceBetweenNodes(nodePosition, oldNodePosition, dimension);
        oldNodePosition = nodePosition;
      }
    }
  );

  std::sort(edgeLengthVector.begin(), edgeLengthVector.end());
  double tolerance = edgeLengthVector[1] * m_tolerance;

  return tolerance;
}

}
}
