/*--------------------------------------------------------------------*/
/*    Copyright 2011 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
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

double distanceBetweenNodes(double * firstNodeCoords, double * secondNodeCoords, const int dimension)
{
    double sum = 0.0;
    for (int idim = 0; idim < dimension; ++idim) {
        sum += std::pow(firstNodeCoords[idim] - secondNodeCoords[idim], 2);
    }
    return std::sqrt(sum);
}

double SecondShortestEdgeFaceSearchTolerance::compute(const stk::mesh::BulkData & mesh,
                                                      const stk::mesh::FieldBase & coordField,
                                                      const stk::mesh::Entity * faceNodes,
                                                      const unsigned numFaceNodes) const
{
    std::vector<double> edgeLengthVector(numFaceNodes);
    const int dimension = mesh.mesh_meta_data().spatial_dimension();
    stk::mesh::Entity lastNode = faceNodes[numFaceNodes-1];
    double * oldNodePosition = static_cast<double*>(stk::mesh::field_data(coordField, lastNode));

    for (size_t inode = 0; inode < numFaceNodes; ++inode) {
      double * nodePosition(static_cast<double*>(stk::mesh::field_data(coordField, faceNodes[inode])));
      edgeLengthVector[inode] = distanceBetweenNodes(nodePosition, oldNodePosition, dimension);
      oldNodePosition = nodePosition;
    }

    std::sort(edgeLengthVector.begin(), edgeLengthVector.end());
    double tolerance = edgeLengthVector[1] * 0.15;

    return tolerance;
}

}
}
