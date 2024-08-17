#ifndef EXAMPLEMESHES_H_
#define EXAMPLEMESHES_H_

#include "stk_util/stk_config.h"
#include <generated/Iogn_DashSurfaceMesh.h>
#include <vector>

namespace unitTestUtils
{

namespace exampleMeshes
{

inline void fillDataForUnitCube(std::vector<double> &coordinates)
{
    // Nodes in Exodus Ordering for hex elements which is same as Flanagan/Belytschko paper
    double *x = coordinates.data();
    double *y = coordinates.data() + 8;
    double *z = coordinates.data() + 16;
    x[4] = 1; x[5] = 1; x[6] = 1; x[7] = 1;
    y[1] = 1; y[2] = 1; y[5] = 1; y[6] = 1;
    z[2] = 1; z[3] = 1; z[6] = 1; z[7] = 1;
}

inline void fillDataForRectangloid(std::vector<double> &coordinates)
{
    // Nodes in Exodus Ordering for hex elements which is same as Flanagan/Belytschko paper
    double *x = coordinates.data();
    double *y = coordinates.data() + 8;
    double *z = coordinates.data() + 16;
    x[4] = 3; x[5] = 3; x[6] = 3; x[7] = 3;
    y[1] = 2; y[2] = 2; y[5] = 2; y[6] = 2;
    z[2] = 0.5; z[3] = 0.5; z[6] = 0.5; z[7] = 0.5;
}

// node numbering in code below
//  ^ y-dir
//  |
//  | 4    3
//  |
//  | 1    2
//  o--------> x-dir
inline
Iogn::ExodusData createExodusDataForDisconnectedHex8s(int numberOfHexes)
{
    int numberOfElementBlocks = numberOfHexes;
    int numberOfNodesInEachElementBlock = 8;
    int globalNumberOfNodes = numberOfElementBlocks * numberOfNodesInEachElementBlock;
    std::vector<double> coordinates(globalNumberOfNodes*3);
    std::vector< std::vector<int> > elementBlockConnectivity(numberOfElementBlocks);
    for(int blockIndex=0; blockIndex < numberOfElementBlocks; blockIndex++)
    {
        elementBlockConnectivity[blockIndex].resize(numberOfNodesInEachElementBlock);
    }
    int numberOfElementsInEachBlock = 1;
    std::vector<int> globalNumberOfElementsInBlock(numberOfElementBlocks, numberOfElementsInEachBlock);
    std::vector<int> localNumberOfElementsInBlock(numberOfElementBlocks, numberOfElementsInEachBlock);
    std::vector<Iogn::Topology> blockTopologicalData(numberOfElementBlocks, Iogn::Hex8);
    int numberOfTotalElements = numberOfElementsInEachBlock * numberOfElementBlocks;
    std::vector<int> globalIdsOfLocalElements(numberOfTotalElements);
    std::vector<int> globalIdsOfLocalNodes(globalNumberOfNodes);
    for(int i=0; i < numberOfTotalElements; i++)
    {
        globalIdsOfLocalElements[i] = i+1;
    }
    for(int i=0; i < globalNumberOfNodes; i++)
    {
        globalIdsOfLocalNodes[i] = i+1;
    }
    for(int blockIndex=0; blockIndex < numberOfElementBlocks; blockIndex++)
    {
        int elementBlockOffset = numberOfNodesInEachElementBlock * blockIndex;
        for(int i=0; i < numberOfNodesInEachElementBlock; i++)
        {
            elementBlockConnectivity[blockIndex][i] = elementBlockOffset+i+1;
        }
    }

    double coords[] = {
        0, 0, 0,
        1, 0, 0,
        1, 1, 0,
        0, 1, 0,
        0, 0, 1,
        1, 0, 1,
        1, 1, 1,
        0, 1, 1,
    };

    coordinates.clear();
    coordinates.resize(globalNumberOfNodes*3);
    for (int i=0;i<numberOfElementBlocks;i++)
    {
        int offset = 3*numberOfNodesInEachElementBlock*i;
        for (int j=0;j<numberOfNodesInEachElementBlock;j++)
        {
            coordinates[offset+3*j + 0] = coords[3*j+0]+i;
            coordinates[offset+3*j + 1] = coords[3*j+1];
            coordinates[offset+3*j + 2] = coords[3*j+2];
        }
    }

    return Iogn::ExodusData(coordinates, elementBlockConnectivity, globalNumberOfElementsInBlock, localNumberOfElementsInBlock,
            blockTopologicalData, globalNumberOfNodes, globalIdsOfLocalElements, globalIdsOfLocalNodes);
}

namespace simple_fields {

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void fillDataForUnitCube(std::vector<double> &coordinates);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void fillDataForRectangloid(std::vector<double> &coordinates);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
Iogn::ExodusData createExodusDataForDisconnectedHex8s(int numberOfHexes);

} // namespace simple_fields

}
}
#endif
