#include <gtest/gtest.h>
#include <unit_tests/UnitTestUtils.hpp>

namespace
{

TEST(NaluPerformance, BoxSphereIntersections)
{

}

void createBoundingBoxesForElementsInElementBlocks(const int procId, const sierra::Mesh &mesh, const std::vector<double> &coordinates, GtkBoxVector& domainBoxes)
{
    size_t numberBoundingBoxes = mesh.getNumberLocalElements();
    domainBoxes.resize(numberBoundingBoxes);

    sierra::Mesh::BlockIdVector blockIds;
    mesh.fillElementBlockIds(blockIds);

    size_t boxCounter = 0;

    std::vector<double> boxCoordinates(6);
    for (size_t elemBlockNum=0;elemBlockNum<mesh.getNumberElementBlocks();elemBlockNum++)
    {
        sierra::Mesh::LocalNodeIdVector connectivity;
        mesh.fillElementToLocalNodeConnectivityForBlock(blockIds[elemBlockNum], connectivity);
        int numNodesPerElement = mesh.getNumberNodesPerElement(blockIds[elemBlockNum]);
        size_t numElementsThisBlock = mesh.getNumberLocalElementsInBlock(blockIds[elemBlockNum]);
        for (size_t elemCounter=0;elemCounter<numElementsThisBlock;elemCounter++)
        {
            createBoundingBoxForElement(&connectivity[numNodesPerElement*elemCounter], numNodesPerElement, coordinates, boxCoordinates);
            Ident domainBoxId(boxCounter, procId);
            Point min();
            Point max();
            domainBoxes[boxCounter] = std::make_pair(GtkBox(boxCoordinates[0], boxCoordinates[1], boxCoordinates[2],
                                                            boxCoordinates[3], boxCoordinates[4], boxCoordinates[5]),
                                                            domainBoxId);
            boxCounter++;
        }
    }
}

}
