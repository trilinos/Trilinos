// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_SETUPREGIONUTILITIES_HPP
#define MUELU_SETUPREGIONUTILITIES_HPP

#include "Teuchos_Array.hpp"

template <class LocalOrdinal>
void findInterface(const int numDimensions, Teuchos::Array<LocalOrdinal> nodesPerDim,
                   const Teuchos::Array<int> boundaryConditions,
                   Teuchos::Array<LocalOrdinal>& interfacesDimensions,
                   Teuchos::Array<LocalOrdinal>& interfacesLIDs) {
  using LO = LocalOrdinal;

  std::cout << "nodesPerDim: " << nodesPerDim << std::endl;

  // The goal of this function is to build a list of LIDs on the skin
  // of the mesh, broken down by logical faces i.e. in 2D four edges,
  // in 3D six faces, while ordering the list lexicographically.

  // Step 1: determine what edges/faces are needed based on BC.
  LO numInterfaces = 0, numTotalLIDs = 0, numInterfaceLIDs;
  interfacesDimensions.resize(18, 1);
  for (LO dim = 0; dim < 3; ++dim) {
    // Check for nodes and no boundary conditions in this direction, otherwise skip it.
    if ((nodesPerDim[dim] == 1) || (boundaryConditions[2 * dim] + boundaryConditions[2 * dim + 1] == 2)) {
      continue;
    }

    // Since we are not skipping this direction
    // we at least need to store data for one surface
    numInterfaceLIDs = 1;
    for (LO dimIdx = 0; dimIdx < 3; ++dimIdx) {
      interfacesDimensions[numInterfaces * 3 + dimIdx] = (dimIdx == dim ? 1 : nodesPerDim[dimIdx]);
      numInterfaceLIDs *= (dimIdx == dim ? 1 : nodesPerDim[dimIdx]);
    }
    numTotalLIDs += numInterfaceLIDs;
    ++numInterfaces;

    // If there are no BC we need to store the surface twice.
    if (boundaryConditions[2 * dim] + boundaryConditions[2 * dim + 1] == 0) {
      for (LO dimIdx = 0; dimIdx < 3; ++dimIdx) {
        interfacesDimensions[numInterfaces * 3 + dimIdx] = interfacesDimensions[(numInterfaces - 1) * 3 + dimIdx];
      }
      numTotalLIDs += numInterfaceLIDs;
      ++numInterfaces;
    }
  }
  interfacesDimensions.resize(3 * numInterfaces);
  interfacesLIDs.resize(numTotalLIDs, -1);

  // Step 2 lazy implementation of all geometrical cases.

  LO nodeOffset = 0;
  if (numDimensions == 2) {
    // left interface
    if (boundaryConditions[0] == 0) {
      for (LO nodeIdx = 0; nodeIdx < nodesPerDim[1]; ++nodeIdx) {
        interfacesLIDs[nodeOffset + nodeIdx] = nodeIdx * nodesPerDim[0];
      }
      nodeOffset += nodesPerDim[1];
    }

    // right interface
    if (boundaryConditions[1] == 0) {
      for (LO nodeIdx = 0; nodeIdx < nodesPerDim[1]; ++nodeIdx) {
        interfacesLIDs[nodeOffset + nodeIdx] = (nodeIdx + 1) * nodesPerDim[0] - 1;
      }
      nodeOffset += nodesPerDim[1];
    }

    // front interface
    if (boundaryConditions[2] == 0) {
      for (LO nodeIdx = 0; nodeIdx < nodesPerDim[0]; ++nodeIdx) {
        interfacesLIDs[nodeOffset + nodeIdx] = nodeIdx;
      }
      nodeOffset += nodesPerDim[0];
    }

    // back interface
    if (boundaryConditions[3] == 0) {
      for (LO nodeIdx = 0; nodeIdx < nodesPerDim[0]; ++nodeIdx) {
        interfacesLIDs[nodeOffset + nodeIdx] = (nodesPerDim[1] - 1) * nodesPerDim[0] + nodeIdx;
      }
      nodeOffset += nodesPerDim[0];
    }
  }

  if (numDimensions == 3) {
    // left interface
    if (boundaryConditions[0] == 0) {
      for (LO k = 0; k < nodesPerDim[2]; ++k) {
        for (LO j = 0; j < nodesPerDim[1]; ++j) {
          interfacesLIDs[nodeOffset + k * nodesPerDim[1] + j] = k * nodesPerDim[1] * nodesPerDim[0] + j * nodesPerDim[0];
        }
      }
      nodeOffset += nodesPerDim[2] * nodesPerDim[1];
    }

    // right interface
    if (boundaryConditions[1] == 0) {
      for (LO k = 0; k < nodesPerDim[2]; ++k) {
        for (LO j = 0; j < nodesPerDim[1]; ++j) {
          interfacesLIDs[nodeOffset + k * nodesPerDim[1] + j] = k * nodesPerDim[1] * nodesPerDim[0] + (j + 1) * nodesPerDim[0] - 1;
        }
      }
      nodeOffset += nodesPerDim[2] * nodesPerDim[1];
    }

    // front interface
    if (boundaryConditions[2] == 0) {
      for (LO k = 0; k < nodesPerDim[2]; ++k) {
        for (LO i = 0; i < nodesPerDim[0]; ++i) {
          interfacesLIDs[nodeOffset + k * nodesPerDim[0] + i] = k * nodesPerDim[1] * nodesPerDim[0] + i;
        }
      }
      nodeOffset += nodesPerDim[2] * nodesPerDim[0];
    }

    // back interface
    if (boundaryConditions[3] == 0) {
      for (LO k = 0; k < nodesPerDim[2]; ++k) {
        for (LO i = 0; i < nodesPerDim[0]; ++i) {
          interfacesLIDs[nodeOffset + k * nodesPerDim[0] + i] = k * nodesPerDim[1] * nodesPerDim[0] + (nodesPerDim[1] - 1) * nodesPerDim[0] + i;
        }
      }
      nodeOffset += nodesPerDim[2] * nodesPerDim[0];
    }

    // bottom interface
    if (boundaryConditions[4] == 0) {
      for (LO j = 0; j < nodesPerDim[1]; ++j) {
        for (LO i = 0; i < nodesPerDim[0]; ++i) {
          interfacesLIDs[nodeOffset + j * nodesPerDim[0] + i] = j * nodesPerDim[0] + i;
        }
      }
      nodeOffset += nodesPerDim[1] * nodesPerDim[0];
    }

    // top interface
    if (boundaryConditions[5] == 0) {
      for (LO j = 0; j < nodesPerDim[1]; ++j) {
        for (LO i = 0; i < nodesPerDim[0]; ++i) {
          interfacesLIDs[nodeOffset + j * nodesPerDim[0] + i] = (nodesPerDim[2] - 1) * nodesPerDim[1] * nodesPerDim[0] + j * nodesPerDim[0] + i;
        }
      }
      nodeOffset += nodesPerDim[1] * nodesPerDim[0];
    }
  }
}  // findInterface

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void createRegionData(const int numDimensions,
                      const bool useUnstructured, const int numDofsPerNode,
                      const Teuchos::ArrayView<GlobalOrdinal> gNodesPerDim,
                      const Teuchos::ArrayView<LocalOrdinal> lNodesPerDim,
                      const Teuchos::ArrayView<GlobalOrdinal> procsPerDim,
                      const Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > nodeMap,
                      const Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > dofMap,
                      int& maxRegPerGID, LocalOrdinal& numLocalRegionNodes,
                      Teuchos::Array<int>& boundaryConditions,
                      Teuchos::Array<GlobalOrdinal>& sendGIDs,  ///< GIDs of nodes
                      Teuchos::Array<int>& sendPIDs,
                      int& numInterfaces,
                      Teuchos::Array<LocalOrdinal>& rNodesPerDim,
                      Teuchos::Array<GlobalOrdinal>& quasiRegionGIDs,
                      Teuchos::Array<GlobalOrdinal>& quasiRegionCoordGIDs,
                      Teuchos::Array<LocalOrdinal>& compositeToRegionLIDs,
                      Teuchos::Array<GlobalOrdinal>& interfaceGIDs,
                      Teuchos::Array<LocalOrdinal>& interfaceLIDsData) {
  using GO = GlobalOrdinal;
  using LO = LocalOrdinal;

  const int myRank = nodeMap->getComm()->getRank();
  Teuchos::Array<GO> startIndices(3);
  Teuchos::Array<GO> endIndices(3);
  const GO startGID = dofMap->getMinGlobalIndex() / numDofsPerNode;
  {
    startIndices[2] = startGID / (gNodesPerDim[1] * gNodesPerDim[0]);
    const GO rem    = startGID % (gNodesPerDim[1] * gNodesPerDim[0]);
    startIndices[1] = rem / gNodesPerDim[0];
    startIndices[0] = rem % gNodesPerDim[0];
    endIndices[0]   = startIndices[0] + lNodesPerDim[0] - 1;
    endIndices[1]   = startIndices[1] + lNodesPerDim[1] - 1;
    endIndices[2]   = startIndices[2] + lNodesPerDim[2] - 1;
  }

  int leftBC = 0, rightBC = 0, frontBC = 0, backBC = 0, bottomBC = 0, topBC = 0;
  if (startIndices[0] == 0) {
    leftBC = 1;
  }
  if (startIndices[1] == 0) {
    frontBC = 1;
  }
  if (startIndices[2] == 0) {
    bottomBC = 1;
  }

  if (endIndices[0] == gNodesPerDim[0] - 1) {
    rightBC = 1;
  }
  if (endIndices[1] == gNodesPerDim[1] - 1) {
    backBC = 1;
  }
  if (endIndices[2] == gNodesPerDim[2] - 1) {
    topBC = 1;
  }

  boundaryConditions.resize(6);
  boundaryConditions[0] = leftBC;
  boundaryConditions[1] = rightBC;
  boundaryConditions[2] = frontBC;
  boundaryConditions[3] = backBC;
  boundaryConditions[4] = bottomBC;
  boundaryConditions[5] = topBC;

  LO numReceive = 0, numSend = 0;
  Teuchos::Array<GO> receiveGIDs;
  Teuchos::Array<int> receivePIDs;
  Teuchos::Array<LO> receiveLIDs, sendLIDs, interfaceLIDs;

  if (numDimensions == 1) {
    maxRegPerGID = 2;
    if (leftBC == 0) {
      numReceive = 1;
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);
      receiveLIDs.resize(numReceive);

      receiveGIDs[0] = startIndices[0] - 1;
      receivePIDs[0] = myRank - 1;
    }
    if (rightBC == 0) {
      numSend = 1;
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);
      sendLIDs.resize(numSend);

      sendGIDs[0] = endIndices[0];
      sendGIDs[0] = myRank + 1;
      sendLIDs[0] = lNodesPerDim[0] - 1;
    }
  } else if (numDimensions == 2) {
    maxRegPerGID = 4;
    // Received nodes
    if (frontBC == 0 && leftBC == 0) {
      numReceive = lNodesPerDim[0] + lNodesPerDim[1] + 1;
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);
      receiveLIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive front-left corner node
      receiveGIDs[countIDs] = startGID - gNodesPerDim[0] - 1;
      receivePIDs[countIDs] = myRank - procsPerDim[0] - 1;
      ++countIDs;
      // Receive front edge nodes
      for (LO i = 0; i < lNodesPerDim[0]; ++i) {
        receiveGIDs[countIDs] = startGID - gNodesPerDim[0] + i;
        receivePIDs[countIDs] = myRank - procsPerDim[0];
        ++countIDs;
      }
      // Receive left edge nodes
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        receiveGIDs[countIDs] = startGID - 1 + j * gNodesPerDim[0];
        receivePIDs[countIDs] = myRank - 1;
        ++countIDs;
      }

      if (useUnstructured) {  // Set parameters for interface aggregation
        numInterfaces = 2;

        interfaceLIDs.resize((lNodesPerDim[0] + 1) + (lNodesPerDim[1] + 1));
        for (LO nodeIdx = 0; nodeIdx < lNodesPerDim[0] + 1; ++nodeIdx) {
          interfaceLIDs[nodeIdx] = nodeIdx;
        }
        for (LO nodeIdx = 0; nodeIdx < lNodesPerDim[1] + 1; ++nodeIdx) {
          interfaceLIDs[lNodesPerDim[0] + 1 + nodeIdx] = nodeIdx * (lNodesPerDim[1] + 1);
        }
      }

    } else if (frontBC == 0) {
      numReceive = lNodesPerDim[0];
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);
      receiveLIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive front edge nodes
      for (LO i = 0; i < lNodesPerDim[0]; ++i) {
        receiveGIDs[countIDs] = startGID - gNodesPerDim[0] + i;
        receivePIDs[countIDs] = myRank - procsPerDim[0];
        ++countIDs;
      }

      if (useUnstructured) {  // Set parameters for interface aggregation
        numInterfaces = 1;

        interfaceLIDs.resize(lNodesPerDim[0] + 1);
        for (LO nodeIdx = 0; nodeIdx < lNodesPerDim[0] + 1; ++nodeIdx) {
          interfaceLIDs[nodeIdx] = nodeIdx;
        }
      }

    } else if (leftBC == 0) {
      numReceive = lNodesPerDim[1];
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);
      receiveLIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive left edge nodes
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        receiveGIDs[countIDs] = startGID - 1 + j * gNodesPerDim[0];
        receivePIDs[countIDs] = myRank - 1;
        ++countIDs;
      }

      if (useUnstructured) {  // Set parameters for interface aggregation
        numInterfaces = 1;

        interfaceLIDs.resize(lNodesPerDim[1] + 1);
        for (LO nodeIdx = 0; nodeIdx < lNodesPerDim[1] + 1; ++nodeIdx) {
          interfaceLIDs[nodeIdx] = nodeIdx * (lNodesPerDim[1] + 1);
        }
      }
    }

    // Sent nodes
    if (rightBC == 0 && backBC == 0) {
      numSend = lNodesPerDim[0] + lNodesPerDim[1] + 1;
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);
      sendLIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of right edge
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        sendGIDs[countIDs] = j * gNodesPerDim[0] + startGID + lNodesPerDim[0] - 1;
        sendPIDs[countIDs] = myRank + 1;
        sendLIDs[countIDs] = (j + 1) * lNodesPerDim[0] - 1;
        ++countIDs;
      }
      // Send nodes of back edge
      for (LO i = 0; i < lNodesPerDim[0]; ++i) {
        sendGIDs[countIDs] = i + startGID + (lNodesPerDim[1] - 1) * gNodesPerDim[0];
        sendPIDs[countIDs] = myRank + procsPerDim[0];
        sendLIDs[countIDs] = (lNodesPerDim[1] - 1) * lNodesPerDim[0] + i;
        ++countIDs;
      }
      // Send node of back-right corner
      sendGIDs[countIDs] = startGID + (lNodesPerDim[1] - 1) * gNodesPerDim[0] + lNodesPerDim[0] - 1;
      sendPIDs[countIDs] = myRank + procsPerDim[0] + 1;
      sendLIDs[countIDs] = lNodesPerDim[1] * lNodesPerDim[0] - 1;
      ++countIDs;
    } else if (backBC == 0) {
      numSend = lNodesPerDim[0];

      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);
      sendLIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of back edge
      for (LO i = 0; i < lNodesPerDim[0]; ++i) {
        sendGIDs[countIDs] = i + startGID + (lNodesPerDim[1] - 1) * gNodesPerDim[0];
        sendPIDs[countIDs] = myRank + procsPerDim[0];
        sendLIDs[countIDs] = (lNodesPerDim[1] - 1) * lNodesPerDim[0] + i;
        ++countIDs;
      }
    } else if (rightBC == 0) {
      numSend = lNodesPerDim[1];
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);
      sendLIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of right edge
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        sendGIDs[countIDs] = j * gNodesPerDim[0] + startGID + lNodesPerDim[0] - 1;
        sendPIDs[countIDs] = myRank + 1;
        sendLIDs[countIDs] = (j + 1) * lNodesPerDim[0] - 1;
        ++countIDs;
      }
    }
  } else if (numDimensions == 3) {
    maxRegPerGID = 8;
    // Received nodes
    if ((bottomBC == 0) && (frontBC == 0) && (leftBC == 0)) {
      numReceive = lNodesPerDim[0] * lNodesPerDim[1]                 // bottom face
                   + lNodesPerDim[0] * (lNodesPerDim[2] + 1)         // front face
                   + (lNodesPerDim[1] + 1) * (lNodesPerDim[2] + 1);  // left face
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);
      receiveLIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive front-left-bottom corner node
      receiveGIDs[countIDs] = startGID - gNodesPerDim[0] - 1 - gNodesPerDim[1] * gNodesPerDim[0];
      receivePIDs[countIDs] = myRank - procsPerDim[0] - 1 - procsPerDim[1] * procsPerDim[0];
      ++countIDs;

      // Receive front-bottom edge nodes
      for (LO i = 0; i < lNodesPerDim[0]; ++i) {
        receiveGIDs[countIDs] = startGID - gNodesPerDim[0] * gNodesPerDim[1] - gNodesPerDim[0] + i;
        receivePIDs[countIDs] = myRank - procsPerDim[0] * procsPerDim[1] - procsPerDim[0];
        ++countIDs;
      }

      // Recieve left-bottom edge nodes
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        receiveGIDs[countIDs] = startGID - gNodesPerDim[0] * gNodesPerDim[1] - 1 + j * gNodesPerDim[0];
        receivePIDs[countIDs] = myRank - procsPerDim[0] * procsPerDim[1] - 1;
        ++countIDs;
        // Recieve bottom face nodes
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          receiveGIDs[countIDs] = startGID - gNodesPerDim[0] * gNodesPerDim[1] + i + j * gNodesPerDim[0];
          receivePIDs[countIDs] = myRank - procsPerDim[0] * procsPerDim[1];
          ++countIDs;
        }
      }

      // Receive front-left edge nodes
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        receiveGIDs[countIDs] = startGID - gNodesPerDim[0] - 1 + k * gNodesPerDim[0] * gNodesPerDim[1];
        receivePIDs[countIDs] = myRank - procsPerDim[0] - 1;
        ++countIDs;
        // Receive front face nodes
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          receiveGIDs[countIDs] = startGID - gNodesPerDim[0] + i + k * (gNodesPerDim[1] * gNodesPerDim[0]);
          receivePIDs[countIDs] = myRank - procsPerDim[0];
          ++countIDs;
        }
        // Receive left face nodes
        for (LO j = 0; j < lNodesPerDim[1]; ++j) {
          receiveGIDs[countIDs] = startGID - 1 + j * gNodesPerDim[0] + k * (gNodesPerDim[1] * gNodesPerDim[0]);
          receivePIDs[countIDs] = myRank - 1;
          ++countIDs;
        }
      }

      if (useUnstructured) {  // Set parameters for interface aggregation
        numInterfaces = 3;

        interfaceLIDs.resize((lNodesPerDim[0] + 1) * (lNodesPerDim[1] + 1) + (lNodesPerDim[0] + 1) * (lNodesPerDim[2] + 1) + (lNodesPerDim[1] + 1) * (lNodesPerDim[2] + 1));

        LO nodeOffset = 0, nodeIdx, nodeLID;
        // Bottom face
        for (nodeIdx = 0; nodeIdx < (lNodesPerDim[0] + 1) * (lNodesPerDim[1] + 1); ++nodeIdx) {
          interfaceLIDs[nodeIdx] = nodeIdx;
        }
        // Front face
        nodeOffset += (lNodesPerDim[0] + 1) * (lNodesPerDim[1] + 1);
        for (LO k = 0; k < lNodesPerDim[2] + 1; ++k) {
          for (LO i = 0; i < lNodesPerDim[0] + 1; ++i) {
            nodeIdx                = k * (lNodesPerDim[0] + 1) + i + nodeOffset;
            nodeLID                = k * (lNodesPerDim[0] + 1) * (lNodesPerDim[1] + 1) + i;
            interfaceLIDs[nodeIdx] = nodeLID;
          }
        }
        // Left face
        nodeOffset += (lNodesPerDim[0] + 1) * (lNodesPerDim[2] + 1);
        for (LO k = 0; k < lNodesPerDim[2] + 1; ++k) {
          for (LO j = 0; j < lNodesPerDim[1] + 1; ++j) {
            nodeIdx                = k * (lNodesPerDim[1] + 1) + j + nodeOffset;
            nodeLID                = k * (lNodesPerDim[0] + 1) * (lNodesPerDim[1] + 1) + j * (lNodesPerDim[0] + 1);
            interfaceLIDs[nodeIdx] = nodeLID;
          }
        }
      }

      // Two faces received
    } else if ((bottomBC == 0) && (frontBC == 0)) {
      numReceive = lNodesPerDim[0] * (lNodesPerDim[1] + lNodesPerDim[2] + 1);
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);
      receiveLIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive front-bottom edge nodes
      for (LO i = 0; i < lNodesPerDim[0]; ++i) {
        receiveGIDs[countIDs] = startGID - gNodesPerDim[0] * gNodesPerDim[1] - gNodesPerDim[0] + i;
        receivePIDs[countIDs] = myRank - procsPerDim[0] * procsPerDim[1] - procsPerDim[0];
        ++countIDs;
      }
      // Receive bottom face nodes
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          receiveGIDs[countIDs] = startGID - gNodesPerDim[0] * gNodesPerDim[1] + i + j * gNodesPerDim[0];
          receivePIDs[countIDs] = myRank - procsPerDim[0] * procsPerDim[1];
          ++countIDs;
        }
      }
      // Receive front face nodes
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          receiveGIDs[countIDs] = startGID - gNodesPerDim[0] + i + k * (gNodesPerDim[1] * gNodesPerDim[0]);
          receivePIDs[countIDs] = myRank - procsPerDim[0];
          ++countIDs;
        }
      }

      if (useUnstructured) {  // Set parameters for interface aggregation
        numInterfaces = 2;

        interfaceLIDs.resize(lNodesPerDim[0] * (lNodesPerDim[1] + 1) + lNodesPerDim[0] * (lNodesPerDim[2] + 1));

        LO nodeOffset = 0, nodeIdx, nodeLID;
        // Bottom face
        for (nodeIdx = 0; nodeIdx < lNodesPerDim[0] * (lNodesPerDim[1] + 1); ++nodeIdx) {
          interfaceLIDs[nodeIdx] = nodeIdx;
        }
        // Front face
        nodeOffset += lNodesPerDim[0] * (lNodesPerDim[1] + 1);
        for (LO k = 0; k < lNodesPerDim[2] + 1; ++k) {
          for (LO i = 0; i < lNodesPerDim[0]; ++i) {
            nodeIdx                = k * lNodesPerDim[0] + i + nodeOffset;
            nodeLID                = k * lNodesPerDim[0] * (lNodesPerDim[1] + 1) + i;
            interfaceLIDs[nodeIdx] = nodeLID;
          }
        }
      }

    } else if ((bottomBC == 0) && (leftBC == 0)) {
      numReceive = lNodesPerDim[1] * (lNodesPerDim[0] + lNodesPerDim[2] + 1);
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);
      receiveLIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive left-bottom edge nodes
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        receiveGIDs[countIDs] = j * gNodesPerDim[0] + startGID - gNodesPerDim[1] * gNodesPerDim[0] - 1;
        receivePIDs[countIDs] = myRank - procsPerDim[1] * procsPerDim[0] - 1;
        ++countIDs;
        // Receive bottom face nodes
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          receiveGIDs[countIDs] = j * gNodesPerDim[0] + i + startGID - gNodesPerDim[1] * gNodesPerDim[0];
          receivePIDs[countIDs] = myRank - procsPerDim[1] * procsPerDim[0];
          ++countIDs;
        }
      }
      // Receive left face nodes
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO j = 0; j < lNodesPerDim[1]; ++j) {
          receiveGIDs[countIDs] = k * gNodesPerDim[1] * gNodesPerDim[0] + j * gNodesPerDim[0] + startGID - 1;
          receivePIDs[countIDs] = myRank - 1;
          ++countIDs;
        }
      }

      if (useUnstructured) {  // Set parameters for interface aggregation
        numInterfaces = 2;

        interfaceLIDs.resize((lNodesPerDim[0] + 1) * lNodesPerDim[1] + lNodesPerDim[1] * (lNodesPerDim[2] + 1));

        LO nodeOffset = 0, nodeIdx, nodeLID;
        // Bottom face
        for (nodeIdx = 0; nodeIdx < (lNodesPerDim[0] + 1) * lNodesPerDim[1]; ++nodeIdx) {
          interfaceLIDs[nodeIdx] = nodeIdx;
        }
        // Left face
        nodeOffset += (lNodesPerDim[0] + 1) * lNodesPerDim[1];
        for (LO k = 0; k < lNodesPerDim[2] + 1; ++k) {
          for (LO j = 0; j < lNodesPerDim[1]; ++j) {
            nodeIdx                = k * lNodesPerDim[1] + j + nodeOffset;
            nodeLID                = k * (lNodesPerDim[0] + 1) * lNodesPerDim[1] + j * (lNodesPerDim[0] + 1);
            interfaceLIDs[nodeIdx] = nodeLID;
          }
        }
      }

    } else if ((frontBC == 0) && (leftBC == 0)) {
      numReceive = lNodesPerDim[2] * (lNodesPerDim[0] + lNodesPerDim[1] + 1);
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);
      receiveLIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive front-left edge nodes
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        receiveGIDs[countIDs] = k * gNodesPerDim[1] * gNodesPerDim[0] + startGID - gNodesPerDim[0] - 1;
        receivePIDs[countIDs] = myRank - procsPerDim[0] - 1;
        ++countIDs;
        // Receive front face nodes
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          receiveGIDs[countIDs] = k * gNodesPerDim[1] * gNodesPerDim[0] + i + startGID - gNodesPerDim[0];
          receivePIDs[countIDs] = myRank - procsPerDim[0];
          ++countIDs;
        }
        // Receive left face nodes
        for (LO j = 0; j < lNodesPerDim[1]; ++j) {
          receiveGIDs[countIDs] = k * gNodesPerDim[1] * gNodesPerDim[0] + j * gNodesPerDim[0] + startGID - 1;
          receivePIDs[countIDs] = myRank - 1;
          ++countIDs;
        }
      }

      if (useUnstructured) {  // Set parameters for interface aggregation
        numInterfaces = 2;

        interfaceLIDs.resize((lNodesPerDim[0] + 1) * lNodesPerDim[2] + (lNodesPerDim[1] + 1) * lNodesPerDim[2]);

        LO nodeOffset = 0, nodeIdx, nodeLID;
        // Front face
        for (LO k = 0; k < lNodesPerDim[2]; ++k) {
          for (LO i = 0; i < lNodesPerDim[0] + 1; ++i) {
            nodeIdx                = k * (lNodesPerDim[0] + 1) + i;
            nodeLID                = k * (lNodesPerDim[0] + 1) * (lNodesPerDim[1] + 1) + i;
            interfaceLIDs[nodeIdx] = nodeLID;
          }
        }
        // Left face
        nodeOffset += (lNodesPerDim[0] + 1) * lNodesPerDim[2];
        for (LO k = 0; k < lNodesPerDim[2]; ++k) {
          for (LO j = 0; j < lNodesPerDim[1] + 1; ++j) {
            nodeIdx                = k * (lNodesPerDim[1] + 1) + j + nodeOffset;
            nodeLID                = k * (lNodesPerDim[0] + 1) * (lNodesPerDim[1] + 1) + j * (lNodesPerDim[0] + 1);
            interfaceLIDs[nodeIdx] = nodeLID;
          }
        }
      }

      // Single face received
    } else if (bottomBC == 0) {
      numReceive = lNodesPerDim[0] * lNodesPerDim[1];
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);
      receiveLIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive bottom face nodes
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          receiveGIDs[countIDs] = j * gNodesPerDim[0] + i + startGID - gNodesPerDim[1] * gNodesPerDim[0];
          receivePIDs[countIDs] = myRank - procsPerDim[1] * procsPerDim[0];
          ++countIDs;
        }
      }

    } else if (frontBC == 0) {
      numReceive = lNodesPerDim[0] * lNodesPerDim[2];
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);
      receiveLIDs.resize(numReceive);

      LO countIDs = 0;
      // Receive front face nodes
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          receiveGIDs[countIDs] = k * gNodesPerDim[1] * gNodesPerDim[0] + i + startGID - gNodesPerDim[0];
          receivePIDs[countIDs] = myRank - procsPerDim[0];
          ++countIDs;
        }
      }

    } else if (leftBC == 0) {
      numReceive = lNodesPerDim[1] * lNodesPerDim[2];
      receiveGIDs.resize(numReceive);
      receivePIDs.resize(numReceive);
      receiveLIDs.resize(numReceive);

      LO countIDs = 0;
      // Recive left face nodes
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO j = 0; j < lNodesPerDim[1]; ++j) {
          receiveGIDs[countIDs] = k * gNodesPerDim[1] * gNodesPerDim[0] + j * gNodesPerDim[0] + startGID - 1;
          receivePIDs[countIDs] = myRank - 1;
          ++countIDs;
        }
      }
    }

    // Sent nodes
    if ((topBC == 0) && (backBC == 0) && (rightBC == 0)) {
      numSend = (lNodesPerDim[0]) * (lNodesPerDim[1]) + (lNodesPerDim[0]) * (lNodesPerDim[2]) + (lNodesPerDim[1]) * (lNodesPerDim[2]) + lNodesPerDim[0] + lNodesPerDim[1] + lNodesPerDim[2] + 1;
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);
      sendLIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of right face
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO j = 0; j < lNodesPerDim[1]; ++j) {
          sendGIDs[countIDs] = k * (gNodesPerDim[1] * gNodesPerDim[0]) + j * gNodesPerDim[0] + startGID + lNodesPerDim[0] - 1;
          sendLIDs[countIDs] = k * (lNodesPerDim[1] * lNodesPerDim[0]) + j * lNodesPerDim[0] + lNodesPerDim[0] - 1;
          sendPIDs[countIDs] = myRank + 1;
          ++countIDs;
        }
      }
      // Send nodes of back face
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          sendGIDs[countIDs] = k * (gNodesPerDim[1] * gNodesPerDim[0]) + i + startGID + (lNodesPerDim[1] - 1) * gNodesPerDim[0];
          sendLIDs[countIDs] = k * (lNodesPerDim[1] * lNodesPerDim[0]) + (lNodesPerDim[1] - 1) * lNodesPerDim[0] + i;
          sendPIDs[countIDs] = myRank + procsPerDim[0];
          ++countIDs;
        }
      }
      // Send nodes of right-back edge
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        sendGIDs[countIDs] = k * (gNodesPerDim[1] * gNodesPerDim[0]) + startGID + (lNodesPerDim[1] - 1) * gNodesPerDim[0] + lNodesPerDim[0] - 1;
        sendLIDs[countIDs] = k * (lNodesPerDim[1] * lNodesPerDim[0]) + lNodesPerDim[1] * lNodesPerDim[0] - 1;
        sendPIDs[countIDs] = myRank + procsPerDim[0] + 1;
        ++countIDs;
      }
      // Send nodes of top face
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          sendGIDs[countIDs] = j * gNodesPerDim[0] + i + startGID + (lNodesPerDim[2] - 1) * gNodesPerDim[1] * gNodesPerDim[0];
          sendLIDs[countIDs] = (lNodesPerDim[2] - 1) * lNodesPerDim[1] * lNodesPerDim[0] + j * lNodesPerDim[0] + i;
          sendPIDs[countIDs] = myRank + procsPerDim[1] * procsPerDim[0];
          ++countIDs;
        }
      }
      // Send nodes of top-right edge
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        sendGIDs[countIDs] = j * gNodesPerDim[0] + startGID + (lNodesPerDim[2] - 1) * gNodesPerDim[1] * gNodesPerDim[0] + lNodesPerDim[0] - 1;
        sendLIDs[countIDs] = (lNodesPerDim[2] - 1) * (lNodesPerDim[1] * lNodesPerDim[0]) + j * lNodesPerDim[0] + lNodesPerDim[0] - 1;
        sendPIDs[countIDs] = myRank + procsPerDim[1] * procsPerDim[0] + 1;
        ++countIDs;
      }
      // Send nodes of top-back edge
      for (LO i = 0; i < lNodesPerDim[0]; ++i) {
        sendGIDs[countIDs] = i + startGID + (lNodesPerDim[2] - 1) * gNodesPerDim[1] * gNodesPerDim[0] + (lNodesPerDim[1] - 1) * gNodesPerDim[0];
        sendLIDs[countIDs] = (lNodesPerDim[2] - 1) * (lNodesPerDim[1] * lNodesPerDim[0]) + (lNodesPerDim[1] - 1) * lNodesPerDim[0] + i;
        sendPIDs[countIDs] = myRank + procsPerDim[1] * procsPerDim[0] + procsPerDim[0];
        ++countIDs;
      }
      // Send node of top-back-right corner
      sendGIDs[countIDs] = startGID + (lNodesPerDim[2] - 1) * gNodesPerDim[1] * gNodesPerDim[0] + (lNodesPerDim[1] - 1) * gNodesPerDim[0] + lNodesPerDim[0] - 1;
      sendLIDs[countIDs] = lNodesPerDim[2] * lNodesPerDim[1] * lNodesPerDim[0] - 1;
      sendPIDs[countIDs] = myRank + procsPerDim[1] * procsPerDim[0] + procsPerDim[0] + 1;
      ++countIDs;

    } else if ((topBC == 0) && (backBC == 0)) {
      numSend = (lNodesPerDim[0] * lNodesPerDim[2])    // back face
                + (lNodesPerDim[0] * lNodesPerDim[1])  // Top face
                + (lNodesPerDim[0]);                   // top-back edge
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);
      sendLIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of back face
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          sendGIDs[countIDs] = k * (gNodesPerDim[1] * gNodesPerDim[0]) + i + startGID + (lNodesPerDim[1] - 1) * gNodesPerDim[0];
          sendLIDs[countIDs] = k * lNodesPerDim[1] * lNodesPerDim[0] + (lNodesPerDim[1] - 1) * lNodesPerDim[0] + i;
          sendPIDs[countIDs] = myRank + procsPerDim[0];
          ++countIDs;
        }
      }
      // Send nodes of top face
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          sendGIDs[countIDs] = j * gNodesPerDim[0] + i + startGID + (lNodesPerDim[2] - 1) * gNodesPerDim[1] * gNodesPerDim[0];
          sendLIDs[countIDs] = (lNodesPerDim[2] - 1) * lNodesPerDim[1] * lNodesPerDim[0] + j * lNodesPerDim[0] + i;
          sendPIDs[countIDs] = myRank + procsPerDim[1] * procsPerDim[0];
          ++countIDs;
        }
      }
      // Send nodes of top-back edge
      for (LO i = 0; i < lNodesPerDim[0]; ++i) {
        sendGIDs[countIDs] = i + startGID + (lNodesPerDim[2] - 1) * gNodesPerDim[1] * gNodesPerDim[0] + (lNodesPerDim[1] - 1) * gNodesPerDim[0];
        sendLIDs[countIDs] = (lNodesPerDim[2] - 1) * lNodesPerDim[1] * lNodesPerDim[0] + (lNodesPerDim[1] - 1) * lNodesPerDim[0] + i;
        sendPIDs[countIDs] = myRank + procsPerDim[1] * procsPerDim[0] + procsPerDim[0];
        ++countIDs;
      }

    } else if ((topBC == 0) && (rightBC == 0)) {
      numSend = (lNodesPerDim[1] * lNodesPerDim[2])    // right face
                + (lNodesPerDim[0] * lNodesPerDim[1])  // Top face
                + (lNodesPerDim[1]);                   // top-right edge
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);
      sendLIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of right face
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO j = 0; j < lNodesPerDim[1]; ++j) {
          sendGIDs[countIDs] = k * (gNodesPerDim[1] * gNodesPerDim[0]) + j * gNodesPerDim[0] + startGID + lNodesPerDim[0] - 1;
          sendLIDs[countIDs] = k * (lNodesPerDim[1] * lNodesPerDim[0]) + j * lNodesPerDim[0] + lNodesPerDim[0] - 1;
          sendPIDs[countIDs] = myRank + 1;
          ++countIDs;
        }
      }
      // Send nodes of top face
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          sendGIDs[countIDs] = j * gNodesPerDim[0] + i + startGID + (lNodesPerDim[2] - 1) * gNodesPerDim[1] * gNodesPerDim[0];
          sendLIDs[countIDs] = (lNodesPerDim[2] - 1) * lNodesPerDim[1] * lNodesPerDim[0] + j * lNodesPerDim[0] + i;
          sendPIDs[countIDs] = myRank + procsPerDim[1] * procsPerDim[0];
          ++countIDs;
        }
      }
      // Send nodes of top-right edge
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        sendGIDs[countIDs] = j * gNodesPerDim[0] + startGID + (lNodesPerDim[2] - 1) * gNodesPerDim[1] * gNodesPerDim[0] + lNodesPerDim[0] - 1;
        sendLIDs[countIDs] = (lNodesPerDim[2] - 1) * lNodesPerDim[1] * lNodesPerDim[0] + j * lNodesPerDim[0] + lNodesPerDim[0] - 1;
        sendPIDs[countIDs] = myRank + procsPerDim[1] * procsPerDim[0] + 1;
        ++countIDs;
      }

    } else if ((backBC == 0) && (rightBC == 0)) {
      numSend = lNodesPerDim[2] * (lNodesPerDim[0] + lNodesPerDim[1] + 1);
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);
      sendLIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of right face
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO j = 0; j < lNodesPerDim[1]; ++j) {
          sendGIDs[countIDs] = k * gNodesPerDim[1] * gNodesPerDim[0] + j * gNodesPerDim[0] + startGID + lNodesPerDim[0] - 1;
          sendLIDs[countIDs] = k * lNodesPerDim[1] * lNodesPerDim[0] + j * lNodesPerDim[0] + lNodesPerDim[0] - 1;
          sendPIDs[countIDs] = myRank + 1;
          ++countIDs;
        }
      }
      // Send nodes of back face
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          sendGIDs[countIDs] = k * gNodesPerDim[1] * gNodesPerDim[0] + i + startGID + (lNodesPerDim[1] - 1) * gNodesPerDim[0];
          sendLIDs[countIDs] = k * lNodesPerDim[1] * lNodesPerDim[0] + (lNodesPerDim[1] - 1) * lNodesPerDim[0] + i;
          sendPIDs[countIDs] = myRank + procsPerDim[0];
          ++countIDs;
        }
      }
      // Send nodes of back-right edge
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        sendGIDs[countIDs] = k * gNodesPerDim[1] * gNodesPerDim[0] + startGID + (lNodesPerDim[1] - 1) * gNodesPerDim[0] + lNodesPerDim[0] - 1;
        sendLIDs[countIDs] = k * (lNodesPerDim[1] * lNodesPerDim[0]) + lNodesPerDim[1] * lNodesPerDim[0] - 1;
        sendPIDs[countIDs] = myRank + procsPerDim[0] + 1;
        ++countIDs;
      }

    } else if (topBC == 0) {
      numSend = lNodesPerDim[0] * lNodesPerDim[1];
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);
      sendLIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of top face
      for (LO j = 0; j < lNodesPerDim[1]; ++j) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          sendGIDs[countIDs] = j * gNodesPerDim[0] + i + startGID + (lNodesPerDim[2] - 1) * gNodesPerDim[0] * gNodesPerDim[1];
          sendLIDs[countIDs] = (lNodesPerDim[2] - 1) * lNodesPerDim[1] * lNodesPerDim[0] + j * lNodesPerDim[0] + i;
          sendPIDs[countIDs] = myRank + procsPerDim[1] * procsPerDim[0];
          ++countIDs;
        }
      }

    } else if (backBC == 0) {
      numSend = lNodesPerDim[0] * lNodesPerDim[2];
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);
      sendLIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of back face
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO i = 0; i < lNodesPerDim[0]; ++i) {
          sendGIDs[countIDs] = k * gNodesPerDim[1] * gNodesPerDim[0] + i + startGID + (lNodesPerDim[1] - 1) * gNodesPerDim[0];
          sendLIDs[countIDs] = k * lNodesPerDim[1] * lNodesPerDim[0] + (lNodesPerDim[1] - 1) * lNodesPerDim[0] + i;
          sendPIDs[countIDs] = myRank + procsPerDim[0];
          ++countIDs;
        }
      }

    } else if (rightBC == 0) {
      numSend = lNodesPerDim[1] * lNodesPerDim[2];
      sendGIDs.resize(numSend);
      sendPIDs.resize(numSend);
      sendLIDs.resize(numSend);

      LO countIDs = 0;
      // Send nodes of right face
      for (LO k = 0; k < lNodesPerDim[2]; ++k) {
        for (LO j = 0; j < lNodesPerDim[1]; ++j) {
          sendGIDs[countIDs] = k * gNodesPerDim[1] * gNodesPerDim[0] + j * gNodesPerDim[0] + startGID + lNodesPerDim[0] - 1;
          sendLIDs[countIDs] = k * lNodesPerDim[1] * lNodesPerDim[0] + j * lNodesPerDim[0] + lNodesPerDim[0] - 1;
          sendPIDs[countIDs] = myRank + 1;
          ++countIDs;
        }
      }
    }
  }

  // The code below could probably be separated and performed in its own function.
  // It pretty much tries to go from the more geometric based information generated
  // above to the quasiRegion map which means finding the GID of nodes-off rank on the interface.
  const LO numLocalCompositeNodes = lNodesPerDim[0] * lNodesPerDim[1] * lNodesPerDim[2];
  numLocalRegionNodes             = numLocalCompositeNodes + numReceive;
  quasiRegionGIDs.resize(numLocalRegionNodes * numDofsPerNode);
  quasiRegionCoordGIDs.resize(numLocalRegionNodes);

  rNodesPerDim[0] = lNodesPerDim[0];
  rNodesPerDim[1] = lNodesPerDim[1];
  rNodesPerDim[2] = lNodesPerDim[2];
  if (leftBC == 0) {
    rNodesPerDim[0] += 1;
  }
  if (frontBC == 0) {
    rNodesPerDim[1] += 1;
  }
  if (bottomBC == 0) {
    rNodesPerDim[2] += 1;
  }

  // Using receiveGIDs, rNodesPerDim and numLocalRegionNodes, build quasi-region row map
  // This will potentially be done by the application or in a MueLu interface but for now
  // let us keep it in this utility function.
  LO interfaceCount = 0, compositeIdx = 0;
  Teuchos::Array<LO> regionIJK(3);
  for (LO nodeRegionIdx = 0; nodeRegionIdx < numLocalRegionNodes; ++nodeRegionIdx) {
    regionIJK[2] = nodeRegionIdx / (rNodesPerDim[1] * rNodesPerDim[0]);
    LO tmp       = nodeRegionIdx % (rNodesPerDim[1] * rNodesPerDim[0]);
    regionIJK[1] = tmp / rNodesPerDim[0];
    regionIJK[0] = tmp % rNodesPerDim[0];

    if ((regionIJK[0] == 0 && leftBC == 0) ||
        (regionIJK[1] == 0 && frontBC == 0) ||
        (regionIJK[2] == 0 && bottomBC == 0)) {
      quasiRegionCoordGIDs[nodeRegionIdx] = receiveGIDs[interfaceCount];
      for (int dof = 0; dof < numDofsPerNode; ++dof) {
        quasiRegionGIDs[nodeRegionIdx * numDofsPerNode + dof] =
            receiveGIDs[interfaceCount] * numDofsPerNode + dof;
      }
      receiveLIDs[interfaceCount] = nodeRegionIdx;
      ++interfaceCount;
    } else {
      compositeIdx                        = (regionIJK[2] + bottomBC - 1) * lNodesPerDim[1] * lNodesPerDim[0] + (regionIJK[1] + frontBC - 1) * lNodesPerDim[0] + (regionIJK[0] + leftBC - 1);
      quasiRegionCoordGIDs[nodeRegionIdx] = nodeMap->getGlobalElement(compositeIdx);
      for (int dof = 0; dof < numDofsPerNode; ++dof) {
        quasiRegionGIDs[nodeRegionIdx * numDofsPerNode + dof]      = dofMap->getGlobalElement(compositeIdx * numDofsPerNode + dof);
        compositeToRegionLIDs[compositeIdx * numDofsPerNode + dof] = nodeRegionIdx * numDofsPerNode + dof;
      }
    }
  }

  // Here we gather the interface GIDs (in composite layout)
  // and the interface LIDs (in region layout) for the local rank
  interfaceLIDsData.resize((sendGIDs.size() + receiveGIDs.size()) * numDofsPerNode);
  interfaceGIDs.resize((sendGIDs.size() + receiveGIDs.size()) * numDofsPerNode);
  using size_type = typename Teuchos::Array<GO>::size_type;
  for (size_type nodeIdx = 0; nodeIdx < sendGIDs.size(); ++nodeIdx) {
    for (int dof = 0; dof < numDofsPerNode; ++dof) {
      LO dofIdx                 = nodeIdx * numDofsPerNode + dof;
      interfaceGIDs[dofIdx]     = sendGIDs[nodeIdx] * numDofsPerNode + dof;
      interfaceLIDsData[dofIdx] = compositeToRegionLIDs[sendLIDs[nodeIdx] * numDofsPerNode + dof];
    }
  }
  for (size_type nodeIdx = 0; nodeIdx < receiveGIDs.size(); ++nodeIdx) {
    for (int dof = 0; dof < numDofsPerNode; ++dof) {
      LO dofIdx                                                    = nodeIdx * numDofsPerNode + dof;
      interfaceGIDs[dofIdx + sendGIDs.size() * numDofsPerNode]     = receiveGIDs[nodeIdx] * numDofsPerNode + dof;
      interfaceLIDsData[dofIdx + sendLIDs.size() * numDofsPerNode] = receiveLIDs[nodeIdx] * numDofsPerNode + dof;
    }
  }

  // Have all the GIDs and LIDs we stort them in place with std::sort()
  // Subsequently we bring unique values to the beginning of the array with
  // std::unique() and delete the duplicates with erase.
  std::sort(interfaceLIDsData.begin(), interfaceLIDsData.end());
  interfaceLIDsData.erase(std::unique(interfaceLIDsData.begin(), interfaceLIDsData.end()),
                          interfaceLIDsData.end());
  std::sort(interfaceGIDs.begin(), interfaceGIDs.end());
  interfaceGIDs.erase(std::unique(interfaceGIDs.begin(), interfaceGIDs.end()),
                      interfaceGIDs.end());

}  // createRegionData

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void MakeRegionPerGIDWithGhosts(const Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >& nodeMap,
                                const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > regionRowMap,
                                const Teuchos::RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> >& rowImport,
                                const int maxRegPerGID,
                                const LocalOrdinal numDofsPerNode,
                                const Teuchos::Array<LocalOrdinal>& lNodesPerDir,
                                const Teuchos::Array<GlobalOrdinal>& sendGIDs,
                                const Teuchos::Array<int>& sendPIDs,
                                const Teuchos::Array<LocalOrdinal>& interfaceRegionLIDs,
                                Teuchos::RCP<Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> >& regionsPerGIDWithGhosts,
                                Teuchos::RCP<Xpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> >& interfaceGIDsMV) {
  using LO = LocalOrdinal;
  using GO = GlobalOrdinal;
  using NO = Node;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;

  const RCP<const Xpetra::Map<LO, GO, NO> > dofMap            = rowImport->getSourceMap();
  const RCP<const Xpetra::Map<LO, GO, NO> > quasiRegionRowMap = rowImport->getTargetMap();
  const int myRank                                            = dofMap->getComm()->getRank();

  RCP<Xpetra::MultiVector<LO, LO, GO, NO> > regionsPerGID =
      Xpetra::MultiVectorFactory<LO, LO, GO, NO>::Build(dofMap, maxRegPerGID, false);
  regionsPerGIDWithGhosts =
      Xpetra::MultiVectorFactory<LO, LO, GO, NO>::Build(quasiRegionRowMap, maxRegPerGID, false);

  {  // Scope for regionsPerGIDView
    Array<ArrayRCP<LO> > regionsPerGIDView(maxRegPerGID);
    for (int regionIdx = 0; regionIdx < maxRegPerGID; ++regionIdx) {
      regionsPerGIDView[regionIdx] = regionsPerGID->getDataNonConst(regionIdx);
    }

    // Initialize all entries to myRank in first column and to -1 in other columns
    for (LO dofIdx = 0; dofIdx < lNodesPerDir[0] * lNodesPerDir[1] * lNodesPerDir[2] * numDofsPerNode; ++dofIdx) {
      regionsPerGIDView[0][dofIdx] = myRank;
      for (int regionIdx = 1; regionIdx < maxRegPerGID; ++regionIdx) {
        regionsPerGIDView[regionIdx][dofIdx] = -1;
      }
    }

    // Now loop over the sendGIDs array to fill entries with values in sendPIDs
    LO nodeIdx = 0;
    for (LO sendIdx = 0; sendIdx < static_cast<LO>(sendPIDs.size()); ++sendIdx) {
      nodeIdx = nodeMap->getLocalElement(sendGIDs[sendIdx]);
      for (int dof = 0; dof < numDofsPerNode; ++dof) {
        LO dofIdx = nodeIdx * numDofsPerNode + dof;
        for (int regionIdx = 1; regionIdx < maxRegPerGID; ++regionIdx) {
          if (regionsPerGIDView[regionIdx][dofIdx] == -1) {
            regionsPerGIDView[regionIdx][dofIdx] = sendPIDs[sendIdx];
            break;
          }
        }
      }
    }
  }

  regionsPerGIDWithGhosts->doImport(*regionsPerGID, *rowImport, Xpetra::INSERT);

  interfaceGIDsMV = Xpetra::MultiVectorFactory<GO, LO, GO, NO>::Build(quasiRegionRowMap, maxRegPerGID, false);
  interfaceGIDsMV->putScalar(Teuchos::OrdinalTraits<GO>::zero());
  const LO numRegionInterfaceLIDs = static_cast<LO>(interfaceRegionLIDs.size());
  {  // Scope for interfaceGIDsPerRegion
    Array<ArrayRCP<LO> > regionsPerGIDWithGhostsData(maxRegPerGID);
    Array<ArrayRCP<GO> > interfaceGIDsMVData(maxRegPerGID);
    for (int regionIdx = 0; regionIdx < maxRegPerGID; ++regionIdx) {
      regionsPerGIDWithGhostsData[regionIdx] = regionsPerGIDWithGhosts->getDataNonConst(regionIdx);
      interfaceGIDsMVData[regionIdx]         = interfaceGIDsMV->getDataNonConst(regionIdx);
      for (LO idx = 0; idx < numRegionInterfaceLIDs; ++idx) {
        LO LID = interfaceRegionLIDs[idx];
        if (regionsPerGIDWithGhostsData[regionIdx][LID] == myRank) {
          interfaceGIDsMVData[regionIdx][LID] = regionRowMap->getGlobalElement(LID);
        }
      }
    }
  }

}  // MakeRegionPerGIDWithGhosts

/*!
\brief Extract list of region GIDs of all interface DOFs from the region row map

Starting from the known list of \c interfaceRegionLIDs, we know, which entries in the \c regionRowMap
refert to interface DOFs, so we can grab them and stick them into the list of \c interfaceRegionGIDs.
*/
template <class LocalOrdinal, class GlobalOrdinal, class Node>
void ExtractListOfInterfaceRegionGIDs(
    Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > regionRowMap,
    const Teuchos::Array<LocalOrdinal>& interfaceRegionLIDs, Teuchos::Array<GlobalOrdinal>& interfaceRegionGIDs) {
  interfaceRegionGIDs.resize(interfaceRegionLIDs.size());
  for (LocalOrdinal interfaceIdx = 0; interfaceIdx < static_cast<LocalOrdinal>(interfaceRegionLIDs.size()); ++interfaceIdx) {
    interfaceRegionGIDs[interfaceIdx] =
        regionRowMap->getGlobalElement(interfaceRegionLIDs[interfaceIdx]);
  }
}  // ExtractListOfInterfaceRegionGIDs

#endif  // MUELU_SETUPREGIONUTILITIES_HPP
