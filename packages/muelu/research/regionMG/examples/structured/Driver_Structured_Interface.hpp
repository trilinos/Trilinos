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
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#include <iostream>

//Teuchos
#include <Teuchos_Array.hpp>

template<class LocalOrdinal>
void findInterface(const int numDimensions, Teuchos::Array<LocalOrdinal> nodesPerDim,
                   Teuchos::Array<int> boundaryConditions,
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
  for(LO dim = 0; dim < 3; ++dim) {
    // Check for nodes and no boundary conditions in this direction, otherwise skip it.
    if((nodesPerDim[dim] == 1) || (boundaryConditions[2*dim] + boundaryConditions[2*dim + 1] == 2)) {
      continue;
    }

    // Since we are not skipping this direction
    // we at least need to store data for one surface
    numInterfaceLIDs = 1;
    for(LO dimIdx = 0; dimIdx < 3; ++dimIdx) {
      interfacesDimensions[numInterfaces*3 + dimIdx] = (dimIdx == dim ? 1 : nodesPerDim[dimIdx]);
      numInterfaceLIDs *= (dimIdx == dim ? 1 : nodesPerDim[dimIdx]);
    }
    numTotalLIDs += numInterfaceLIDs;
    ++numInterfaces;

    // If there are no BC we need to store the surface twice.
    if(boundaryConditions[2*dim] + boundaryConditions[2*dim + 1] == 0) {
      for(LO dimIdx = 0; dimIdx < 3; ++dimIdx) {
        interfacesDimensions[numInterfaces*3 + dimIdx]
          = interfacesDimensions[(numInterfaces - 1)*3 + dimIdx];
      }
      numTotalLIDs += numInterfaceLIDs;
      ++numInterfaces;
    }
  }
  interfacesDimensions.resize(3*numInterfaces);
  interfacesLIDs.resize(numTotalLIDs, -1);

  // Step 2 lazy implementation of all geometrical cases.

  LO nodeOffset = 0;
  if(numDimensions == 2) {
    // left interface
    if(boundaryConditions[0] == 0) {
      for(LO nodeIdx = 0; nodeIdx < nodesPerDim[1]; ++nodeIdx) {
        interfacesLIDs[nodeOffset + nodeIdx] = nodeIdx*nodesPerDim[0];
      }
      nodeOffset += nodesPerDim[1];
    }

    // right interface
    if(boundaryConditions[1] == 0) {
      for(LO nodeIdx = 0; nodeIdx < nodesPerDim[1]; ++nodeIdx) {
        interfacesLIDs[nodeOffset + nodeIdx] = (nodeIdx + 1)*nodesPerDim[0] - 1;
      }
      nodeOffset += nodesPerDim[1];
    }

    // front interface
    if(boundaryConditions[2] == 0) {
      for(LO nodeIdx = 0; nodeIdx < nodesPerDim[0]; ++nodeIdx) {
        interfacesLIDs[nodeOffset + nodeIdx] = nodeIdx;
      }
      nodeOffset += nodesPerDim[0];
    }

    // back interface
    if(boundaryConditions[3] == 0) {
      for(LO nodeIdx = 0; nodeIdx < nodesPerDim[0]; ++nodeIdx) {
        interfacesLIDs[nodeOffset + nodeIdx] = (nodesPerDim[1] - 1)*nodesPerDim[0] + nodeIdx;
      }
      nodeOffset += nodesPerDim[0];
    }
  }

  if(numDimensions == 3) {
    // left interface
    if(boundaryConditions[0] == 0) {
      for(LO k = 0; k < nodesPerDim[2]; ++k) {
        for(LO j = 0; j < nodesPerDim[1]; ++j) {
          interfacesLIDs[nodeOffset + k*nodesPerDim[1] + j]
            = k*nodesPerDim[1]*nodesPerDim[0] + j*nodesPerDim[0];
        }
      }
      nodeOffset += nodesPerDim[2]*nodesPerDim[1];
    }

    // right interface
    if(boundaryConditions[1] == 0) {
      for(LO k = 0; k < nodesPerDim[2]; ++k) {
        for(LO j = 0; j < nodesPerDim[1]; ++j) {
          interfacesLIDs[nodeOffset + k*nodesPerDim[1] + j]
            = k*nodesPerDim[1]*nodesPerDim[0] + (j + 1)*nodesPerDim[0] - 1;
        }
      }
      nodeOffset += nodesPerDim[2]*nodesPerDim[1];
    }

    // front interface
    if(boundaryConditions[2] == 0) {
      for(LO k = 0; k < nodesPerDim[2]; ++k) {
        for(LO i = 0; i < nodesPerDim[0]; ++i) {
          interfacesLIDs[nodeOffset + k*nodesPerDim[0] + i]
            = k*nodesPerDim[1]*nodesPerDim[0] + i;
        }
      }
      nodeOffset += nodesPerDim[2]*nodesPerDim[0];
    }

    // back interface
    if(boundaryConditions[3] == 0) {
      for(LO k = 0; k < nodesPerDim[2]; ++k) {
        for(LO i = 0; i < nodesPerDim[0]; ++i) {
          interfacesLIDs[nodeOffset + k*nodesPerDim[0] + i]
            = k*nodesPerDim[1]*nodesPerDim[0] + (nodesPerDim[1] - 1)*nodesPerDim[0] + i;
        }
      }
      nodeOffset += nodesPerDim[2]*nodesPerDim[0];
    }

    // bottom interface
    if(boundaryConditions[4] == 0) {
      for(LO j = 0; j < nodesPerDim[1]; ++j) {
        for(LO i = 0; i < nodesPerDim[0]; ++i) {
          interfacesLIDs[nodeOffset + j*nodesPerDim[0] + i]
            = j*nodesPerDim[0] + i;
        }
      }
      nodeOffset += nodesPerDim[1]*nodesPerDim[0];
    }

    //top interface
    if(boundaryConditions[5] == 0) {
      for(LO j = 0; j < nodesPerDim[1]; ++j) {
        for(LO i = 0; i < nodesPerDim[0]; ++i) {
          interfacesLIDs[nodeOffset + j*nodesPerDim[0] + i]
            = (nodesPerDim[2] - 1)*nodesPerDim[1]*nodesPerDim[0] + j*nodesPerDim[0] + i;
        }
      }
      nodeOffset += nodesPerDim[1]*nodesPerDim[0];
    }
  }
}
