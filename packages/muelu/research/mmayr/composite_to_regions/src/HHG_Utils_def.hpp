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
#ifndef MUELU_HHG_UTILS_DEF_HPP
#define MUELU_HHG_UTILS_DEF_HPP

#include <Kokkos_DefaultNode.hpp>

#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MultiVector.hpp>

#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_RCP.hpp>

// Input data is read into a generic vector.
// Use these enums to access entries in this vector.
enum InputDataIndices
{
  inpData_isStructured,
  inpData_ownedX,
  inpData_ownedY,
  inpData_regionX,
  inpData_regionY,
  inpData_cornerX,
  inpData_cornerY,
  inpData_nGhosts,
  inpData_firstLIDsOfGhosts
};

// this little widget handles application specific data
// used to implement LIDregion() in the HHG driver
struct widget {

  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef KokkosClassic::DefaultNode::DefaultNodeType Node;

  int maxRegPerProc;
  int *minGIDComp;
  int *maxGIDComp;
  int *myRegions;
  Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> *colMap;
  int maxRegPerGID;
  Teuchos::RCP<Xpetra::MultiVector<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node>> regionsPerGIDWithGhosts;
  int *gDim, *lDim, *lowInd;
  int *trueCornerx; // global coordinates of region
  int *trueCornery; // corner within entire 2D mesh
  int *relcornerx;  // coordinates of corner relative
  int *relcornery;  // to region corner
  int *lDimx;
  int *lDimy;
  int nx;
  int myRank;
};

//! Print an object in regional layout to screen
template <class T>
void printRegionalObject(const std::string objName, ///< string to be used for screen output
    const std::vector<Teuchos::RCP<T> > regObj, ///< regional object to be printed to screen
    const int myRank, ///< rank of calling proc
    Teuchos::FancyOStream& outstream ///< output stream
    )
{
  for (int j = 0; j < (int) regObj.size(); j++) {
    outstream << myRank << ": " << objName << " " << j << std::endl;
    regObj[j]->describe(outstream, Teuchos::VERB_EXTREME);
  }
}

#endif // MUELU_HHG_UTILS_DEF_HPP
