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
#ifndef MUELU_STRUCTUREDAGGREGATIONFACTORY_DECL_HPP_
#define MUELU_STRUCTUREDAGGREGATIONFACTORY_DECL_HPP_


#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_StructuredAggregationFactory_fwd.hpp"

#include "MueLu_AggregationAlgorithmBase.hpp"

#include "MueLu_Level_fwd.hpp"
//#include "MueLu_Graph_fwd.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_LocalLexicographicIndexManager.hpp"
#include "MueLu_GlobalLexicographicIndexManager.hpp"

namespace MueLu {

/*!
    @class StructuredAggregationFactory class.
    @brief Factory for building aggregates on structured grids.

    Factory for creating aggregates from grid structure of the problem. The structured aggregation
    method can return an aggregate structure or a geometric used by prolongator factories.

    Internally, each node has a status which can be one of the following:

    Node status | Meaning
    ------------|---------
    READY       | Node is not aggregated and can be used for building a new aggregate or can be added to an existing aggregate.
    AGGREGATED  | Node is aggregated.
    IGNORED     | Node is not considered for aggregation (it may have been dropped or put into a singleton aggregate)
    BOUNDARY    | Node is a Dirichlet boundary node (with one or more Dirichlet boundary conditions).
    ONEPT       | The user forces the aggregation algorithm to treat the node as a singleton. Important: Do not forget to set aggregation: allow user-specified singletons to true! Otherwise Phase3 will just handle the ONEPT nodes and probably not build singletons

    @ingroup Aggregation

    ## Input/output of StructuredAggregationFactory ##

    ### User parameters of StructuredAggregationFactory ###
    Parameter | type | default | master.xml | validated | requested | description
    ----------|------|---------|:----------:|:---------:|:---------:|------------
     DofsPerNode        | Factory | null |   | * | * | Generating factory for variable 'DofsPerNode', usually the same as for 'Graph'
     OnePt aggregate map name  | string |  | | * | * | Name of input map for single node aggregates (default=''). Makes only sense if the parameter 'aggregation: allow user-specified singletons' is set to true.
     OnePt aggregate map factory | Factory | null |   | * | * | Generating factory of (DOF) map for single node aggregates.  Makes only sense if the parameter 'aggregation: allow user-specified singletons' is set to true.


    The * in the @c master.xml column denotes that the parameter is defined in the @c master.xml file.<br>
    The * in the @c validated column means that the parameter is declared in the list of valid input parameters (see StructuredAggregationFactory::GetValidParameters).<br>
    The * in the @c requested column states that the data is requested as input with all dependencies (see StructuredAggregationFactory::DeclareInput).

    ### Variables provided by StructuredAggregationFactory ###

    After StructuredAggregationFactory::Build the following data is available (if requested)

    Parameter | generated by | description
    ----------|--------------|------------
    | Aggregates   | StructuredAggregationFactory   | Container class with aggregation information. See also Aggregates.
*/

  template <class Scalar = double,
            class LocalOrdinal = int,
            class GlobalOrdinal = LocalOrdinal,
            class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class StructuredAggregationFactory : public SingleLevelFactoryBase {
#undef MUELU_STRUCTUREDAGGREGATIONFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    StructuredAggregationFactory();

    //! Destructor.
    virtual ~StructuredAggregationFactory() { }

    RCP<const ParameterList> GetValidParameterList() const;

    //@}

    //! @name Set/get methods.
    //@{
    // set information about 1-node aggregates (map name and generating factory)
    void SetOnePtMapName(const std::string name, Teuchos::RCP<const FactoryBase> mapFact) {
      SetParameter("OnePt aggregate map name", ParameterEntry(std::string(name))); // revalidate
      SetFactory("OnePt aggregate map factory",mapFact);
    }

    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const;

    //@}

    //! @name Build methods.
    //@{

    /*! @brief Build aggregates. */
    void Build(Level &currentLevel) const;

    //@}

private:

    //! boolean flag: definition phase
    //! if true, the aggregation algorithms still can be set and changed.
    //! if false, no change in aggregation algorithms is possible any more
    mutable bool bDefinitionPhase_;

    struct GeometricData {
      // Geometric algorithm require a copious amount of data to be passed around so this struct
      // will reduce the amount of input/output parameters of methods in the class. Additionally
      // the struct can be rewritten to accomodate constraints of Kokkos/CUDA data types

      std::string meshLayout = "Global Lexicographic";
      int numDimensions = -1, myRank = -1, numRanks = -1;
      LO lNumFineNodes = -1, lNumCoarseNodes = -1, lNumGhostNodes = -1,lNumGhostedNodes  = -1;
      LO myBlock = -1, numBlocks = -1, lNumFineNodes10 = -1, lNumCoarseNodes10 = -1;
      LO numGhostedCoarseNodes = -1, numGhostedCoarseNodes10 = -1, pi = -1, pj = -1, pk = -1;
      LO myRankIndex = -1;
      GO gNumFineNodes = -1, gNumCoarseNodes = -1, gNumFineNodes10 = -1, gNumCoarseNodes10 = -1;
      GO minGlobalIndex = -1;
      Array<int> coarseRate, endRate, rankIndices;
      Array<LO> lFineNodesPerDir, lCoarseNodesPerDir, offsets, ghostedCoarseNodesPerDir;
      Array<GO> startIndices, gFineNodesPerDir, gCoarseNodesPerDir, startGhostedCoarseNode;
      std::vector<std::vector<GO> > meshData; // These are sorted later so they are in std::vector
      std::vector<std::vector<GO> > coarseMeshData; // These are sorted while being computed
      bool ghostInterface[6] = {false}, ghostedDir[6] = {false};
      typename std::vector<std::vector<GO> >::iterator myBlockStart, myBlockEnd;

      GeometricData() {
        coarseRate.resize(3);
        endRate.resize(3);
        lFineNodesPerDir.resize(3);
        lCoarseNodesPerDir.resize(3);
        offsets.resize(6);
        ghostedCoarseNodesPerDir.resize(3);
        startIndices.resize(6);
        gFineNodesPerDir.resize(3);
        gCoarseNodesPerDir.resize(3);
        startGhostedCoarseNode.resize(3);
      }

      void computeGeometricParameters(const std::string MeshLayout, const int NumDimensions,
                                      const int MyRank, const int NumRanks,
                                      const Array<GO> GFineNodesPerDir,
                                      const Array<LO> LFineNodesPerDir,
                                      const Array<LO> CoarseRate, const Array<GO> MeshData,
                                      const GO MinGlobalIndex) {

        // First copy over the input data
        meshLayout       = MeshLayout;
        numDimensions    = NumDimensions;
        myRank           = MyRank;
        numRanks         = NumRanks;
        gFineNodesPerDir = GFineNodesPerDir;
        lFineNodesPerDir = LFineNodesPerDir;

        // Load coarse rate, being careful about formating
        for(int dim = 0; dim < 3; ++dim) {
          if(dim < numDimensions) {
            if(CoarseRate.size() == 1) {
              coarseRate[dim] = CoarseRate[0];
            } else if(CoarseRate.size() == numDimensions) {
              coarseRate[dim] = CoarseRate[dim];
            }
          } else {
            coarseRate[dim] = 1;
          }
        }

        // Load meshData for local lexicographic case
        if(meshLayout == "Local Lexicographic") {
          meshData.resize(numRanks);
          for(int rank = 0; rank < numRanks; ++rank) {
            meshData[rank].resize(10);
            for(int entry = 0; entry < 10; ++entry) {
              meshData[rank][entry] = MeshData[10*rank + entry];
            }
          }
        }

        // Start simple parameter calculation
        gNumFineNodes10 = gFineNodesPerDir[1]*gFineNodesPerDir[0];
        gNumFineNodes   = gFineNodesPerDir[2]*gNumFineNodes10;
        lNumFineNodes10 = lFineNodesPerDir[1]*lFineNodesPerDir[0];
        lNumFineNodes   = lFineNodesPerDir[2]*lNumFineNodes10;

        if(meshLayout == "Global Lexicographic") {
          GO tmp = 0;
          startIndices[2] = MinGlobalIndex / (gFineNodesPerDir[1]*gFineNodesPerDir[0]);
          tmp             = MinGlobalIndex % (gFineNodesPerDir[1]*gFineNodesPerDir[0]);
          startIndices[1] = tmp / gFineNodesPerDir[0];
          startIndices[0] = tmp % gFineNodesPerDir[0];

          for(int dim = 0; dim < 3; ++dim) {
            startIndices[dim + 3] = startIndices[dim] + lFineNodesPerDir[dim] - 1;
          }

        } else if(meshLayout == "Local Lexicographic") {
          myBlock = meshData[myRank][2];
          for(int i = 0; i < 3; ++i) {
            startIndices[i]     = meshData[myRank][2*i + 3];
            startIndices[i + 3] = meshData[myRank][2*i + 4];
          }
          sortLocalLexicographicData();
        }

        for(int dim = 0; dim < 3; ++dim) {
          if(dim < numDimensions) {
            offsets[dim]     = Teuchos::as<LO>(startIndices[dim]) % coarseRate[dim];
            offsets[dim + 3] = Teuchos::as<LO>(startIndices[dim]) % coarseRate[dim];

            if(startIndices[dim] % coarseRate[dim] != 0 ||
               startIndices[dim] == gFineNodesPerDir[dim]-1) {
              ghostInterface[2*dim] = true;
            }
            if((startIndices[dim + 3] != gFineNodesPerDir[dim] - 1) &&
               ((lFineNodesPerDir[dim] == 1) || (startIndices[dim + 3] % coarseRate[dim] != 0))) {
              ghostInterface[2*dim+1] = true;
            }
          }
        }


        // Here one element can represent either the degenerate case of one node or the more general
        // case of two nodes, i.e. x---x is a 1D element with two nodes and x is a 1D element with
        // one node. This helps generating a 3D space from tensorial products...
        // A good way to handle this would be to generalize the algorithm to take into account the
        // discretization order used in each direction, at least in the FEM sense, since a 0 degree
        // discretization will have a unique node per element. This way 1D discretization can be
        // viewed as a 3D problem with one 0 degree element in the y direction and one 0 degre
        // element in the z direction.
        // !!! Operations below are aftecting both local and global values that have two         !!!
        // different orientations. Orientations can be interchanged using mapDirG2L and mapDirL2G.
        // coarseRate, endRate and offsets are in the global basis, as well as all the variables
        // starting with a g.
        // !!! while the variables starting with an l are in the local basis.                    !!!
        for(int dim = 0; dim < 3; ++dim) {
          if(dim < numDimensions) {
            // This array is passed to the RAPFactory and eventually becomes gFineNodePerDir on the next
            // level.
            gCoarseNodesPerDir[dim] = (gFineNodesPerDir[dim] - 1) / coarseRate[dim];
            endRate[dim] = Teuchos::as<LO>((gFineNodesPerDir[dim] - 1) % coarseRate[dim]);
            if(endRate[dim] == 0) {
              endRate[dim] = coarseRate[dim];
              ++gCoarseNodesPerDir[dim];
            } else {
              gCoarseNodesPerDir[dim] += 2;
            }
          } else {
            endRate[dim] = 1;
            gCoarseNodesPerDir[dim] = 1;
          }
        }

        gNumCoarseNodes10 = gCoarseNodesPerDir[0]*gCoarseNodesPerDir[1];
        gNumCoarseNodes   = gNumCoarseNodes10*gCoarseNodesPerDir[2];

        for(LO dim = 0; dim < 3; ++dim) {
          if(dim < numDimensions) {
            // Check whether the partition includes the "end" of the mesh which means that endRate
            // will apply. Also make sure that endRate is not 0 which means that the mesh does not
            // require a particular treatment at the boundaries.
            if( (startIndices[dim] + lFineNodesPerDir[dim]) == gFineNodesPerDir[dim] ) {
              lCoarseNodesPerDir[dim] = (lFineNodesPerDir[dim] - endRate[dim] + offsets[dim] - 1)
                / coarseRate[dim] + 1;
              if(offsets[dim] == 0) {++lCoarseNodesPerDir[dim];}
            } else {
              lCoarseNodesPerDir[dim] = (lFineNodesPerDir[dim] + offsets[dim] - 1) /coarseRate[dim];
              if(offsets[dim] == 0) {++lCoarseNodesPerDir[dim];}
            }
          } else {
            lCoarseNodesPerDir[dim] = 1;
          }
          // This would happen if the rank does not own any nodes but in that case a subcommunicator
          // should be used so this should really not be a concern.
          if(lFineNodesPerDir[dim] < 1) {lCoarseNodesPerDir[dim] = 0;}
        }

        // Assuming linear interpolation, each fine point has contribution from 8 coarse points
        // and each coarse point value gets injected.
        // For systems of PDEs we assume that all dofs have the same P operator.
        lNumCoarseNodes10 = lCoarseNodesPerDir[0]*lCoarseNodesPerDir[1];
        lNumCoarseNodes   = lNumCoarseNodes10*lCoarseNodesPerDir[2];

        // For each direction, determine how many points (including ghosts) are required.
        for(int dim = 0; dim < 3; ++dim) {
          // The first branch of this if-statement will be used if the rank contains only one layer
          // of nodes in direction i, that layer must also coincide with the boundary of the mesh
          // and coarseRate[i] == endRate[i]...
          if(dim < numDimensions) {
            if((startIndices[dim] == gFineNodesPerDir[dim] - 1) &&
               (startIndices[dim] % coarseRate[dim] == 0)) {
              startGhostedCoarseNode[dim] = startIndices[dim] / coarseRate[dim] - 1;
            } else {
              startGhostedCoarseNode[dim] = startIndices[dim] / coarseRate[dim];
            }
          }
          ghostedCoarseNodesPerDir[dim] = lCoarseNodesPerDir[dim];
          // Check whether face *low needs ghost nodes
          if(ghostInterface[2*dim]) {ghostedCoarseNodesPerDir[dim] += 1;}
          // Check whether face *hi needs ghost nodes
          if(ghostInterface[2*dim + 1]) {ghostedCoarseNodesPerDir[dim] += 1;}
        }
        numGhostedCoarseNodes10 = ghostedCoarseNodesPerDir[1]*ghostedCoarseNodesPerDir[0];
        numGhostedCoarseNodes   = numGhostedCoarseNodes10*ghostedCoarseNodesPerDir[2];
        lNumGhostNodes = numGhostedCoarseNodes - lNumCoarseNodes;
      }

      void sortLocalLexicographicData() {
        std::sort(meshData.begin(), meshData.end(),
                  [](const std::vector<GO>& a, const std::vector<GO>& b)->bool {
                    // The below function sorts ranks by blockID, kmin, jmin and imin
                    if(a[2] < b[2]) {
                      return true;
                    } else if(a[2] == b[2]) {
                      if(a[7] < b[7]) {
                        return true;
                      } else if(a[7] == b[7]) {
                        if(a[5] < b[5]) {
                          return true;
                        } else if(a[5] == b[5]) {
                          if(a[3] < b[3]) {return true;}
                        }
                      }
                    }
                    return false;
                  });

        numBlocks = meshData[numRanks - 1][2] + 1;
        // Find the range of the current block
        myBlockStart = std::lower_bound(meshData.begin(), meshData.end(), myBlock - 1,
                                        [] (const std::vector<GO>& vec, const GO val)->bool {
                                          return (vec[2] < val) ? true : false;
                                        });
        myBlockEnd = std::upper_bound(meshData.begin(), meshData.end(), myBlock,
                                      [] (const GO val, const std::vector<GO>& vec)->bool {
                                        return (val < vec[2]) ? true : false;
                                      });
        // Assuming that i,j,k and ranges are split in pi, pj and pk processors
        // we search for these numbers as they will allow us to find quickly the PID of processors
        // owning ghost nodes.
        auto myKEnd = std::upper_bound(myBlockStart, myBlockEnd, (*myBlockStart)[3],
                                       [] (const GO val, const std::vector<GO>& vec)->bool {
                                         return (val < vec[7]) ? true : false;
                                       });
        auto myJEnd = std::upper_bound(myBlockStart, myKEnd, (*myBlockStart)[3],
                                       [] (const GO val, const std::vector<GO>& vec)->bool {
                                         return (val < vec[5]) ? true : false;
                                       });
        pi = std::distance(myBlockStart, myJEnd);
        pj = std::distance(myBlockStart, myKEnd) / pi;
        pk = std::distance(myBlockStart, myBlockEnd) / (pj*pi);

        // We also look for the index of the local rank in the current block.
        const int MyRank = myRank;
        myRankIndex = std::distance(meshData.begin(),
                                    std::find_if(myBlockStart, myBlockEnd,
                                                 [MyRank] (const std::vector<GO>& vec)->bool {
                                                   return (vec[0] == MyRank) ? true : false;
                                                 })
                                    );
        // We also construct a mapping of rank to rankIndex in the meshData vector,
        // this will allow us to access data quickly later on.
        rankIndices.resize(numRanks);
        for(int rankIndex = 0; rankIndex < numRanks; ++rankIndex) {
          rankIndices[meshData[rankIndex][0]] = rankIndex;
        }
      }

      void computeCoarseLocalLexicographicData() {
        coarseMeshData.resize(numRanks);
        Array<LO> rankOffset(3);
        for(int rank = 0; rank < numRanks; ++rank) {
          coarseMeshData[rank].resize(10);
          coarseMeshData[rank][0] = meshData[rank][0];
          coarseMeshData[rank][1] = meshData[rank][1];
          coarseMeshData[rank][2] = meshData[rank][2];
          for(int dim = 0; dim < 3; ++dim) {
            coarseMeshData[rank][3 + 2*dim] = meshData[rank][3 + 2*dim] / coarseRate[dim];
            if(meshData[rank][3 + 2*dim] % coarseRate[dim] > 0) {
              ++coarseMeshData[rank][3 + 2*dim];
            }
            coarseMeshData[rank][3 + 2*dim + 1] = meshData[rank][3 + 2*dim + 1] / coarseRate[dim];
            if(meshData[rank][3 + 2*dim + 1] == gFineNodesPerDir[dim] - 1 &&
               endRate[dim] < coarseRate[dim]) {
              ++coarseMeshData[rank][3 + 2*dim + 1];
            }
          }
          if(rank > 0) {
            coarseMeshData[rank][9] = coarseMeshData[rank - 1][9]
              + (coarseMeshData[rank - 1][8] - coarseMeshData[rank - 1][7] + 1)
              * (coarseMeshData[rank - 1][6] - coarseMeshData[rank - 1][5] + 1)
              * (coarseMeshData[rank - 1][4] - coarseMeshData[rank - 1][3] + 1);
          }
        }
      }

      void getFineNodeGlobalTuple(const GO myGID, GO& i, GO& j, GO& k) const {
        GO tmp;
        k   = myGID / gNumFineNodes10;
        tmp = myGID % gNumFineNodes10;
        j   = tmp / gFineNodesPerDir[0];
        i   = tmp % gFineNodesPerDir[0];
      }

      void getFineNodeLocalTuple(const LO myLID, LO& i, LO& j, LO& k) const {
        LO tmp;
        k   = myLID / lNumFineNodes10;
        tmp = myLID % lNumFineNodes10;
        j   = tmp / lFineNodesPerDir[0];
        i   = tmp % lFineNodesPerDir[0];
      }

      void getFineNodeGhostedTuple(const LO myLID, LO& i, LO& j, LO& k) const {
        LO tmp;
        k   = myLID / lNumFineNodes10;
        tmp = myLID % lNumFineNodes10;
        j   = tmp / lFineNodesPerDir[0];
        i   = tmp % lFineNodesPerDir[0];

        k += offsets[2];
        j += offsets[1];
        i += offsets[0];
      }

      void getFineNodeGID(const GO i, const GO j, const GO k, GO& myGID) const {
        myGID = k*gNumFineNodes10 + j*gFineNodesPerDir[0] + i;
      }

      void getFineNodeLID(const LO i, const LO j, const LO k, LO& myLID) const {
        myLID = k*lNumFineNodes10 + j*lFineNodesPerDir[0] + i;
      }

      void getCoarseNodeGlobalTuple(const GO myGID, GO& i, GO& j, GO& k) const {
        GO tmp;
        k   = myGID / gNumCoarseNodes10;
        tmp = myGID % gNumCoarseNodes10;
        j   = tmp / gCoarseNodesPerDir[0];
        i   = tmp % gCoarseNodesPerDir[0];
      }

      void getCoarseNodeLocalTuple(const LO myLID, LO& i, LO& j, LO& k) const {
        LO tmp;
        k   = myLID / lNumCoarseNodes10;
        tmp = myLID % lNumCoarseNodes10;
        j   = tmp / lCoarseNodesPerDir[0];
        i   = tmp % lCoarseNodesPerDir[0];
      }

      void getCoarseNodeGID(const GO i, const GO j, const GO k, GO& myGID) const {
        myGID = k*gNumCoarseNodes10 + j*gCoarseNodesPerDir[0] + i;
      }

      void getCoarseNodeLID(const LO i, const LO j, const LO k, LO& myLID) const {
        myLID = k*lNumCoarseNodes10 + j*lCoarseNodesPerDir[0] + i;
      }

      void getCoarseNodeGhostedLID(const LO i, const LO j, const LO k, LO& myLID) const {
        myLID = k*numGhostedCoarseNodes10 + j*ghostedCoarseNodesPerDir[0] + i;
      }

      void getCoarseNodeFineLID(const LO i, const LO j, const LO k, LO& myLID) const {
        myLID = 0;
        if(k*coarseRate[2] < lFineNodesPerDir[2]) {
          myLID += k*coarseRate[2]*lNumCoarseNodes10;
        } else {
          myLID += (lFineNodesPerDir[2] - 1)*lNumCoarseNodes10;
        }

        if(j*coarseRate[1] < lFineNodesPerDir[1]) {
          myLID += j*coarseRate[1]*lFineNodesPerDir[0];
        } else {
          myLID += (lFineNodesPerDir[1] - 1)*lFineNodesPerDir[1];
        }

        if(i*coarseRate[0] < lFineNodesPerDir[0]) {
          myLID += i*coarseRate[0];
        } else {
          myLID += lFineNodesPerDir[0] - 1;
        }
      }

      void getGhostedNodeFineLID(const LO i, const LO j, const LO k, LO& myLID) const {
        LO itmp = i - (offsets[0] > 0 ? 1 : 0);
        LO jtmp = j - (offsets[1] > 0 ? 1 : 0);
        LO ktmp = k - (offsets[2] > 0 ? 1 : 0);
        myLID = 0;
        if(ktmp*coarseRate[2] < lFineNodesPerDir[2]) {
          myLID += ktmp*coarseRate[2]*lNumCoarseNodes10;
        } else {
          myLID += (lFineNodesPerDir[2] - 1)*lNumCoarseNodes10;
        }

        if(jtmp*coarseRate[1] < lFineNodesPerDir[1]) {
          myLID += jtmp*coarseRate[1]*lFineNodesPerDir[0];
        } else {
          myLID += (lFineNodesPerDir[1] - 1)*lFineNodesPerDir[1];
        }

        if(itmp*coarseRate[0] < lFineNodesPerDir[0]) {
          myLID += itmp*coarseRate[0];
        } else {
          myLID += lFineNodesPerDir[0] - 1;
        }
      }

      void getGhostedNodeCoarseLID(const LO i, const LO j, const LO k, LO& myLID) const {
        LO itmp = i - (offsets[0] > 0 ? 1 : 0);
        LO jtmp = j - (offsets[1] > 0 ? 1 : 0);
        LO ktmp = k - (offsets[2] > 0 ? 1 : 0);
        myLID = ktmp*lNumCoarseNodes10 + jtmp*lCoarseNodesPerDir[0] + itmp;
      }

      void GetGIDLocalLexicographic(const LO iGhosted, const LO jGhosted, const LO kGhosted,
                                    const Array<LO> coarseNodeFineIndices,
                                    GO& myGID, LO& myPID, LO& myLID) const {

        LO ni = -1, nj = -1, li = -1, lj = -1, lk = -1;
        LO myRankGuess = myRankIndex;
        // We try to make a logical guess as to which PID owns the current coarse node
        if(iGhosted == 0 && ghostInterface[0]) {
          --myRankGuess;
        } else if((iGhosted == ghostedCoarseNodesPerDir[0] - 1) && ghostInterface[1]) {
          ++myRankGuess;
        }
        if(jGhosted == 0 && ghostInterface[2]) {
          myRankGuess -= pi;
        } else if((jGhosted == ghostedCoarseNodesPerDir[1] - 1) && ghostInterface[3]) {
          myRankGuess += pi;
        }
        if(kGhosted == 0 && ghostInterface[4]) {
          myRankGuess -= pj*pi;
        } else if((kGhosted == ghostedCoarseNodesPerDir[2] - 1) && ghostInterface[5]) {
          myRankGuess += pj*pi;
        }
        if(coarseNodeFineIndices[0] >= meshData[myRankGuess][3]
           && coarseNodeFineIndices[0] <= meshData[myRankGuess][4]
           && coarseNodeFineIndices[1] >= meshData[myRankGuess][5]
           && coarseNodeFineIndices[1] <= meshData[myRankGuess][6]
           && coarseNodeFineIndices[2] >= meshData[myRankGuess][7]
           && coarseNodeFineIndices[2] <= meshData[myRankGuess][8]
           && myRankGuess < numRanks - 1) {
          myPID = meshData[myRankGuess][0];
          ni = meshData[myRankGuess][4] - meshData[myRankGuess][3] + 1;
          nj = meshData[myRankGuess][6] - meshData[myRankGuess][5] + 1;
          li = coarseNodeFineIndices[0] - meshData[myRankGuess][3];
          lj = coarseNodeFineIndices[1] - meshData[myRankGuess][5];
          lk = coarseNodeFineIndices[2] - meshData[myRankGuess][7];
          myLID = lk*nj*ni + lj*ni + li;
          myGID = meshData[myRankGuess][9] + myLID;
        } else { // The guess failed, let us use the heavy artilery: std::find_if()
          // It could be interesting to monitor how many times this branch of the code gets
          // used as it is far more expensive than the above one...
          auto nodeRank = std::find_if(myBlockStart, myBlockEnd,
                                       [coarseNodeFineIndices](const std::vector<GO>& vec){
                                         if(coarseNodeFineIndices[0] >= vec[3]
                                            && coarseNodeFineIndices[0] <= vec[4]
                                            && coarseNodeFineIndices[1] >= vec[5]
                                            && coarseNodeFineIndices[1] <= vec[6]
                                            && coarseNodeFineIndices[2] >= vec[7]
                                            && coarseNodeFineIndices[2] <= vec[8]) {
                                           return true;
                                         } else {
                                           return false;
                                         }
                                       });
          myPID = (*nodeRank)[0];
          ni = (*nodeRank)[4] - (*nodeRank)[3] + 1;
          nj = (*nodeRank)[6] - (*nodeRank)[5] + 1;
          li = coarseNodeFineIndices[0] - (*nodeRank)[3];
          lj = coarseNodeFineIndices[1] - (*nodeRank)[5];
          lk = coarseNodeFineIndices[2] - (*nodeRank)[7];
          myLID = lk*nj*ni + lj*ni + li;
          myGID = (*nodeRank)[9] + myLID;
        }
      } // End GetGIDLocalLexicographic

    };

}; // class StructuredAggregationFactory

}

#define MUELU_STRUCTUREDAGGREGATIONFACTORY_SHORT
#endif /* MUELU_STRUCTUREDAGGREGATIONFACTORY_DECL_HPP_ */
