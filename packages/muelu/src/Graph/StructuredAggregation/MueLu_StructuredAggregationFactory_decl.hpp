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
      int numDimensions;
      LO lNumFineNodes = -1, lNumCoarseNodes = -1, lNumGhostNodes = -1,lNumGhostedNodes  = -1;
      LO myBlock = -1, numBlocks = -1, lNumFineNodes10 = -1, lNumCoarseNodes10 = -1;
      LO numGhostedCoarseNodes = -1, numGhostedCoarseNodes10 = -1;
      GO gNumFineNodes = -1, gNumCoarseNodes = -1, gNumFineNodes10 = -1, gNumCoarseNodes10 = -1;
      GO minGlobalIndex = -1;
      Array<int> coarseRate, endRate;
      Array<LO> lFineNodesPerDir, lCoarseNodesPerDir, offsets, ghostedCoarseNodesPerDir;
      Array<GO> startIndices, gFineNodesPerDir, gCoarseNodesPerDir, startGhostedCoarseNode;
      std::vector<std::vector<GO> > meshData; // These are sorted later so they are in std::vector
      bool ghostInterface[6] = {false};

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
        LO itmp = i - offsets[0];
        LO jtmp = j - offsets[1];
        LO ktmp = k - offsets[2];
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
        LO itmp = i - offsets[0];
        LO jtmp = j - offsets[1];
        LO ktmp = k - offsets[2];
        myLID = k*lNumCoarseNodes10 + j*lCoarseNodesPerDir[0] + i;
      }

    };

    //! @name mesh layout handling methods.
    //@{

    /*! @brief assumes global lexicographic layout of the mesh to build aggregates */
    void GlobalLexicographicLayout(const RCP<const Map> coordMap, RCP<GeometricData> geoData,
                                   RCP<Aggregates> aggregates, std::vector<unsigned>& aggStat,
                                   LO& numNonAggregatedNodes) const;

    /*! @brief assumes local lexicographic layout of the mesh to build aggregates */
    void LoxalLexicographicLayout(const RCP<const Map> coordMap, RCP<GeometricData> geoData,
                                  RCP<Aggregates> aggregates) const;

    //@}

}; // class StructuredAggregationFactory

}

#define MUELU_STRUCTUREDAGGREGATIONFACTORY_SHORT
#endif /* MUELU_STRUCTUREDAGGREGATIONFACTORY_DECL_HPP_ */
