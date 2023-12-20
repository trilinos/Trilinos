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
#ifndef MUELU_GENERALGEOMETRICPFACTORY_DECL_HPP
#define MUELU_GENERALGEOMETRICPFACTORY_DECL_HPP

#include <Teuchos_SerialDenseVector.hpp>

#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_GeneralGeometricPFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"

namespace MueLuTests {
// Forward declaration of friend tester class used to UnitTest GeneralGeometricPFactory
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class GeneralGeometricPFactoryTester;
}  // namespace MueLuTests

namespace MueLu {

/*!
  @class GeneralGeometricPFactory
  @ingroup MueLuTransferClasses
  @brief Prolongator factory performing geometric coarsening.

  The geometric algorithm assumes the underlying mesh is reasonably structured.
  Any rate of coarsening can be applied, and the rate is automatically decrease along
  an edge if the number of element is not divisible by the coarsening rate.
  The coarsening rate is allowed to be different in all direction which means that
  semi-coarsening can be achieved within this algorithm in 1 or 2 directions.
  The main difficulty is to obtain the number of elements/nodes in each directions
  to identify coarse nodes and fine nodes.

  ## Input/output of GeneralGeometricPFactory ##

  ### User parameters of SemiCoarsenPFactory ###
  | Parameter   | type    | default   | master.xml | validated | requested | description                                                                      |
  | ------------|---------|-----------|:----------:|:---------:|:---------:|----------------------------------------------------------------------------------|
  | Coarsen     | string  | -1        |            |           |           | A string that specify the coarsening rate, if it is a single character, it will  |
  |             |         |           |            |           |           | indicate a unique coarsening rate in each direction, if it is longer, it will be |
  |             |         |           |            |           |           | processed as a vector with 3 entries, one for each spatial direction             |
  | A           | Factory | null      |            | *         | *         | Generating factory of the matrix A used during the prolongator smoothing process |
  | Nullspace   | Factory | null      |            | *         | *         | Generating factory of the nullspace. The GeneralGeometricPFactory provides       |
  |             |         |           |            |           |           | a coarse version of the given Nullspace.                                         |
  | Coordinates | Factory | NoFactory |            | *         | *         | Generating factory for coorindates. The coordinates are expected to be provided  |
  |             |         |           |            |           |           | on the finest level using the NoFactory mechanism. The coordinates are used to   |
  |             |         |           |            |           |           | compute the coarsening stencil and coarse coordinates are generated for the next |
  |             |         |           |            |           |           | level.                                                                           |



  The * in the @c master.xml column denotes that the parameter is defined in the @c master.xml file.<br>
  The * in the @c validated column means that the parameter is declared in the list of valid input parameters (see GeneralGeometricCoarsenPFactory::GetValidParameters).<br>
  The * in the @c requested column states that the data is requested as input with all dependencies (see GeneralGeometricCoarsenPFactory::DeclareInput).

  ### Variables provided by GeneralGeometricPFactory ###

  After GeneralGeometricPFactory::Build the following data is available (if requested)

  | Parameter         | generated by             | description
  |-------------------|--------------------------|------------------------------------------------------------------------------------------------------------------|
  | P                 | GeneralGeometricPFactory | Prolongator                                                                                                      |
  | Nullspace         | GeneralGeometricPFactory | Coarse nullspace (the fine level nullspace information is coarsened using P to generate a coarse version         |
  |                   |                          | of the nullspace. No scaling is applied.                                                                         |
  | coarseCoordinates | NoFactory                | Coarse coordinates that will be used on the next fine level to compute the coarsening stencils                   |

*/
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class GeneralGeometricPFactory : public PFactory {
#undef MUELU_GENERALGEOMETRICPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  friend class MueLuTests::GeneralGeometricPFactoryTester<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  //! @name Constructors/Destructors.
  //@{

  //! Constructor
  GeneralGeometricPFactory() {}

  //! Destructor.
  virtual ~GeneralGeometricPFactory() {}
  //@}

  RCP<const ParameterList> GetValidParameterList() const;

  //! Input
  //@{

  void DeclareInput(Level& fineLevel, Level& coarseLevel) const;

  //@}

  //! @name Build methods.
  //@{

  void Build(Level& fineLevel, Level& coarseLevel) const;
  void BuildP(Level& fineLevel, Level& coarseLevel) const;

  //@}

 private:
  struct GeometricData {
    // Geometric algorithm require a copious amount of data to be passed around so this struct
    // will reduce the amount of input/output parameters of methods in the class. Additionally
    // the struct can be rewritten to accomodate constraints of Kokkos/CUDA data types

    std::string meshLayout = "Global Lexicographic";
    int numDimensions;
    LO lNumFineNodes = -1, lNumCoarseNodes = -1, lNumGhostNodes = -1, lNumGhostedNodes = -1;
    LO myBlock = -1, numBlocks = -1, lNumFineNodes10 = -1;
    GO gNumFineNodes = -1, gNumCoarseNodes = -1, gNumFineNodes10 = -1, minGlobalIndex = -1;
    Array<int> coarseRate, endRate;
    Array<LO> lFineNodesPerDir, lCoarseNodesPerDir, offsets, ghostedCoarseNodesPerDir;
    Array<GO> startIndices, gFineNodesPerDir, gCoarseNodesPerDir, startGhostedCoarseNode;
    std::vector<std::vector<GO> > meshData;  // These are sorted later so they are in std::vector
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
  };

  struct NodesIDs {
    // This small struct just carries basic data associated with coarse nodes that is needed
    // to compute colMapP and to fillComplete P,

    Array<GO> GIDs, coarseGIDs;
    Array<int> PIDs;
    Array<LO> LIDs;
    std::vector<GO> colInds;
  };

  struct NodeID {
    // This small struct is similar to the one above but only for one node.
    // It is used to create a vector of NodeID that can easily be sorted

    GO GID;
    int PID;
    LO LID, lexiInd;
  };

  void MeshLayoutInterface(const int interpolationOrder, const LO blkSize,
                           RCP<const Map> fineCoordsMap, RCP<GeometricData> myGeometry,
                           RCP<NodesIDs> ghostedCoarseNodes,
                           Array<Array<GO> >& lCoarseNodesGIDs) const;

  void GetCoarsePoints(const int interpolationOrder, const LO blkSize,
                       RCP<const Map> fineCoordsMap, RCP<GeometricData> myGeometry,
                       RCP<NodesIDs> ghostedCoarseNodes,
                       Array<Array<GO> >& lCoarseNodesGIDs) const;

  void MakeGeneralGeometricP(RCP<GeometricData> myGeo,
                             const RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LO, GO, NO> >& fCoords,
                             const LO nnzP, const LO dofsPerNode,
                             RCP<const Map>& stridedDomainMapP,
                             RCP<Matrix>& Amat, RCP<Matrix>& P,
                             RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LO, GO, NO> >& cCoords,
                             RCP<NodesIDs> ghostedCoarseNodes, Array<Array<GO> > coarseNodesGIDs,
                             int interpolationOrder) const;

  void ComputeStencil(const LO numDimension, const Array<GO> currentNodeIndices,
                      const Array<GO> coarseNodeIndices, const LO rate[3],
                      const Array<Array<typename Teuchos::ScalarTraits<Scalar>::coordinateType> > coord, const int interpolationOrder,
                      std::vector<double>& stencil) const;

  void ComputeConstantInterpolationStencil(const LO numDimension,
                                           const Array<GO> currentNodeIndices,
                                           const Array<GO> coarseNodeIndices,
                                           const LO rate[3], std::vector<double>& stencil) const;

  void ComputeLinearInterpolationStencil(const LO numDimension, const Array<Array<typename Teuchos::ScalarTraits<Scalar>::coordinateType> > coord,
                                         std::vector<double>& stencil) const;
  void GetInterpolationFunctions(const LO numDimension,
                                 const Teuchos::SerialDenseVector<LO, double> parameters,
                                 double functions[4][8]) const;

  void sh_sort_permute(
      const typename Teuchos::Array<LocalOrdinal>::iterator& first1,
      const typename Teuchos::Array<LocalOrdinal>::iterator& last1,
      const typename Teuchos::Array<LocalOrdinal>::iterator& first2,
      const typename Teuchos::Array<LocalOrdinal>::iterator& last2) const;

  void sh_sort2(
      const typename Teuchos::Array<LocalOrdinal>::iterator& first1,
      const typename Teuchos::Array<LocalOrdinal>::iterator& last1,
      const typename Teuchos::Array<LocalOrdinal>::iterator& first2,
      const typename Teuchos::Array<LocalOrdinal>::iterator& last2) const;

  void GetGIDLocalLexicographic(const GO i, const GO j, const GO k,
                                const Array<LO> coarseNodeFineIndices,
                                const RCP<GeometricData> myGeo, const LO myRankIndex, const LO pi,
                                const LO pj, const LO pk,
                                const typename std::vector<std::vector<GO> >::iterator blockStart,
                                const typename std::vector<std::vector<GO> >::iterator blockEnd,
                                GO& myGID, LO& myPID, LO& myLID) const;

};  // class GeneralGeometricPFactory

}  // namespace MueLu

#define MUELU_GENERALGEOMETRICPFACTORY_SHORT
#endif  // MUELU_GENERALGEOMETRICPFACTORY_DECL_HPP
