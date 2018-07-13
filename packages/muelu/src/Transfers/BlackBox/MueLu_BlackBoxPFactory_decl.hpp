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
#ifndef MUELU_BLACKBOXPFACTORY_DECL_HPP
#define MUELU_BLACKBOXPFACTORY_DECL_HPP

#include <Teuchos_SerialDenseVector.hpp>

#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_BlackBoxPFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"

namespace MueLuTests {
  // Forward declaration of friend tester class used to UnitTest BlackBoxPFactory
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class BlackBoxPFactoryTester;
}

namespace MueLu {

/*!
  @class BlackBoxPFactory
  @ingroup MueLuTransferClasses
  @brief Prolongator factory performing geometric coarsening.

  The geometric algorithm assumes the underlying mesh is reasonably structured.
  Any rate of coarsening can be applied, and the rate is automatically decrease along
  an edge if the number of element is not divisible by the coarsening rate.
  The coarsening rate is allowed to be different in all direction which means that
  semi-coarsening can be achieved within this algorithm in 1 or 2 directions.
  The main difficulty is to obtain the number of elements/nodes in each directions
  to identify coarse nodes and fine nodes.

  ## Input/output of BlackBoxPFactory ##

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
  The * in the @c validated column means that the parameter is declared in the list of valid input parameters (see BlackBoxPFactory::GetValidParameters).<br>
  The * in the @c requested column states that the data is requested as input with all dependencies (see BlackBoxPFactory::DeclareInput).

  ### Variables provided by BlackBoxPFactory ###

  After BlackBoxPFactory::Build the following data is available (if requested)

  | Parameter         | generated by             | description
  |-------------------|--------------------------|------------------------------------------------------------------------------------------------------------------|
  | P                 | BlackBoxPFactory         | Prolongator                                                                                                      |
  | Nullspace         | BlackBoxPFactory         | Coarse nullspace (the fine level nullspace information is coarsened using P to generate a coarse version         |
  |                   |                          | of the nullspace. No scaling is applied.                                                                         |
  | coarseCoordinates | BlackBoxPFactory         | coarseCoordinates are to be picked up on the coarse level by the coordinates transfer factory and renamed        |
  |                   |                          | coordinates so that on coarser levels coordinates are available in case another factory needs them.              |

*/
  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class BlackBoxPFactory : public PFactory {
#undef MUELU_BLACKBOXPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    friend class MueLuTests::BlackBoxPFactoryTester<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

    //! @name Constructors/Destructors.
    //@{

    //! Constructor
    BlackBoxPFactory() { }

    //! Destructor.
    virtual ~BlackBoxPFactory() { }
    //@}

    RCP<const ParameterList> GetValidParameterList() const;

    //! Input
    //@{

    void DeclareInput(Level& fineLevel, Level& coarseLevel) const;

    //@}

    //! @name Build methods.
    //@{

    void Build (Level& fineLevel, Level& coarseLevel) const;
    void BuildP(Level& fineLevel, Level& coarseLevel) const;

    //@}

  private:

    struct NodesIDs {
      // This small struct just carries basic data associated with coarse nodes that is needed
      // to compute colMapP and to fillComplete P,

      Array<GO>  GIDs, coarseGIDs;
      Array<int> PIDs;
      Array<LO>  LIDs;
      std::vector<GO> colInds;
    };

    struct NodeID {
      // This small struct is similar to the one above but only for one node.
      // It is used to create a vector of NodeID that can easily be sorted

      GO  GID;
      int PID;
      LO  LID, lexiInd;
    };

    void GetGeometricData(RCP<Xpetra::MultiVector<double,LO,GO,NO> >& coordinates,
                          const Array<LO> coarseRate, const Array<GO> gFineNodesPerDir,
                          const Array<LO> lFineNodesPerDir, const LO BlkSize, Array<GO>& gIndices,
                          Array<LO>& myOffset, Array<bool>& ghostInterface, Array<LO>& endRate,
                          Array<GO>& gCoarseNodesPerDir, Array<LO>& lCoarseNodesPerDir,
                          Array<LO>& glCoarseNodesPerDir, Array<GO>& ghostGIDs,
                          Array<GO>& coarseNodesGIDs, Array<GO>& colGIDs, GO& gNumCoarseNodes,
                          LO& lNumCoarseNodes, ArrayRCP<Array<double> > coarseNodes,
                          Array<int>& boundaryFlags, RCP<NodesIDs> ghostedCoarseNodes) const;

    void ComputeLocalEntries(const RCP<const Matrix>& Aghost, const Array<LO> coarseRate,
                             const Array<LO> endRate, const LO BlkSize, const Array<LO> elemInds,
                             const Array<LO> lCoarseElementsPerDir,
                             const LO numDimensions, const Array<LO> lFineNodesPerDir,
                             const Array<GO> gFineNodesPerDir, const Array<GO> gIndices,
                             const Array<LO> lCoarseNodesPerDir, const Array<bool> ghostInterface,
                             const Array<int> elementFlags, const std::string stencilType,
                             const std::string blockStrategy, const Array<LO> elementNodesPerDir,
                             const LO numNodesInElement, const Array<GO> colGIDs,
                             Teuchos::SerialDenseMatrix<LO,SC>& Pi,
                             Teuchos::SerialDenseMatrix<LO,SC>& Pf,
                             Teuchos::SerialDenseMatrix<LO,SC>& Pe,
                             Array<LO>& dofType, Array<LO>& lDofInd) const;

    void CollapseStencil(const int type, const int orientation, const int collapseFlags[3],
                         Array<SC>& stencil) const ;

    void FormatStencil(const LO BlkSize, const Array<bool> ghostInterface, const LO ie,
                       const LO je, const LO ke, const ArrayView<const SC> rowValues,
                       const Array<LO> elementNodesPerDir, const int collapseFlags[3],
                       const std::string stencilType, Array<SC>& stencil) const;

    void GetNodeInfo(const LO ie, const LO je, const LO ke, const Array<LO> elementNodesPerDir,
                     int* type, LO& ind, int* orientation) const;

    void sh_sort_permute(
                const typename Teuchos::Array<LocalOrdinal>::iterator& first1,
                const typename Teuchos::Array<LocalOrdinal>::iterator& last1,
                const typename Teuchos::Array<LocalOrdinal>::iterator& first2,
                const typename Teuchos::Array<LocalOrdinal>::iterator& last2) const;

  }; //class BlackBoxPFactory

} //namespace MueLu

#define MUELU_BLACKBOXPFACTORY_SHORT
#endif // MUELU_BLACKBOXPFACTORY_DECL_HPP
