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

#ifndef MUELU_SUBBLOCKAFACTORY_DECL_HPP_
#define MUELU_SUBBLOCKAFACTORY_DECL_HPP_

#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_MapExtractor_fwd.hpp>
#include <Xpetra_StridedMap_fwd.hpp>
#include <Xpetra_StridedMapFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_SubBlockAFactory_fwd.hpp"

namespace MueLu {

/*!
  @class SubBlockAFactory class.
  @brief Factory for building a thresholded operator.

  This is a very simple class to access a single matrix block in a blocked operator A.

  Example
  \code{.cpp}
  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> > bOp = Teuchos::rcp(new Xpetra::BlockedCrsMatrix<Scalar,LO,GO>(mapExtractor,mapExtractor,10));
  // ... let bOp be a 2x2 blocked operator ...
  bOp->fillComplete();

  // define factory for accessing block (0,0) in blocked operator A (assuming that the blocked operator is stored in Level class with NoFactory as generating factory)
  RCP<SubBlockAFactory> A11Fact = Teuchos::rcp(new SubBlockAFactory());
  A11Fact->SetFactory("A", MueLu::NoFactory::getRCP());
  A11Fact->SetParameter("block row", 0);
  A11Fact->SetParameter("block col", 0);

  // define factory for accessing block (1,1) in blocked operator A
  RCP<SubBlockAFactory> A22Fact = Teuchos::rcp(new SubBlockAFactory());
  A22Fact->SetFactory("A", MueLu::NoFactory::getRCP());
  A22Fact->SetParameter("block row", 1);
  A22Fact->SetParameter("block col", 1);

  RCP<Matrix> A11 = level.Get<RCP<Matrix> >("A", A11Fact); // extract (0,0) block from blocked operator A
  RCP<Matrix> A22 = level.Get<RCP<Matrix> >("A", A22Fact); // extract (1,1) block from blocked operator A
  \endcode
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class SubBlockAFactory : public SingleLevelFactoryBase {
#undef MUELU_SUBBLOCKAFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  SubBlockAFactory() = default;

  //! Destructor.
  virtual ~SubBlockAFactory() = default;
  //@}

  //! Input
  //@{

  RCP<const ParameterList> GetValidParameterList() const override;

  void DeclareInput(Level& currentLevel) const override;

  //@}

  //@{
  //! @name Build methods.

  /*! @brief Build an object with this factory.
   *
   * Extract sub block matrix from a given blocked crs operator.
   * Strided or block information is extracted in the following way:
   *   1) first check whether the corresponding sub maps are strided
   *      If yes, use the fixed block size and strided block id
   *   2) If no, get the full map of the map extractors.
   *      Check whether the full map is strided.
   *      If yes, use the strided information of the full map and build
   *      partial (strided) maps with it
   *      If no, throw an exception
   *
   * For blocked operators with block maps one should use the striding
   * information from the sub maps. For strided operators, the striding
   * information of the full map is the best choice.
   */
  void Build(Level& currentLevel) const override;

  //@}

 private:
  bool CheckForUserSpecifiedBlockInfo(bool bRange, std::vector<size_t>& stridingInfo, LocalOrdinal& stridedBlockId) const;

};  // class SubBlockAFactory

}  // namespace MueLu

#define MUELU_SUBBLOCKAFACTORY_SHORT
#endif /* MUELU_SUBBLOCKAFACTORY_DECL_HPP_ */
