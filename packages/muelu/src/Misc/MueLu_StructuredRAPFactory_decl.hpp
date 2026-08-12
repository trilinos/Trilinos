// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef MUELU_STRUCTUREDRAPFACTORY_DECL_HPP
#define MUELU_STRUCTUREDRAPFACTORY_DECL_HPP

#include <string>
#include <vector>

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_CrsMatrix_fwd.hpp>
#include <Xpetra_MatrixFactory_fwd.hpp>
#include <Xpetra_MatrixUtils_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_StructuredRAPFactory_fwd.hpp"
#include "MueLu_RAPFactory_fwd.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_PerfUtils_fwd.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace MueLu {
/*!
  @class StructuredRAPFactory
  @brief Factory for building coarse matrices.
*/
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class StructuredRAPFactory : public TwoLevelFactoryBase {
#undef MUELU_STRUCTUREDRAPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  StructuredRAPFactory();

  virtual ~StructuredRAPFactory();

  //@}

  //! @name Input
  //@{

  RCP<const ParameterList> GetValidParameterList() const;

  void DeclareInput(Level& fineLevel, Level& coarseLevel) const;

  //@}

  //! @name Build methods.
  //@{
  void Build(Level& fineLevel, Level& coarseLevel) const;
  //@}

  //@{
  /*! @brief Add transfer factory in the end of list of transfer factories in RepartitionAcFactory.

  Transfer factories are derived from TwoLevelFactoryBase and project some data from the fine level to
  the next coarser level.
  */
  void AddTransferFactory(const RCP<const FactoryBase>& factory);

  // TODO add a function to remove a specific transfer factory?

  //! Returns number of transfer factories.
  size_t NumTransferFactories() const { return transferFacts_.size(); }

  //@}

 private:
  //@{

  mutable bool hasDeclaredInput_;

  //@}

  //@{

  //! list of user-defined transfer Factories
  std::vector<RCP<const FactoryBase> > transferFacts_;

  //@}

  //@{

  struct StencilOffset {
    int x;
    int y;
    int z;
  };

  struct StructuredGraphSpec {
    int numDimensions;
    LocalOrdinal dofsPerNode;
    std::vector<StencilOffset> stencilOffsets;
    std::string description;
  };

  StructuredGraphSpec GetStructuredGraphSpec(const std::string& matrixType, int interpolationOrder) const;

  void GetStructuredGraph(RCP<Matrix>& Ac, const RCP<Matrix> P,
                          const Teuchos::Array<LocalOrdinal>& lCoarseNodesPerDim,
                          const StructuredGraphSpec& graphSpec) const;

  void ConfigureRAPFactoryDelegate() const;

  mutable RCP<MueLu::RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>> rapFactoryDelegate_;

  //}

};  // class StructuredRAPFactory

}  // namespace MueLu

#define MUELU_STRUCTUREDRAPFACTORY_SHORT
#endif  // MUELU_STRUCTUREDRAPFACTORY_DECL_HPP
