// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_REBALANCEBLOCKRESTRICTIONFACTORY_DECL_HPP_
#define MUELU_REBALANCEBLOCKRESTRICTIONFACTORY_DECL_HPP_

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MatrixFactory_fwd.hpp>
#include <Xpetra_MapExtractor_fwd.hpp>
#include <Xpetra_MapExtractorFactory_fwd.hpp>
#include "Xpetra_MultiVector_fwd.hpp"
#include "Xpetra_MultiVectorFactory_fwd.hpp"
#include "Xpetra_Vector_fwd.hpp"
#include "Xpetra_VectorFactory_fwd.hpp"
#include "Xpetra_Import_fwd.hpp"

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_PerfUtils_fwd.hpp"
#include "MueLu_RebalanceBlockRestrictionFactory_fwd.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_Types.hpp"

namespace MueLu {

/*!
  @class RebalanceBlockRestrictionFactory class.
  @brief Applies permutation to restriction operators.
  @ingroup MueLuTransferClasses
*/
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class RebalanceBlockRestrictionFactory : public TwoLevelFactoryBase {
#undef MUELU_REBALANCEBLOCKRESTRICTIONFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  RebalanceBlockRestrictionFactory()
    : UseSingleSourceImporters_(false) {}

  //! Destructor.
  virtual ~RebalanceBlockRestrictionFactory() {}

  RCP<const ParameterList> GetValidParameterList() const;

  //@}

  //! @name Input
  //@{

  /*! @brief Specifies the data that this class needs, and the factories that generate that data.

      If the Build method of this class requires some data, but the generating factory is not specified in DeclareInput, then this class
      will fall back to the settings in FactoryManager.
  */
  void DeclareInput(Level &fineLevel, Level &coarseLevel) const;

  //! Add a factory manager
  void AddFactoryManager(RCP<const FactoryManagerBase> FactManager);

  //! Get the importers from a single source
  void SetSingleSourceImporters(bool value) { UseSingleSourceImporters_ = value; }

  //@}

  //! @name Build methods.
  //@{

  //! Build an object with this factory.
  void Build(Level &fineLevel, Level &coarseLevel) const;

  //@}

 private:
  //! Input factories
  std::vector<Teuchos::RCP<const FactoryManagerBase> > FactManager_;

  bool UseSingleSourceImporters_;

};  // class RebalanceBlockTransferFactory

}  // namespace MueLu

#define MUELU_REBALANCEBLOCKRESTRICTIONFACTORY_SHORT
#endif /* MUELU_REBALANCEBLOCKRESTRICTIONFACTORY_DECL_HPP_ */
