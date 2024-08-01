// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_RAPFACTORY_DECL_HPP
#define MUELU_RAPFACTORY_DECL_HPP

#include <string>

#include <Xpetra_MatrixFactory_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MatrixUtils_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_RAPFactory_fwd.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_PerfUtils_fwd.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace MueLu {
/*!
  @class RAPFactory
  @brief Factory for building coarse matrices.
*/
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class RAPFactory : public TwoLevelFactoryBase {
#undef MUELU_RAPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  RAPFactory();

  virtual ~RAPFactory();

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

};  // class RAPFactory

}  // namespace MueLu

#define MUELU_RAPFACTORY_SHORT
#endif  // MUELU_RAPFACTORY_DECL_HPP
