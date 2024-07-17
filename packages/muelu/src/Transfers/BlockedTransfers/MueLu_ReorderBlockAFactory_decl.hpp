// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_REORDERBLOCKAFACTORY_DECL_HPP_
#define MUELU_REORDERBLOCKAFACTORY_DECL_HPP_

#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_MapExtractor_fwd.hpp>
#include <Xpetra_StridedMap_fwd.hpp>
#include <Xpetra_StridedMapFactory_fwd.hpp>
#include "Xpetra_ReorderedBlockedCrsMatrix_fwd.hpp"

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"

namespace MueLu {

/*!
  @class ReorderBlockAFactory class.
  @brief Factory for building a reordered (nested) block operator

  Example
  \code
  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> > bOp = ...;
  bOp->fillComplete();

  TODO
  \endcode
*/

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class ReorderBlockAFactory : public SingleLevelFactoryBase {
#undef MUELU_REORDERBLOCKAFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! Input
  //@{

  RCP<const ParameterList> GetValidParameterList() const;

  void DeclareInput(Level &currentLevel) const;

  //@}

  //@{
  //! @name Build methods.

  /*! @brief Build an object with this factory.
   */
  void Build(Level &currentLevel) const;

  //@}
};  // class ReorderBlockAFactory

}  // namespace MueLu

#define MUELU_REORDERBLOCKAFACTORY_SHORT
#endif /* MUELU_REORDERBLOCKAFACTORY_DECL_HPP_ */
