// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_MERGEDBLOCKEDMATRIXFACTORY_DECL_HPP_
#define MUELU_MERGEDBLOCKEDMATRIXFACTORY_DECL_HPP_

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_CrsMatrix_fwd.hpp>
#include <Xpetra_BlockedCrsMatrix_fwd.hpp>
#include <Xpetra_MatrixFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"

namespace MueLu {
/*!
  @class MergedBlockedMatrix
  @brief Factory provides a merged version of a blocked matrix
*/
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class MergedBlockedMatrixFactory : public SingleLevelFactoryBase {
#undef MUELU_MERGEDBLOCKEDMATRIXFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  MergedBlockedMatrixFactory();

  virtual ~MergedBlockedMatrixFactory() {}
  //@}

  //! @name Input
  //@{

  RCP<const ParameterList> GetValidParameterList() const;

  void DeclareInput(Level &currentLevel) const;

  //@}

  //! @name Build methods.
  //@{
  void Build(Level &currentLevel) const;
  //@}

 private:
};  // class MergedBlockedMatrixFactory

}  // namespace MueLu

#define MUELU_MERGEDBLOCKEDMATRIXFACTORY_SHORT

#endif /* MUELU_MERGEDBLOCKEDMATRIXFACTORY_DECL_HPP_ */
