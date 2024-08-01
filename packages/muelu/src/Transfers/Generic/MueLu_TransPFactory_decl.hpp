// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_TRANSPFACTORY_DECL_HPP
#define MUELU_TRANSPFACTORY_DECL_HPP

#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TransPFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_PerfUtils_fwd.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace MueLu {

/*!
  @class TransPFactory class.
  @brief Factory for building restriction operators.

  This factory currently depends on an underlying matrix-matrix multiply with the identity
  matrix to do the transpose.  This should probably be fixed at some point.
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class TransPFactory : public TwoLevelFactoryBase {
#undef MUELU_TRANSPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  TransPFactory() {}

  //! Destructor.
  virtual ~TransPFactory() {}

  RCP<const ParameterList> GetValidParameterList() const;

  //@}

  //! Input
  //@{

  void DeclareInput(Level &fineLevel, Level &coarseLevel) const;

  //@}

  //! @name Build methods.
  //@{

  void Build(Level &fineLevel, Level &coarseLevel) const;

  //@}

};  // class TransPFactory

}  // namespace MueLu

#define MUELU_TRANSPFACTORY_SHORT
#endif  // MUELU_TRANSPFACTORY_DECL_HPP
