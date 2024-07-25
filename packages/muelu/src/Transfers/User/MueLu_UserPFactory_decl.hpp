// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_USERPFACTORY_DECL_HPP
#define MUELU_USERPFACTORY_DECL_HPP

#include <Xpetra_MultiVectorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_PerfUtils_fwd.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_UserPFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"

namespace MueLu {

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class UserPFactory : public PFactory {
#undef MUELU_USERPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor
  UserPFactory() {}

  //! Destructor.
  virtual ~UserPFactory() {}
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

};  // class UserPFactory

}  // namespace MueLu

#define MUELU_USERPFACTORY_SHORT
#endif  // MUELU_USERPFACTORY_DECL_HPP
