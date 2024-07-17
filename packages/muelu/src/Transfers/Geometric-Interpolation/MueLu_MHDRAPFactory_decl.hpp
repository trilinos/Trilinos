// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_MHDRAPFACTORY_DECL_HPP
#define MUELU_MHDRAPFACTORY_DECL_HPP

#include <string>

#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
//#include "MueLu_MHDRAPFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"

namespace MueLu {
/*!
  @class MHDRAPFactory
  @brief Factory for building coarse matrices.
*/
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class MHDRAPFactory : public TwoLevelFactoryBase {
#undef MUELU_MHDRAPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  MHDRAPFactory();

  virtual ~MHDRAPFactory() {}

  //@}

  //! @name Input
  //@{

  void DeclareInput(Level &fineLevel, Level &coarseLevel) const;

  //@}

  //! @name Build methods.
  //@{
  void Build(Level &fineLevel, Level &coarseLevel) const;
  //@}

  //! @name Handling of user-defined transfer factories
  //@{

  //! Indicate that the restriction operator action should be implicitly defined by the transpose of the prolongator.
  void SetImplicitTranspose(bool const &implicit) {
    implicitTranspose_ = implicit;
  }

  //@}

  //! @name internal print methods.
  //    static std::string PerfUtils::PrintMatrixInfo(const Matrix & Ac, const std::string & msgTag);

  static std::string PrintLoadBalancingInfo(const Matrix &Ac, const std::string &msgTag);

 private:
  //! If true, the action of the restriction operator action is implicitly defined by the transpose of the prolongator.
  bool implicitTranspose_;

};  // class MHDRAPFactory

}  // namespace MueLu

#define MUELU_MHDRAPFACTORY_SHORT
#endif  // MUELU_MHDRAPFACTORY_DECL_HPP
