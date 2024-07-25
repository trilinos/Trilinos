// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_RIGIDBODYMODEFACTORY_DECL_HPP
#define MUELU_RIGIDBODYMODEFACTORY_DECL_HPP

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_RigidBodyModeFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"

namespace MueLu {

/*!
  @class RigidBodyModeFactory class.
  @brief Nullspace Factory for building rigid body modes.
  @ingroup MueLuTransferClasses
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class RigidBodyModeFactory : public SingleLevelFactoryBase {
#undef MUELU_RIGIDBODYMODEFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor
  RigidBodyModeFactory(const int numPDEs)
    : nspName_("Nullspace")
    , numPDEs_(numPDEs) {}
  //! Constructor
  RigidBodyModeFactory(const std::string &nspName = "Nullspace")
    : nspName_(nspName)
    , numPDEs_(3) {}

  //! Destructor.
  virtual ~RigidBodyModeFactory();

  //@}

  //! @name Input
  //@{

  /*! @brief Specifies the data that this class needs, and the factories that generate that data.

      If the Build method of this class requires some data, but the generating factory is not specified in DeclareInput, then this class
      will fall back to the settings in FactoryManager.
  */
  void DeclareInput(Level &currentLevel) const;

  //@}

  //! @name Build methods.
  //@{

  //! Build an object with this factory.
  void Build(Level &currentLevel) const;

  //@}
  void setNumPDEs(int numPDEs) {
    numPDEs_ = numPDEs;
  }

 private:
  //! name of nullspace vector on finest level
  std::string nspName_;

  int numPDEs_;

};  // class RigidBodyModeFactory

}  // namespace MueLu

#define MUELU_RIGIDBODYMODEFACTORY_SHORT
#endif  // MUELU_RIGIDBODYMODEFACTORY_DECL_HPP
