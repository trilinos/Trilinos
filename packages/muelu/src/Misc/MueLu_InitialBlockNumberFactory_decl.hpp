// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_INITIALBLOCKNUMBER_FACTORY_DECL_HPP
#define MUELU_INITIALBLOCKNUMBER_FACTORY_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "Xpetra_Matrix_fwd.hpp"
#include "Xpetra_Vector_fwd.hpp"
#include "Xpetra_VectorFactory_fwd.hpp"

#include "MueLu_InitialBlockNumberFactory_fwd.hpp"

namespace MueLu {

/*!
  @class InitialBlockNumberFactory class.
  @brief Class for generating an initial LocalOrdinal-type BlockNumber vector, based on an input paraemter for interleaved dofs.


*/
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class InitialBlockNumberFactory : public SingleLevelFactoryBase {
#undef MUELU_INITIALBLOCKNUMBERFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.

  /*! @brief Constructor.
   */
  InitialBlockNumberFactory() {}

  //! Destructor.
  virtual ~InitialBlockNumberFactory() {}

  RCP<const ParameterList> GetValidParameterList() const;

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

 private:
};  // class InitialBlockNumberFactory

}  // namespace MueLu

#define MUELU_INITIALBLOCKNUMBERFACTORY_SHORT
#endif  // MUELU_INITIALBLOCKNUMBER_FACTORY_DECL_HPP
