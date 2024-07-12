// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PACKAGES_MUELU_SRC_MISC_MUELU_FINELEVELINPUTDATAFACTORY_DECL_HPP_
#define PACKAGES_MUELU_SRC_MISC_MUELU_FINELEVELINPUTDATAFACTORY_DECL_HPP_

#include <string>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_FineLevelInputDataFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"

#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_SmootherPrototype_fwd.hpp"
#include "MueLu_SmootherBase_fwd.hpp"

namespace MueLuTests {
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class FineLevelInputDataFactoryTester;
}

namespace MueLu {

/*!
  @class FineLevelInputData class.
  @brief Factory for piping in input data from the finest level into the MueLu data dependency system
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class FineLevelInputDataFactory : public SingleLevelFactoryBase {
  friend class MueLuTests::FineLevelInputDataFactoryTester<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
#undef MUELU_FINELEVELINPUTDATAFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  FineLevelInputDataFactory() {}

  //! Destructor.
  virtual ~FineLevelInputDataFactory() {}

  RCP<const ParameterList> GetValidParameterList() const;

  //@}

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const;

  //@}

  //! @name Build methods.
  //@{

  /*!
    @brief Build method.
    */
  void Build(Level& currentLevel) const;

  //@}
 private:
  void test() const { std::cout << "TEST" << std::endl; }

};  // class FineLevelInputDataFactory

}  // namespace MueLu

#define MUELU_FINELEVELINPUTDATAFACTORY_SHORT

#endif /* PACKAGES_MUELU_SRC_MISC_MUELU_FINELEVELINPUTDATAFACTORY_DECL_HPP_ */
