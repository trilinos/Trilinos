// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_MUELU_FACADECLASSFACTORY_DECL_HPP_
#define PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_MUELU_FACADECLASSFACTORY_DECL_HPP_

#include <Teuchos_ParameterList.hpp>

#include "MueLu_BaseClass.hpp"
#include "MueLu_FacadeClassBase.hpp"

#include "MueLu_ConfigDefs.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class FacadeClassFactory
  : public virtual BaseClass {
#undef MUELU_FACADECLASSFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors
  //@{

  //! Constructor.
  FacadeClassFactory();

  //! Destructor.
  virtual ~FacadeClassFactory() {}

  //@}

  /*! @brief Set parameter list for FacadeClassFactory interpreter.

     @param[in] paramList: ParameterList containing the MueLu parameters.
  */
  Teuchos::RCP<Teuchos::ParameterList> SetParameterList(const Teuchos::ParameterList& paramList);

  /*! @brief Register new facade class
   *
   * Register new externally provided facade class in FacadeClassFactory
   *
   * @param[in] name: name that is used to access Facade class
   * @param[in] facadeclass: RCP pointer to facade class instance
   */
  void RegisterFacadeClass(std::string name, Teuchos::RCP<FacadeClassBase> facadeclass) {
    facadeClasses_[name] = facadeclass;
  }

 private:
  std::map<std::string, Teuchos::RCP<FacadeClassBase> > facadeClasses_;
};

}  // namespace MueLu

#define MUELU_FACADECLASSFACTORY_SHORT

#endif /* PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_MUELU_FACADECLASSFACTORY_DECL_HPP_ */
