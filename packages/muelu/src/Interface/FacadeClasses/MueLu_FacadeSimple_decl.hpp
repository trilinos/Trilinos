// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_Simple_DECL_HPP_
#define PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_Simple_DECL_HPP_

#include <Teuchos_ParameterList.hpp>

#include "MueLu_FacadeClassBase.hpp"

#include "MueLu_ConfigDefs.hpp"

namespace MueLu {

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class FacadeSimple : public FacadeClassBase<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors
  //@{

  //! Constructor.
  FacadeSimple();

  //! Destructor.
  virtual ~FacadeSimple() {}

  //@}

  /*! @brief Set parameter list for FacadeClass interpreter.

      @param[in] paramList: ParameterList containing the MueLu parameters for chosen facade class.
  */
  Teuchos::RCP<Teuchos::ParameterList> SetParameterList(const Teuchos::ParameterList& paramList);

 private:
};

}  // namespace MueLu

#endif /* PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_Simple_DECL_HPP_ */
