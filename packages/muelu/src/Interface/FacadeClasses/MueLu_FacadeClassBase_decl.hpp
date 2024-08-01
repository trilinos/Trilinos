// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_MUELU_FACADECLASSBASE_DECL_HPP_
#define PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_MUELU_FACADECLASSBASE_DECL_HPP_

#include "MueLu_BaseClass.hpp"

namespace MueLu {

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class FacadeClassBase
  : public virtual BaseClass {
#undef MUELU_FACADECLASSBASE_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors
  //@{

  //! Constructor.
  FacadeClassBase();

  //! Destructor.
  virtual ~FacadeClassBase() {}

  //@}

  /*! @brief Set parameter list for FacadeClass (abstract member).

      @param[in] paramList: ParameterList containing the MueLu parameters.
  */
  virtual Teuchos::RCP<Teuchos::ParameterList> SetParameterList(const Teuchos::ParameterList& paramList) = 0;

 protected:
  /*! @brief Replace all occurrences of search string "search" by the string in "replace" given the string "subject"
   */
  std::string ReplaceString(std::string& subject, const std::string& search, const std::string& replace) {
    size_t pos = 0;
    while ((pos = subject.find(search, pos)) != std::string::npos) {
      subject.replace(pos, search.length(), replace);
      pos += replace.length();
    }
    return subject;
  }
};

}  // namespace MueLu

#define MUELU_FACADECLASSBASE_SHORT

#endif /* PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_MUELU_FACADECLASSBASE_DECL_HPP_ */
