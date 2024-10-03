// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_DISABLEMULTIPLECALLCHECK_HPP
#define MUELU_DISABLEMULTIPLECALLCHECK_HPP

#include <Teuchos_RCP.hpp>

#include "MueLu_TwoLevelFactoryBase.hpp"

namespace MueLu {

//! An exception safe way to call the method TwoLevelFactoryBase::DisableMultipleCallCheck
class DisableMultipleCallCheck {
 public:
  DisableMultipleCallCheck(const RCP<const TwoLevelFactoryBase>& fact)
    : fact_(fact) { fact_->DisableMultipleCallCheck(); }
  ~DisableMultipleCallCheck() { fact_->EnableMultipleCallCheck(); }

 private:
  const RCP<const TwoLevelFactoryBase> fact_;
};

}  // namespace MueLu

#endif  // MUELU_DISABLEMULTIPLECALLCHECK_HPP
