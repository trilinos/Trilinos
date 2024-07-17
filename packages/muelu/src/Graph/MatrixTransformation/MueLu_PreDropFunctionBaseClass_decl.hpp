// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_PREDROPFUNCTIONBASECLASS_DECL_HPP
#define MUELU_PREDROPFUNCTIONBASECLASS_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_PreDropFunctionBaseClass_fwd.hpp"

namespace MueLu {

/*!
 * Base class you can derive from to allow user defined dropping
 *
 */
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class PreDropFunctionBaseClass : public BaseClass {
#undef MUELU_PREDROPFUNCTIONBASECLASS_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! Destructor
  virtual ~PreDropFunctionBaseClass() {}

  //! Drop
  virtual bool Drop(size_t lrow, GlobalOrdinal grow, size_t k, LocalOrdinal lcid, GlobalOrdinal gcid, const Teuchos::ArrayView<const LocalOrdinal>& indices, const Teuchos::ArrayView<const Scalar>& vals) = 0;
};
}  // namespace MueLu

#define MUELU_PREDROPFUNCTIONBASECLASS_SHORT
#endif  // MUELU_PREDROPFUNCTIONBASECLASS_DECL_HPP
