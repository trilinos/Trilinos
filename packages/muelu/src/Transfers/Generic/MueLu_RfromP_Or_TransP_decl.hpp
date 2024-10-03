// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_RFROMP_OR_TRANSP_DECL_HPP
#define MUELU_RFROMP_OR_TRANSP_DECL_HPP

#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_RfromP_Or_TransP_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_PerfUtils_fwd.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace MueLu {

/*!
  @class RfromP_Or_TransP class.
  @brief Factory for building restriction operators.

  This factory currently depends on an underlying matrix-matrix multiply with the identity
  matrix to do the transpose.  This should probably be fixed at some point.
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class RfromP_Or_TransP : public TwoLevelFactoryBase {
#undef MUELU_RFROMP_OR_TRANSP_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  RfromP_Or_TransP() {}

  //! Destructor.
  virtual ~RfromP_Or_TransP() {}

  RCP<const ParameterList> GetValidParameterList() const;

  //@}

  //! Input
  //@{

  void DeclareInput(Level &fineLevel, Level &coarseLevel) const;

  //@}

  //! @name Build methods.
  //@{

  void Build(Level &fineLevel, Level &coarseLevel) const;

  //@}

};  // class RfromP_Or_TransP

}  // namespace MueLu

#define MUELU_RFROMP_OR_TRANSP_SHORT
#endif  // MUELU_RFROMP_OR_TRANSP_DECL_HPP
