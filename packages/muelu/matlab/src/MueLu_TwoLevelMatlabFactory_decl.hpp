// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_TWOLEVELMATLABFACTORY_DECL_HPP
#define MUELU_TWOLEVELMATLABFACTORY_DECL_HPP

#include <string>

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_CrsMatrix_fwd.hpp>
#include <Xpetra_MatrixFactory_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelMatlabFactory_fwd.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_PerfUtils_fwd.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_Utilities_fwd.hpp"

#ifdef HAVE_MUELU_MATLAB
#include "mex.h"

namespace MueLu {
/*!
  @class TwoLevelMatlabFactory
  @ingroup MueMexClasses
  @brief Factory for interacting with Matlab
*/
template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class TwoLevelMatlabFactory : public TwoLevelFactoryBase {
#undef MUELU_TWOLEVELMATLABFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  TwoLevelMatlabFactory();

  virtual ~TwoLevelMatlabFactory() {}

  //@}

  //! @name Input
  //@{
  RCP<const ParameterList> GetValidParameterList() const;

  void DeclareInput(Level& fineLevel, Level& coarseLevel) const;

  //@}

  //! @name Build methods.
  //@{
  void Build(Level& fineLevel, Level& coarseLevel) const;
  //@}

  //! @ name Description
  //@{
  std::string description() const;
  //@}
 private:
  //@{

  mutable bool hasDeclaredInput_;

  //@}

  //@{

  //! List of arguments to the MATLAB function, in order.  These args must correspond to MueLu "Needs" objects for the fine level.  These must be listed before coarse needs.
  mutable std::vector<std::string> needsFine_;

  //! List of arguments to the MATLAB function, in order.  These args must correspond to MueLu "Needs" objects for the coarse level.  These must be listed after fine needs.
  mutable std::vector<std::string> needsCoarse_;

  //@}

};  // class TwoLevelMatlabFactory

}  // namespace MueLu

#define MUELU_TWOLEVELMATLABFACTORY_SHORT

#endif  // HAVE_MUELU_MATLAB
#endif  // MUELU TWOLEVELMATLABFACTORY_DECL_HPP
