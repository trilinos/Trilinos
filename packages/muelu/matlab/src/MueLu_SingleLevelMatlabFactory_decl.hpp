// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SINGLELEVELMATLABFACTORY_DECL_HPP
#define MUELU_SINGLELEVELMATLABFACTORY_DECL_HPP

#include <string>

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_CrsMatrix_fwd.hpp>
#include <Xpetra_MatrixFactory_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelMatlabFactory_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_PerfUtils_fwd.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_Utilities_fwd.hpp"

#ifdef HAVE_MUELU_MATLAB
#include "mex.h"

namespace MueLu {
/*!
  @class SingleLevelMatlabFactory
  @ingroup MueMexClasses
  @brief Factory for interacting with Matlab
*/
template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class SingleLevelMatlabFactory : public SingleLevelFactoryBase {
#undef MUELU_SINGLELEVELMATLABFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  SingleLevelMatlabFactory();

  virtual ~SingleLevelMatlabFactory() {}

  //@}

  //! @name Input
  //@{
  RCP<const ParameterList> GetValidParameterList() const;

  void DeclareInput(Level& currentLevel) const;

  //@}

  //! @name Build methods.
  //@{
  void Build(Level& currentLevel) const;
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

  //! List of arguments to the MATLAB function, in order.  These args must correspond to MueLu "Needs" objects.
  mutable std::vector<std::string> needs_;

  //@}

};  // class SingleLevelMatlabFactory

}  // namespace MueLu

#define MUELU_SINGLELEVELMATLABFACTORY_SHORT

#endif  // HAVE_MUELU_MATLAB
#endif  // MUELU SINGLELEVELMATLABFACTORY_DECL_HPP
