// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 * MueLu_PermutationFactory_decl.hpp
 *
 *  Created on: Nov 28, 2012
 *      Author: wiesner
 */

#ifndef MUELU_PERMUTATIONFACTORY_DECL_HPP_
#define MUELU_PERMUTATIONFACTORY_DECL_HPP_

#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"

#include "MueLu_AlgebraicPermutationStrategy_fwd.hpp"
#include "MueLu_LocalPermutationStrategy_fwd.hpp"

namespace MueLu {

/*!
  @class PermutationFactory class.
  @brief factory generates a row- and column permutation operators P and Q such that
  P*A*Q^T is a (hopefully) diagonal-dominant matrix.
  It's meant to be used with PermutingSmoother.
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class PermutationFactory : public SingleLevelFactoryBase {
#undef MUELU_PERMUTATIONFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  PermutationFactory();

  //! Destructor.
  virtual ~PermutationFactory();

  //@}

  RCP<const ParameterList> GetValidParameterList() const;

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
};  // class PermutationFactory

}  // namespace MueLu

#define MUELU_PERMUTATIONFACTORY_SHORT

#endif /* MUELU_PERMUTATIONFACTORY_DECL_HPP_ */
