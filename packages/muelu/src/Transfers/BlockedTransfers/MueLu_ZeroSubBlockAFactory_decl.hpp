// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_ZEROSUBBLOCKAFACTORY_DECL_HPP_
#define MUELU_ZEROSUBBLOCKAFACTORY_DECL_HPP_

#include <Xpetra_Map_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_ZeroSubBlockAFactory_fwd.hpp"

namespace MueLu {

/*!
  @class ZeroSubBlockAFactory class.
  @brief Factory for extracting a zero block from a BlockedCrsMatrix.

  This is a very simple class to access a single matrix block in a blocked operator A,
  where the matrix block of interest is an actual zero block (as e.g. in saddle point systems).

*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class ZeroSubBlockAFactory : public SingleLevelFactoryBase {
#undef MUELU_ZEROSUBBLOCKAFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! Input
  //@{

  RCP<const ParameterList> GetValidParameterList() const override;

  void DeclareInput(Level &currentLevel) const override;

  //@}

  //@{
  //! @name Build methods

  /*!
    @brief Build a zero sub-block object with this factory

    Create a zero sub block matrix, that fits into a given blocked crs operator.

    Strided or block information is extracted in the following way:
      1) first check whether the corresponding sub maps are strided
         If yes, use the fixed block size and strided block id
      2) If no, get the full map of the map extractors.
         Check whether the full map is strided.
         If yes, use the strided information of the full map and build
         partial (strided) maps with it
         If no, throw an exception

    For blocked operators with block maps one should use the striding
    information from the sub maps. For strided operators, the striding
    information of the full map is the best choice.
  */
  void Build(Level &currentLevel) const override;

  //@}

};  // class ZeroSubBlockAFactory

}  // namespace MueLu

#define MUELU_ZEROSUBBLOCKAFACTORY_SHORT
#endif /* MUELU_ZEROSUBBLOCKAFACTORY_DECL_HPP_ */
