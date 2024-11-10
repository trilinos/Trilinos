// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_THRESHOLDAFILTERFACTORY_DECL_HPP
#define MUELU_THRESHOLDAFILTERFACTORY_DECL_HPP

/*
 * MueLu_ThresholdAFilterFactory.hpp
 *
 *  Created on: 14.10.2011
 *      Author: tobias
 */

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_ThresholdAFilterFactory_fwd.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace MueLu {

/*!
  @class ThresholdAFilterFactory class.
  @brief Factory for building a thresholded operator.

*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class ThresholdAFilterFactory : public SingleLevelFactoryBase {
#undef MUELU_THRESHOLDAFILTERFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  using magnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  ThresholdAFilterFactory(const std::string& ename, const magnitudeType threshold, const bool keepDiagonal = true, const GlobalOrdinal expectedNNZperRow = -1);

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const;

  //@}

  //@{
  //! @name Build methods.

  //! Build an object with this factory.
  void Build(Level& currentLevel) const;

  //@}

 private:
  std::string varName_;            ///< name of input and output variable
  const magnitudeType threshold_;  ///< threshold parameter
  const bool keepDiagonal_;
  const GlobalOrdinal expectedNNZperRow_;

};  // class ThresholdAFilterFactory

}  // namespace MueLu

#define MUELU_THRESHOLDAFILTERFACTORY_SHORT
#endif  // MUELU_THRESHOLDAFILTERFACTORY_DECL_HPP
