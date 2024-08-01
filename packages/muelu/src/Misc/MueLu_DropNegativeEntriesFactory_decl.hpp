// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_DROPNEGATIVEENTRIESFACTORY_DECL_HPP
#define MUELU_DROPNEGATIVEENTRIESFACTORY_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_DropNegativeEntriesFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"

namespace MueLu {

/*!
  @class DropNegativeEntriesFactory class.
  @brief Application-specific filtering for A. Can be used in context of graph coarsening and aggregation.

  This factory drops all negative entries (or entries with a magnitude < 0). Only weak positive connections are kept.
  Do not use this kind of filtering for regular PDEs unless you have very good reasons.
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class DropNegativeEntriesFactory : public SingleLevelFactoryBase {
#undef MUELU_DROPNEGATIVEENTRIESFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  DropNegativeEntriesFactory() {}

  //! Destructor.
  virtual ~DropNegativeEntriesFactory() {}

  RCP<const ParameterList> GetValidParameterList() const;

  //@}

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const;

  //@}

  //! @name Build methods.
  //@{

  /*!
    @brief Build method.

    Builds filtered matrix and returns it in <tt>currentLevel</tt>.
    */
  void Build(Level& currentLevel) const;

  //@}

};  // class DropNegativeEntriesFactory

}  // namespace MueLu

#define MUELU_DROPNEGATIVEENTRIESFACTORY_SHORT
#endif  // MUELU_DROPNEGATIVEENTRIESFACTORY_DECL_HPP
