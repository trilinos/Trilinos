// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_TWOLEVELFACTORY_HPP
#define MUELU_TWOLEVELFACTORY_HPP

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_Factory.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_TimeMonitor.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

/*!
  @class TwoLevelFactoryBase class.
  @brief Base class for factories that use two levels (fineLevel and coarseLevel).

  Examples of such factories are R, P, and A_coarse.

  @ingroup MueLuBaseClasses
*/

class TwoLevelFactoryBase : public Factory {
 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  TwoLevelFactoryBase();

  //! Destructor.
  virtual ~TwoLevelFactoryBase();

  //@}

  //! Input
  //@{

  /*! @brief Specifies the data that this class needs, and the factories that generate that data.

      If the Build method of this class requires some data, but the generating factory is not specified in DeclareInput, then this class
      will fall back to the settings in FactoryManager.
  */
  virtual void DeclareInput(Level& fineLevel, Level& coarseLevel) const = 0;

  //!
  virtual void CallDeclareInput(Level& requestedLevel) const;

  //@}

  //! @name Build methods.
  //@{

  //! Build an object with this factory.
  virtual void Build(Level& fineLevel, Level& coarseLevel) const = 0;

  //!
  virtual void CallBuild(Level& requestedLevel) const;

  //@}

};  // class TwoLevelFactoryBase

}  // namespace MueLu

#define MUELU_TWOLEVELFACTORY_SHORT
#endif  // ifndef MUELU_TWOLEVELFACTORY_HPP
