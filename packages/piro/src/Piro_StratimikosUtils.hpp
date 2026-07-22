// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_STRATIMIKOSUTILS_HPP
#define PIRO_STRATIMIKOSUTILS_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

namespace Piro {

  //! \brief Extracts the Stratimikos sublist from the Piro Solver parameter list
  Teuchos::RCP<Teuchos::ParameterList>
    extractStratimikosParams(const Teuchos::RCP<Teuchos::ParameterList> &piroParams);
  
  //! \brief Rename the preconditioner and parameter list
  void
    renamePreconditionerParamList(const Teuchos::RCP<Teuchos::ParameterList> &stratParams, 
              const std::string &oldname, const std::string &newname);

} // namespace Piro

#endif /* PIRO_STRATIMIKOSUTILS_HPP */
