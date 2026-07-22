// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_PARAMETER_LIBRARY_ACCEPTOR_HPP
#define PANZER_PARAMETER_LIBRARY_ACCEPTOR_HPP

#include "Panzer_ParameterLibrary.hpp"
#include "Teuchos_RCP.hpp"

namespace panzer {

  /** \brief Pure Virtual base class for accepting the parameter library

      This class is used to retrieve the parameter library from an
      object.
  */
  class ParameterLibraryAcceptor {

  public:

    virtual ~ParameterLibraryAcceptor() {}

    virtual Teuchos::RCP<panzer::ParamLib> getParameterLibrary() const = 0;

  };

}

#endif
