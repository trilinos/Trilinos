// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RBGEN_CREATE_PARAMS_H
#define RBGEN_CREATE_PARAMS_H

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

namespace RBGen {

  //! Create a Teuchos::ParameterList from an XML file.
  /*!
   */
  Teuchos::RCP<Teuchos::ParameterList> createParams( const std::string& filename );

  //! Extract the filename list from a Teuchos::ParameterList.
  /*!
   */
  Teuchos::RCP<std::vector<std::string> > genFileList( const Teuchos::ParameterList& params );

}
#endif
