// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_STK_LOCAL_MESH_UTILITIES_HPP
#define PANZER_STK_LOCAL_MESH_UTILITIES_HPP

#include "Teuchos_RCP.hpp"

namespace panzer
{
  struct LocalMeshInfo;
}

namespace panzer_stk
{
  class STK_Interface;

  /**
   * \brief Create a structure containing information about the local portion of a given element block
   *
   * \param[in] mesh Reference to STK mesh interface
   *
   * \returns Structure containing local mesh information
   */
  Teuchos::RCP<panzer::LocalMeshInfo>
  generateLocalMeshInfo(const panzer_stk::STK_Interface & mesh);
}

#endif
