// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_STK_UTILITIES_BUILD_MESH_HPP_
#define PANZER_STK_UTILITIES_BUILD_MESH_HPP_

#include "Teuchos_RCP.hpp"

#include <vector>

namespace Teuchos
{
class ParameterList;
}

namespace panzer_stk
{

class STK_Interface;

Teuchos::RCP<panzer_stk::STK_Interface>
buildMesh(const std::vector<int> & N,
          const std::vector<int> & B,
          const std::vector<double> & L);

Teuchos::RCP<panzer_stk::STK_Interface>
buildParallelMesh(const std::vector<int> & N,
                  const std::vector<int> & B,
                  const std::vector<int> & P,
                  const std::vector<double> & L,
                  const std::vector<int> & periodic_dims);

}

#endif
