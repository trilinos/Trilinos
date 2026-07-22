// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_STK_STK_INTERFACE_GENERATOR_HPP_
#define PANZER_STK_STK_INTERFACE_GENERATOR_HPP_

#include "Teuchos_RCP.hpp"

namespace Teuchos
{
class ParameterList;
template<typename T>
class Comm;
}

namespace panzer_stk
{

class STK_Interface;

/**\brief Generalized mesh interface
 *
 * Generates a mesh given a comm and a parameter list
 *
 * \param[in] comm Communicator for domain
 * \param[in] parameter_list The mesh description should be in this
 *
 * \return A mesh
 */
Teuchos::RCP<panzer_stk::STK_Interface>
generateMesh(const Teuchos::ParameterList & parameter_list);

}

#endif
