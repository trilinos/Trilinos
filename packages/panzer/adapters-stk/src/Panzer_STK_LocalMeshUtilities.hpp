// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_STK_LOCAL_MESH_UTILITIES_HPP
#define PANZER_STK_LOCAL_MESH_UTILITIES_HPP

#include "Panzer_LocalMeshInfo.hpp"
#include <string>

namespace panzer_stk
{
  class STK_Interface;

/** Create a structure containing information about the local portion of a given element block
  *
  * \param[in] mesh Reference to STK mesh interface
  * \param[in] element_block_name Name of the element block of interest
  *
  * \returns Structure containing local mesh information
  */
template <typename LO, typename GO>
panzer::LocalMeshInfo<LO,GO>
generateLocalMeshInfo(const panzer_stk::STK_Interface & mesh,
                      const std::string & element_block_name);


/** Create a structure containing information about the local portion of a given element block's sideset
  *
  * \param[in] mesh Reference to STK mesh interface
  * \param[in] element_block_name Name of the element block of interest
  * \param[in] sideset_name Name of the sideset on the element block of interest
  *
  * \returns Structure containing local mesh information
  */
template <typename LO, typename GO>
panzer::LocalMeshInfo<LO,GO>
generateLocalSidesetInfo(const panzer_stk::STK_Interface & mesh,
                      const std::string & element_block_name,
                      const std::string & sideset_name);

}
#endif
