// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_CHECK_BC_CONSISTENCY
#define PANZER_CHECK_BC_CONSISTENCY

#include <vector>
#include <string>

namespace panzer {

  class BC;

  /** \brief Validates that each boundary condition's element block and
    * sideset exist in the mesh, throwing std::runtime_error on the
    * first one that does not.
    *
    * \param[in] element_block_names the element block IDs present in the mesh.
    * \param[in] sideset_names the sideset IDs present in the mesh.
    * \param[in] bcs the boundary conditions to validate.
    */
  void checkBCConsistency(const std::vector<std::string>& element_block_names,
                          const std::vector<std::string>& sideset_names,
                          const std::vector<panzer::BC>& bcs);

}

#endif
