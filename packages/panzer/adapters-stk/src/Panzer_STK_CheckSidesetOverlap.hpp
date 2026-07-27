// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_CHECK_SIDESET_OVERLAP_HPP
#define PANZER_CHECK_SIDESET_OVERLAP_HPP

#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include <string>

namespace panzer_stk {

  /** \brief Checks whether two sidesets share any nodes.
    *
    * \param[in] side_a_name the name of the first sideset.
    * \param[in] side_b_name the name of the second sideset.
    * \param[in] mesh the mesh containing both sidesets.
    * \return true if the two sidesets overlap (share at least one node).
    */
  bool checkSidesetOverlap(const std::string& side_a_name,
                           const std::string& side_b_name,
                           const panzer_stk::STK_Interface& mesh);
}

#endif
