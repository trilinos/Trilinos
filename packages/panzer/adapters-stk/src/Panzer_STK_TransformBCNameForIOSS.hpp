// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_TRANSFORM_BC_NAME_FOR_IOSS
#define PANZER_TRANSFORM_BC_NAME_FOR_IOSS

#include "PanzerAdaptersSTK_config.hpp"
#include <string>

namespace panzer_stk {
  /** Takes a sideset or nodeset name and transforms it for IOSS
      restrictions. Cubit replaces whitespace in sideset and nodeset
      names with underscores when writing to exodus. When IOSS reads
      exodus files into STK, it replaces all capital letters with
      lowercase letters. This function allows you to use the name
      labels used in Cubit and tie them to the correct nodeset or
      sideset in a STK mesh database.

      @param bc_name String to be converted - this object is transformed in-place
      @param make_lower_case If set to true, convert to lower case.
      @return Returns the transformed name

      NOTE: On 2026.05.30, IOSS was changed to allow for upper case in
      sideset names. We now only have to account for the spaces being
      changed to underscores. This behavior can be reversed with IOSS
      properties.
  */
  std::string transformBCNameForIOSS(std::string& bc_name, const bool make_lower_case=false);
}

#endif
