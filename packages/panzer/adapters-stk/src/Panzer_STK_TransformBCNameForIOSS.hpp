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
  */
  std::string transformBCNameForIOSS(std::string& bc_name);
}

#endif
