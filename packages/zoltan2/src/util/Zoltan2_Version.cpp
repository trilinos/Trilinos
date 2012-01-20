//@HEADER
// ************************************************************************
//         copyright
// ************************************************************************
//@HEADER

#include <Zoltan2_Version.hpp>

/*! \file Zoltan2_Version.cpp
    \brief Simple function for returning the current version number
*/

namespace Zoltan2 {

  std::string Zoltan2_Version() { 
    return("Zoltan2 in Trilinos " TRILINOS_VERSION_STRING);
  }

} // namespace Zoltan2 
