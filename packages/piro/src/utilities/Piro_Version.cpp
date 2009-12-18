// @HEADER
// @HEADER

#include "Piro_Version.hpp"
#include "Trilinos_version.h"

std::string Piro::Piro_Version()
{ 
  return("Piro in Trilinos " TRILINOS_VERSION_STRING); 
}
