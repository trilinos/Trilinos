// @HEADER
// *****************************************************************************
//                            Sphynx
//
// Copyright 2020 NTESS and the Sphynx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Trilinos_version.h"
#include <string> 

namespace Zoltan2 {

   std::string Sphynx_Version() { 
     return("Sphynx in Trilinos " TRILINOS_VERSION_STRING); 
   }

} // namespace Zoltan2
