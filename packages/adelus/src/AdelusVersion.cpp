/*
//@HEADER
// *****************************************************************************
//                        Adelus
//
// Copyright 2020 NTESS and the Adelus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER
*/

/*! \file AdelusVersion.cpp
    \brief Simple function for returning the current version number [necessary for portability]
*/

#include "AdelusVersion.hpp"
namespace Adelus {

  std::string Adelus_Version() {
  	return("Adelus Version 1.0 - 09/07/????");
  }

}// namespace Adelus
