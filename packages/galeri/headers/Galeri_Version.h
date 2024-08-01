// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_VERSION_H
#define GALERI_VERSION_H

#include "Galeri_ConfigDefs.h"
#include "Trilinos_version.h"

std::string Galeri_Version() { 
	return("Galeri in Trilinos " TRILINOS_VERSION_STRING); 
};

#endif /* GALERI_VERSION_H */
