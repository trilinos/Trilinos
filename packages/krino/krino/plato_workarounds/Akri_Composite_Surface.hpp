// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Compatibility shim: krino::Composite_Surface was removed in favor of
// krino::CSG_Surface / krino::CSGTree.  This header is retained only so external
// callers that still #include <Akri_Composite_Surface.hpp> continue to compile;
// the minimal CSG-backed replacement lives in Akri_PlatoWorkarounds.hpp.
#ifndef Akri_Composite_Surface_h
#define Akri_Composite_Surface_h

#include <Akri_PlatoWorkarounds.hpp>

#endif // Akri_Composite_Surface_h
