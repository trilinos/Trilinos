// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_PARSER_AKRI_YAML_HPP_
#define KRINO_KRINO_PARSER_AKRI_YAML_HPP_
#include <Akri_config.hpp>

#ifdef KRINO_HAVE_YAML

#ifdef __INTEL_COMPILER
#include <yaml-cpp/yaml.h>
#else
//YAML has shadowed variables
//this disables the checking on GCC only
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <yaml-cpp/yaml.h>
#pragma GCC diagnostic pop
#endif


#endif

#endif /* KRINO_KRINO_PARSER_AKRI_YAML_HPP_ */
