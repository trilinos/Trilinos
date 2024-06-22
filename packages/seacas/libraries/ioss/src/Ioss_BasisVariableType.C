// Copyright(C) 2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "Ioss_BasisVariableType.h"
#include "Ioss_ComposedVariableType.h"
#include "Ioss_CompositeVariableType.h"
#include "Ioss_ConstructedVariableType.h"
#include "Ioss_NamedSuffixVariableType.h"
#include "Ioss_QuadratureVariableType.h"
#include "Ioss_Utils.h"
#include "Ioss_VariableType.h"
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "Ioss_CodeTypes.h"
