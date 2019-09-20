// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/Stacktrace.hpp>

namespace percept {

#if DO_STACKTRACE_POS
#if DO_STACKTRACE_POS_USE_INT
  int Stacktrace::s_position = 0;
#else
  std::string Stacktrace::s_position = "";
#endif
#endif
}
