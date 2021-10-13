// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Unit_LogRedirecter.hpp>

#include <Akri_DiagWriter.hpp>

namespace krino {

LogRedirecter::LogRedirecter()
: myOriginalBuffer(krinolog.getStream().rdbuf())
{
  krinolog.getStream().rdbuf(&myBuffer);
}

LogRedirecter::~LogRedirecter()
{
  krinolog.getStream().rdbuf(myOriginalBuffer);
}

}

