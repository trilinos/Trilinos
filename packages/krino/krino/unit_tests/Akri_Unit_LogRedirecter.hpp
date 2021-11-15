// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef AKRI_UNIT_LOGREDIRECTER_H_
#define AKRI_UNIT_LOGREDIRECTER_H_

#include <sstream>

namespace krino {

class LogRedirecter
{
public:
  LogRedirecter();
  ~LogRedirecter();

  void clear() { myBuffer.str(""); }
  std::string get_log() const { return myBuffer.str(); }
private:
  std::stringbuf myBuffer;
  std::streambuf * myOriginalBuffer;
};

}


#endif
