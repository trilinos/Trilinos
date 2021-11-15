// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Startup_h
#define Akri_Startup_h

#include <fstream>
#include <streambuf>

namespace krino{
  
class Startup
{
 public:

  Startup(int argc, char ** argv);
  ~Startup();

  bool exit_early() const { return my_flag_exit_early; }
  void handle_exception(const char * what, const bool is_parsing);

  static void setup_commandline_options();
  static void report_handler(const char *message, int type);

private:
  int my_flag_exit_early;
};

} // namespace krino

#endif // Akri_Startup_h
