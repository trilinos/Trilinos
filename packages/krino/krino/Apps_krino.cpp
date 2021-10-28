// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Parser.hpp>
#include <Akri_Startup.hpp>
#include <Akri_Simulation.hpp>
#include <stk_util/environment/Trace.hpp>

int main( int argc, char ** argv )
{
  krino::Startup startup(argc, argv);

  if (startup.exit_early()) {
    return 0;
  }

  bool is_parsing = true;

  try {
    krino::Parser::parse();

    is_parsing = false;
    krino::Simulation & simulation = krino::Simulation::get();
    simulation.commit();
    simulation.initialize();
    simulation.execute();
  }
  catch (std::exception &x) {
    stk::diag::Trace::Preserve preserve__;
    startup.handle_exception(x.what(),is_parsing);
  }
  catch (...) {
    stk::diag::Trace::Preserve preserve__;
    startup.handle_exception("Unknown",is_parsing);
  }

  // call Simulation::reset() manually to make sure the associated static objects are destroyed before MPI_FINALIZE is called by ~Startup
  krino::Simulation::reset();

  // all done
  return 0;
}
