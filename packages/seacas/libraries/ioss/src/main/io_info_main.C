/*
 * Copyright(C) 1999-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "io_info.h"

// ========================================================================

namespace {
  std::string codename;
  std::string version = "1.0";
} // namespace

int main(int argc, char *argv[])
{
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  Info::Interface interface;
  interface.parse_options(argc, argv);

  Ioss::Init::Initializer io;

  if (interface.show_config()) {
    Ioss::IOFactory::show_configuration();
    exit(EXIT_SUCCESS);
  }

  codename   = argv[0];
  size_t ind = codename.find_last_of('/', codename.size());
  if (ind != std::string::npos) {
    codename = codename.substr(ind + 1, codename.size());
  }

  OUTPUT << "Input:    '" << interface.filename() << "', Type: " << interface.type() << '\n';
  OUTPUT << '\n';

  if (interface.list_groups()) {
    Ioss::io_info_group_info(interface);
  }
  else {
    Ioss::io_info_file_info(interface);
  }

  OUTPUT << "\n" << codename << " execution successful.\n";
#ifdef SEACAS_HAVE_MPI
  MPI_Finalize();
#endif
  return EXIT_SUCCESS;
}
