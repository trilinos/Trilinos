/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "fmt/format.h"
#include "io_info.h"
#include <Ioss_ScopeGuard.h>

// ========================================================================

namespace {
  std::string codename;
  std::string version = "1.05";
} // namespace

int main(int argc, char *argv[])
{
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
  ON_BLOCK_EXIT(MPI_Finalize);
#endif

  Info::Interface interFace;
  interFace.parse_options(argc, argv);

  Ioss::Init::Initializer io;

  if (interFace.show_config()) {
    Ioss::IOFactory::show_configuration();
    exit(EXIT_SUCCESS);
  }

  codename   = argv[0];
  size_t ind = codename.find_last_of('/', codename.size());
  if (ind != std::string::npos) {
    codename = codename.substr(ind + 1, codename.size());
  }

  if (interFace.list_groups()) {
    Ioss::io_info_group_info(interFace);
  }
  else {
    Ioss::io_info_file_info(interFace);
  }

  fmt::print("\n{} execution successful.\n", codename);
  return EXIT_SUCCESS;
}
