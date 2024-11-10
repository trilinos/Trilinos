/*
 * Copyright(C) 1999-2020, 2023, 2024 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <fmt/core.h>
#include <ostream>
#include <stdlib.h>
#include <string>

#include "Ionit_Initializer.h"
#include "Ioss_IOFactory.h"
#include "Ioss_ScopeGuard.h"
#include "Ioss_Utils.h"
#include "info_interface.h"
#include "io_info.h"

// ========================================================================

namespace {
  std::string codename;
  std::string version = "1.06";

#ifdef SEACAS_HAVE_MPI
  void mpi_finalize()
  {
    MPI_Comm parentcomm;
    MPI_Comm_get_parent(&parentcomm);
    if (parentcomm != MPI_COMM_NULL) {
      int istatus = EXIT_SUCCESS;
      MPI_Send(&istatus, 1, MPI_INT, 0, 0, parentcomm);
    }
    MPI_Finalize();
  }
#endif
} // namespace

int main(int argc, char *argv[])
{
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
  ON_BLOCK_EXIT(mpi_finalize);
#endif

  Info::Interface interFace;
  interFace.parse_options(argc, argv);

  Ioss::Init::Initializer io;

  if (interFace.show_config()) {
    Ioss::OUTPUT() << "\n" << Ioss::IOFactory::show_configuration();
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
