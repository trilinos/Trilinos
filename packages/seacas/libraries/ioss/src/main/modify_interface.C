/*
 * Copyright(C) 1999-2023 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <cstdlib> // for exit, EXIT_SUCCESS, getenv
#include <fmt/core.h>
#include <iostream> // for operator<<, basic_ostream, etc
#include <stdio.h>
#include <string> // for char_traits, string

#include "Ioss_GetLongOpt.h"
#include "Ioss_Utils.h"
#include "modify_interface.h"

Modify::Interface::Interface() { enroll_options(); }

void Modify::Interface::enroll_options()
{
  options_.usage("[options] basename");

  options_.enroll("help", Ioss::GetLongOption::NoValue, "Print this summary and exit", nullptr);

  options_.enroll("db_type", Ioss::GetLongOption::MandatoryValue,
                  "Database Type: generated"
#if defined(SEACAS_HAVE_PAMGEN)
                  ", pamgen"
#endif
#if defined(SEACAS_HAVE_EXODUS)
                  ", exodus"
#endif
#if defined(SEACAS_HAVE_CGNS)
                  ", cgns"
#endif
#if defined(SEACAS_HAVE_DATAWAREHOUSE)
                  ", data_warehouse"
#endif
                  ".",
                  "unknown");
  options_.enroll("in_type", Ioss::GetLongOption::MandatoryValue, "(alias for db_type)", nullptr,
                  nullptr, true);

  options_.enroll("allow_modifications", Ioss::GetLongOption::NoValue,
                  "By default, io_modify will only allow creation of new assemblies.\n"
                  "\t\tIf this option is specified, then can modify assemblies that already exist "
                  "in database.\n"
                  "\t\tThis will cause the database to be rewritten. Without this option, it is "
                  "updated in place.",
                  nullptr, nullptr, true);
  options_.enroll("version", Ioss::GetLongOption::NoValue, "Print version and exit", nullptr);

  options_.enroll("copyright", Ioss::GetLongOption::NoValue, "Show copyright and license data.",
                  nullptr);
}

bool Modify::Interface::parse_options(int argc, char **argv)
{
  // Get options from environment variable also...
  char *options = getenv("IO_MODIFY_OPTIONS");
  if (options != nullptr) {
    fmt::print(
        stderr,
        "\nThe following options were specified via the IO_MODIFY_OPTIONS environment variable:\n"
        "\t{}\n\n",
        options);
    options_.parse(options, Ioss::GetLongOption::basename(*argv));
  }

  int option_index = options_.parse(argc, argv);
  if (option_index < 1) {
    return false;
  }

  if (options_.retrieve("help") != nullptr) {
    options_.usage(std::cerr);
    fmt::print(stderr, "\n\tCan also set options via IO_MODIFY_OPTIONS environment variable.\n\n");
    fmt::print(stderr,
               "\tDocumentation: "
               "https://sandialabs.github.io/seacas-docs/sphinx/html/index.html#io-modify\n\n");
    fmt::print(stderr, "\t->->-> Send email to gdsjaar@sandia.gov for {} support.<-<-<-\n",
               options_.program_name());
    exit(EXIT_SUCCESS);
  }

  if (options_.retrieve("version") != nullptr) {
    // Version is printed up front, just exit...
    exit(0);
  }

  allowModification_ = options_.retrieve("allow_modifications") != nullptr;
  filetype_          = options_.get_option_value("db_type", filetype_);

  if (options_.retrieve("copyright") != nullptr) {
    Ioss::Utils::copyright(std::cerr, "2020-2022");
    exit(EXIT_SUCCESS);
  }

  // Parse remaining options as directory paths.
  if (option_index < argc) {
    filename_ = argv[option_index];
  }
  else {
    fmt::print(stderr, "\nERROR: filename not specified\n\n");
    return false;
  }

  if (filetype_ == "unknown") {
    filetype_ = Ioss::Utils::get_type_from_file(filename_);
  }

  return true;
}
