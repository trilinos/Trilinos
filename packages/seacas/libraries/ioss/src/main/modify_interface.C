/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "Ioss_CodeTypes.h"
#include "Ioss_FileInfo.h"
#include "Ioss_GetLongOpt.h" // for GetLongOption, etc
#include "fmt/ostream.h"
#include "modify_interface.h"
#include <cstddef>  // for nullptr
#include <cstdlib>  // for exit, EXIT_SUCCESS, getenv
#include <iostream> // for operator<<, basic_ostream, etc
#include <string>   // for char_traits, string

namespace {
  std::string get_type_from_file(const std::string &filename)
  {
    Ioss::FileInfo file(filename);
    auto           extension = file.extension();
    if (extension == "e" || extension == "g" || extension == "gen" || extension == "exo") {
      return "exodus";
    }
    else if (extension == "cgns") {
      return "cgns";
    }
    else {
      // "exodus" is default...
      return "exodus";
    }
  }
} // namespace

Modify::Interface::Interface() { enroll_options(); }

Modify::Interface::~Interface() = default;

void Modify::Interface::enroll_options()
{
  options_.usage("[options] basename");

  options_.enroll("help", Ioss::GetLongOption::NoValue, "Print this summary and exit", nullptr);

  options_.enroll("version", Ioss::GetLongOption::NoValue, "Print version and exit", nullptr);

  options_.enroll("allow_modifications", Ioss::GetLongOption::NoValue,
                  "By default, io_modify will only allow creation of new assemblies.\n"
                  "\t\tIf this option is specified, then can modify assemblies that already exist "
                  "in database.\n"
                  "\t\tThis will cause the database to be rewritten. Without this option, it is "
                  "updated in place.",
                  nullptr);
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
  options_.enroll("in_type", Ioss::GetLongOption::MandatoryValue, "(alias for db_type)", nullptr);
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
    fmt::print(stderr, "\t->->-> Send email to gdsjaar@sandia.gov for {} support.<-<-<-\n",
               options_.program_name());
    exit(EXIT_SUCCESS);
  }

  if (options_.retrieve("version") != nullptr) {
    // Version is printed up front, just exit...
    exit(0);
  }

  if (options_.retrieve("allow_modifications") != nullptr) {
    allowModification_ = true;
  }

  {
    const char *temp = options_.retrieve("db_type");
    if (temp != nullptr) {
      filetype_ = temp;
    }
  }

  if (options_.retrieve("copyright") != nullptr) {
    fmt::print(stderr, "\n"
                       "Copyright(C) 2020 National Technology & Engineering Solutions\n"
                       "of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with\n"
                       "NTESS, the U.S. Government retains certain rights in this software.\n\n"
                       "Redistribution and use in source and binary forms, with or without\n"
                       "modification, are permitted provided that the following conditions are\n"
                       "met:\n\n "
                       "    * Redistributions of source code must retain the above copyright\n"
                       "      notice, this list of conditions and the following disclaimer.\n\n"
                       "    * Redistributions in binary form must reproduce the above\n"
                       "      copyright notice, this list of conditions and the following\n"
                       "      disclaimer in the documentation and/or other materials provided\n"
                       "      with the distribution.\n\n"
                       "    * Neither the name of NTESS nor the names of its\n"
                       "      contributors may be used to endorse or promote products derived\n"
                       "      from this software without specific prior written permission.\n\n"
                       "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS\n"
                       "\" AS IS \" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT\n"
                       "LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR\n"
                       "A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT\n"
                       "OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,\n"
                       "SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT\n"
                       "LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,\n"
                       "DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY\n"
                       "THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT\n"
                       "(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE\n"
                       "OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n\n");
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
    filetype_ = get_type_from_file(filename_);
  }

  return true;
}
