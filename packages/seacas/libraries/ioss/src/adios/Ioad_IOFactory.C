// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "Ioss_DBUsage.h"          // for DatabaseUsage
#include "Ioss_IOFactory.h"        // for IOFactory
#include <adios/Ioad_DatabaseIO.h> // for DatabaseIO
#include <adios/Ioad_IOFactory.h>
#include <stddef.h> // for nullptr
#include <string>   // for string

#include <adios2/common/ADIOSConfig.h>
#include <fmt/ostream.h>

namespace Ioss {
  class PropertyManager;
}

namespace Ioad {

  const IOFactory *IOFactory::factory()
  {
    static IOFactory registerThis;
    return &registerThis;
  }

  IOFactory::IOFactory() : Ioss::IOFactory("adios") { Ioss::IOFactory::alias("adios", "adios2"); }

  Ioss::DatabaseIO *IOFactory::make_IO(const std::string &filename, Ioss::DatabaseUsage db_usage,
                                       MPI_Comm                     communicator,
                                       const Ioss::PropertyManager &properties) const
  {
    return new DatabaseIO(nullptr, filename, db_usage, communicator, properties);
  }

  void IOFactory::show_config() const
  {
    fmt::print(Ioss::OUTPUT(), "\tADIOS2 Library Version: {}.{}.{}\n", ADIOS2_VERSION_MAJOR,
               ADIOS2_VERSION_MINOR, ADIOS2_VERSION_PATCH);
#if defined(ADIOS2_HAVE_BZIP2)
    fmt::print(Ioss::OUTPUT(), "\t\tBZip2 (http://www.bzip.org/) compression enabled\n");
#endif

#if defined(ADIOS2_HAVE_ZFP)
    fmt::print(Ioss::OUTPUT(), "\t\tZFP (https://github.com/LLNL/zfp) compression enabled\n");
#endif

#if defined(ADIOS2_HAVE_SZ)
    fmt::print(Ioss::OUTPUT(), "\t\tSZ compression enabled\n");
#endif

#if defined(ADIOS2_HAVE_MPI)
    fmt::print(Ioss::OUTPUT(), "\t\tParallel (MPI) enabled\n");
#else
    fmt::print(Ioss::OUTPUT(), "\t\tParallel *NOT* enabled\n");
#endif

#if defined(ADIOS2_HAVE_SST)
    fmt::print(Ioss::OUTPUT(), "\t\tStaging engine enabled\n");
#else
    fmt::print(Ioss::OUTPUT(), "\t\tStaging engine *NOT* enabled\n");
#endif

#if defined(ADIOS2_HAVE_HDF5)
    fmt::print(Ioss::OUTPUT(), "\t\tHDF5 (https://www.hdfgroup.org) engine enabled\n\n");
#else
    fmt::print(Ioss::OUTPUT(), "\t\tHDF5 engine *NOT* enabled\n\n");
#endif

    /* #if defined(ADIOS2_HAVE_ZEROMQ) */
    /* #if defined(ADIOS2_HAVE_WDM) */
    /* #if defined(ADIOS2_HAVE_DATAMAN) */
    /* #if defined(ADIOS2_HAVE_MGARD) */
    /* #undef ADIOS2_HAVE_PYTHON */
    /* #define ADIOS2_HAVE_FORTRAN */
    /* #define ADIOS2_HAVE_SYSVSHMEM */
    /* #undef ADIOS2_HAVE_ENDIAN_REVERSE */
  }
} // namespace Ioad
