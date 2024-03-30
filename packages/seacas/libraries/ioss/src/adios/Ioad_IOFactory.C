// Copyright(C) 1999-2020, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "Ioss_DBUsage.h"          // for DatabaseUsage
#include "Ioss_IOFactory.h"        // for IOFactory
#include "adios/Ioad_DatabaseIO.h" // for DatabaseIO
#include "adios/Ioad_IOFactory.h"
#include <cstddef> // for nullptr
#include <string>  // for string

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
                                       Ioss_MPI_Comm                communicator,
                                       const Ioss::PropertyManager &properties) const
  {
    return new DatabaseIO(nullptr, filename, db_usage, communicator, properties);
  }

  std::string IOFactory::show_config() const
  {
    std::stringstream config;
    fmt::print(config, "\tADIOS2 Library Version: {}.{}.{}\n", ADIOS2_VERSION_MAJOR,
               ADIOS2_VERSION_MINOR, ADIOS2_VERSION_PATCH);
#if defined(ADIOS2_HAVE_BZIP2)
    fmt::print(config, "\t\tBZip2 (http://www.bzip.org/) compression enabled\n");
#endif

#if defined(ADIOS2_HAVE_ZFP)
    fmt::print(config, "\t\tZFP (https://github.com/LLNL/zfp) compression enabled\n");
#endif

#if defined(ADIOS2_HAVE_SZ)
    fmt::print(config, "\t\tSZ compression enabled\n");
#endif

#if defined(ADIOS2_HAVE_MPI)
    fmt::print(config, "\t\tParallel (MPI) enabled\n");
#else
    fmt::print(config, "\t\tParallel *NOT* enabled\n");
#endif

#if defined(ADIOS2_HAVE_SST)
    fmt::print(config, "\t\tStaging engine enabled\n");
#else
    fmt::print(config, "\t\tStaging engine *NOT* enabled\n");
#endif

#if defined(ADIOS2_HAVE_HDF5)
    fmt::print(config, "\t\tHDF5 (https://www.hdfgroup.org) engine enabled\n\n");
#else
    fmt::print(config, "\t\tHDF5 engine *NOT* enabled\n\n");
#endif

    /* #if defined(ADIOS2_HAVE_ZEROMQ) */
    /* #if defined(ADIOS2_HAVE_WDM) */
    /* #if defined(ADIOS2_HAVE_DATAMAN) */
    /* #if defined(ADIOS2_HAVE_MGARD) */
    /* #undef ADIOS2_HAVE_PYTHON */
    /* #define ADIOS2_HAVE_FORTRAN */
    /* #define ADIOS2_HAVE_SYSVSHMEM */
    /* #undef ADIOS2_HAVE_ENDIAN_REVERSE */
    return config.str();
  }
} // namespace Ioad
