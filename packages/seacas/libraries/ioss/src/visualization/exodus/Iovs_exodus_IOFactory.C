// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "Ioss_DBUsage.h"   // for DatabaseUsage
#include "Ioss_IOFactory.h" // for IOFactory
#include <cstddef>          // for nullptr
#include <string>           // for string
#include <visualization/exodus/Iovs_exodus_DatabaseIO.h>
#include <visualization/exodus/Iovs_exodus_IOFactory.h>

namespace Ioss {
  class PropertyManager;
} // namespace Ioss

namespace Iovs_exodus {

  const IOFactory *IOFactory::factory()
  {
    static IOFactory registerThis;
    return &registerThis;
  }

  IOFactory::IOFactory() : Ioss::IOFactory("catalyst_exodus")
  {
    // Tell the database to register itself with sierra's product registry.
    // XXX exodus doesn't do this, do we need to?
    // register_library_versions();
  }

  Ioss::DatabaseIO *IOFactory::make_IO(const std::string &filename, Ioss::DatabaseUsage db_usage,
                                       Ioss_MPI_Comm                communicator,
                                       const Ioss::PropertyManager &properties) const
  {

    return new DatabaseIO(nullptr, filename, db_usage, communicator, properties);
  }

  /**
   * Call the sierra product registry and register all dependent third-party libraries
   */
  void IOFactory::register_library_versions() const
  {
    // Internals::register_library_versions();
  }
} // namespace Iovs_exodus
