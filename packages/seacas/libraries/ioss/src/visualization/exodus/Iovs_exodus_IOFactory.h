// Copyright(C) 1999-2021, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef IOSS_Iovs_exodus_IOFactory_h
#define IOSS_Iovs_exodus_IOFactory_h

#include "iovs_export.h"

#include "Ioss_CodeTypes.h"
#include "Ioss_DBUsage.h"    // for DatabaseUsage
#include "Ioss_DatabaseIO.h" // for DatabaseIO
#include "Ioss_IOFactory.h"  // for IOFactory
#include <string>            // for string
namespace Ioss {
  class PropertyManager;
} // namespace Ioss

namespace Iovs_exodus {

  class IOVS_EXPORT IOFactory : public Ioss::IOFactory
  {
  public:
    static const IOFactory *factory();

  private:
    IOFactory();
    Ioss::DatabaseIO *make_IO(const std::string &filename, Ioss::DatabaseUsage db_usage,
                              Ioss_MPI_Comm                communicator,
                              const Ioss::PropertyManager &properties) const override;

    /**
     * Call the sierra product registry and register all dependent third-party libraries
     */
    void register_library_versions() const;
  };
} // namespace Iovs_exodus
#endif // IOSS_Iovs_exodus_IOFactory_h
