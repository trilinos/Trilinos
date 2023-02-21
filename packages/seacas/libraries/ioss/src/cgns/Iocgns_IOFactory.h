// Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include "iocgns_export.h"

#include "Ioss_DatabaseIO.h" // for DatabaseIO
#include <Ioss_CodeTypes.h>
#include <Ioss_DBUsage.h>   // for DatabaseUsage
#include <Ioss_IOFactory.h> // for IOFactory
#include <string>           // for string

namespace Ioss {
  class PropertyManager;
} // namespace Ioss

namespace Iocgns {

  class IOCGNS_EXPORT IOFactory : public Ioss::IOFactory
  {
  public:
    static const IOFactory *factory();
    static void             finalize();

  private:
    IOFactory();
    Ioss::DatabaseIO *make_IO(const std::string &filename, Ioss::DatabaseUsage db_usage,
                              Ioss_MPI_Comm                communicator,
                              const Ioss::PropertyManager &properties) const override;
    std::string       show_config() const override;
  };
} // namespace Iocgns
