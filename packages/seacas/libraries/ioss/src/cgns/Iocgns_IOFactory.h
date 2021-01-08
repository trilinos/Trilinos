// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef IOSS_Iocgns_IOFactory_h
#define IOSS_Iocgns_IOFactory_h

#include "Ioss_DatabaseIO.h" // for DatabaseIO
#include <Ioss_CodeTypes.h>
#include <Ioss_DBUsage.h>   // for DatabaseUsage
#include <Ioss_IOFactory.h> // for IOFactory
#include <string>           // for string

namespace Ioss {
  class PropertyManager;
} // namespace Ioss

namespace Iocgns {

  class IOFactory : public Ioss::IOFactory
  {
  public:
    static const IOFactory *factory();
    static void             finalize();

  private:
    IOFactory();
    Ioss::DatabaseIO *make_IO(const std::string &filename, Ioss::DatabaseUsage db_usage,
                              MPI_Comm                     communicator,
                              const Ioss::PropertyManager &properties) const override;
    std::string       show_config() const override;
  };
} // namespace Iocgns
#endif // IOSS_Iocgns_IOFactory_h
