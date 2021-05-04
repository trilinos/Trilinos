// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef Iovs_cgns_IOFactory_h
#define Iovs_cgns_IOFactory_h

#include "Ioss_DatabaseIO.h" // for DatabaseIO
#include <Ioss_CodeTypes.h>
#include <Ioss_DBUsage.h>   // for DatabaseUsage
#include <Ioss_IOFactory.h> // for IOFactory
#include <string>           // for string
namespace Ioss {
  class PropertyManager;
} // namespace Ioss

namespace Iovs_cgns {

  class IOFactory : public Ioss::IOFactory
  {
  public:
    static const IOFactory *factory();

  private:
    IOFactory();
    Ioss::DatabaseIO *make_IO(const std::string &filename, Ioss::DatabaseUsage db_usage,
                              MPI_Comm                     communicator,
                              const Ioss::PropertyManager &properties) const override;
  };
} // namespace Iovs_cgns
#endif // Iovs_cgns_IOFactory_h
