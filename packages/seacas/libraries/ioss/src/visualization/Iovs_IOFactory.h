/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*--------------------------------------------------------------------*/
/*    Copyright 2010 NTESS.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef IOSS_Iovs_IOFactory_h
#define IOSS_Iovs_IOFactory_h

#include "Ioss_DatabaseIO.h" // for DatabaseIO
#include <Ioss_CodeTypes.h>
#include <Ioss_DBUsage.h>   // for DatabaseUsage
#include <Ioss_IOFactory.h> // for IOFactory
#include <string>           // for string
namespace Ioss {
  class PropertyManager;
} // namespace Ioss

namespace Iovs {

  class IOFactory : public Ioss::IOFactory
  {
  public:
    static const IOFactory *factory();

  private:
    IOFactory();
    Ioss::DatabaseIO *make_IO(const std::string &filename, Ioss::DatabaseUsage db_usage,
                              MPI_Comm                     communicator,
                              const Ioss::PropertyManager &properties) const override;

    /**
     * Call the sierra product registry and register all dependent third-party libraries
     */
    void register_library_versions() const;
  };
} // namespace Iovs
#endif // IOSS_Iovs_IOFactory_h
