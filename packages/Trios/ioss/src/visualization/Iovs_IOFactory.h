/*--------------------------------------------------------------------*/
/*    Copyright 2010 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef SIERRA_Iovs_IOFactory_h
#define SIERRA_Iovs_IOFactory_h

#include <Ioss_IOFactory.h>
#include <Ioss_DBUsage.h>
#include <string>

namespace Iovs {

  class IOFactory : public Ioss::IOFactory
    {
    public:
      static const IOFactory* factory();


    private:
      IOFactory();
      Ioss::DatabaseIO* make_IO(const std::string& filename,
				Ioss::DatabaseUsage db_usage,
				MPI_Comm communicator) const;

      /**
       * Call the sierra product registry and register all dependent third-party libraries
       */
      void register_library_versions() const;
    };
}
#endif // SIERRA_Iovs_IOFactory_h
