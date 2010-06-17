/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef IOSS_Ioxf_IOFactory_h
#define IOSS_Ioxf_IOFactory_h

#include <Ioss_IOFactory.h>
#include <Ioss_DBUsage.h>
#include <string>

namespace Ioxf {

  class IOFactory : public Ioss::IOFactory
    {
    public:
      static const IOFactory* factory();
      static void finalize();

    private:
      IOFactory();
      Ioss::DatabaseIO* make_IO(const std::string& filename,
				Ioss::DatabaseUsage db_usage,
				MPI_Comm communicator) const;

    };
}
#endif // IOSS_Ioxf_IOFactory_h
