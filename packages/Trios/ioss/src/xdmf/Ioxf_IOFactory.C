/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <xdmf/Ioxf_IOFactory.h>
#include <xdmf/Ioxf_DatabaseIO.h>
#include <string>

namespace Ioxf {

  const IOFactory* IOFactory::factory()
  {
    static IOFactory registerThis;
    return &registerThis;
  }

  IOFactory::IOFactory()
    : Ioss::IOFactory("xdmf")
  {
    Ioss::IOFactory::alias("xdmf", "xdmf_output");
  }

  Ioss::DatabaseIO* IOFactory::make_IO(const std::string& filename,
				       Ioss::DatabaseUsage db_usage,
				       MPI_Comm communicator) const
  { return new DatabaseIO(NULL, filename, db_usage, communicator); }

  void IOFactory::finalize()
  {
    DatabaseIO::finalize();
  }


}
