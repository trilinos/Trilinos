/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <exodusII/Ioex_IOFactory.h>
#include <exodusII/Ioex_DatabaseIO.h>
#include <exodusII/Ioex_Internals.h>
#include <string>

namespace Ioex {

  const IOFactory* IOFactory::factory()
  {
    static IOFactory registerThis;
    return &registerThis;
  }

  IOFactory::IOFactory()
    : Ioss::IOFactory("exodusII")
  {
    Ioss::IOFactory::alias("exodusII", "exodusii");
    Ioss::IOFactory::alias("exodusII", "exodus");
    Ioss::IOFactory::alias("exodusII", "genesis");
    Ioss::IOFactory::alias("exodusII", "exodusII_input");
    Ioss::IOFactory::alias("exodusII", "exodusII_output");
  }

  Ioss::DatabaseIO* IOFactory::make_IO(const std::string& filename,
				       Ioss::DatabaseUsage db_usage,
				       MPI_Comm communicator) const
  { return new DatabaseIO(NULL, filename, db_usage, communicator); }

}
