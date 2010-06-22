/*--------------------------------------------------------------------*/
/*    Copyright 2000, 2008, 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_IOFactory.h>

#include <Ioss_Utils.h>
#include <string>

namespace {
  const char *
  get_product_name()
  {
    return "I/O System";
  }
}

Ioss::DatabaseIO* Ioss::IOFactory::create(const std::string& type,
					  const std::string& filename,
					  Ioss::DatabaseUsage db_usage,
					  MPI_Comm communicator)
{
  Ioss::DatabaseIO *db = NULL;
  Ioss::IOFactoryMap::iterator iter = registry()->find(type);
  if (iter == registry()->end()) {
    if (registry()->size() == 0) {
      std::ostringstream errmsg;
      errmsg << "FATAL: No database types have been registered.\n"
	     << "       Was Ioss::Init::Initializer() called?\n\n";
      IOSS_ERROR(errmsg);
    } else {
      std::ostringstream errmsg;
      errmsg << "FATAL: The database type '" << type
	     << "' is not supported.\n";
      IOSS_ERROR(errmsg);
    }
  } else {
    Ioss::IOFactory* factory = (*iter).second;
    db = factory->make_IO(filename, db_usage, communicator);
  }
  return db;
}

int Ioss::IOFactory::describe(NameList *names)
{
  int count = 0;
  Ioss::IOFactoryMap::const_iterator I;
  for (I = registry()->begin(); I != registry()->end(); ++I) {
    names->push_back((*I).first);
    ++count;
  }
  return count;
}

Ioss::IOFactory::IOFactory(const std::string& type)
{
  registry()->insert(IOFactoryValuePair(type, this));
}

void Ioss::IOFactory::alias(const std::string& base, const std::string& syn)
{
  Ioss::IOFactory* factory = (*registry()->find(base)).second;
  registry()->insert(IOFactoryValuePair(syn, factory));
}

Ioss::IOFactoryMap* Ioss::IOFactory::registry()
{
  static IOFactoryMap registry_;
  return &registry_;
}

void Ioss::IOFactory::clean()
{
}
