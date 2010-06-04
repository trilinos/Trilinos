/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_Transform.h>
#include <Ioss_Utils.h>
#include <string>

namespace Iotr {

Ioss::Transform* Factory::create(const std::string& type)
{
  Ioss::Transform *transform = NULL;
  FactoryMap::iterator iter = registry()->find(type);
  if (iter == registry()->end()) {
    if (registry()->size() == 0) {
      std::ostringstream errmsg;
      errmsg << "FATAL: No transformations have been registered.\n"
	     << "       Was Iotr::Initializer::initialize() called?\n\n";
      IOSS_ERROR(errmsg);
    } else {
      std::ostringstream errmsg;
      errmsg << "FATAL: The transform named '" << type
	     << "' is not supported.\n";
      IOSS_ERROR(errmsg);
    }
  } else {
    Factory* factory = (*iter).second;
    transform = factory->make(type);
  }
  return transform;
}

int Factory::describe(NameList *names)
{
  int count = 0;
  FactoryMap::const_iterator I;
  for (I = registry()->begin(); I != registry()->end(); ++I) {
    names->push_back((*I).first);
    count++;
  }
  return count;
}

Factory::Factory(const std::string& type)
{
  registry()->insert(FactoryValuePair(type, this));
}

void Factory::alias(const std::string& base, const std::string& syn)
{
  Factory* factory = (*registry()->find(base)).second;
  registry()->insert(FactoryValuePair(syn, factory));
}

FactoryMap* Factory::registry()
{
  static FactoryMap registry_;
  return &registry_;
}

}
