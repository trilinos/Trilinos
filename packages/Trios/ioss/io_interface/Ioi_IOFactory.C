/*--------------------------------------------------------------------*/
/*    Copyright 2007 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioi_IOFactory.h>
#include <assert.h>

Ioss::Interface::IOFactory* Ioss::Interface::IOFactory::factory_ = NULL;

Ioss::Interface::IOBroker* Ioss::Interface::IOFactory::create(const std::string& region_name)
{
  assert(factory_ != NULL);
  Ioss::Interface::IOBroker *db = factory_->make_IO(region_name);
  return db;
}

Ioss::Interface::IOFactory::IOFactory()
{
  factory_ = this;
}
