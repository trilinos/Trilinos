/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <transform/Iotr_Scale.h>
#include <Ioss_Field.h>
#include <Ioss_VariableType.h>
#include <string>

#include <math.h>
#include <assert.h>

namespace Iotr {

  const Scale_Factory* Scale_Factory::factory()
  {
    static Scale_Factory registerThis;
    return &registerThis;
  }

  Scale_Factory::Scale_Factory()
    : Factory("scale")
  {
    Factory::alias("scale", "multiply");
  }

  Ioss::Transform* Scale_Factory::make(const std::string&) const
  { return new Scale(); }

  Scale::Scale()
    : intMultiplier(1), realMultiplier(1.0)
  {}

  void Scale::set_property(const std::string&, int value)
  { intMultiplier = value; }

  void Scale::set_property(const std::string&, double value)
  { realMultiplier = value; }

  const Ioss::VariableType
  *Scale::output_storage(const Ioss::VariableType *in) const
  {
    return in;
  }

  int Scale::output_count(int in) const
  {
    // Does not modify the entity count...
    return in;
  }

  bool Scale::internal_execute(const Ioss::Field &field, void *data)
  {
    size_t count = field.transformed_count();
    int components = field.transformed_storage()->component_count();

    if (field.get_type() == Ioss::Field::REAL) {
      double *rdata = static_cast<double*>(data);

      for (size_t i = 0; i < count*components; i++) {
	rdata[i] = rdata[i] * realMultiplier;
      }
    } else if (field.get_type() == Ioss::Field::INTEGER) {
      int *idata = static_cast<int*>(data);

      for (size_t i = 0; i < count*components; i++) {
	idata[i] = idata[i] * intMultiplier;
      }
    } else {
    }
    return true;
  }
}
