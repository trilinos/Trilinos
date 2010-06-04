/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <transform/Iotr_Offset.h>
#include <Ioss_Field.h>
#include <Ioss_VariableType.h>
#include <string>

#include <math.h>
#include <assert.h>

namespace Iotr {

  const Offset_Factory* Offset_Factory::factory()
  {
    static Offset_Factory registerThis;
    return &registerThis;
  }

  Offset_Factory::Offset_Factory()
    : Factory("offset")
  {
    Factory::alias("offset", "add");
  }

  Ioss::Transform* Offset_Factory::make(const std::string&) const
  { return new Offset(); }

  Offset::Offset()
    : intOffset(0), realOffset(0.0)
  {}

  void Offset::set_property(const std::string&, int value)
  { intOffset = value; }

  void Offset::set_property(const std::string&, double value)
  { realOffset = value; }

  const Ioss::VariableType
  *Offset::output_storage(const Ioss::VariableType *in) const
  {
    return in;
  }

  int Offset::output_count(int in) const
  {
    // Does not modify the entity count...
    return in;
  }

  bool Offset::internal_execute(const Ioss::Field &field, void *data)
  {
    size_t count = field.transformed_count();
    int components = field.transformed_storage()->component_count();

    if (field.get_type() == Ioss::Field::REAL) {
      double *rdata = static_cast<double*>(data);

      for (size_t i = 0; i < count*components; i++) {
	rdata[i] = rdata[i] + realOffset;
      }
    } else if (field.get_type() == Ioss::Field::INTEGER) {
      int *idata = static_cast<int*>(data);

      for (size_t i = 0; i < count*components; i++) {
	idata[i] = idata[i] + intOffset;
      }
    } else {
    }
    return true;
  }
}
