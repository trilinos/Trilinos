/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <transform/Iotr_Offset3D.h>
#include <Ioss_Field.h>
#include <Ioss_VariableType.h>
#include <string>

#include <assert.h>

namespace Iotr {

  const Offset3D_Factory* Offset3D_Factory::factory()
  {
    static Offset3D_Factory registerThis;
    return &registerThis;
  }

  Offset3D_Factory::Offset3D_Factory()
    : Factory("offset3D")
  {
    Factory::alias("offset3D", "add3D");
  }

  Ioss::Transform* Offset3D_Factory::make(const std::string&) const
  { return new Offset3D(); }

  Offset3D::Offset3D()
  {
    intOffset[0] = intOffset[1] = intOffset[2] = 0;
    realOffset[0] = realOffset[1] = realOffset[2] = 0.0;
  }

  void Offset3D::set_properties(const std::string&, const std::vector<int> &values)
  {
    assert(values.size() == 3);
    intOffset[0] = values[0];
    intOffset[1] = values[1];
    intOffset[2] = values[2];
  }

  void Offset3D::set_properties(const std::string&, const std::vector<double> &values)
  {
    assert(values.size() == 3);
    realOffset[0] = values[0];
    realOffset[1] = values[1];
    realOffset[2] = values[2];
  }

  const Ioss::VariableType
  *Offset3D::output_storage(const Ioss::VariableType *in) const
  {
    return in;
  }

  int Offset3D::output_count(int in) const
  {
    // Does not modify the entity count...
    return in;
  }

  bool Offset3D::internal_execute(const Ioss::Field &field, void *data)
  {
    size_t count = field.transformed_count();
    assert(field.transformed_storage()->component_count() == 3);
    
    if (field.get_type() == Ioss::Field::REAL) {
      double *rdata = static_cast<double*>(data);

      for (size_t i = 0; i < count*3; i+=3) {
	rdata[i+0] += realOffset[0];
	rdata[i+1] += realOffset[1];
	rdata[i+2] += realOffset[2];
      }
    } else if (field.get_type() == Ioss::Field::INTEGER) {
      int *idata = static_cast<int*>(data);

      for (size_t i = 0; i < count*3; i+=3) {
	idata[i+0] += intOffset[0];
	idata[i+1] += intOffset[1];
	idata[i+2] += intOffset[2];
      }
    } else {
    }
    return true;
  }
}
