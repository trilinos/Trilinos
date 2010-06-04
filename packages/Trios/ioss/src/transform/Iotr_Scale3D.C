/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <transform/Iotr_Scale3D.h>
#include <Ioss_Field.h>
#include <Ioss_VariableType.h>
#include <string>

#include <assert.h>

namespace Iotr {

  const Scale3D_Factory* Scale3D_Factory::factory()
  {
    static Scale3D_Factory registerThis;
    return &registerThis;
  }

  Scale3D_Factory::Scale3D_Factory()
    : Factory("scale3D")
  {
    Factory::alias("scale3D", "multiply3D");
  }

  Ioss::Transform* Scale3D_Factory::make(const std::string&) const
  { return new Scale3D(); }

  Scale3D::Scale3D()
  {
    intScale[0] = intScale[1] = intScale[2] = 1;
    realScale[0] = realScale[1] = realScale[2] = 1.0;
  }

  void Scale3D::set_properties(const std::string&, const std::vector<int> &values)
  {
    assert(values.size() == 3);
    intScale[0] = values[0];
    intScale[1] = values[1];
    intScale[2] = values[2];
  }

  void Scale3D::set_properties(const std::string&, const std::vector<double> &values)
  {
    assert(values.size() == 3);
    realScale[0] = values[0];
    realScale[1] = values[1];
    realScale[2] = values[2];
  }

  const Ioss::VariableType
  *Scale3D::output_storage(const Ioss::VariableType *in) const
  {
    return in;
  }

  int Scale3D::output_count(int in) const
  {
    // Does not modify the entity count...
    return in;
  }

  bool Scale3D::internal_execute(const Ioss::Field &field, void *data)
  {
    size_t count = field.transformed_count();
    assert(field.transformed_storage()->component_count() == 3);
    
    if (field.get_type() == Ioss::Field::REAL) {
      double *rdata = static_cast<double*>(data);

      for (size_t i = 0; i < count*3; i+=3) {
	rdata[i+0] *= realScale[0];
	rdata[i+1] *= realScale[1];
	rdata[i+2] *= realScale[2];
      }
    } else if (field.get_type() == Ioss::Field::INTEGER) {
      int *idata = static_cast<int*>(data);

      for (size_t i = 0; i < count*3; i+=3) {
	idata[i+0] *= intScale[0];
	idata[i+1] *= intScale[1];
	idata[i+2] *= intScale[2];
      }
    } else {
    }
    return true;
  }
}
