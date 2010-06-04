/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <transform/Iotr_VectorMagnitude.h>
#include <Ioss_Field.h>
#include <Ioss_VariableType.h>
#include <string>

#include <math.h>
#include <assert.h>

namespace Iotr {

  const VM_Factory* VM_Factory::factory()
  {
    static VM_Factory registerThis;
    return &registerThis;
  }

  VM_Factory::VM_Factory()
    : Factory("vector magnitude")
  {
    Factory::alias("vector magnitude", "length");
  }

  Ioss::Transform* VM_Factory::make(const std::string&) const
  { return new VectorMagnitude(); }

  VectorMagnitude::VectorMagnitude() {}

  const Ioss::VariableType
  *VectorMagnitude::output_storage(const Ioss::VariableType *in) const
  {
    static const Ioss::VariableType *v2d = Ioss::VariableType::factory("vector_2d");
    static const Ioss::VariableType *v3d = Ioss::VariableType::factory("vector_3d");
    static const Ioss::VariableType *sca = Ioss::VariableType::factory("scalar");
    if (in == v2d || in == v3d) {
      return sca;
    } else {
      return NULL;
    }
  }

  int VectorMagnitude::output_count(int in) const
  {
    // Does not modify the entity count...
    return in;
  }

  bool VectorMagnitude::internal_execute(const Ioss::Field &field,
					 void *data)
  {
    double *rdata = static_cast<double*>(data);

    size_t count = field.transformed_count();
    if (field.transformed_storage()->component_count() == 3) {
      int j = 0;
      for (size_t i = 0; i < count; i++) {
	rdata[i] = sqrt(rdata[j]*rdata[j] + rdata[j+1]*rdata[j+1] +
			rdata[j+2]*rdata[j+2]);
	j+=3;
      }
    } else {
      int j = 0;
      for (size_t i = 0; i < count; i++) {
	rdata[i] = sqrt(rdata[j]*rdata[j] + rdata[j+1]*rdata[j+1]);
	j+=2;
      }
    }
    return true;
  }
}
