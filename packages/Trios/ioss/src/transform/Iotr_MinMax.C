/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <transform/Iotr_MinMax.h>
#include <Ioss_Field.h>
#include <Ioss_VariableType.h>
#include <string>

#include <algorithm>

#include <math.h>
#include <assert.h>

namespace {
  bool IntAbsLess(int elem1, int elem2)
  {return abs(elem1) < abs(elem2);}
  bool doubleAbsLess(double elem1, double elem2)
  {return fabs(elem1) < fabs(elem2);}
}

namespace Iotr {

  const MinMax_Factory* MinMax_Factory::factory()
  {
    static MinMax_Factory registerThis;
    return &registerThis;
  }

  MinMax_Factory::MinMax_Factory()
    : Factory("generic_minmax")
  {
    Factory::alias("generic_minmax", "minimum");
    Factory::alias("generic_minmax", "maximum");
    Factory::alias("generic_minmax", "absolute_minimum");
    Factory::alias("generic_minmax", "absolute_maximum");
  }

  Ioss::Transform* MinMax_Factory::make(const std::string& type) const
  { return new MinMax(type); }

  MinMax::MinMax(const std::string &type)
  {
    if (type == "minimum") {
      doMin = true;
      doAbs = false;
    } else if (type == "maximum") {
      doMin = false;
      doAbs = false;
    } else if (type == "absolute_minimum") {
      doMin = true;
      doAbs = true;
    } else if (type == "absolute_maximum") {
      doMin = false;
      doAbs = true;
    }
  }

  const Ioss::VariableType
  *MinMax::output_storage(const Ioss::VariableType *in) const
  {
    // Only operates on scalars...
    static const Ioss::VariableType *sca = Ioss::VariableType::factory("scalar");
    if (in == sca) {
      return sca;
    } else {
      return NULL;
    }
  }

  int MinMax::output_count(int /* in */) const
  {
    // Returns a single value...
    return 1;
  }

  bool MinMax::internal_execute(const Ioss::Field &field, void *data)
  {
    size_t count = field.transformed_count();
    size_t components = field.transformed_storage()->component_count();
    size_t n = count * components;
    if (field.get_type() == Ioss::Field::REAL) {
      double *rdata = static_cast<double*>(data);
      double value;
      if (doMin) {
	if (doAbs) {
	  value = *std::min_element(&rdata[0], &rdata[n], doubleAbsLess);
	} else {
	  value = *std::min_element(&rdata[0], &rdata[n]);
	}
      } else { // doMax
	if (doAbs) {
	  value = *std::max_element(&rdata[0], &rdata[n], doubleAbsLess);
	} else {
	  value = *std::max_element(&rdata[0], &rdata[n]);
	}
      }
      rdata[0] = value;
    } else if (field.get_type() == Ioss::Field::INTEGER) {
      int *idata = static_cast<int*>(data);
      int value;
      if (doMin) {
	if (doAbs) {
	  value = *std::min_element(&idata[0], &idata[n], IntAbsLess);
	} else {
	  value = *std::min_element(&idata[0], &idata[n]);
	}
      } else { // doMax
	if (doAbs) {
	  value = *std::max_element(&idata[0], &idata[n], IntAbsLess);
	} else {
	  value = *std::max_element(&idata[0], &idata[n]);
	}
      }
      idata[0] = value;
    }
    return true;
  }
}
