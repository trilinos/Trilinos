/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef SIERRA_Iotr_Tensor_h
#define SIERRA_Iotr_Tensor_h

#include <Ioss_Transform.h>
#include <string>

namespace Ioss {
  class Field;
}

namespace Iotr {

  class Tensor_Factory : public Factory
    {
    public:
      static const Tensor_Factory* factory();
    private:
      Tensor_Factory();
      Ioss::Transform* make(const std::string& type) const;
    };

  class Tensor: public Ioss::Transform
    {
      friend class Tensor_Factory;
      enum TranType {TRACE, SPHERICAL, DEVIATOR, MAGNITUDE,
		     INVARIANTS, INVARIANT1, INVARIANT2, INVARIANT3};
    public:
      const Ioss::VariableType *output_storage(const Ioss::VariableType *in) const;
      int output_count(int in) const;

    protected:
      explicit Tensor(const std::string& type);

      bool internal_execute(const Ioss::Field &field, void *data);

    private:
      TranType type_;
    };
}

#endif // SIERRA_Iotr_Tensor_h
