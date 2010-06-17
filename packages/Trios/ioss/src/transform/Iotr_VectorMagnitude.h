/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef IOSS_Iotr_VectorMagnitude_h
#define IOSS_Iotr_VectorMagnitude_h

#include <Ioss_Transform.h>
#include <string>

namespace Ioss {
  class Field;
}

namespace Iotr {

  class VM_Factory : public Factory
    {
    public:
      static const VM_Factory* factory();
    private:
      VM_Factory();
      Ioss::Transform* make(const std::string&) const;
    };

  class VectorMagnitude: public Ioss::Transform
    {
      friend class VM_Factory;

    public:
      const Ioss::VariableType *output_storage(const Ioss::VariableType *in) const;
      int output_count(int in) const;

    protected:
      VectorMagnitude();

      bool internal_execute(const Ioss::Field &field, void *data);
    };
}

#endif // IOSS_Iotr_VectorMagnitude_h
