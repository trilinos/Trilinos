/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef SIERRA_Iotr_Scale_h
#define SIERRA_Iotr_Scale_h

#include <Ioss_Transform.h>
#include <string>

namespace Ioss {
  class Field;
}

namespace Iotr {

  class Scale_Factory : public Factory
    {
    public:
      static const Scale_Factory* factory();
    private:
      Scale_Factory();
      Ioss::Transform* make(const std::string&) const;
    };

  class Scale: public Ioss::Transform
    {
      friend class Scale_Factory;

    public:
      const Ioss::VariableType *output_storage(const Ioss::VariableType *in) const;
      int output_count(int in) const;

      void set_property(const std::string &name, int value);
      void set_property(const std::string &name, double value);

    protected:
      Scale();

      bool internal_execute(const Ioss::Field &field, void *data);

    private:
      int    intMultiplier;
      double realMultiplier;
    };
}

#endif // SIERRA_Iotr_Scale_h
