/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef SIERRA_Iotr_MinMax_h
#define SIERRA_Iotr_MinMax_h

#include <Ioss_Transform.h>
#include <string>

namespace Ioss {
  class Field;
}

namespace Iotr {

  class MinMax_Factory : public Factory
    {
    public:
      static const MinMax_Factory* factory();
    private:
      MinMax_Factory();
      Ioss::Transform* make(const std::string& type) const;
    };

  class MinMax: public Ioss::Transform
    {
      friend class MinMax_Factory;

    public:
      const Ioss::VariableType *output_storage(const Ioss::VariableType *in) const;
      int output_count(int in) const;

    protected:
      explicit MinMax(const std::string& type);

      bool internal_execute(const Ioss::Field &field, void *data);

    private:
      bool doMin;
      bool doAbs;
    };
}

#endif // SIERRA_Iotr_MinMax_h
