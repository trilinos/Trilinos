/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef IOSS_Iotr_Offset3D_h
#define IOSS_Iotr_Offset3D_h

#include <Ioss_Transform.h>
#include <string>

namespace Ioss {
  class Field;
}

namespace Iotr {

  class Offset3D_Factory : public Factory
    {
    public:
      static const Offset3D_Factory* factory();
    private:
      Offset3D_Factory();
      Ioss::Transform* make(const std::string&) const;
    };

  class Offset3D: public Ioss::Transform
    {
      friend class Offset3D_Factory;

    public:
      const Ioss::VariableType *output_storage(const Ioss::VariableType *in) const;
      int output_count(int in) const;

      void set_properties(const std::string &name, const std::vector<int>    &values);
      void set_properties(const std::string &name, const std::vector<double> &values);

    protected:
      Offset3D();

      bool internal_execute(const Ioss::Field &field, void *data);

    private:
      int    intOffset[3];
      double realOffset[3];
    };
}

#endif // IOSS_Iotr_Offset3D_h
