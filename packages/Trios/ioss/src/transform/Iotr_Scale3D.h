/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef SIERRA_Iotr_Scale3D_h
#define SIERRA_Iotr_Scale3D_h

#include <Ioss_Transform.h>
#include <string>

namespace Ioss {
  class Field;
}

namespace Iotr {

  class Scale3D_Factory : public Factory
    {
    public:
      static const Scale3D_Factory* factory();
    private:
      Scale3D_Factory();
      Ioss::Transform* make(const std::string&) const;
    };

  class Scale3D: public Ioss::Transform
    {
      friend class Scale3D_Factory;

    public:
      const Ioss::VariableType *output_storage(const Ioss::VariableType *in) const;
      int output_count(int in) const;

      void set_properties(const std::string &name, const std::vector<int>    &values);
      void set_properties(const std::string &name, const std::vector<double> &values);

    protected:
      Scale3D();

      bool internal_execute(const Ioss::Field &field, void *data);

    private:
      int    intScale[3];
      double realScale[3];
    };
}

#endif // SIERRA_Iotr_Scale3D_h
