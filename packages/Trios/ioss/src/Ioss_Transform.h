/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef IOSS_Ioss_Transform_h
#define IOSS_Ioss_Transform_h

#include <string>
#include <vector>
#include <map>

namespace Ioss {
  class Field;
  class VariableType;
}
namespace Ioss {
  class Transform
    {
    public:
      virtual ~Transform();
      virtual const
	Ioss::VariableType *output_storage(const Ioss::VariableType *in) const = 0;
      virtual int output_count(int in) const = 0;

      bool execute(const Ioss::Field &field, void *data);

      virtual void set_property(const std::string &name, int value);
      virtual void set_property(const std::string &name, double value);
      virtual void set_properties(const std::string &name,
				  const std::vector<int> &values);
      virtual void set_properties(const std::string &name,
				  const std::vector<double> &values);
    protected:
      Transform();

      virtual bool internal_execute(const Ioss::Field &field, void *data) = 0;
    };
}

namespace Iotr {

  class Factory;

  typedef std::vector<std::string> NameList;
  typedef std::map<std::string, Factory*, std::less<std::string> > FactoryMap;
  typedef FactoryMap::value_type FactoryValuePair;

  class Factory {
  public:
    virtual ~Factory() {};
    static Ioss::Transform* create(const std::string& type);

    static int describe(NameList *names);

  protected:
    explicit Factory(const std::string& type);
    virtual Ioss::Transform* make(const std::string&) const = 0;
    static void alias(const std::string& base, const std::string& syn);

  private:
    static FactoryMap* registry();
  };
}

#endif // IOSS_Ioss_Transform_h
