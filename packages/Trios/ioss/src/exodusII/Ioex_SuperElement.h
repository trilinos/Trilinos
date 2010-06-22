#ifndef IOSS_Ioex_SuperElement_h
#define IOSS_Ioex_SuperElement_h

#include <Ioss_CodeTypes.h>
#include <Ioss_GroupingEntity.h>
#include <string>

namespace Ioss {
  class Property;
  class Field;
}

namespace Ioex {
  class SuperElement : public Ioss::GroupingEntity
  {
  public:
    SuperElement(const std::string &filename, const std::string& name);
    ~SuperElement();

    std::string type_string() const {return "SuperElement";}
    Ioss::EntityType type() const {return Ioss::SUPERELEMENT;}
      
    // Handle implicit properties -- These are calcuated from data stored
    // in the grouping entity instead of having an explicit value assigned.
    // An example would be 'element_block_count' for a region.
    Ioss::Property get_implicit_property(const std::string& name) const;

  protected:
    int internal_get_field_data(const Ioss::Field& field,
				void *data, size_t data_size) const;

    int internal_put_field_data(const Ioss::Field& field,
				void *data, size_t data_size) const;

  private:
    std::string fileName;
    size_t numDOF;
    size_t numEIG;
    int filePtr;
  };
}
#endif
