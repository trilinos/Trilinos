/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef IOSS_Ioss_Property_h
#define IOSS_Ioss_Property_h

#include <Ioss_CodeTypes.h>
#include <string>

namespace Ioss {
  class GroupingEntity;

  class Property {
  public:
    enum BasicType {INVALID = -1, REAL, INTEGER, POINTER, STRING};
    enum VariableType { UNKNOWN_VAR_TYPE = -1, SCALAR};

    Property();
    Property(const std::string &name, const BasicType type,
		  const VariableType storage, void *data,
		  bool is_implicit = false);
    Property(const std::string &name, const int    value,
		  bool is_implicit = false);
    Property(const std::string &name, const double   value,
		  bool is_implicit = false);
    Property(const std::string &name, const std::string &value,
		  bool is_implicit = false);

    // To set implicit property
    Property(const GroupingEntity* ge,
		  const std::string &name, const BasicType type);

    Property(const Property&);

    bool operator<(const Property& other) const;

    ~Property();

    bool get_value(int    *value) const;
    bool get_value(double   *value) const;
    bool get_value(std::string *value) const;
    bool get_value(void   *value) const;

    std::string get_string()  const;
    int    get_int()     const;
    double   get_real()    const;
    void*  get_pointer() const;

    bool is_implicit() const {return isImplicit_;}
    bool is_explicit() const {return !isImplicit_;}
    bool is_valid()    const {return type_ != INVALID;}
    bool is_invalid()  const {return type_ == INVALID;}

    std::string get_name() const {return name_;}
    BasicType get_type() const {return type_;}

  private:
    Property& operator=(const Property&); // Do not implement
    std::string  name_;
    BasicType       type_;
    VariableType    storage_;

    // True if property is calculated rather than stored.
    // False if property is stored in 'data_'
    bool            isImplicit_;

    // The actual value of the property.  Use 'type_' and 'storage_' to
    // discriminate the actual type of the property.
    union Data {
      std::string* sval;
      void*   pval;
      const GroupingEntity *ge;
      double    rval;
      int     ival;
    };
    Data data_;
  };
}
#endif
