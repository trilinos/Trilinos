// Copyright(C) 1999-2010
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
//         
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef IOSS_Ioss_Property_h
#define IOSS_Ioss_Property_h

#include <Ioss_CodeTypes.h>
#include <string>
#include <stdint.h>

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
    Property(const std::string &name, int64_t  value,
		  bool is_implicit = false);
    Property(const std::string &name, int      value,
		  bool is_implicit = false);
    Property(const std::string &name, double   value,
		  bool is_implicit = false);
    Property(const std::string &name, const std::string &value,
		  bool is_implicit = false);
    Property(const std::string &name, void *value,
		  bool is_implicit = false);

    // To set implicit property
    Property(const GroupingEntity* ge,
		  const std::string &name, const BasicType type);

    Property(const Property&);

    bool operator<(const Property& other) const;

    ~Property();

    std::string get_string()  const;
    int64_t  get_int()     const;
    double   get_real()    const;
    void*    get_pointer() const;

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

    bool get_value(int64_t *value) const;
    bool get_value(double  *value) const;
    bool get_value(std::string *value) const;
    bool get_value(void   *&value) const;

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
      int64_t   ival;
    };
    Data data_;
  };
}
#endif
