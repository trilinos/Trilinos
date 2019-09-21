// Copyright(C) 1999-2017 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
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
//
//     * Neither the name of NTESS nor the names of its
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

#include <cstdint> // for int64_t
#include <string>  // for string
namespace Ioss {
  class GroupingEntity;
} // namespace Ioss

namespace Ioss {

  /** \brief A named value that has a known type.
   *
   */
  class Property
  {
  public:
    enum BasicType { INVALID = -1, REAL, INTEGER, POINTER, STRING };

    Property();
    Property(std::string name, int64_t value);
    Property(std::string name, int value);
    Property(std::string name, double value);
    Property(std::string name, const std::string &value);
    Property(std::string name, const char *value);
    Property(std::string name, void *value);

    // To set implicit property
    Property(const GroupingEntity *ge, std::string name, BasicType type);

    Property(const Property &from);
    Property &operator=(Property rhs);

    ~Property();

    std::string get_string() const;
    int64_t     get_int() const;
    double      get_real() const;
    void *      get_pointer() const;

    /** \brief Tells whether the property is calculated, rather than stored.
     *
     * \returns True if property is calculated; False if it is stored directly.
     */
    bool is_implicit() const { return isImplicit_; }

    /** \brief Tells whether the property is stored, rather than calculated.
     *
     * \returns True if property is stored directly; False if it is calculated.
     */
    bool is_explicit() const { return !isImplicit_; }

    /** Tells whether the property has a valid type (currently REAL, INTEGER, POINTER, or STRING)
     *
     *  \returns True if the property type is valid.
     */
    bool is_valid() const { return type_ != INVALID; }

    /** Tells whether the property has an invalid type (currently not one of REAL, INTEGER, POINTER,
     * or STRING)
     *
     *  \returns True if the property type is invalid.
     */
    bool is_invalid() const { return type_ == INVALID; }

    /** \brief Get the property name.
     *
     *  \returns The property name.
     */
    std::string get_name() const { return name_; }

    /** \brief Get the property type.
     *
     *  \returns The property type.
     */
    BasicType get_type() const { return type_; }

  private:
    std::string name_;
    BasicType   type_;

    bool get_value(int64_t *value) const;
    bool get_value(double *value) const;
    bool get_value(std::string *value) const;
    bool get_value(void *&value) const;

    /// True if property is calculated rather than stored.
    /// False if property is stored in 'data_'
    bool isImplicit_{false};

    /// The actual value of the property.  Use 'type_' and 'storage_' to
    /// discriminate the actual type of the property.
    union Data {
      std::string *         sval;
      void *                pval;
      const GroupingEntity *ge;
      double                rval;
      int64_t               ival;
    };
    Data data_{};
  };
} // namespace Ioss
#endif
