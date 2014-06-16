#ifndef SAMBA_SAMBA_UTILITY_MAKE_PRIMITIVE_HPP
#define SAMBA_SAMBA_UTILITY_MAKE_PRIMITIVE_HPP

#include <samba/utility/integer_max.hpp>
#include <samba/utility/primitive_comparison.hpp>
#include <samba/utility/primitive_arithmetic_operators.hpp>
#include <ostream>

//
// Use these macros to define primitive types. Sometimes, this won't be
// possible if your primitive requires some specialized behavior; in
// that case, just try to make the struct as similar to the macro-defined
// primitives as possible.
//
// There are a few key concepts for primitive types:
//   1) They should be plain-old data (POD)
//   2) They should have a printable tag
//   3) They should have single integer type of type value_type
//   3) They should be assignable/covertable to their value_type
//   4) They should be serializable
//   5) They should be creatable via a static create factory method
//   6) They should have a well-known invalid value
//

// Use this macro to make standard primitives
#define SAMBA_MAKE_PRIMITIVE(primitive,output,value_t)              \
  struct primitive                                                  \
  {                                                                 \
    struct tag                                                      \
    {                                                               \
      friend inline std::ostream& operator<<(std::ostream& out,tag) \
      { return out << #output; }                                    \
    };                                                              \
                                                                    \
    typedef value_t value_type;                                     \
                                                                    \
    static const primitive invalid()                                \
    {                                                               \
      primitive d = {samba::integer_max<value_type>::value};        \
      return d;                                                     \
    }                                                               \
                                                                    \
    static const primitive create(value_type v)                     \
    {                                                               \
      primitive d = {v};                                            \
      return d;                                                     \
    }                                                               \
                                                                    \
    value_type operator()() const { return m_value; }               \
                                                                    \
    SAMBA_PRIMITIVE_COMPARISONS(primitive,value_t)                  \
    SAMBA_ARITHMETIC_OPERATORS(primitive,value_t)                   \
                                                                    \
    primitive& operator=(value_type v)                              \
    { m_value = v; return *this; }                                  \
                                                                    \
    template <typename Archive>                                     \
    void serialize(Archive &ar, const unsigned version)             \
    { ar & m_value; }                                               \
                                                                    \
    value_type m_value;                                             \
  }

// Use this macro to make bit-restricted primitives; these are primitives
// that do not use all bits in their value_type;
#define SAMBA_MAKE_BIT_RESTRICTED_PRIMITIVE(primitive,output,value_t,n_bits)  \
  struct primitive                                                      \
  {                                                                     \
    struct tag                                                          \
    {                                                                   \
      friend inline std::ostream& operator<<(std::ostream& out,tag)     \
      { return out << #output; }                                        \
    };                                                                  \
                                                                        \
    typedef value_t value_type;                                         \
    static const int num_bits = n_bits;                                 \
                                                                        \
    static const primitive invalid()                                    \
    {                                                                   \
      primitive d = {(1<<num_bits)-1};                                  \
      return d;                                                         \
    }                                                                   \
                                                                        \
    static const primitive create(value_type v)                         \
    {                                                                   \
      primitive d = {v};                                                \
      primitive i = invalid();                                          \
      return d() < i() ? d : i;                                         \
    }                                                                   \
                                                                        \
    value_type operator()() const { return m_value; }                   \
                                                                        \
    SAMBA_PRIMITIVE_COMPARISONS(primitive,value_t)                      \
    SAMBA_ARITHMETIC_OPERATORS(primitive,value_t)                       \
                                                                        \
    primitive& operator=(value_type v)                                  \
    {                                                                   \
      m_value = v < invalid()() ? v : invalid()();                      \
      return *this;                                                     \
    }                                                                   \
                                                                        \
    template <typename Archive>                                         \
    void serialize(Archive &ar, const unsigned version)                 \
    { ar & m_value; }                                                   \
                                                                        \
    value_type m_value;                                                 \
  }

// Use this macro to make primitives that can be used in unordered maps.
// This can be useful when using multiple primitives in the same maps e.g.
// we don't want entity_rank(8) to hash to the same value as entity_topology(8).
#define SAMBA_MAKE_PRIMITIVE_WITH_HASHABLE_TAG(primitive,output,value_t,tag_value) \
  struct primitive                                                                 \
  {                                                                                \
    struct tag                                                                     \
    {                                                                              \
      typedef int value_type;                                                      \
      static const int value = tag_value;                                          \
      friend inline std::ostream& operator<<(std::ostream& out,tag)                \
      { return out << #output; }                                                   \
    };                                                                             \
                                                                                   \
    typedef value_t value_type;                                                    \
                                                                                   \
    static const primitive invalid()                                               \
    {                                                                              \
      primitive d = {samba::integer_max<value_type>::value};                       \
      return d;                                                                    \
    }                                                                              \
                                                                                   \
    static const primitive create(value_type v)                                    \
    {                                                                              \
      primitive d = {v};                                                           \
      return d;                                                                    \
    }                                                                              \
                                                                                   \
    value_type operator()() const { return m_value; }                              \
                                                                                   \
    SAMBA_PRIMITIVE_COMPARISONS(primitive,value_t)                                 \
    SAMBA_ARITHMETIC_OPERATORS(primitive,value_t)                                  \
                                                                                   \
    primitive& operator=(value_type v)                                             \
    { m_value = v; return *this; }                                                 \
                                                                                   \
    template <typename Archive>                                                    \
    void serialize(Archive &ar, const unsigned version)                            \
    { ar & m_value; }                                                              \
                                                                                   \
    value_type m_value;                                                            \
  }

#endif //SAMBA_SAMBA_UTILITY_MAKE_PRIMITIVE_HPP
