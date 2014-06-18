#ifndef SAMBA_SAMBA_UTILITY_PRIMITIVE_ARITHMETIC_OPERATORS_HPP
#define SAMBA_SAMBA_UTILITY_PRIMITIVE_ARITHMETIC_OPERATORS_HPP

#define SAMBA_ARITHMETIC_OPERATORS(primitive,value_type)            \
  /* ++p */                                                         \
 friend inline primitive & operator++(primitive & p)                \
 { ++p.m_value; return p; }                                         \
  /* p++ */                                                         \
 friend inline primitive  operator++(primitive & p, int)            \
 { primitive tmp = p; ++p.m_value; return tmp; }                    \
  /* --p */                                                         \
 friend inline primitive & operator--(primitive & p)                \
 { --p.m_value; return p; }                                         \
  /* p-- */                                                         \
 friend inline primitive  operator--(primitive & p, int)            \
 { primitive tmp = p; --p.m_value; return tmp; }                    \
  /* p += q */                                                      \
 friend inline primitive & operator+=(primitive & p, primitive q)   \
 { p.m_value += q.m_value; return p; }                              \
  /* p += v */                                                      \
 friend inline primitive & operator+=(primitive & p, value_type v)  \
 { p.m_value += v; return p; }                                      \
  /* p -= q */                                                      \
 friend inline primitive & operator-=(primitive & p, primitive q)   \
 { p.m_value -= q.m_value; return p; }                              \
  /* p -= v */                                                      \
 friend inline primitive & operator-=(primitive & p, value_type v)  \
 { p.m_value -= v; return p; }                                      \
  /* p /= q */                                                      \
 friend inline primitive & operator/=(primitive & p, primitive q)   \
 { p.m_value /= q.m_value; return p; }                              \
  /* p /= v */                                                      \
 friend inline primitive & operator/=(primitive & p, value_type v)  \
 { p.m_value /= v; return p; }                                      \
  /* p *= q */                                                      \
 friend inline primitive & operator*=(primitive & p, primitive q)   \
 { p.m_value *= q.m_value; return p; }                              \
  /* p *= v */                                                      \
 friend inline primitive & operator*=(primitive & p, value_type v)  \
 { p.m_value *= v; return p; }                                      \
  /* p %= q */                                                      \
 friend inline primitive & operator%=(primitive & p, primitive q)   \
 { p.m_value %= q.m_value; return p; }                              \
  /* p %= v */                                                      \
 friend inline primitive & operator%=(primitive & p, value_type v)  \
 { p.m_value %= v; return p; }                                      \
  /* p + q */                                                       \
 friend inline primitive operator+(primitive p, primitive q)        \
 { p += q; return p; }                                              \
  /* p + v */                                                       \
 friend inline primitive operator+(primitive p, value_type v)       \
 { p += v; return p; }                                              \
  /* v + p */                                                       \
 friend inline primitive operator+(value_type v, primitive p)       \
 { p += v; return p; }                                              \
  /* p - q */                                                       \
 friend inline primitive operator-(primitive p, primitive q)        \
 { p -= q; return p; }                                              \
  /* p - v */                                                       \
 friend inline primitive operator-(primitive p, value_type v)       \
 { p -= v; return p; }                                              \
  /* v - p */                                                       \
 friend inline primitive operator-(value_type v, primitive p)       \
 { p -= v; return p; }                                              \
  /* p * q */                                                       \
 friend inline primitive operator*(primitive p, primitive q)        \
 { p *= q; return p; }                                              \
  /* p * v */                                                       \
 friend inline primitive operator*(primitive p, value_type v)       \
 { p *= v; return p; }                                              \
  /* v * p */                                                       \
 friend inline primitive operator*(value_type v, primitive p)       \
 { p *= v; return p; }                                              \
  /* p / q */                                                       \
 friend inline primitive operator/(primitive p, primitive q)        \
 { p /= q; return p; }                                              \
  /* p / v */                                                       \
 friend inline primitive operator/(primitive p, value_type v)       \
 { p /= v; return p; }                                              \
  /* v / p */                                                       \
 friend inline primitive operator/(value_type v, primitive p)       \
 { p /= v; return p; }                                              \
  /* p % q */                                                       \
 friend inline primitive operator%(primitive p, primitive q)        \
 { p %= q; return p; }                                              \
  /* p % v */                                                       \
 friend inline primitive operator%(primitive p, value_type v)       \
 { p %= v; return p; }                                              \
  /* v % p */                                                       \
 friend inline primitive operator%(value_type v, primitive p)       \
 { p %= v; return p; }


#endif //SAMBA_SAMBA_UTILITY_PRIMITIVE_ARITHMETIC_OPERATORS_HPP

