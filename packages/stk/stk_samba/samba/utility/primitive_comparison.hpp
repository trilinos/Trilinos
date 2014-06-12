#ifndef SAMBA_SAMBA_UTILITY_PRIMITIVE_COMPARISON_HPP
#define SAMBA_SAMBA_UTILITY_PRIMITIVE_COMPARISON_HPP

#define SAMBA_PRIMITIVE_COMPARISONS(primitive,value_type)     \
  /* a < b */                                                 \
  friend inline bool operator < (primitive a, primitive b)    \
  { return a() < b(); }                                       \
  friend inline bool operator < (primitive a, value_type b)   \
  { return a() < b; }                                         \
  friend inline bool operator < (value_type a, primitive b)   \
  { return a < b(); }                                         \
  /* a <= b */                                                \
  friend inline bool operator <= (primitive a, primitive b)   \
  { return a() <= b(); }                                      \
  friend inline bool operator <= (primitive a, value_type b)  \
  { return a() <= b; }                                        \
  friend inline bool operator <= (value_type a, primitive b)  \
  { return a <= b(); }                                        \
  /* a > b */                                                 \
  friend inline bool operator > (primitive a, primitive b)    \
  { return a() > b(); }                                       \
  friend inline bool operator > (primitive a, value_type b)   \
  { return a() > b; }                                         \
  friend inline bool operator > (value_type a, primitive b)   \
  { return a > b(); }                                         \
  /* a >= b */                                                \
  friend inline bool operator >= (primitive a, primitive b)   \
  { return a() >= b(); }                                      \
  friend inline bool operator >= (primitive a, value_type b)  \
  { return a() >= b; }                                        \
  friend inline bool operator >= (value_type a, primitive b)  \
  { return a >= b(); }                                        \
  /* a == b */                                                \
  friend inline bool operator == (primitive a, primitive b)   \
  { return a() == b(); }                                      \
  friend inline bool operator == (primitive a, value_type b)  \
  { return a() == b; }                                        \
  friend inline bool operator == (value_type a, primitive b)  \
  { return a == b(); }                                        \
  /* a != b */                                                \
  friend inline bool operator != (primitive a, primitive b)   \
  { return a() != b(); }                                      \
  friend inline bool operator != (primitive a, value_type b)  \
  { return a() != b; }                                        \
  friend inline bool operator != (value_type a, primitive b)  \
  { return a != b(); }


#endif //SAMBA_SAMBA_UTILITY_PRIMITIVE_COMPARISON_HPP


