#ifndef SAMBA_SAMBA_UTILITY_VALUE_OF_HPP
#define SAMBA_SAMBA_UTILITY_VALUE_OF_HPP

namespace samba {

template <typename T>
struct value_of
{ typedef typename T::value_type type; };

} //namespace samba

#endif //SAMBA_SAMBA_UTILITY_VALUE_OF_HPP

