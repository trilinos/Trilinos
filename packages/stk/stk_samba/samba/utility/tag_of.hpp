#ifndef SAMBA_SAMBA_UTILITY_TAG_OF_HPP
#define SAMBA_SAMBA_UTILITY_TAG_OF_HPP

namespace samba {

template <typename T>
struct tag_of
{ typedef typename T::tag type; };

} //namespace samba

#endif //SAMBA_SAMBA_UTILITY_TAG_OF_HPP

