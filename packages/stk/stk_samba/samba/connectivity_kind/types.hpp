#ifndef SAMBA_SAMBA_CONNECTIVITY_KIND_TYPES_HPP
#define SAMBA_SAMBA_CONNECTIVITY_KIND_TYPES_HPP

namespace samba {

template <connectivity_kind::value_type Kind>
inline std::ostream& operator<<(std::ostream& out, connectivity_kind::kind_type<Kind>)
{ return out << static_cast<int>(Kind); }

inline std::ostream& operator<<(std::ostream& out, connectivity_kind::fixed_type)
{ return out << "fixed"; }

inline std::ostream& operator<<(std::ostream& out, connectivity_kind::dynamic_type)
{ return out << "dynamic"; }

inline std::ostream& operator<<(std::ostream& out, connectivity_kind::invalid_type)
{ return out << "invalid"; }

} //namespace samba

#endif //SAMBA_SAMBA_CONNECTIVITY_KIND_TYPES_HPP

