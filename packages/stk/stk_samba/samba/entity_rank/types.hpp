#ifndef SAMBA_SAMBA_ENTITY_RANK_TYPES_HPP
#define SAMBA_SAMBA_ENTITY_RANK_TYPES_HPP

namespace samba {

template <entity_rank::value_type Rank>
inline std::ostream& operator<<(std::ostream& out, entity_rank::rank_type<Rank>)
{ return out << static_cast<int>(Rank); }

inline std::ostream& operator<<(std::ostream& out, entity_rank::node_type)
{ return out << "node"; }

inline std::ostream& operator<<(std::ostream& out, entity_rank::edge_type)
{ return out << "edge"; }

inline std::ostream& operator<<(std::ostream& out, entity_rank::face_type)
{ return out << "face"; }

inline std::ostream& operator<<(std::ostream& out, entity_rank::element_type)
{ return out << "element"; }

inline std::ostream& operator<<(std::ostream& out, entity_rank::invalid_type)
{ return out << "invalid"; }

} //namespace samba

#endif //SAMBA_SAMBA_ENTITY_RANK_TYPES_HPP

