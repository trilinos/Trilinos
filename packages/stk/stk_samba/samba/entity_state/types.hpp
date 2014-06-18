#ifndef SAMBA_SAMBA_ENTITY_STATE_TYPES_HPP
#define SAMBA_SAMBA_ENTITY_STATE_TYPES_HPP

namespace samba {

template <entity_state::value_type State>
inline std::ostream& operator<<(std::ostream& out, entity_state::state_type<State>)
{ return out << static_cast<int>(State); }

inline std::ostream& operator<<(std::ostream& out, entity_state::universe_type)
{ return out << "universe"; }

inline std::ostream& operator<<(std::ostream& out, entity_state::modified_type)
{ return out << "modified"; }

inline std::ostream& operator<<(std::ostream& out, entity_state::owned_type)
{ return out << "owned"; }

inline std::ostream& operator<<(std::ostream& out, entity_state::shared_type)
{ return out << "shared"; }

inline std::ostream& operator<<(std::ostream& out, entity_state::ghosted_type)
{ return out << "ghosted"; }

inline std::ostream& operator<<(std::ostream& out, entity_state::invalid_type)
{ return out << "invalid"; }

} //namespace samba

#endif //SAMBA_SAMBA_ENTITY_STATE_TYPES_HPP

