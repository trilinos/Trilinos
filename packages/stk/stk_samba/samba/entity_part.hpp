#ifndef SAMBA_SAMBA_ENTITY_PART_HPP
#define SAMBA_SAMBA_ENTITY_PART_HPP

#include <samba/utility/hash_value.hpp>
#include <samba/utility/debug_message.hpp>

#include <samba/entity_rank.hpp>
#include <samba/entity_topology.hpp>
#include <samba/entity_state.hpp>
#include <samba/entity_block_key.hpp>

#include <boost/assert.hpp>

namespace samba {

/**
 * Used internally by partitions to manage what parts it's in and why.
 *
 * TODO: If clients aren't supposed to create these, why is the header
 * in the public area?
 */
class entity_part
{

 public:
  enum which_type{ Rank, Topology, State, Set, Invalid };

  entity_part()
    : m_which(Invalid)
    , m_value(integer_max<uint64_t>::value)
    , m_induced(false)
  {}

  entity_part(entity_rank r)
    : m_which(Rank)
    , m_rank(r)
    , m_induced(false)
  {}

  entity_part(entity_topology t, bool induced = false)
    : m_which(Topology)
    , m_topology(t)
    , m_induced(induced)
  {}

  entity_part(entity_state s)
    : m_which(State)
    , m_state(s)
    , m_induced(false)
  {}

  entity_part(entity_block_key s, bool induced = false)
    : m_which(Set)
    , m_set(s)
    , m_induced(induced)
  {}

  which_type which() const { return m_which; }

  bool induced() const { return m_induced; }

  friend inline bool equivalent(entity_part l, entity_part r)
  {
    if (l.m_which != r.m_which) return false;
    switch(l.m_which)
    {
      case Rank:     return l.m_rank == r.m_rank;
      case Topology: return l.m_topology == r.m_topology;
      case State:    return l.m_state == r.m_state;
      case Set:      return l.m_set == r.m_set;
      case Invalid:  return true;
      default:       BOOST_ASSERT_MSG(false, debug_message() << "Unhandled part-type: " << l.m_which);
    }
    // will never get here, but makes compiler happy
    return false;
  }

  friend inline std::ostream& operator<<(std::ostream& out, entity_part t)
  {
    if(t.induced())  out << "{induced->";

    switch(t.m_which)
    {
      case Rank:     out << t.m_rank; break;
      case Topology: out << t.m_topology; break;
      case State:    out << t.m_state; break;
      case Set:      out << t.m_set; break;
      case Invalid:  out << "{entity_part:invalid}"; break;
      default:       BOOST_ASSERT_MSG(false, debug_message() << "Unhandled part-type: " << t.m_which);
    }

    if (t.induced()) out << "}";

    return out;
  }

  friend inline bool operator< (entity_part l, entity_part r )
  {
    if (l.induced() == r.induced()) {
      if (l.m_which == r.m_which) {
        switch(l.m_which)
        {
          case Rank:     return l.m_rank < r.m_rank;
          case Topology: return l.m_topology < r.m_topology;
          case State:    return l.m_state < r.m_state;
          case Set:      return l.m_set < r.m_set;
          case Invalid:  return false;
          default:       BOOST_ASSERT_MSG(false, debug_message() << "Unhandled part-type: " << l.m_which);
        }
      }
      return l.m_which < r.m_which;
    }
    return l.induced() < r.induced();
  }

  friend inline bool operator== (entity_part l, entity_part r )
  {
    if (l.induced() == r.induced())
      return equivalent(l,r);
    return false;
  }

  friend inline bool operator<=(entity_part l, entity_part r )
  { return (l < r) || (l == r); };

  friend inline bool operator> (entity_part  l, entity_part r )
  { return (r > l); }

  friend inline bool operator>=(entity_part const& l, entity_part const& r )
  { return (r < l) || (l == r); };

  friend inline bool operator!=(entity_part const& l, entity_part const& r )
  { return !(l == r); }

  friend inline size_t hash_value(entity_part part)
  {
    size_t seed = 0;
    boost::hash_combine(seed,part.induced());
    switch(part.m_which)
    {
      case Rank:     boost::hash_combine(seed, part.m_rank);     break;
      case Topology: boost::hash_combine(seed, part.m_topology); break;
      case State:    boost::hash_combine(seed, part.m_state);    break;
      case Set:      boost::hash_combine(seed, part.m_set);      break;
      case Invalid:  seed = 0;                                   break;
      default:       BOOST_ASSERT_MSG(false, debug_message() << "Unhandled part-type: " << part.m_which);
    }
    return seed;
  }

  entity_rank rank() const
  {
    BOOST_ASSERT_MSG(m_which == Rank,
                     debug_message() << "Can only call rank() on entity_parts of type Rank; "
                     << "this entity_part has type: " << m_which);
    return m_rank;
  }

  entity_topology topology() const
  {
    BOOST_ASSERT_MSG(m_which == Topology,
                     debug_message() << "Can only call topology() on entity_parts of type Topology; "
                     << "this entity_part has type: " << m_which);
    return m_topology;
  }

  entity_state state() const
  {
    BOOST_ASSERT_MSG(m_which == State,
                     debug_message() << "Can only call state() on entity_parts of type State; "
                     << "this entity_part has type: " << m_which);
    return m_state;
  }

  entity_block_key set() const
  {
    BOOST_ASSERT_MSG(m_which == Set,
                     debug_message() << "Can only call set() on entity_parts of type Set; "
                     << "this entity_part has type: " << m_which);
    return m_set;
  }

 private:

  which_type m_which;
  union {
    uint64_t         m_value;
    entity_rank      m_rank;
    entity_topology  m_topology;
    entity_state     m_state;
    entity_block_key m_set;
  };

  // This will be set to true if the partition has this part because its entities
  // were induced into the part. It is critical that this information be maintained
  // for part-ordering and so we can know if entities were induced into a part or
  // manually added.
  bool m_induced;
};

namespace detail {

template <typename T> struct convertable_to_entity_part
  : public boost::mpl::false_ {};

template <> struct convertable_to_entity_part<entity_part>
  : public boost::mpl::true_ {};

template <> struct convertable_to_entity_part<entity_rank>
  : public boost::mpl::true_ {};

template <> struct convertable_to_entity_part<entity_topology>
  : public boost::mpl::true_ {};

template <> struct convertable_to_entity_part<entity_state>
  : public boost::mpl::true_ {};

template <> struct convertable_to_entity_part<entity_block_key>
  : public boost::mpl::true_ {};

} //namespace detail

} //namespace samba

#endif  //SAMBA_SAMBA_ENTITY_PART_HPP
