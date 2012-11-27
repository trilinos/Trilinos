#ifndef STKTOPOLOGY_PRETTY_PRINT_HPP
#define STKTOPOLOGY_PRETTY_PRINT_HPP

#include <stk_topology/topology.hpp>
#include <iosfwd>

namespace stk {

std::ostream & operator<<(std::ostream &out, topology::rank_t r);
std::ostream & operator<<(std::ostream &out, topology t);

void verbose_print_topology(std::ostream &out, topology t);

} //namespace stk

#endif //STKTOPOLOGY_PRETTY_PRINT_HPP
