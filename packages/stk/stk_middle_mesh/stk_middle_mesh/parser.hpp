#ifndef PARSER_H
#define PARSER_H

#include "mesh_input.hpp"
#include <string>
#include <vector>

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

// splits the input string at the semicolons (removing the semicolons)
// If there are no semicolons, the entire string is returned
std::vector<std::string> split_semicolons(const std::string& names);

// splits each string at the colon (removing the colon from output strings)
// Throws an error if there is not exactly one colon in each string
std::vector<mesh::impl::MeshInput::NamePair> split_colons(const std::vector<std::string>& names);

// removes leading and trailing whitespace from the string (as defined by
// std::isspace())
std::string trim_whitespace(const std::string& str);

// parses a string containing all the sideset names into a vector of
// pairs of sideset names.  The requirements are:
//   1. pairs of sideset names are separated a semicolons
//   2. the two names within the pair are separated by a colon
//
//  Whitespace around each sideset name will be removed.
//
//  Examples:
//    1.  "foo:bar" parses to the std::vector<NamePair>{ {"foo", "bar"}}
//    2.  "foo : bar" is the same as the above
//    3.  "foo : bar;abc : def" parses to std::vector<NamePair>{ {"foo", "bar"}, {"abc", def"}}
//    4. "foo : bar ; abc : def" is the same as the abovea
//
//  Errors are thrown if sideset names are empty
std::vector<mesh::impl::MeshInput::NamePair> parse_sideset_names(const std::string& names);

} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
#endif
