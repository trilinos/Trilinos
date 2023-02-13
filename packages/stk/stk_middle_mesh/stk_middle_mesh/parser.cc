#include "parser.h"
#include <algorithm>
#include <cctype>
#include <iostream>
#include <stdexcept>

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

std::vector<std::string> split_semicolons(const std::string& names)
{
  if (names.size() == 0)
    throw std::runtime_error("names is empty");

  std::vector<std::string> namePairs;
  unsigned int prevIdx = 0;
  for (unsigned int i = 0; i < names.size(); ++i)
    if (names[i] == ';')
    {
      if (i - prevIdx == 0)
        throw std::runtime_error("detected empty name pair");

      namePairs.push_back(names.substr(prevIdx, i - prevIdx));
      prevIdx = i + 1;
    }

  // allow the last string to not end in a semicolon
  if (prevIdx != names.size())
  {
    namePairs.push_back(names.substr(prevIdx, names.size()));
  }

  return namePairs;
}

std::string trim_whitespace(const std::string& str)
{
  std::cout << std::boolalpha;
  bool foundNonspace    = false;
  unsigned int idxStart = 0, idxEnd = str.size() - 1;
  for (unsigned int i = 0; i < str.size(); ++i)
  {
    if (!std::isspace(static_cast<unsigned char>(str[i])))
    {
      foundNonspace = true;
      idxStart      = i;
      break;
    }
  }

  const unsigned int endVal = -1;
  for (unsigned int i = str.size() - 1; i != endVal; --i)
  {
    if (!std::isspace(static_cast<unsigned char>(str[i])))
    {
      idxEnd = i;
      break;
    }
  }

  if (!foundNonspace)
    throw std::runtime_error("sideset name is entirely whitespace");

  return str.substr(idxStart, idxEnd - idxStart + 1);
}

std::vector<mesh::impl::MeshInput::NamePair> split_colons(const std::vector<std::string>& names)
{
  std::vector<mesh::impl::MeshInput::NamePair> namePairs;

  for (auto& name : names)
  {
    auto n = std::count(name.begin(), name.end(), ':');
    if (n > 1)
      throw std::runtime_error("found more than one colon separator in name pair: " + name);
    else if (n < 1)
      throw std::runtime_error("found no colon separator in name pair: " + name);

    auto idx   = name.find(':');
    auto name1 = trim_whitespace(name.substr(0, idx));
    auto name2 = trim_whitespace(name.substr(idx + 1, name.size()));
    namePairs.emplace_back(name1, name2);
  }

  return namePairs;
}

std::vector<mesh::impl::MeshInput::NamePair> parse_sideset_names(const std::string& names)
{
  auto namesSplit = split_semicolons(names);
  return split_colons(namesSplit);
}

} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
