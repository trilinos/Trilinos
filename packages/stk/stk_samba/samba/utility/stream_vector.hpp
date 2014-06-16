#ifndef SAMBA_SAMBA_UTILITY_STREAM_VECTOR_HPP
#define SAMBA_SAMBA_UTILITY_STREAM_VECTOR_HPP

#include <iostream>

namespace samba {

template <typename DataType>
std::ostream &streamit(std::ostream &os, const std::vector<DataType> &vec, size_t max_elts = 100)
{
  bool truncated = false;
  size_t num_elts = vec.size();
  size_t num_to_output = num_elts;
  if (num_to_output > max_elts)
  {
    truncated = true;
    num_to_output = max_elts;
  }
  os << "{";
  for (size_t i = 0; i < num_to_output; ++i)
  {
    os << " " << vec[i];
  }
  if (truncated)
  {
    os << "\n VECTOR OUTPUT TRUNCATED AT " << max_elts << " ELEMENTS\n";
  }
  os << "}";

  return os;
}

} // namespace samba

#endif
