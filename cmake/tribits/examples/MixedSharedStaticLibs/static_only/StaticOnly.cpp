#include "StaticOnly.hpp"
#include "SharedOnly.hpp"

std::string MixedSharedStaticLibs::staticPassThrough(const std::string &str)
{
  return MixedSharedStaticLibs::sharedPassThrough(str);
}
