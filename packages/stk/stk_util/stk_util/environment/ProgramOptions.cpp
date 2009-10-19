#include <stk_util/environment/ProgramOptions.hpp>

namespace stk {

boost::program_options::options_description &
get_options_description()
{
  static boost::program_options::options_description s_optionsDescription;

  return s_optionsDescription;
}

boost::program_options::variables_map &
get_variables_map()
{
  static boost::program_options::variables_map s_variablesMap;

  return s_variablesMap;
}

} // namespace stk
