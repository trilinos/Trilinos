/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/environment/ProgramOptions.hpp>

namespace stk_classic {

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

} // namespace stk_classic
