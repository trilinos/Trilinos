/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_UTIL_ENVIRONMENT_PROGRAMOPTIONS_HPP
#define STK_UTIL_ENVIRONMENT_PROGRAMOPTIONS_HPP

#include <boost/program_options.hpp>

namespace stk {

///
/// @addtogroup command_line_options_detail
/// @{
///

/**
 * @brief Function <b>get_options_description</b> is a singleton used to store the command
 * line option descriptions for the boost::program_options library.  This option
 * descriptions should be populated with options by each module using <b>Bootstrap</b>
 * object callback functions.  This allows modules to populate these prior to main's
 * execution of the boost::parse_command_line() functions.
 *
 * @return	        an <b>boost::program_options::options_description</b> reference to the
 *                      program options to be used for all command line option descriptions.
 */
boost::program_options::options_description &get_options_description();

/**
 * @brief Function <b>get_variabel_map</b> is a singleton used to store the variables parsed from
 * the line option descriptions.
 *
 * @return	        an <b>boost::program_options::variable_map</b> reference to the
 *                      program options to be used for all command line option descriptions.
 */
boost::program_options::variables_map &get_variables_map();

///
/// @}
///

} // namespace stk

#endif // STK_UTIL_ENVIRONMENT_PROGRAMOPTIONS_HPP
