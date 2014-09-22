// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef STK_UTIL_ENVIRONMENT_PROGRAMOPTIONS_HPP
#define STK_UTIL_ENVIRONMENT_PROGRAMOPTIONS_HPP

#include <boost/program_options.hpp>

namespace stk {

///
/// @addtogroup command_line_options_detail
/// @{
///

/**
 * @brief Function <b>get_options_description</b> accesses a singleton used to store the command
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
 * @brief Function <b>get_variables_map</b> accesses a singleton used to store the variables parsed from
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
