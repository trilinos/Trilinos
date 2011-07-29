// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_PARAMETERS_HPP_
#define _ZOLTAN2_PARAMETERS_HPP_

/*! \file Zoltan2_Parameters.hpp

  Macros describing parameter values.
*/

// Definitions for ERROR_CHECK_LEVEL parameter.

/*!  We should always check basic assertions.
*/
#define Z2_BASIC_ASSERTION      0

/*!  Extra, more expensive level of checking.
 *
 * A parameter will state whether "extra" checking should be
 * done.  An example of extra checking is checking that an
 * input graph is valid.  
 */
#define Z2_COMPLEX_ASSERTION    1

/*!  Even more extensive checking.
 *
 * This is extra checking we would do when debugging
 * a problem.
 *
 */
#define Z2_DEBUG_MODE_ASSERTION  2

#define Z2_MAX_CHECK_LEVEL Z2_DEBUG_MODE_ASSERTION  

#endif
