// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_EXCEPTIONS_HPP_
#define _ZOLTAN2_EXCEPTIONS_HPP_

/*! \file Zoltan2_Exceptions.hpp

  Exception handling macros.  Zoltan2 throws 3 types of error:

   \li \c  std::runtime_error   for an error in input
   \li \c  std::bad_alloc       for failure to allocate memory
   \li \c  std::logic_error     for an apparent bug in the code
*/

#include <stdexcept>

/*!  \brief Throw an error returned from outside the Zoltan2 library.
 *
 *   A \c try block that calls another library should be followed
 *   by this \c catch block.  The error message if any is printed and the
 *   exception is then thrown.
 */
#define Z2_THROW_OUTSIDE_ERROR(env, e) \
  catch (std::exception &e) {          \
    std::cerr<<(env).myRank_<<" "<<__FILE__<<","<<__LINE__<<","<<e.what()<<std::endl; \
    throw e; \
  }

/*! \brief  Forward an exception back through call stack.
 *
 *  A \c try block that calls another Zoltan2 function should
 *  be followed by this \c catch series.  It ensures that the
 *  original runtime, logic or bad_alloc exception gets passed up.
 *  For example, if a runtime_error is caught as a std::exception,
 *  it gets passed up as a std::exception and the \c what() message
 *  is lost. 
 */

#define Z2_FORWARD_EXCEPTIONS \
  catch (std::runtime_error &e) { throw e; } \
  catch (std::logic_error   &e) { throw e; } \
  catch (std::bad_alloc     &e) { throw e; } \
  catch (std::exception     &e) { throw e; } 

#endif
