// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#ifndef _ZOLTAN2_EXCEPTIONS_HPP_
#define _ZOLTAN2_EXCEPTIONS_HPP_

/*! \file Zoltan2_Exceptions.hpp
    \brief Defines exception handling macros.
*/

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>

/*!  \brief Throw an error returned from outside the Zoltan2 library.
 *
 *   A \c try block that calls another library should be followed
 *   by this \c catch block.  The error message if any is printed and the
 *   exception is then thrown.
 */
#define Z2_THROW_OUTSIDE_ERROR(env) \
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

/*! \brief Throw an error when experimental code is requested but not compiled.
 *
 *  Experimental code must be enabled with CMAKE Option
 *  -D Zoltan2_ENABLE_Experimental:BOOL=ON
 *  If it is not enabled but it is called, throw an exception.
 *  The input string mystr is a use-specific message included in the throw
 *  message.
 */

#define Z2_THROW_EXPERIMENTAL(mystr) \
  { \
  std::ostringstream oss; \
  oss << (mystr) << std::endl \
      << "To experiment with this software, configure with " \
      << "-D Zoltan2_ENABLE_Experimental:BOOL=ON " \
      << std::endl; \
  throw std::runtime_error(oss.str()); \
  }

#endif
