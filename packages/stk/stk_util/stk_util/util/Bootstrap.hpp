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

#ifndef STK_UTIL_UTIL_BOOTSTRAP_HPP
#define STK_UTIL_UTIL_BOOTSTRAP_HPP

namespace stk {

///
/// @addtogroup bootstrap_detail
/// @{
///

/**
 * @brief Class <b>Bootstrap</b> serves as a bootstrapping mechanism for products in the
 * sierra toolkit and elsewhere.
 *
 * Often, it is convenient to have a product perform registrations and other operations
 * when linked into an application.  One method of accomplishing this is to utilize a
 * static object whose constructor perform these operations.  However, static
 * constructions are executed prior to main and in non-deterministc order.  The later is
 * particularly problematic if the constructor results in the usage of another static
 * object which may not have been constructed yet.
 *
 * So, the Bootstrap class creates a stack of callback functions that are executed, by
 * main(), when Bootstrap::bootstrap() is called.  These functions are still executed in a
 * non-deterministic order, but all static constructions have occurred.
 *
 */
class Bootstrap
{
public:

  typedef void (*FunctionPtr)();

  /**
   * @brief Member function <b>bootstrap</b> runs through the stored bootstrap function
   * pointers and executes each function.
   *
   */
  static void bootstrap();
  
  /**
   * @brief Creates a new <b>Bootstrap</b> instance.
   *
   * The instance serves only to insert the specified function pointer to the bootstrap
   * function list.  If the bootstrapper has already been executed, the function is added
   * to the list and executed immediately.
   *
   */
  Bootstrap(void (*f)());

private:
  Bootstrap(const Bootstrap &);
  Bootstrap &operator=(const Bootstrap &);

public:
  ~Bootstrap()
  {}

private:
  static Bootstrap *    s_front;                        ///< Front of bootstrap registry
  static bool           s_bootstrapped;                 ///< Bootstart::bootstrap has been called
  
  Bootstrap *           m_next;                         ///< Next bootstrap object
  FunctionPtr           m_f;                            ///< Function to call during bootstrap()
};

///
/// @}
///

} // namespace stk

#endif // STK_UTIL_UTIL_BOOTSTRAP_HPP
