// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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

#ifndef stk_expreval_function_hpp
#define stk_expreval_function_hpp

#include <string>
#include <limits>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <stdexcept>
#include <cctype>

#include <stk_util/util/string_case_compare.hpp>

namespace stk {
namespace expreval {

/**
 * @brief Class <b>CFunctionBase</b> is a base class for calling a function from the
 * expression evaluator.  Classes which inherit from this function must implement a
 * <b>operator()</b> function which accepts an argument count and an array of length
 * that count of pointer to the argument values.  This function should perform whatever
 * argument manipulation is required to call the actual implementation of the function.
 *
 */
class CFunctionBase
{
public:
  /**
   * Creates a new <b>CFunctionBase</b> instance.
   *
   * @param arg_count		an <b>int</b> value of the number of arguments
   *				expected when <b>operator()</b> is called.
   */
  explicit CFunctionBase(int arg_count)
    : m_argCount(arg_count)
  {}

  virtual ~CFunctionBase()
  {}

  /**
   * @brief Member function <b>operator()</b> is the pure virtual interface to the
   * function.  This function is overloaded to implement the calculation to be performed
   * when the function is called.
   *
   * @param argc		an <b>int</b> value of the number of arguments.
   *
   * @param argv		a <b>double</b> pointer to an array of double
   *				pointers to the actual function arguments.
   *
   * @return			a <b>double</b> value of the result of the execution
   *				of the function.
   */
  virtual double operator()(int argc, const double *argv) = 0;

  /**
   * @brief Member function <b>getArgCount</b> returns the argument count for this
   * function.
   *
   * @return			an <b>int</b> value of the argument count for this
   *				function.
   */
  int getArgCount() const {
    return m_argCount;
  }

private:
  int		m_argCount;			///< Argument count for this function
};


/**
 * @brief Class <b>CFunction</b> declares a base template for template instantiation
 * of C-like called functions.
 *
 */
template <class S>
class CFunction;

/**
 * @brief Class <b>CFunctionMap</b> maps (function names, num arguments) to c-style function pointers via
 * the <b>CFunctionBase</b> base class.  The mapping is case insensitive.
 *
 */
class CFunctionMap : public std::multimap< std::string, CFunctionBase *, LessCase>
{
public:
  /**
   * Creates a new <b>CFunctionMap</b> instance.
   *
   * The <b>CFunctionMap</b> is populated with common C functions on construction.
   *
   */
  CFunctionMap();

  /**
   * Destroys a <b>CFunctionMap</b> instance.
   *
   */
  ~CFunctionMap();
};

/**
 * @brief Member function <b>getCFunctionMap</b> returns a reference to the
 * C-style function map.
 *
 * @return			a <b>CFunctionMap</b> reference to the C-style
 *				function map.
 */
CFunctionMap &getCFunctionMap();

/**
 * @brief Member function <b>addFunction</b> adds a C-style function to the
 * C-style function map.
 *
 * @param name		a <b>std::string</b> const reference to the function
 *				name.
 *
 * @param function		a <b>CFunctionBase</b> pointer to the C-style
 *				function.
 */
inline
void addFunction(const std::string &name, CFunctionBase *function) {
  getCFunctionMap().insert(std::make_pair(name, function));
}

} // namespace expreval
} // namespace stk

#endif // stk_expreval_function_hpp
