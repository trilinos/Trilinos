// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

// Callback class basics from Herb Sutter, http://www.gotw.ca, GOTW #83.

#ifndef STK_UTIL_UTIL_Callback_hpp
#define STK_UTIL_UTIL_Callback_hpp

/**
 * @file
 *
 * Callback implements a functor to call a member function of an object.
 *
 */

namespace sierra {

template<class T>
class Callback;

/**
 * @brief Class <b>Callback</b> ...
 *
 */
template<>
class Callback<void>
{
public:
  /**
   * @brief Member function <b>operator()</b> calls the member function on the
   * object.
   *
   */
  virtual void operator()() const = 0;

  /**
   * Destroys a <b>Callback</b> instance.
   *
   */
  virtual ~Callback()
  {}
};

typedef Callback<void> CallbackBase;			///< Shorthand for Callback<void>

/**
 * @brief Class <b>Callback</b> ...
 *
 */
template<typename T>
class Callback : public Callback<void>
{
public:
  typedef void (T::*F)();				///< Member function signature

  /**
   * Creates a new <b>Callback</b> instance.
   *
   * @param t			a <b>T</b> reference to the object.
   *
   * @param f			a <b>F</b> pointer to member function to call.
   *
   */
  Callback(T &t, F f)
    : m_t(&t),
      m_f(f)
  {}

  /**
   * Destroys a <b>Callback</b> instance.
   *
   */
  virtual ~Callback()
  {}

  /**
   * @brief Member function <b>operator()</b> calls the member function on the
   * object.
   *
   */
  virtual void operator()() const override {
    (m_t->*m_f)();
  }

private:
  T *		m_t;					///< Pointer to object
  F		m_f;					///< Member function to call
};

/**
 * @brief Member function <b>create_callback</b> creates a new callback object which
 * calls the member function <b>f</b> on the object <b>t</b>.
 *
 * @param t			a <b>T</b> reference to the object that wishes to be
 *				called back.
 *
 * @return			a <b>CallbackBase</b> pointer to the newly created
 *				callback object.
 */
template<typename T>
CallbackBase *
create_callback(
  T &		t,
  void		(T::*f)())
{
  return new Callback<T>(t, f);
}

} // namespace sierra

#endif // STK_UTIL_UTIL_Callback_hpp
