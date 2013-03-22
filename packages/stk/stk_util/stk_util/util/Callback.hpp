/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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
  virtual void operator()() const {
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
