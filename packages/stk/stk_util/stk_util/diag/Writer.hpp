/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_UTIL_DIAG_WRITER_HPP
#define STK_UTIL_DIAG_WRITER_HPP

#include <ostream>
#include <string>
#include <vector>
#include <utility>
#include <stdint.h>

#include <stk_util/diag/Writer_fwd.hpp>

namespace stk {
namespace diag {

///
/// @addtogroup diag_writer_detail
/// @{
///

class WriterThrowSafe
{
public:
  explicit WriterThrowSafe(Writer &writer);

  ~WriterThrowSafe();

private:
  Writer &      m_writer;
  int           m_depth;

  WriterThrowSafe(const WriterThrowSafe &);
  void operator = (const WriterThrowSafe &);
};


/**
 * @brief Class Writer implements a runtime selectable diagnostic output writer to aid in the
 * development and diagnostics of massively parallel applications.
 *
 */
class Writer
{
public:
  /**
   * @brief Enumeration <b>Flags</b>.
   *
   * <UL>
   * <LI><b>DISABLED</b> disables Writer output<BR>
   * <LI><b>ENABLED</b> enabled Writer output<BR>
   * </UL>
   */
  enum Flags {
    DISABLED    =  0x00,
    ENABLED    =  0x01
  };

private:
  /**
   * @brief Class <b>LineMaskStack</b> implements a special stack for maintaining the
   * depth and line mask for the diagnostic stream.  It has the following special
   * properties:
   *
   * <UL>
   * <LI>Never pops the top element off.
   * <LI>Pushing a new depth increments the top depth and replicates the top print mask.
   * <LI>Pushing a line mask replicates the top depth and adds the new line mask.
   * <LI>Popping does not allow the last depth/line mask to be removed
   * <LI>Reseting the depth, pops from the stack leaving a the entry on the top which is
   *     greater in depth than the previous entry's depth.
   * </UL>
   *
   */
  struct LineMaskStack : public std::vector<std::pair<int, PrintMask> >
  {
    /**
     * @brief Creates a new <b>LineMaskStack</b> instance, inserting the sentinel indentation
     * level/line mask value.
     *
     */
    LineMaskStack() {
      push_back(std::pair<int, PrintMask>(0, LOG_ALWAYS));
    }

    /**
     * @brief Member function <b>pushDepth</b> pushes the a new indentation level in the stack,
     * replicating the line mask value.
     *
     * @return      a <b>LineMaskStack</b> reference to the depth stack.
     */
    LineMaskStack &pushDepth() {
      push_back(std::make_pair(back().first + 1, back().second));
      return *this;
    }

    /**
     * @brief Member function <b>push</b> pushes a new line mask on the stack, replicating the
     * current indentation level.
     *
     * @param line_mask    a <b>PrintMask</b> value of the line mask.
     *
     * @return      a <b>LineMaskStack</b> reference to the depth stack.
     */
    LineMaskStack &push(PrintMask line_mask) {
      push_back(std::make_pair(back().first, line_mask));
      return *this;
    }

    /**
     * @brief Member function <b>pop</b> pops the top entry from the stack, leaving the sentinel
     * value if near the bottom.
     *
     * @return      a <b>LineMaskStack</b> reference to the depth stack.
     */
    LineMaskStack &pop() {
      if (size() > 1)
        pop_back();
      return *this;
    }

    /**
     * @brief Member function <b>popLineMask</b> removes the top line mask only if
     * it is not the last one of the current depth.
     *
     * @return      a <b>LineMaskStack</b> reference to the depth stack.
     */
    LineMaskStack &popLineMask() {
      if (size() > 1 && getNextDepth() == getDepth())
        pop_back();
      return *this;
    }

    /**
     * @brief Member function <b>getDepth</b> returns the depth from the top entry on the stack.
     *
     * @return      an <b>int</b> value of the depth from the top entry
     *        on the stack.
     */
    int getDepth() const {
      return back().first;
    }

    /**
     * @brief Member function <b>getLineMask</b> returns the line mask from the top entry on the
     * stack.
     *
     * @return      an <b>int</b> value of the line mask from the top
     *        entry on the stack.
     */
    int getLineMask() const {
      return back().second;
    }

    /**
     * @brief Member function <b>getNextDepth</b> returns the depth from the entry immediately prior
     * to the top entry on the stack.
     *
     * @return      an <b>int</b> returns the depth from the entry
     *        immediately prior to the top entry on the stack.
     */
    int getNextDepth() const {
      return (end() - 2)->first;
    }

    /**
     * @brief Member function <b>resetDepth</b> resets the stack so that the depth of the top entry
     * on the stack exceeds the depth of the entry prior to the top entry on the stack.
     *
     * @return      a <b>LineMaskStack</b> reference to the depth stack.
     */
    LineMaskStack &resetDepth() {
      while (size() > 1 && getNextDepth() == getDepth())
        pop_back();
      return *this;
    }
  };

public:
  /**
   * @brief Creates a new <b>Writer</b> instance with the specified print mask and output
   * flags.
   *
   * @param print_mask    a <b>PrintMask</b> value of the print mask to apply
   *        to this writer.
   *
   * @param flags    a <b>Flags</b> value of the selected output flags.
   *
   */
  explicit Writer(std::streambuf *streambuf, PrintMask print_mask = static_cast<PrintMask>(LOG_MEMBERS), Flags flags = static_cast<Flags>(ENABLED));

  /**
   * @brief Destroys a <b>Writer</b> instance.
   *
   */
  ~Writer();

  /**
   * @brief Member function <b>getStream</b> returns the output stream.
   *
   * @return      a <b>std::ostream</b> reference to the output
   *        stream.
   */
  std::ostream &getStream() {
    return m_writerStream;
  }

  /**
   * @brief Member function <b>setFlags</b> sets the flags bitmask which describes
   * the output line prefix content.
   *
   * @param flags    an <b>int</b> of the bitmask of flags.
   *
   * @return      a <b>Writer</b> reference to this object
   */
  Writer &setFlags(int flags) {
    m_flags = (Flags) flags;
    return *this;
  }

  /**
   * @brief Member function <b>getFlags</b> returns the flags bitmask.
   *
   * @return      an <b>int</b> of the flags bitmask.
   */
  int getFlags() {
    return m_flags;
  }

  int getDepth() const {
    return m_lineMaskStack.getDepth();
  }

  Writer &restoreDepth(int depth) {
    while (m_lineMaskStack.getDepth() > depth)
      pop();
    return *this;
  }

  /**
   * @brief Member function <b>setPrintMask</b> sets the print output mask.
   *
   * @param mask    an <b>PrintMask</b> value of the new print
   *        mask.
   *
   * @return      a <b>Writer</b> reference to this diagnostic
   *        writer.
   */
  Writer &setPrintMask(PrintMask mask = 0) {
    m_printMask = mask;
    return *this;
  }

  /**
   * @brief Member function <b>setPrintMask</b> sets the print output mask.
   *
   * @param mask_string    an <b>PrintMask</b> value of the new print
   *        mask.
   *
   * @return      a <b>Writer</b> reference to this diagnostic
   *        writer.
   */
  //Writer &setPrintMask(const char *mask_string);

  /**
   * @brief Member function <b>setLineMask</b> sets the line mask of this line.
   *
   * @param line_mask    an <b>PrintMask</b> of the mask for this
   *        line.
   *
   * @return      a <b>Writer</b> reference to this diagnostic
   *        writer.
   */
  Writer &setLineMask(PrintMask line_mask) {
    m_lineMaskStack.push(line_mask);

    return *this;
  }

  /**
   * @brief Member function <b>m</b> sets the line mask of this line.
   *
   * @param line_mask    an <b>PrintMask</b> of the mask for this
   *        line.
   *
   * @return      a <b>Writer</b> reference to this object.
   */
  Writer &m(PrintMask line_mask) {
    setLineMask(line_mask);

    return *this;
  }

  /**
   * @brief Member function <b>m</b> sets the line mask of this line.
   *
   * @param line_mask    an <b>PrintMask</b> of the mask for this
   *        line.
   *
   * @return      a <b>Writer</b> reference to this object.
   */
  Writer &w(bool on, PrintMask line_mask) {
    setLineMask(on ? line_mask : 0x80000000);

    return *this;
  }

  /**
   * @brief Member function <b>t</b> sets the line mask of this line to <i>line_make</i> bitwise
   * or'ed with LOG_TRACE.
   *
   * @param line_mask    an <b>PrintMask</b> of the mask for this
   *        line.
   *
   * @return      a <b>Writer</b> reference to this object.
   */
  Writer &t(PrintMask line_mask = 0) {
    setLineMask(line_mask | stk::LOG_TRACE);

    return *this;
  }

  /**
   * @brief Member function <b>getLineMask</b> returns the current line mask.
   *
   * @return      an <b>int</b> value of the current line mask.
   */
  PrintMask getPrintMask() {
    return m_printMask;
  }

  /**
   * @brief Member function <b>isEnabled</b> returns true if the ENABLED bit is set in the flags
   * bitmask.
   *
   * @return      a <b>bool</b> of true if the ENABLED bit is set in
   *        the flags bitmask.
   */
  bool isEnabled() {
    return (m_flags & ENABLED) != 0;
  }

  /**
   * @brief Member function <b>isLoggable</b> returns true if any corresponding bit in the line mask
   * matches a bit in the print mask, except LOG_TRACE which also requires isTracing() to be true.
   *
   * @return      a <b>bool</b> of true if any corresponding bit in
   *        the line mask matches a bit inthe print mask.
   */
  bool isLoggable(PrintMask line_mask) {
    return line_mask == 0            // Always
      || ((line_mask & m_printMask & stk::LOG_TRACE)      // LOG_TRACE?
          ? isTracing()              // Yes, must be tracing
//    : (line_mask & m_printMask) != 0);        // No, any matching bits
          : (line_mask & m_printMask) == line_mask);      // No, all matching bits
  }

  /**
   * @brief Member function <b>shouldPrint</b> returns true if the line should print.
   *
   * @return      a <b>bool</b> of true if this line should be printed.
   */
  bool shouldPrint() {
    return shouldPrint(m_lineMaskStack.getLineMask());
  }

  /**
   * @brief Member function <b>shouldPrint</b> returns true if the line should print.
   *
   * @param line_mask    a <b>PrintMask</b> value of the line mask.
   *
   * @return      a <b>bool</b> of true if this line should be printed.
   */
  bool shouldPrint(PrintMask line_mask) {
    return isEnabled() && isLoggable(line_mask);
  }

  /**
   * @brief Member function <b>shouldTrace</b> returns true if any corresponding bit in the line
   * mask matches a bit in the print mask, except LOG_TRACE which also requires isTracing() to be
   * true.
   *
   * @return      a <b>bool</b> of true if any corresponding bit in
   *        the line mask matches a bit in the print mask.
   */
  bool shouldTrace(int line_mask) {
    return line_mask == 0            // Always
      || (line_mask & m_printMask) != 0;        // Any set
//      || (line_mask & m_printMask) == line_mask;      // All set
  }

  /**
   * @brief Member function <b>dflush</b> flushes the output stream.
   *
   * @return      a <b>Writer</b> reference to this object
   */
  Writer &dflush();

  /**
   * @brief Member function <b>dendl</b> is a manipulator which sets the output
   * stream to a new line.<P>
   *
   * The std::endl manipulator is sent to the output stream.
   *
   * @return      a <b>Writer</b> reference to this object
   */
  Writer &dendl();

  /**
   * @brief Member function <b>push</b> is a manipulator which increases the line
   * mask depth by one.
   *
   * @return      a <b>Writer</b> reference to this object
   */
  Writer &push();

  /**
   * @brief Member function <b>pop</b> is a manipulator which decreases the line mask depth by one,
   * but not less than zero(0).
   *
   * @return      a <b>Writer</b> reference to this object
   */
  Writer &pop();

  /**
   * @brief Member function <b>pop</b> is a manipulator which decreases the line mask depth by one,
   * but not less than zero(0).
   *
   * @return      a <b>Writer</b> reference to this object
   */
  Writer &resetLineMask();

#ifndef SWIG
  /**
   * @brief Member function <b>operator<<</b> is the manipulator instantiation function
   *
   * @return      a <b>Writer</b> reference to this object
   */
  Writer& operator<<(Writer& (*f)(Writer&));

  /**
   * @brief Member function <b>operator<<</b> passes the ios_base manipulator function to the output
   * stream.
   *
   * @return      a <b>Writer</b> reference to this object
   */
  Writer& operator<<(std::ios_base& (*f)(std::ios_base&));

  /**
   * @brief Member function <b>operator<<</b> passes the iostream manipulator function to the output
   * stream.
   *
   * @return      a <b>Writer</b> reference to this object
   */
  Writer& operator<<(std::ostream& (*f)(std::ostream&));
#endif // SWIG
  /**
   * @brief Member function <b>incTraceDepth</b> increments the tracing count.
   *
   * @return      an <b>int</b> value of the new tracing count.
   */
  int incTraceDepth() {
    return ++m_traceDepth;
  }

  /**
   * @brief Member function <b>decTraceDepth</b> decrements the tracing count.
   *
   * @return      an <b>int</b> value of the new tracing count.
   */
  int decTraceDepth() {
    return --m_traceDepth;
  }

  /**
   * @brief Member function <b>isTracing</b> returns true of the trace depth is greater than zero.
   * The value of -1 is initially stored in the depth as a flag that the trace counters have never
   * been called. (This may be and obsolete requirement).
   *
   * @return      a <b>bool</b> value of true of the trace depth is greater
   *        than zero
   */
  bool isTracing() {
    return m_traceDepth <= 0 ? false
      : (m_traceDepth == 1 || (m_traceDepth > 1 && (m_printMask & stk::LOG_TRACE_SUB_CALLS)));
  }

  /**
   * @brief Member function <b>isTraceable</b> returns true if currently tracing or
   * tracing is enabled.
   *
   * @return      a <b>bool</b> of true if tracing is enabled and
   *        active, or if tracing is disabled.
   */
  bool isTraceable() {
    return isTracing() || (m_printMask & stk::LOG_TRACE) != 0; // Currently in a trace or tracing bit set
  }

private:
  Flags        m_flags;    ///< Describes the output and line prefix information to be printed.
  PrintMask      m_printMask;    ///< Print mask that the line mask must the match to print
  LineMaskStack      m_lineMaskStack;  ///< Stack of pushed line masks
  int        m_traceDepth;    ///< Trace depth
  std::ostream                  m_writerStream;
};

/**
 * @brief Writer function <b>dendl</b> calls the Writer::dendl manipulator.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to dendl.
 *
 * @return    a <b>Writer</b> reference to this object
 */
inline Writer &dendl(Writer &dout) {
  return dout.dendl();
}

/**
 * @brief Writer function <b>dflush</b> calls the Writer::dflush
 * manipulator.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to flush.
 *
 * @return    a <b>Writer</b> reference to this object
 */
inline Writer &dflush(Writer &dout) {
  return dout.dflush();
}

/**
 * @brief Function <b>push</b> calls the Writer::push manipulator.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      push.
 *
 * @return    a <b>Writer</b> reference to this object
 */
inline Writer &push(Writer &dout) {
  return dout.push();
}

/**
 * @brief Member function <b>pop</b> calls the Writer::pop manipulator.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      pop.
 *
 * @return    a <b>Writer</b> reference to this object
 */
inline Writer &pop(Writer &dout) {
  return dout.pop();
}


/**
 * @brief Class <b>_setlinemask</b> is the line mask manipulator.
 *
 */
struct _setlinemask
{
  /**
   * @brief Creates a new <b>setlinemask</b> instance.
   *
   * @param line_mask  an <b>PrintMask</b> value of the new line mask.
   */
  _setlinemask(PrintMask line_mask)
    : m_lineMask(line_mask)
  {}

  PrintMask    m_lineMask;
};

/**
 * @brief Function <code>setlinemask</code> sets the active line mask bits as a manipulator.
 *
 * @param line_mask     a <b>PrintMask</b> value of the bits to set.
 *
 */
inline _setlinemask setlinemask(PrintMask line_mask) {
  return _setlinemask(line_mask);
}

/**
 * @brief Function <b>operator<<</b> class the Writer::setLineMask manipulator.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      set the line mask.
 *
 * @param set_line_mask a <b>_setlinemask</b> value of the line mask to set.
 *
 * @return    a <b>Writer</b> reference to this object
 */
#ifndef SWIG
inline Writer &operator<<(Writer &dout, _setlinemask set_line_mask) {
  return dout.setLineMask(set_line_mask.m_lineMask);
}
#endif // SWIG

/**
 * @brief Function <b>resetlinemask</b> calls the Writer::resetLineMask manipulator.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      dendl.
 *
 * @return    a <b>Writer</b> reference to this object
 */
inline Writer &resetlinemask(Writer &dout) {
  return dout.resetLineMask();
}

/**
 * @brief Function <b>operator<<</b> writes the c sytle string to the output stream.
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to write the
 *      c style string to.
 *
 * @param c_str    a <b>char</b> const pointer to the start of the c style
 *      string.
 *
 * @return    a <b>Writer</b> reference to this object
 */
#ifndef SWIG
Writer &operator<<(Writer &dout, const char *c_str);
Writer &operator<<(Writer &dout, const std::string &str);
Writer &operator<<(Writer &dout, const void *ptr);
Writer &operator<<(Writer &dout, const float &x);
Writer &operator<<(Writer &dout, const double &x);
Writer &operator<<(Writer &dout, const long double &x);
Writer &operator<<(Writer &dout, const int &x);
Writer &operator<<(Writer &dout, const unsigned int &x);
Writer &operator<<(Writer &dout, const long &x);
Writer &operator<<(Writer &dout, const unsigned long &x);
Writer &operator<<(Writer &dout, const short &x);
Writer &operator<<(Writer &dout, const unsigned short &x);
Writer &operator<<(Writer &dout, const long long &x);
Writer &operator<<(Writer &dout, const unsigned long long &x);
#endif // SWIG
/**
 * @brief Class <b>c_ptr_</b> simply stores a pointer to an object of type T.  This
 * allows pointers which want to be deferenced if they are not null to be output using
 * operator<< on a c_ptr function.
 *
 */
template <class T>
class c_ptr_
{
public:
  /**
   * Creates a new <b>c_ptr_</b> instance.
   *
   * @param t    a <b>T</b> pointer to object
   */
  explicit c_ptr_(const T *t)
    : m_t(t)
  {}

public:
  const T *  m_t;      ///< Pointer to object
};

/**
 * Member function <b>c_ptr</b> creates a c_ptr_ object of type T ala std::make_pair.
 *
 * @param t    a <b>T</b> pointer to an object that is to be dereferenced.
 *
 * @return    a <b>c_ptr_</b> object which contains the pointer t.
 */
template <class T>
c_ptr_<T> c_ptr(const T *t) {
  return c_ptr_<T>(t);
}

/**
 * @brief Class <b>c_ptr_func_</b> simply stores a pointer to an object of type T.
 * This allows pointers which want to call the specified member function if they are not
 * null to be output using operator<< on a c_ptr_func function.
 *
 */
template <class T, typename R>
class c_ptr_func_
{
public:
  /**
   * Creates a new <b>c_ptr_func_</b> instance.
   *
   * @param t    a <b>T</b> pointer to object
   *
   * @param pmf    a <b>T::*</b> member function pointer to call
   *
   */
#ifndef SWIG
  explicit c_ptr_func_(const T *t, R (T::*pmf)() const)
    : m_t(t),
      m_pmf(pmf)
  {}
public:
  const T *  m_t;      ///< Pointer to object
  R (T::*m_pmf)() const;    ///< Function to call for dump
#endif // SWIG
};

/**
 * @brief Template function <b>c_ptr</b> creates a c_ptr_func_ object of type T ala
 * std::make_pair.  This T must implement a member function which takes no arguments and
 * returns a value of type R.
 *
 * @param t    a <b>T</b> pointer to an object that is call the specified
 *      member function.
 *
 *
 * @param pmf    a <b>T::*</b> member function pointer to call
 *
 * @return    a <b>c_ptr_</b> object which contains the pointer t and a
 *      member function whch takes no arguments.
 */
#ifndef SWIG
template <class T, typename R>
c_ptr_func_<T, R> c_ptr_func(const T *t, R (T::*pmf)() const) {
  return c_ptr_func_<T, R>(t, pmf);
}

/**
 * @brief Template function <b>operator<<</b> dereferences the <b>c_ptr_</b> object's member m_t if
 * it is not null and writes that to the diagnostic writer.  If the object's member is null, it
 * writes "<not created>".
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      write the <T> object to if the pointer to it is not null.
 *
 * @param c    a <b>c_ptr_</b> reference with a member to dereference and
 *      write to ethe diagnostic writer if not null.
 *
 * @return    a <b>Writer</b> reference to this object
 */
template <class T>
Writer &operator<<(Writer &dout, const c_ptr_<T> &c) {
  dout << "(pointer " << (void *) c.m_t << "), ";

  if (c.m_t)
    dout << *c.m_t;
  else
    dout << "<not created>";

  return dout;
}

/**
 * @brief Template function <b>operator<<</b> dereferences the <b>c_ptr_func_</b>
 * object's member m_t if it is not null and calls the m_pmf member function and writes
 * the result of that to the diagnostic writer.  If the object's member is null, it writes
 * "<not created>".
 *
 * @param dout    a <b>Writer</b> reference to the diagnostic writer to
 *      write the <b>T</b> object to if the pointer to it is not null.
 *
 * @param c    a <b>c_ptr_func_</b> reference with a member to dereference
 *      and call the member function m_pmt if m_t is not null.
 *
 * @return    a <b>Writer</b> reference to this object
 */
template <class T, typename R>
Writer &operator<<(Writer &dout, const c_ptr_func_<T, R> &c) {
  if (c.m_t)
    dout << "(pointer), " << (c.m_t->*c.m_pmf)();
  else
    dout << "(pointer), <not created>";

  return dout;
}
#endif // SWIG

///
/// @}
///

} // namespace diag
} // namespace stk

#include <stk_util/diag/WriterManip.hpp>

namespace sierra {


using stk::diag::push;
using stk::diag::pop;
using stk::diag::dendl;
using stk::diag::dflush;

namespace Diag {

using stk::diag::push;
using stk::diag::pop;
using stk::diag::dendl;
using stk::diag::dflush;
using stk::diag::setlinemask;
using stk::diag::resetlinemask;
using stk::diag::c_ptr;
using stk::diag::c_ptr_func;
using stk::diag::c_ptr_func_;

}

} // namespace sierra

#include <stk_util/diag/WriterExt.hpp>

#endif // STK_UTIL_DIAG_WRITER_HPP
