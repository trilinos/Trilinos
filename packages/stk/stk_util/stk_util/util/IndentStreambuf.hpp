#ifndef STK_UTIL_UTIL_INDENTSTREAMBUF_HPP
#define STK_UTIL_UTIL_INDENTSTREAMBUF_HPP

#include <streambuf>
#include <utility>

namespace stk {

static const char PUSH = '\016';                                ///< Meta-character to increase indentation
static const char POP = '\017';                                 ///< Meta-character to decrease indentation
static const char LEFT = '\021';                                ///< Meta-character to force left justification

/**
 * @brief Class <b>basic_indent_streambuf</b> implements a output streambuf that performs
 * indentation, blank line removal and outline bracing, sending the result character stream to
 * another output stream buffer.
 *
 * When, the meta-characters PUSH, POP, and LEFT are inserted into the stream buffer, the characters
 * is discarded and the appropriate operation occurs.  In the case of PUSH, the next line will
 * indented by an additional <B>m_indentSize</B> spaces and the currently line will end with an
 * open brace ({) is BRACES are enabled.  In the case of POP, the next line will be indented by
 * <B>m_indentSize</B> spaces, it may be preceeded by a line containing a properly indented close
 * brace (}).  In the case of LEFT, the next line (current line if at the start), will not be
 * indented.
 *
 * The stream bufer can be created with BRACES and BLANK_LINES enabled or disabled.  When BRACE is
 * enabled, then indentation will be produces braces which allows editors to traverse blocks.  Blank
 * lines are normally eliminated, but may be allowed if BLANK_LINES is specified.
 *
 */
template<class Ch, class Tr = std::char_traits<Ch> >
class basic_indent_streambuf : public std::basic_streambuf<Ch, Tr>
{
public:
  enum {
    MAX_INDENT_LEVEL = 50                                       ///< Maximum indentation level
  };

  /**
   * @brief Class <b>Flags</b> ...
   *
   */
  enum Flags {
    NO_BRACES      = 0x00,                                      ///< No braces on indentation shift
    NO_BLANK_LINES = 0x00,                                      ///< No blank line are written
    BRACES         = 0x01,                                      ///< Braces on indentation shift
    BLANK_LINES    = 0x02                                       ///< Blank line are written
  };

  /**
   * Creates a new <b>basic_indent_streambuf</b> instance.
   *
   * @param indent_size		a <b>size_t</b> value of the number of spaces for each indentation
   *                            level. 
   *
   * @param flags		an <b>unsigned int</b> value of Flags to enable or disable
   *                            BLANK_LINES and BRACES.
   *
   */
  explicit basic_indent_streambuf(std::basic_streambuf<Ch, Tr> *stream_buffer, size_t indent_size = 2, unsigned flags = BRACES)
    : m_streamBuffer(stream_buffer),
      m_atLineBegin(true),
      m_leftJustify(false),
      m_indentLevel(0),
      m_nextIndentLevel(0),
      m_indentSize(indent_size),
      m_flags((Flags) flags),
      m_indentString(0)
  {
    set_indent_size(indent_size);
  }

  /**
   * Destroys a <b>basic_indent_streambuf</b> instance.
   *
   */
  virtual ~basic_indent_streambuf() {
    delete[] m_indentString;
  }

  /**
   * @brief Member function <b>redirect</b> sets the destination output stream buffer.
   *
   */
  void redirect(std::basic_streambuf<Ch, Tr> *stream_buffer) {
    m_streamBuffer = stream_buffer;
  }

  /**
   * @brief Member function <b>get_stream_buffer</b> returns the current destination output stream
   * buffer. 
   *
   * @return			a <b>std::streambuf</b> pointer to the current destination output
   *                            stream buffer.
   */
  std::streambuf *get_stream_buffer() {
    return m_streamBuffer;
  }
  
  /**
   * @brief Member function <b>set_indent_size</b> set the number of spaces to write for each
   * indentation level.
   *
   * @param indent_size		a <b>size_t</b> value of the number of spaces for each indentation
   *                            level. 
   *
   */
  void set_indent_size(size_t indent_size) {
    m_indentSize = indent_size;

    delete[] m_indentString;
    m_indentString = new Ch[MAX_INDENT_LEVEL*m_indentSize];
    std::fill(m_indentString, m_indentString + MAX_INDENT_LEVEL*m_indentSize, static_cast<Ch>(' '));
  }

  /**
   * @brief Member function <b>set_flags</b> enables or disables the BLANK_LINES and BRACES written
   * to the destination stream.
   *
   * @param flags		an <b>unsigned int</b> value of the Flags to enable/disable.
   *
   */
  void set_flags(unsigned flags) {
    m_flags = (Flags) flags;
  }


private:
  /**
   * @brief Member function <b>indent_level</b> returns the current indentation level.
   *
   * @return			a <b>size_t</b> value of the current indentation level.
   */
  size_t indent_level() {
    return std::min(m_indentLevel*m_indentSize, (size_t) MAX_INDENT_LEVEL*m_indentSize);
  }

  /**
   * @brief Member function <b>prefix</b> prints the indentation spaces only at the beginning of the
   * line and no left justification.
   *
   */
  void prefix() {
    if (m_atLineBegin) {
      if (!m_leftJustify)
        m_streamBuffer->sputn(m_indentString, indent_level());
      m_leftJustify = false;
      m_atLineBegin = false;
    }
  }

  /**
   * @brief Member function <b>next_line</b> prints the braces for indentation.
   *
   */
  void next_line() {
    if (m_nextIndentLevel > m_indentLevel) {
      if (m_flags & BRACES) 
	m_streamBuffer->sputn(" {", 2);
      m_streamBuffer->sputc(Tr::to_int_type('\n'));
    }
    else if (m_nextIndentLevel < m_indentLevel) {
      m_indentLevel = m_nextIndentLevel;
      if (!m_atLineBegin)
        m_streamBuffer->sputc(Tr::to_int_type('\n'));
      if (m_flags & BRACES) {
        m_streamBuffer->sputn(m_indentString, indent_level());
        m_streamBuffer->sputc(Tr::to_int_type('}'));
        m_streamBuffer->sputc(Tr::to_int_type('\n'));
      }
    }
    else if (!m_atLineBegin || (m_flags & BLANK_LINES))
      m_streamBuffer->sputc(Tr::to_int_type('\n'));

    m_indentLevel = m_nextIndentLevel;
    m_atLineBegin = true;
  }

public:
  /**
   * @brief Member function <b>overflow</b> interprets a meta-character or writes the specified
   * character to the destination output stream buffer.  If the character is a meta-character, the
   * aproprate action is performed.
   *
   * @param c			an <b>int</b> value of the character to write.
   *
   * @return			an <b>int</b> value of the character written.
   */
  virtual typename std::basic_streambuf<Ch, Tr>::int_type overflow(typename std::basic_streambuf<Ch, Tr>::int_type c) {
    if (c == Tr::to_int_type('\n'))
      next_line();
    else if (c == Tr::to_int_type(POP)) {
      if (m_nextIndentLevel != m_indentLevel)
	next_line();
      if (m_indentLevel > 0)
	m_nextIndentLevel = m_indentLevel - 1;
    }
    else if (c == Tr::to_int_type(PUSH)) {
      if (m_nextIndentLevel != m_indentLevel)
	next_line();
      m_nextIndentLevel = m_indentLevel + 1;
    }
    else if (c == Tr::to_int_type(LEFT)) {
      m_leftJustify = true;
    }
    else {
      prefix();
      m_streamBuffer->sputc(c);
    }

    return c;
  }

  /**
   * @brief Member function <b>xsputn</b> interprets the meta-characters or writes the specified
   * characters to the destination output stream buffer.  If a character is a meta-character, the
   * aproprate action is performed.
   *
   * @param p			a <b>Ch</b> const pointer to the character string to write.
   *
   * @param n			a <b>std::streamsize</b> value of the number of characters in the
   *                            string. 
   *
   * @return			a <b>std::streamsize</b> value of the number of characters os the
   *                            string which were interpreted or written.
   */
  virtual std::streamsize xsputn(const Ch *p,  std::streamsize n) {
    const Ch *p_end = p + n;
    for (const Ch *q = p; q != p_end; ++q) {

// If at start of line, PUSH and POP have immediate effect
      if (p == q && m_atLineBegin) {
	if (Tr::to_int_type(*p) == Tr::to_int_type('\n')) {
	  next_line();
	  ++p;
	}
	else if (Tr::to_int_type(*p) == Tr::to_int_type(POP)) {
	  ++p;
	  if (m_nextIndentLevel != m_indentLevel)
	    next_line();
	  if (m_indentLevel > 0)
	    m_nextIndentLevel = m_indentLevel - 1;
	}
	else if (Tr::to_int_type(*p) == Tr::to_int_type(PUSH)) {
	  ++p;
	  if (m_nextIndentLevel != m_indentLevel)
	    next_line();
	  m_nextIndentLevel = m_indentLevel + 1;
	}
        else if (Tr::to_int_type(*p) == Tr::to_int_type(LEFT)) {
          ++p;
          m_leftJustify = true;
        }
      }

// If not at start, PUSH and POP are for the next line.
      else {
	if (Tr::to_int_type(*q) == Tr::to_int_type('\n')) {
	  prefix();
	  m_streamBuffer->sputn(p, q - p);
	  next_line();
	  p = q + 1;
	}
	else if (Tr::to_int_type(*q) == Tr::to_int_type(POP)) {
	  prefix();
	  m_streamBuffer->sputn(p, q - p);
	  p = q + 1;
	  if (m_nextIndentLevel != m_indentLevel)
	    next_line();
	  if (m_indentLevel > 0)
	    m_nextIndentLevel = m_indentLevel - 1;
	}
	else if (Tr::to_int_type(*q) == Tr::to_int_type(PUSH)) {
	  prefix();
	  m_streamBuffer->sputn(p, q - p);
	  p = q + 1;
	  if (m_nextIndentLevel != m_indentLevel)
	    next_line();
	  m_nextIndentLevel = m_indentLevel + 1;
	}
        else if (Tr::to_int_type(*q) == Tr::to_int_type(LEFT)) {
          m_leftJustify = true;
	  p = q + 1;
        }
      }
    }
    if (p != p_end) {
      prefix();
      m_streamBuffer->sputn(p, p_end - p);
      m_atLineBegin = false;
    }

    return n;
  }

  /**
   * @brief Member function <b>sync</b> syncs the destination output stream buffer. 
   *
   * @return			an <b>int</b> value result of the pub sync operation.
   */
  virtual int sync() {
    return m_streamBuffer->pubsync();
  }

private:
  basic_indent_streambuf(const basic_indent_streambuf &);
  basic_indent_streambuf &operator=(const basic_indent_streambuf &);
  
private:
  std::streambuf *	m_streamBuffer;         ///< Pointer to destination output stream buffer
  bool			m_atLineBegin;          ///< Flag indicating at beginning of line
  bool                  m_leftJustify;          ///< Flag indicating next (current) line is to be left justified
  size_t		m_indentLevel;          ///< Current indentation level
  size_t		m_nextIndentLevel;      ///< Next line indentation level
  size_t		m_indentSize;           ///< Number of characters to print for each indentation level
  Flags			m_flags;                ///< BLANK_LINES and BRACES flags
  Ch *			m_indentString;         ///< Character string used to write indentation spaces
};

template<class Ch, class Tr>
std::basic_ostream<Ch, Tr> &push(std::basic_ostream<Ch, Tr> &os) {
  os.put(PUSH);
//  os.put('\n');
  os.flush();
  return os;
}

template<class Ch, class Tr>
std::basic_ostream<Ch, Tr> &pop(std::basic_ostream<Ch, Tr> &os) {
  os.put(POP);
//  os.put('\n');
  os.flush();
  return os;
}


template<class Ch, class Tr>
class indent_streambuf_throwsafe
{
  explicit indent_streambuf_throwsafe(basic_indent_streambuf<Ch, Tr> &sb)
    : m_indentStreambuf(sb),
      m_indentLevel(sb.indent_level())
  {}

  ~indent_streambuf_throwsafe() {
    while (m_indentStreambuf.indent_level() > m_indentLevel)
      m_indentStreambuf.pop();
  }

private:
  basic_indent_streambuf<Ch, Tr> &      m_indentStreambuf;
  size_t                                m_indentLevel;
};


struct IndentFlags {
  int           m_flags;
};

inline IndentFlags indent_flags(int flags) {
  IndentFlags f;
  f.m_flags = flags;
  return f;
}

template<class Ch, class Tr>
std::basic_ostream<Ch, Tr> &
operator<<(std::basic_ostream<Ch, Tr> &os, IndentFlags indent_flags) {
  basic_indent_streambuf<Ch, Tr> *osb = dynamic_cast<basic_indent_streambuf<Ch, Tr> *>(os.rdbuf());
  if (osb)
    osb->set_flags(indent_flags.m_flags);
      
  return os;
};

typedef stk::basic_indent_streambuf<char> indent_streambuf;

} // namespace stk

#endif //  STK_UTIL_UTIL_INDENTSTREAMBUF_HPP
