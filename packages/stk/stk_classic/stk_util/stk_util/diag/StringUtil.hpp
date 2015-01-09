#ifndef STK_UTIL_DIAG_StringUtil_h
#define STK_UTIL_DIAG_StringUtil_h

#include <algorithm>
#include <functional>
#include <string>
#include <cctype>
#include <limits>

#include <stk_util/util/FeatureTest.hpp>
#include <stk_util/diag/String.hpp>

namespace sierra {

///
/// @addtogroup StringUtilDetail
/// @{
///

/**
 * @brief Function <b>case_strcmp</b> compares two null terminated strings case
 * insenstively.  It returns zero if they are case insensitively equal, a negative value
 * if <b>c1</b> is case insensitively less than <b>c2</b>, or a positive value if
 * <b>c1</b> is case insensitively greater than <b>c2</b>
 *
 * @param c1		a <b>char</b> const pointer to the first string.
 *
 * @param c2		a <b>char</b> const pointer to the second string.
 *
 * @return		an <b>int</b> value of zero if they are equal, negative if
 *			<b>c1</b> is less than <b>c2</b>, or positive if <b>c1</b> is
 *			greater than <b>c2</b>, without regard to case.
 */
int case_strcmp(const char *c1, const char *c2);

int case_strncmp(const char *c1, const char *c2, size_t n);

/**
 * @brief Function <b>case_strcmp</b> compares two null terminated strings case
 * insenstively.  It returns zero if they are case insensitively equal, a negative value
 * if <b>s1</b> is case insensitively less than <b>s2</b>, or a positive value if
 * <b>s1</b> is case insensitively greater than <b>s2</b>
 *
 * @param s1		a <b>std::string</b> reference to the first string.
 *
 * @param s2		a <b>std::string</b> reference to the second string.
 *
 * @return		an <b>int</b> value of zero if they are equal, negative if
 *			<b>s1</b> is less than <b>s2</b>, or positive if <b>s1</b> is
 *			greater than <b>s2</b>, without regard to case.
 */
inline int case_strcmp(const std::string &s1, const std::string &s2) {
  return case_strcmp(s1.c_str(), s2.c_str());
}

/**
 * @brief Function <b>case_strcmp</b> compares two null terminated strings case
 * insenstively.  It returns zero if they are case insensitively equal, a negative value
 * if <b>s1</b> is case insensitively less than <b>s2</b>, or a positive value if
 * <b>s1</b> is case insensitively greater than <b>s2</b>
 *
 * @param s1		a <b>String</b> reference to the first string.
 *
 * @param s2		a <b>String</b> reference to the second string.
 *
 * @return		an <b>int</b> value of zero if they are equal, negative if
 *			<b>s1</b> is less than <b>s2</b>, or positive if <b>s1</b> is
 *			greater than <b>s2</b>, without regard to case.
 */
inline int case_strcmp(const String &s1, const String &s2) {
  return case_strcmp(s1.c_str(), s2.c_str());
}

const char *case_strstr(const char *t, const char *find);

inline const char *case_strstr(const std::string &s1, const std::string &s2) {
  return case_strstr(s1.c_str(), s2.c_str());
}

inline const char *case_strstr(const String &s1, const String &s2) {
  return case_strstr(s1.c_str(), s2.c_str());
}


/*
 * @brief Member function <b>undef_if_empty</b> returns the string value
 * "<undefined>" if the string <b>s</b> is empty, otherwise it returns the string.
 *
 * @param s		a <b>String</b> reference to the string to check if empty.
 *
 * @return		a <b>String</b> value of "<undefined>" if the string
 *			<b>s</b> is empty, otherwise the string.
 */
inline String undef_if_empty(const String &s) {
  return s.empty() ? "<undefined>" : s;
}

/**
 * @brief Function <b>to_upper</b> returns the value converted to upper case.
 *
 * @param c		an <b>int</b> value to be converted.
 *
 * @return		an <b>int</b> value converted to upper case.
 */
inline int to_upper(int c) {
  return std::toupper(c);
}

/**
 * @brief Function <b>to_lower</b> returns the value converted to lower case.
 *
 * @param c		an <b>int</b> value to be converted
 *
 * @return		an <b>int</b> value converted to lower case.
 */
inline int to_lower(int c) {
  return std::tolower(c);
}

/**
 * @brief Function <b>to_identifier_upper</b> returns the value converted to underscore if
 * it is a space or to upper case.
 *
 * @param c		an <b>int</b> value to be converted
 *
 * @return		an <b>int</b> value converted to an underscore or upper case.
 */
inline int to_identifier_upper(int c) {
  return std::isspace(c) ? '_' : std::toupper(c);
}

/**
 * @brief Function <b>to_identifier_lower</b> returns the value converted to underscore if
 * it is a space or to lower case.
 *
 * @param c		an <b>int</b> value to be converted
 *
 * @return		an <b>int</b> value converted to an underscore or lower case.
 */
inline int to_identifier_lower(int c) {
  return std::isspace(c) ? '_' : std::tolower(c);
}

/**
 * @brief Function <b>make_upper</b> converts strgin <b>name</b> to upper case.  The
 * conversion happens in place.
 *
 * @param name		a <b>T</b> reference to convert to a upper case.
 *
 * @return		an <b>T</b> reference to the converted upper case string.
 */
template <class T>
inline T &make_upper(T &name) {
  std::transform(name.begin(), name.end(), name.begin(), to_upper);

  return name;
}

/**
 * @brief Function <b>make_lower</b> converts string <b>name</b> to lower case.  The
 * convertion happens in place.
 *
 * @param name		a <b>T</b> reference to convert to an lower.
 *
 * @return		an <b>T</b> reference to the converted lower case string.
 */
template <class T>
inline T &make_lower(T &name) {
  std::transform(name.begin(), name.end(), name.begin(), to_lower);

  return name;
}

/**
 * @brief Function <b>make_identifier</b> converts string <b>name</b> to an identifier case.  The
 * convertion happens in place.
 *
 * @param name		a <b>T</b> reference to convert to an identifier.
 *
 * @return		an <b>T</b> reference to the converted identifier string.
 */
template <class T>
inline T &make_identifier(T &name) {
  std::transform(name.begin(), name.end(), name.begin(), to_identifier_upper);

  return name;
}

/**
 * @brief Function <b>make_lower_identifier</b> converts string <b>name</b> to a lower
 * case identifier case.  The convertion happens in place.
 *
 * @param name		a <b>T</b> reference to convert to a lower case identifier.
 *
 * @return		an <b>T</b> reference to the converted lower case identifier
 *			string.
 */
template <class T>
inline T &make_lower_identifier(T &name) {
  std::transform(name.begin(), name.end(), name.begin(), to_identifier_lower);

  return name;
}

/**
 * @brief Function <b>make_identifier</b> converts string <b>name</b> to an identifier
 * case.  The convertion happens in place.
 *
 * @param name		a <b>char</b> pointer to convert to an identifier.
 *
 * @return		an <b>char</b> pointer to the converted identifier string.
 */
inline char *make_identifier(char *name) {
  for (char *c = name; *c != 0; ++c )
    *c = to_identifier_upper(*c);

  return name;
}

/**
 * @brief Function <b>make_lower_identifier</b> converts string <b>name</b> to a lower
 * case identifier case.  The convertion happens in place.
 *
 * @param name		a <b>char</b> pointer to convert to a lower case identifier.
 *
 * @return		an <b>char</b> pointer to the converted lower case identifier
 *			string.
 */
inline char *make_lower_identifier(char *name) {
  for (char *c = name; *c != 0; ++c )
    *c = to_identifier_lower(*c);

  return name;
}

/**
 * @brief Function <b>trim</b> trims preceeding and trailing spaces from <b>name</b>.  The
 * convertion happens in place.
 *
 * @param name		a <b>T</b> reference to trim.
 *
 * @return		an <b>T</b> reference to the trimmed string.
 */
template <class T>
inline T &
trim(
  T &			name)
{
  typename T::iterator it0 = name.begin();
  while (it0 != name.end() && std::isspace(*it0))
    ++it0;

  typename T::iterator it1 = name.end();
  while (it1 != it0 && std::isspace(*(it1 - 1)))
    --it1;
  return name = T(it0, it1);
}

/**
 * @brief Function <b>title</b> returns a first letter of each word capitalized of the
 * string.
 *
 * @param s		a <b>std::string</b> const reference to be word capitalized.
 *
 * @return		a <b>std::string</b> value of the word capitalized string.
 */
std::string title(const std::string &s);

inline std::string title(const String &s) {
  return title(std::string(s.c_str()));
}


/**
 * @brief Function <b>lower</b> returns a lower case version of the string.
 *
 * @param s		a <b>String</b> const reference to be lower cased.
 *
 * @return		a <b>String</b> value of the lower cased string.
 */
inline String lower(const String &s)
{
  String t = s;
  return make_lower(t);
}

/**
 * @brief Function <b>to_string</b> returns a string representation of <b>r</b>.
 *
 * @param r		a <b>double</b> value to be stringized.
 *
 * @param precision	an <b>int</b> value of the precision.
 *
 * @return		a <b>std::string</b> value of the string representation of
 *			<b>r</b>.
 */
std::string to_string(const double &r, int precision = 4);

/**
 * @brief Function <b>to_string</b> returns a string representation of <b>r</b>.
 *
 * @param r		a <b>float</b> value to be stringized.
 *
 * @param precision	an <b>int</b> value of the precision.
 *
 * @return		a <b>std::string</b> value of the string representation of
 *			<b>r</b>.
 */
std::string to_string(const float &r, int precision = 4);

/**
 * Template function <b>to_string</b> returns a string representation of
 * <b>r</b>.
 *
 * @param t		a <b>T</b> value to be stringized.
 *
 * @return		a <b>std::string</b> value of the string representation of
 *			<b>r</b>.
 */
template <typename T>
std::string to_string(const T &t);

/**
 * @brief Function <b>demangle</b> returns the demangled C++ symbol from the mangled
 * C++ symbol.  The mangled named is obtained from the <b>type_info</b>
 * <b>name()</b> function.  From some compilers, the name is already demangled.
 *
 * @param symbol	a <b>char</b> const pointer to the symbol.
 *
 * @return		a <b>std::string</b> value of the demangled name.
 */
#ifdef SIERRA_USE_PLATFORM_DEMANGLER
// Implement this is Slib_EnvPlatform.C
std::string demangle(const char *symbol);
#else
const char *demangle(const char *symbol);
#endif

/**
 * @brief Function <b>format_time</b> encodes the time using the format specified.
 * The format is described in <b>stdftime</b>.
 *
 * @param t		a <b>time_t</b> value of the time to format.
 *
 * @param format	a <b>char</b> const pointer to the format.
 *
 * @return		a <b>String</b> value of the encoded time.
 */
std::string format_time(double t, const char *format = "%b %e %Y %H:%M:%S");

/**
 * @brief Function <b>word_wrap</b> reformats a string into multiple lines, none longer
 * that <b>line_length</b>, the first line prefixed with <b>prefix_first_line</b> and the
 * remaining lines prefixed with <b>prefix</b>.
 *
 * @param s			a <b>std::string</b> const reference to the string
 *				to word wrap.
 *
 * @param line_length		an <b>unsigned int</b> value of the line length.
 *
 * @param prefix		a <b>std::string</b> const reference to a line
 *				prefix string.
 *
 * @param prefix_first_line	a <b>std::string</b> const reference to a line
 *				prefix string for the first line.
 *
 * @return			a <b>std::string</b> value of the word wrapped string.
 */
std::string word_wrap(const std::string &s, unsigned int line_length,
		      const std::string &prefix, const std::string &prefix_first_line);

/**
 * @brief Function <b>word_wrap</b> reformats a string into multiple lines, none longer
 * that <b>line_length</b> and each line prefixed with <b>prefix</b>.
 *
 * @param s			a <b>std::string</b> const reference to the string
 *				to word wrap.
 *
 * @param line_length		an <b>unsigned int</b> value of the line length.
 *
 * @param prefix		a <b>std::string</b> const reference to a line
 *				prefix string.
 *
 * @return			a <b>std::string</b> value of the word wrapped string.
 */
inline std::string word_wrap(const std::string &s, unsigned int line_length = 72, const std::string &prefix = "") {
  return word_wrap(s, line_length, prefix, prefix);
}


/**
 * @brief Class <b>object_phrase</b> makes a pretty string for those annoying plural or
 * singular noun/verb phrases.
 *
 * The output is <i>plural</i> no <i>noun</i>s, when <b>n</b> is zero, <br>
 * <i>singular</i> 1 <i>noun</i>, when <b>n</b> is 1, or<br>
 * <i>plural</i> <i>n</i> <i>noun</i>s, when <b>n</b> is greater than 1.
 *
 */
class object_phrase
{
public:
  /**
   * Creates a new <b>object_phrase</b> instance.
   *
   * @param n		an <b>int</b> value of the quantity of objects
   *
   * @param noun	a <b>char</b> const pointer to the name of the object.
   *
   * @param singlar	a <b>char</b> singular form of the verb acting on the object.
   *
   * @param plural	a <b>char</b> plural form of the verb acting on the object.
   *
   */
  object_phrase(int n, const char *noun, const char *singlar = "is", const char *plural = "are")
    : m_n(n),
      m_noun(noun),
      m_singular(singlar),
      m_plural(plural)
  {}

  /**
   * @brief Member function <b>print</b> writes the object phrase to the output stream.
   *
   * @param os		a <b>std::ostream</b> reference to the output stream to write to.
   *
   * @return		a <b>std::ostream</b> reference to the output stream.
   */
  std::ostream &print(std::ostream &os) const;

  /**
   * @brief Member function <b>operator std::string</b> returns a string of the object
   * phrase.
   *
   * @return		a <b>std::string</b> value of the object string.
   */
  operator std::string () const;

private:
  int			m_n;			///< Object count
  const char *		m_noun;			///< Object name
  const char *		m_singular;		///< Singular verb
  const char *		m_plural;		///< Plural verb
};

/**
 * @brief Function <b>operator<<</b> writes an object phrase to the output stream.
 *
 * @param os		a <b>std::ostream</b> reference to the output stream to write to.
 *
 * @param phrase	an <b>object_phrase</b> const reference to the object phrase to be
 *			written.
 *
 * @return		a <b>std::ostream</b> reference to the output stream.
 */
inline std::ostream &operator<<(std::ostream &os, const object_phrase &phrase) {
  return phrase.print(os);
}


/**
 * @brief Function <b>getline</b> returns a string from the input stream which has been
 * terminated by the newline character.
 *
 * @param is		a <b>std::istream</b> reference to the input stream.
 *
 * @param s		a <b>sierra::String</b> reference to received the line.
 *
 * @param eol		a <b>char</b> value to use of the newline character.
 *
 * @return		a <b>std::istream</b> reference to the input stream.
 */
std::istream &getline(std::istream &is, sierra::String &s, char eol = '\n');

/**
 * @brief Class <b>less_nocase</b> implements a case insensitive compare functor.
 *
 */
template <class T>
struct less_nocase : public std::binary_function<T, T, bool>
{
  /**
   * @brief Member function <b>operator()</b> returns true if the <b>lhs</b> is less than
   * <b>rhs</b>.
   *
   * @param lhs		a <b>T</b> const reference to the left hand side value.
   *
   * @param rhs		a <b>T</b> cosnt reference to the right hand side value.
   *
   * @return		a <b>bool</b> value of true if the <b>lhs</b> is less than
   *			<b>rhs</b>
   */
  bool operator()(const T &lhs, const T &rhs) const {
    typename T::const_iterator lhs_it = lhs.begin();
    typename T::const_iterator rhs_it = rhs.begin();
    typename T::const_iterator lhs_it_end = lhs.end();
    typename T::const_iterator rhs_it_end = rhs.end();

    for (; lhs_it != lhs_it_end && rhs_it != rhs_it_end; ++lhs_it, ++rhs_it) {
      int i =  std::tolower(*lhs_it) - std::tolower(*rhs_it);
      if (i != 0)
	return i < 0;
    }

    if (lhs_it == lhs_it_end)
      return rhs_it != rhs_it_end;
    else
      return false;
  }
};


/**
 * @brief Class specialization <b>less_nocase</b> for <b>String</b>.
 *
 */
template <>
struct less_nocase<String> : public std::binary_function<String, String, bool>
{
  /**
   * @brief Member function <b>operator()</b> returns true if the <b>lhs</b> is less than
   * <b>rhs</b>.
   *
   * @param lhs		a <b>String</b> const reference to the left hand side value.
   *
   * @param rhs		a <b>String</b> cosnt reference to the right hand side value.
   *
   * @return		a <b>bool</b> value of true if the <b>lhs</b> is less than
   *			<b>rhs</b>
   */
  bool operator()(const String &lhs, const String &rhs) const {
    const char * lhs_it = lhs.c_str();
    const char * rhs_it = rhs.c_str();
    const char * lhs_it_end = lhs_it + lhs.length();
    const char * rhs_it_end = rhs_it + rhs.length();

    for (; lhs_it != lhs_it_end && rhs_it != rhs_it_end; ++lhs_it, ++rhs_it) {
      int i =  std::tolower(*lhs_it) - std::tolower(*rhs_it);
      if (i != 0)
	return i < 0;
    }

    if (lhs_it == lhs_it_end)
      return rhs_it != rhs_it_end;
    else
      return false;
  }
};


/**
 * @brief Class specialization <b>less_nocase</b> for <b>std::string</b>.
 *
 */
template <>
struct less_nocase<std::string> : public std::binary_function<std::string, std::string, bool>
{
  /**
   * @brief Member function <b>operator()</b> returns true if the <b>lhs</b> is less than
   * <b>rhs</b>.
   *
   * @param lhs		a <b>std::string</b> const reference to the left hand side value.
   *
   * @param rhs		a <b>std::string</b> cosnt reference to the right hand side value.
   *
   * @return		a <b>bool</b> value of true if the <b>lhs</b> is less than
   *			<b>rhs</b>
   */
  bool operator()(const std::string &lhs, const std::string &rhs) const {
    const char * lhs_it = lhs.c_str();
    const char * rhs_it = rhs.c_str();
    const char * lhs_it_end = lhs_it + lhs.length();
    const char * rhs_it_end = rhs_it + rhs.length();

    for (; lhs_it != lhs_it_end && rhs_it != rhs_it_end; ++lhs_it, ++rhs_it) {
      int i =  std::tolower(*lhs_it) - std::tolower(*rhs_it);
      if (i != 0)
	return i < 0;
    }

    if (lhs_it == lhs_it_end)
      return rhs_it != rhs_it_end;
    else
      return false;
  }
};


/**
 * @brief Class specialization <b>less_nocase</b> for <b>char const pointer</b>.
 *
 */
template<>
struct less_nocase<const char *> : public std::binary_function<const char *, const char *, bool>
{
  /**
   * @brief Member function <b>operator()</b> returns true if the <b>lhs</b> is less than
   * <b>rhs</b>.
   *
   * @param lhs		a <b>char</b> const pointer to the left hand side value.
   *
   * @param rhs		a <b>char</b> const pointer to the right hand side value.
   *
   * @return		a <b>bool</b> value of true if the <b>lhs</b> is less than
   *			<b>rhs</b>
   */
  bool operator()(const char *lhs , const char *rhs) const {
    bool result ;
    if (NULL == lhs )
      result = NULL != rhs;
    else {
      for (; *lhs && *rhs && std::tolower(*lhs) == std::tolower(*rhs); ++lhs, ++rhs)
	;
      result = std::tolower(*lhs) < std::tolower(*rhs);
    }
    return result ;
  }
};


/**
 * @brief Class <b>equal_nocase</b> implements a case insensitive compare functor.
 *
 */
template <class T>
struct equal_nocase : public std::binary_function<T, T, bool>
{
  /**
   * @brief Member function <b>operator()</b> returns true if the <b>lhs</b> is equal to
   * <b>rhs</b>.
   *
   * @param lhs		a <b>T</b> const reference to the left hand side value.
   *
   * @param rhs		a <b>T</b> cosnt reference to the right hand side value.
   *
   * @return		a <b>bool</b> value of true if the <b>lhs</b> is less than
   *			<b>rhs</b>
   */
  bool operator()(const T &lhs, const T &rhs) const {
    typename T::const_iterator lhs_it = lhs.begin();
    typename T::const_iterator rhs_it = rhs.begin();
    typename T::const_iterator lhs_it_end = lhs.end();
    typename T::const_iterator rhs_it_end = rhs.end();

    for (; lhs_it != lhs_it_end && rhs_it != rhs_it_end; ++lhs_it, ++rhs_it) {
      if (std::tolower(*lhs_it) != std::tolower(*rhs_it))
	return false;
    }

    return lhs_it == lhs_it_end && rhs_it == rhs_it_end;
  }
};

/**
 * @brief Class specialization <b>equal_nocase</b> for <b>String</b>.
 *
 */
template <>
struct equal_nocase<String> : public std::binary_function<String, String, bool>
{
  /**
   * @brief Member function <b>operator()</b> returns true if the <b>lhs</b> is equal to
   * <b>rhs</b>.
   *
   * @param lhs		a <b>String</b> const reference to the left hand side value.
   *
   * @param rhs		a <b>String</b> cosnt reference to the right hand side value.
   *
   * @return		a <b>bool</b> value of true if the <b>lhs</b> is less than
   *			<b>rhs</b>
   */
  bool operator()(const String &lhs, const String &rhs) const {
    const char * lhs_it = lhs.c_str();
    const char * rhs_it = rhs.c_str();
    const char * lhs_it_end = lhs_it + lhs.length();
    const char * rhs_it_end = rhs_it + rhs.length();

    for (; lhs_it != lhs_it_end && rhs_it != rhs_it_end; ++lhs_it, ++rhs_it) {
      if (std::tolower(*lhs_it) != std::tolower(*rhs_it))
	return false;
    }

    return lhs_it == lhs_it_end && rhs_it == rhs_it_end;
  }
};

/**
 * @brief Class specialization <b>equal_nocase</b> for <b>std::string</b>.
 *
 */
template <>
struct equal_nocase<std::string> : public std::binary_function<std::string, std::string, bool>
{
  /**
   * @brief Member function <b>operator()</b> returns true if the <b>lhs</b> is equal to
   * <b>rhs</b>.
   *
   * @param lhs		a <b>std::string</b> const reference to the left hand side value.
   *
   * @param rhs		a <b>std::string</b> cosnt reference to the right hand side value.
   *
   * @return		a <b>bool</b> value of true if the <b>lhs</b> is less than
   *			<b>rhs</b>
   */
  bool operator()(const std::string &lhs, const std::string &rhs) const {
    const char * lhs_it = lhs.c_str();
    const char * rhs_it = rhs.c_str();
    const char * lhs_it_end = lhs_it + lhs.length();
    const char * rhs_it_end = rhs_it + rhs.length();

    for (; lhs_it != lhs_it_end && rhs_it != rhs_it_end; ++lhs_it, ++rhs_it) {
      if (std::tolower(*lhs_it) != std::tolower(*rhs_it))
	return false;
    }

    return lhs_it == lhs_it_end && rhs_it == rhs_it_end;
  }
};

/**
 * @brief Class specialization <b>equal_nocase</b> for <b>char const pointer</b>.
 *
 */
template <>
struct equal_nocase<const char *> : public std::binary_function<const char *, const char *, bool> {
  /**
   * @brief Member function <b>operator()</b> returns true if the <b>lhs</b> is equal to
   * <b>rhs</b>.
   *
   * @param lhs		a <b>char</b> const pointer to the left hand side value.
   *
   * @param rhs		a <b>char</b> cosnt pointer to the right hand side value.
   *
   * @return		a <b>bool</b> value of true if the <b>lhs</b> is less than
   *			<b>rhs</b>
   */
  bool operator()(const char *lhs , const char *rhs) const {
    bool result = lhs == rhs ;
    if ( ! result && NULL != lhs && NULL != rhs ) {
      for (; *lhs && *rhs && std::tolower(*lhs) == std::tolower(*rhs) ; ++lhs, ++rhs);
      result = 0 == *lhs && 0 == *rhs ;
    }
    return result ;
  }
};

/**
 * @brief Class <b>hash_nocase</b> is a traits class for hash functions.
 *
 */
template <class _Key> struct hash_nocase {};

/**
 * @brief Function <b>hash_string_nocase</b> hash a character string case insensitively.
 *
 * @param p		a <b>char</b> const pointer to a character string.
 *
 * @return		a <b>size_t</b> value of the hashed string.
 */
inline
size_t hash_string_nocase(
  const char *		p)
{
  size_t h = 0;
  const size_t sr = std::numeric_limits<unsigned char>::digits *  sizeof(size_t) - 8;
  const size_t mask = ((size_t) 0xF) << (sr + 4);
  while (*p) {
    h = (h << 4) + std::tolower(*p++);
    size_t g = h & mask;
    h ^= g | (g >> sr);
  }
  return h;
}

/**
 * @brief Class specialization <b>hash_nocase</b> for <b>std::string</b>.
 *
 */
template<>
struct hash_nocase<std::string>
{
  /**
   * @brief Member function <b>operator()</b> returns the hash value of the
   * specified string.
   *
   * @param __s		a <b>std::string</const reference to the string to hash.
   *
   * @return		a <b>size_t</b> value of the hash value of the specified
   *			string.
   */
  size_t operator()(const std::string & __s) const { return hash_string_nocase(__s.c_str()); }
};

/**
 * @brief Class specialization <b>hash_nocase</b> for <b>String</b>.
 *
 */
template<>
struct hash_nocase<String>
{
  /**
   * @brief Member function <b>operator()</b> returns the hash value of the
   * specified string.
   *
   * @param __s		a <b>std::string</const reference to the string to hash.
   *
   * @return		a <b>size_t</b> value of the hash value of the specified
   *			string.
   */
  size_t operator()(const String & __s) const { return hash_string_nocase(__s.c_str()); }
};

// /**
//  * @brief Class specialization <b>hash_nocase</b> for <b>std::string</b>.
//  *
//  */
// template<>
// struct hash_nocase<const std::string>
// {
//   /**
//    * @brief Member function <b>operator()</b> returns the hash value of the
//    * specified string.
//    *
//    * @param __s		a <b>std::string</const reference to the string to hash.
//    *
//    * @return		a <b>size_t</b> value of the hash value of the specified
//    *			string.
//    */
//   size_t operator()(const std::string & __s) const { return hash_string_nocase(__s.c_str()); }
// };

// /**
//  * @brief Class specialization <b>hash_nocase</b> for <b>std::string</b>.
//  *
//  */
// template<>
// struct hash_nocase<const String>
// {
//   /**
//    * @brief Member function <b>operator()</b> returns the hash value of the
//    * specified string.
//    *
//    * @param __s		a <b>std::string</const reference to the string to hash.
//    *
//    * @return		a <b>size_t</b> value of the hash value of the specified
//    *			string.
//    */
//   size_t operator()(const String & __s) const { return hash_string_nocase(__s.c_str()); }
// };

///
/// @}
///

template <class T>
T convert_cast(const String &s);

} // namespace sierra

#endif // STK_UTIL_DIAG_StringUtil_h
