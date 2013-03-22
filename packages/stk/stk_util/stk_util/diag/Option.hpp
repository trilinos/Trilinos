#ifndef STK_UTIL_DIAG_Option_h
#define STK_UTIL_DIAG_Option_h

#include <map>

namespace stk {
namespace diag {



typedef unsigned long OptionMask;

/**
 * @brief Class <b>OptionDescription</b> is an interface class for describing a
 * command line option or option value.
 *
 */
struct OptionDescription
{
  /**
   * Destroys a <b>OptionDescription</b> instance.
   *
   */
  virtual ~OptionDescription()
  {}

  /**
   * @brief Member function <b>describe</b> prints a description of the option to
   * the stream.
   *
   * @param os      a <b>std::ostream</b> reference to print to
   *        description to.
   *
   * @return      a <b>std::ostream</b> reference to the output
   *        stream.
   */
  virtual std::ostream &describe(std::ostream &os) const = 0;
};


/**
 * @brief Class <b>Option</b> holds the command line name, environment variable name
 * and the current value of an option.  It implements the <b>OptionDescription</b>
 * interface so that a help description may be generated.
 *
 */
struct Option : public OptionDescription
{
  /**
   * Creates a new <b>Option</b> instance.
   *
   */
  Option()
    : m_name(""),
      m_envName(""),
      m_description(""),
      m_value(""),
      m_subOptionDescription(0)
  {}

  /**
   * Creates a new <b>Option</b> instance.
   *
   * @param name      a <b>std::string</b> const reference to the name
   *          of the option.  This is used for the command line
   *          option argument.
   *
   * @param env_name      a <b>std::string</b> const reference to the
   *          environment variable name.
   *
   * @param value      a <b>std::string</b> const reference to the
   *          default/initial value of the option.
   *
   * @param description      a <b>std::string</b> const reference to the
   *          description of the option.  This is printed when
   *          the -h option is parsed.
   *
   * @param sub_option_description  an <b>OptionDescription</b> const pointer to sub
   *          options which are available for the option.
   */
  Option(const std::string &name, const std::string &env_name, const std::string &value = std::string(),
         const std::string &description = "No description available",
         const OptionDescription *sub_option_description = 0)
    : m_name(name),
      m_envName(env_name),
      m_description(description),
      m_value(value),
      m_subOptionDescription(sub_option_description)
  {}

  /**
   * Creates a new <b>Option</b> instance.
   *
   * @param option    an <b>Option</b> const reference to the Option to
   *        copy.
   *
   */
  Option(const Option &option)
    : m_name(option.m_name),
      m_envName(option.m_envName),
      m_description(option.m_description),
      m_value(option.m_value),
      m_subOptionDescription(option.m_subOptionDescription)
  {}


  /**
   * @brief Member function <b>operator=</b> assigns an option from another option.
   *
   * @param option    an <b>Option</b> const reference to the rhs option.
   *
   * @return      an <b>Option</b> reference to the lhs option.
   */
  Option &operator=(const Option &option)
  {
    m_name = option.m_name;
    m_envName = option.m_envName;
    m_description = option.m_description;
    m_value = option.m_value;
    m_subOptionDescription = option.m_subOptionDescription;
    return *this;
  }

  /**
   * Destroys a <b>Option</b> instance.
   *
   */
  virtual ~Option()
  {}

  const std::string &getName() const {
    return m_name;
  }

  const std::string &getValue() const {
    return m_value;
  }

  /**
   * @brief Member function <b>describe</b> prints a description of the option to
   * the stream.
   *
   * @param os      a <b>std::ostream</b> reference to print to
   *        description to.
   *
   * @return      a <b>std::ostream</b> reference to the output
   *        stream.
   */
  virtual std::ostream &describe(std::ostream &os) const;

  operator std::string &() {
    return m_value;
  }

  std::string      m_name;      ///< Name/Command line option name
  std::string      m_envName;    ///< Environment variable name
  std::string      m_description;    ///< Brief '-h' description
  std::string      m_value;    ///< Value of option
  const OptionDescription *  m_subOptionDescription;  ///< Suboptions (used for '-h' parsing)
};


/**
 * @brief Class <b>OptionMaskName</b> ...
 *
 */
struct OptionMaskName
{
  /**
   * Creates a new <b>OptionMaskName</b> instance.
   *
   */
  OptionMaskName()
    : m_name(""),
      m_mask(0),
      m_description("")
  {}

  /**
   * Creates a new <b>OptionMaskName</b> instance.
   *
   * @param name    a <b>std::string</b> variable ...
   * @param mask    an <b>OptionMask</b> variable ...
   * @param description    a <b>std::string</b> variable ...
   */
  OptionMaskName(const std::string &name, const OptionMask &mask, const std::string &description = "No description available")
    : m_name(name),
      m_mask(mask),
      m_description(description)
  {}

  /**
   * Destroys a <b>OptionMaskName</b> instance.
   *
   */
  virtual ~OptionMaskName()
  {}

  /**
   * @brief Member function <b><</b> ...
   *
   * @param o      an <b>OptionMaskName</b> variable ...
   * @return      a <b>bool</b> ...
   */
  bool operator<(const OptionMaskName &o) const {
    return m_name < o.m_name;
  }

  /**
   * @brief Member function <b>describe</b> ...
   *
   * @param os      a <b>std::ostream</b> variable ...
   * @return      a <b>std::ostream</b> ...
   */
  virtual std::ostream &describe(std::ostream &os) const;

  std::string    m_name;
  OptionMask    m_mask;
  std::string    m_description;
};


/**
 * @brief Class <b>OptionMaskNameMap</b> ...
 *
 */
class OptionMaskNameMap: public std::map<std::string, OptionMaskName>
{
public:
  /**
   * @brief Member function <b>mask</b> ...
   *
   * @param name    a <b>std::string</b> variable ...
   * @param mask    an <b>OptionMask</b> variable ...
   * @param description  a <b>std::string</b> variable ...
   */
  void mask(const std::string &name, const OptionMask l_mask, const std::string &description) {
    iterator it = find(name);
    if (it == end())
      insert(std::make_pair(name, OptionMaskName(name, l_mask, description)));
    else {
      (*it).second.m_mask = l_mask;
      (*it).second.m_description = description;
    }
  }
};


/**
 * Class <b>OptionMaskParser</b> defines a mapping between strings and bit masks and
 * output streams.
 *
 * After populating a Parser object, parse() will parse the input string.  The
 * getMask() and virtual getOutputStream() functions will return the print parsed
 * print mask and the selected output stream.
 *
 */
class OptionMaskParser : public OptionDescription
{
public:
  typedef OptionMask Mask;    ///< Mask for this option

public:
  /**
   * Creates a new <b>OptionMaskParser</b> instance.
   *
   */
  OptionMaskParser()
    : m_optionMaskNameMap(),
      m_optionMask(0),
      m_status(true)
  {}

  /**
   * Destroys a <b>OptionMaskParser</b> instance.
   *
   */
  virtual ~OptionMaskParser()
  {}

  const OptionMaskNameMap &getOptionMaskNameMap() const {
    return m_optionMaskNameMap;
  }

  /**
   * Member function <b>parse</b> parses the string
   *
   * @param mask    a <b>char</b> const pointer of the string to parse.
   *
   * @return      a <b>Mask</b> value of the parsed bitmask.
   */
  virtual Mask parse(const char *mask) const;

  /**
   * Member function <b>parseArg</b> parses the argument and its argument values.
   *
   * @param name    a <b>std::string</b> const reference to the argument
   *        name.
   *
   * @param arg      a <b>std::string</b> const reference to the argument
   *        values.
   */
  virtual void parseArg(const std::string &name, const std::string &arg) const;

  /**
   * Member function <b>operator[]</b> returns the print mask with the specified
   * name.  If the name is not found, a new entry is added.
   *
   * @param name    a <b>std::string</b> const reference of the name of
   *        the mask.
   *
   * @return      a <b>Mask</b> reference of the print
   *        mask associated with the print mask name.
   */
  Mask &operator[](const std::string &name) {
    return m_optionMaskNameMap[name].m_mask;
  }

  /**
   * Member function <b>mask</b> adds a named mask to the parser.  The mask can also
   * be given a description which is displayed using the <b>describe()</b> function.
   *
   * @param name    a <b>std::string</b> const reference of the name of
   *        the mask.
   *
   * @param mask    a <b>Mask</b> value to associate with the name
   *
   * @param description    a <b>std::string</b> const reference which describes
   *        the mask.
   *
   */
  void mask(const std::string &name, const Mask l_mask, const std::string &description) {
    m_optionMaskNameMap.mask(name, l_mask, description);
  }

  /**
   * Member function <b>describe</b> ...
   *
   * @param os    a <b>std::ostream</b> variable ...
   * @return a <b>std::ostream</b> ...
   */
  std::ostream &describe(std::ostream &os) const;

protected:
  OptionMaskNameMap    m_optionMaskNameMap;  ///< Mask name vector
  mutable OptionMask    m_optionMask;    ///< Most recently parsed mask
  mutable bool      m_status;    ///< Result of most recent parse
};

} // namespace diag
} // namespace stk

namespace sierra {
typedef stk::diag::OptionMask OptionMask;
typedef stk::diag::OptionDescription OptionDescription;
typedef stk::diag::Option Option;
typedef stk::diag::OptionMaskName OptionMaskName;
typedef stk::diag::OptionMaskNameMap OptionMaskNameMap;
typedef stk::diag::OptionMaskParser OptionMaskParser;
} // namespace sierra

#endif // STK_UTIL_DIAG_Option_h
