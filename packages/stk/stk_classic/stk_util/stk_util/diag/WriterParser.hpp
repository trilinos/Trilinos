#ifndef STK_UTIL_Diag_WriterParser_h
#define STK_UTIL_Diag_WriterParser_h

#include <stk_util/diag/Option.hpp>

namespace stk {
namespace diag {

///
/// @addtogroup DiagWriterDetail
/// @{
///

/**
 * @brief Class <b>WriterParser</b> implements a parser a Writer PrintMask string.
 *
 */
class WriterParser : public OptionMaskParser
{
public:
  /**
   * @brief Typedef <b>Mask</b> bring the OptionMaskParser Mask definition into this
   * namespace.
   *
   */
  typedef OptionMaskParser::Mask  Mask;

public:
  /**
   * @brief Creates a new <b>WriterParser</b> instance containing the lowerest
   * level PrintMask names.
   *
   */
  WriterParser();

  /**
   * @brief Member function <b>parse</b> returns the mask which results from parsing the
   * <b>mask_string</b>.
   *
   * @param mask_string    a <b>std::string</b> const reference to the string to be
   *        parsed.
   *
   * @return      a <b>Mask</b> value of the result from parsing the mask
   *        string.
   */
  Mask parse(const char *mask_string) const;

  /**
   * @brief Member function <b>parseArg</b> parses the argument and its argument
   * values.
   *
   * @param name    a <b>std::string</b> const reference to the argument
   *        name.
   *
   * @param arg      a <b>std::string</b> const reference to the argument
   *        values.
   */
  virtual void parseArg(const std::string &name, const std::string &arg) const;
};

///
/// @}
///

} // namespace diag
} // namespace stk

namespace sierra {
namespace Diag {

typedef stk::diag::WriterParser WriterParser;

} // namespace Diag
} // namespace sierra

#endif // STK_UTIL_Diag_WriterParser_h
