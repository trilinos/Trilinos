#ifndef STK_UTIL_DIAG_WriterRegsitry_h
#define STK_UTIL_DIAG_WriterRegsitry_h

#include <map>
#include <string>
#include <vector>

#include <stk_util/util/string_case_compare.hpp>
#include <stk_util/diag/Writer.hpp>

#include <stk_util/util/FeatureTest.hpp>
#include <stk_util/diag/Writer_fwd.hpp>
#include <stk_util/diag/WriterParser.hpp>

namespace sierra {
namespace Diag {

///
/// @addtogroup DiagWriterDetail
/// @{
///

/**
 * @brief Typedef <b>WriterRegistry</b> is a mapping from name to diagnostic
 * writer.
 *
 */
class WriterRegistry : public std::map<std::string, std::pair<stk::diag::Writer *, OptionMaskParser *>, stk::LessCase>
{
public:  
  WriterRegistry();
  
  ~WriterRegistry();
};

class WriterThrowSafe 
{
public:
  WriterThrowSafe();

  ~WriterThrowSafe();

private:
  std::vector<stk::diag::WriterThrowSafe *>     m_writerVector;
};
  
/**
 * @brief Function <b>getWriterRegistry</b> returns a reference to the diagnostic
 * writer registry.
 *
 * @return		a <b>WriterRegistry</b> reference to the diagnostic writer
 *			registry.
 */
WriterRegistry &getWriterRegistry();

/**
 * @brief Function <b>registerWriter</b> registers a diagnostic writer with the
 * diagnostic writer registry.
 *
 * @param name		a <b>std::string</b> const reference to the name to use for the
 *			diagnostic writer.
 *
 * @param diag_writer	a <b>Writer</b> reference to the diagnostic writer.
 *
 */
void registerWriter(const std::string &name, Writer &diag_writer, OptionMaskParser &option_parser);

/**
 * @brief Member function <b>unregisterWriter</b> unregisters a diagnostic writer
 * from the diagnostic writer registry.
 *
 * @param name		a <b>std::string</b> const reference to the name to use for the
 *			diagnostic writer.
 *
 * @param diag_writer	a <b>Writer</b> reference to the diagnostic writer.
 *
 */
void unregisterWriter(const std::string &name, Writer &diag_writer);

///
/// @}
///

} // namespace Diag
} // namespace sierra

#include <stk_util/diag/WriterExt.hpp>

#endif // STK_UTIL_DIAG_WriterRegsitry_h

