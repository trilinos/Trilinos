#include <stk_util/diag/WriterRegistry.hpp>
#include <stk_util/diag/Option.hpp>

namespace sierra {
namespace Diag {

WriterThrowSafe::WriterThrowSafe()
{
  for (Diag::WriterRegistry::iterator it = Diag::getWriterRegistry().begin(); it != Diag::getWriterRegistry().end(); ++it)
    m_writerVector.push_back(new stk::diag::WriterThrowSafe(*(*it).second.first));
}


WriterThrowSafe::~WriterThrowSafe()
{
  for (std::vector<stk::diag::WriterThrowSafe *>::iterator it = m_writerVector.begin(); it != m_writerVector.end(); ++it)
    delete (*it);
}


WriterRegistry::WriterRegistry()
{}


WriterRegistry::~WriterRegistry()
{}


WriterRegistry &
getWriterRegistry()
{
  static WriterRegistry s_writerRegistry;

  return s_writerRegistry;
}


void
registerWriter(
  const std::string &	name,
  Writer &		diag_writer,
  OptionMaskParser &    option_parser)
{
  getWriterRegistry().insert(std::make_pair(name, std::make_pair(&diag_writer, &option_parser)));
}


void
unregisterWriter(
  const std::string &	name,
  Writer &		writer)
{
  WriterRegistry::iterator it = getWriterRegistry().find(name);
  if (it != getWriterRegistry().end() && (*it).second.first == &writer)
    getWriterRegistry().erase(it);
}

} // namespace Diag
} // namespace sierra
