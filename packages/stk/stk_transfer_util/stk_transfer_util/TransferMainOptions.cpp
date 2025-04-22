/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "TransferMainOptions.hpp"
#include <stk_util/command_line/CommandLineParserUtils.hpp>

namespace stk {
namespace transfer_util {

TransferMainOptions::TransferMainOptions(MPI_Comm comm, int argc, const char** argv)
: m_optionNames(),
  m_cmdLineParser(comm),
  m_parseState(stk::CommandLineParser::ParseError)
{
  stk::CommandLineOption fromMeshFilename{m_optionNames.fromMesh, "f",
          "Mesh to read from, i.e., source mesh."};
  stk::CommandLineOption toMeshFilename{m_optionNames.toMesh, "t",
          "Mesh to write to, i.e., destination or output mesh."};

  m_cmdLineParser.add_optional<std::string>(fromMeshFilename);
  m_cmdLineParser.add_optional<std::string>(toMeshFilename);

  m_parseState = m_cmdLineParser.parse(argc, argv);
}

std::string TransferMainOptions::get_fromMesh_filename() const
{
  return m_cmdLineParser.is_option_parsed(m_optionNames.fromMesh) ?
    m_cmdLineParser.get_option_value<std::string>(m_optionNames.fromMesh) : "";
}

std::string TransferMainOptions::get_toMesh_filename() const
{
  return m_cmdLineParser.is_option_parsed(m_optionNames.toMesh) ?
    m_cmdLineParser.get_option_value<std::string>(m_optionNames.toMesh) : "";
}

stk::CommandLineParser::ParseState TransferMainOptions::get_parse_state() const
{
  return m_parseState;
}

void TransferMainOptions::print_usage_help(std::ostream& os) const
{
  os << m_cmdLineParser.get_usage() << std::endl;
}

} // namespace transfer_util
} // namespace stk

