/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "TransferMainParser.hpp"
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <stk_transfer_util/LogMessage.hpp>
#include <stk_util/Version.hpp>

namespace stk {
namespace transfer_util {

TransferMainParser::TransferMainParser(MPI_Comm comm)
  : m_comm(comm),
    m_cmdLineParser(comm),
    m_quickExample("quick example not provided yet"),
    m_longExamples("long examples not provided yet")
{
  add_options_to_parser();
}

void TransferMainParser::add_options_to_parser() 
{

  stk::CommandLineOption fromMeshFilename{m_optionNames.fromMesh, "f",
          "Mesh to read from, i.e., source mesh. (required)"};
  stk::CommandLineOption toMeshFilename{m_optionNames.toMesh, "t",
          "Mesh to write to, i.e., destination or output mesh. (required)"};
  stk::CommandLineOption fieldList{m_optionNames.fieldList, "l",
          "Comma-separated list of fields to transfer. Syntax example: field_1, field_2. If this option is not specified, all fields will be transferred. "};
  stk::CommandLineOption ExtrapolateOption{m_optionNames.ExtrapolateOption, "e",
          "Extrapolation option to be used. If this option is not specified, it will default to IGNORE."};
  m_cmdLineParser.add_required<std::string>(fromMeshFilename);
  m_cmdLineParser.add_required<std::string>(toMeshFilename);
  m_cmdLineParser.add_optional<std::string>(fieldList);
  m_cmdLineParser.add_optional<std::string>(ExtrapolateOption);

  m_cmdLineParser.disallow_unrecognized();

}

void TransferMainParser::parse_command_line_options(int argc, const char** argv)
{
  TransferMainSettings settings;
  stk::CommandLineParser::ParseState parseState = m_cmdLineParser.parse(argc, argv);

  if (parseState == stk::CommandLineParser::ParseComplete) {
    m_transferParserStatus = TransferParserStatus::SUCCESS;
  }
  else {
    switch(parseState) {
      case stk::CommandLineParser::ParseError:
        m_transferParserStatus = TransferParserStatus::PARSE_ERROR;
        stk::outputP0() << m_quickExample << std::endl;
        stk::outputP0() << "Use stk_transfer --help' for more information." << std::endl;
        break;
      case stk::CommandLineParser::ParseHelpOnly:
        m_transferParserStatus = TransferParserStatus::PARSE_ONLY;
        stk::outputP0() << m_quickExample << std::endl;
        stk::outputP0() << m_cmdLineParser.get_usage() << std::endl;
        stk::outputP0() << m_longExamples << std::endl;
        break;
      case stk::CommandLineParser::ParseVersionOnly:
        m_transferParserStatus = TransferParserStatus::PARSE_ONLY;
        stk::outputP0() << "STK Version: " << stk::version_string() << std::endl;
        break;
      default: break;
    }
  }
}

TransferMainSettings TransferMainParser::generate_transfer_settings()
{
  TransferMainSettings settings;

  if (m_transferParserStatus == TransferParserStatus::SUCCESS) {
    set_num_input_processors(settings);
    set_num_output_processors(settings);
    set_fromMesh_filename(settings);
    set_toMesh_filename(settings);
    set_transfer_fields(settings);
    set_extrapolate_option(settings);
  }

  return settings;
}

void TransferMainParser::set_num_input_processors(TransferMainSettings& settings)
{
  settings.set_num_input_processors(1);
}

void TransferMainParser::set_num_output_processors(TransferMainSettings& settings)
{
  settings.set_num_output_processors(1);
}

void TransferMainParser::set_fromMesh_filename(TransferMainSettings& settings)
{
  const std::string fromMesh = m_cmdLineParser.get_option_value<std::string>(m_optionNames.fromMesh);
  settings.set_fromMesh_filename(fromMesh);
}

void TransferMainParser::set_toMesh_filename(TransferMainSettings& settings)
{
  const std::string toMesh = m_cmdLineParser.get_option_value<std::string>(m_optionNames.toMesh);
  settings.set_toMesh_filename(toMesh);
}

void TransferMainParser::set_transfer_fields(TransferMainSettings& settings)
{
  if (m_cmdLineParser.is_option_parsed(m_optionNames.fieldList))
  {
    const std::string fieldList = m_cmdLineParser.get_option_value<std::string>(m_optionNames.fieldList);
    std::vector<std::string> fieldNames = stk::split_csv_string(fieldList);
    for (const std::string & fieldName : fieldNames) {
      settings.set_transfer_field(fieldName);  
    }
  }
}

void TransferMainParser::set_extrapolate_option(TransferMainSettings& settings)
{
  const std::string policy = m_cmdLineParser.is_option_parsed(m_optionNames.ExtrapolateOption) ?
    m_cmdLineParser.get_option_value<std::string>(m_optionNames.ExtrapolateOption) : "IGNORE";

  bool successfulPolicy = settings.set_extrapolate_option(policy);  
  
  if (!successfulPolicy) {
    m_transferParserStatus = TransferParserStatus::PARSE_ERROR;
  }  
}

} // namespace transfer_util
} // namespace stk

