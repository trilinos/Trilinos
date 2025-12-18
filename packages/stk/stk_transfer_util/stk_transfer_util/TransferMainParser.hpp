/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_TRANSFER_UTIL_TRANSFERMAINPARSER_HPP
#define STK_TRANSFER_UTIL_TRANSFERMAINPARSER_HPP

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/command_line/CommandLineParserParallel.hpp>
#include <string>
#include <ostream>
#include <stk_transfer_util/TransferMainSettings.hpp>

namespace stk {
namespace transfer_util {

struct OptionNames
{
  const std::string sendMesh = "send-mesh";
  const std::string recvMesh = "recv-mesh";
  const std::string fieldList = "field-list";
  const std::string ExtrapolateOption = "extrapolate-option";
};

enum class TransferParserStatus {
  SUCCESS                    = 0,
  PARSE_ERROR                = 1,
  PARSE_ONLY                 = 2
};

class TransferMainParser
{
public:
  TransferMainParser(MPI_Comm comm);

  virtual ~TransferMainParser() = default;

  void add_options_to_parser(); 
  void parse_command_line_options(int argc, const char** argv);
  TransferParserStatus get_parser_status() const { return m_transferParserStatus; }
  TransferMainSettings generate_transfer_settings();

protected:

  void set_num_input_processors(TransferMainSettings& settings);
  void set_num_output_processors(TransferMainSettings& settings);
  void set_sendMesh_filename(TransferMainSettings& settings);
  void set_recvMesh_filename(TransferMainSettings& settings);
  void set_transfer_fields(TransferMainSettings& settings);
  void set_extrapolate_option(TransferMainSettings& settings);

  const MPI_Comm m_comm;
  OptionNames m_optionNames;
  stk::CommandLineParserParallel m_cmdLineParser;
  std::string m_quickExample;
  std::string m_longExamples;
  std::string m_quickError;
  TransferParserStatus m_transferParserStatus;

};

} // namespace transfer_util
} // namespace stk

#endif // STK_TRANSFER_UTIL_TRANSFERMAINPARSER_HPP
