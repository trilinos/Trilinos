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

#include <optional>
#include <ostream>
#include <stk_transfer_util/TransferMainInputFile.hpp>
#include <stk_transfer_util/TransferMainOptionNames.hpp>
#include <stk_transfer_util/TransferMainSettings.hpp>
#include <stk_util/command_line/CommandLineParserParallel.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <string>

namespace stk {
namespace transfer_util {

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

  // Parse the contents of a YAML input file (see TransferMainInputFile.hpp for the syntax).
  // Options specified on the command line take precedence over options in the input file.
  // Lets callers (e.g., unit tests) provide input-file text directly instead of via
  // --input-file. Must be called BEFORE parse_command_line_options (enforced; throws
  // otherwise) so that required options can be satisfied by the input file. If
  // --input-file is also given, the file replaces this text. Throws an exception derived
  // from std::exception on syntax errors.
  void set_input_file_contents(const std::string& yamlText);

#ifdef STK_BUILT_FOR_SIERRA
  static inline const std::string default_master_element_name = "SIERRA";
#else
  static inline const std::string default_master_element_name = "INTREPID2";
#endif

  bool input_file_specified() const;
  std::string get_input_file_name() const;

protected:
 void read_input_file_if_specified();
 void require_mesh_specified(const std::string& optionName, const std::optional<std::string>& fileValue) const;
 std::optional<std::string> get_string_option(
     const std::string& optionName, const std::optional<std::string>& fileValue) const;

 void set_num_input_processors(TransferMainSettings& settings);
 void set_num_output_processors(TransferMainSettings& settings);
 void set_sendMesh_filename(TransferMainSettings& settings);
 void set_recvMesh_filename(TransferMainSettings& settings);
 void set_outputMesh_filename(TransferMainSettings& settings);
 void set_transfer_fields(TransferMainSettings& settings);
 void set_transfer_parts(TransferMainSettings& settings);
 void set_extrapolate_option(TransferMainSettings& settings);
 void set_master_elements_name(TransferMainSettings& settings);
 void set_recv_type(TransferMainSettings& settings);
 void set_time_steps_spec(TransferMainSettings& settings);
 void set_transfer_type(TransferMainSettings& settings);
 void set_2d_to_3d_mapping_type(TransferMainSettings& settings);
 void set_2d_to_3d_axisymmetric_transfer_params(TransferMainSettings& settings);
 void set_3d_to_3d_axisymmetric_transfer_params(TransferMainSettings& settings);
 void set_coord_transf_z_expr(TransferMainSettings& settings);

 const MPI_Comm m_comm;
 OptionNames m_optionNames;
 stk::CommandLineParserParallel m_cmdLineParser;
 std::string m_quickExample;
 std::string m_longExamples;
 std::string m_quickError;
 TransferParserStatus m_transferParserStatus;
 TransferInputFileOptions m_inputFileOptions;
 bool m_commandLineParsed = false;
};

} // namespace transfer_util
} // namespace stk

#endif // STK_TRANSFER_UTIL_TRANSFERMAINPARSER_HPP
