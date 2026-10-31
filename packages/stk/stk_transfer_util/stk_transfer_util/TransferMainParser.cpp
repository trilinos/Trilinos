/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "TransferMainParser.hpp"

#include <climits>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <stk_transfer_util/LogMessage.hpp>
#include <stk_util/Version.hpp>
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <stk_util/diag/StringUtil.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/parallel/mpi_filebuf.hpp>
#include <utility>

#include "stk_transfer_util/TransferMainSettings.hpp"

namespace stk {
namespace transfer_util {

TransferMainParser::TransferMainParser(MPI_Comm comm)
    : m_comm(comm),
      m_cmdLineParser(comm),
      m_quickExample("Basic example: mpirun --np 1 stk_transfer -s send.exo -r recv.exo -o xfer_result.exo"),
      m_longExamples(
          "Examples:\nmpirun --np 1 stk_transfer -s send.exo -r recv.exo -o xfer_result.exo --xfer-type INTERP "
          "--recv-type FACE_CENTROID --part-list block_1:surface_1\n"
          "mpirun --np 1 stk_transfer -s send.exo -r recv.exo -o xfer_result.exo --xfer-type PATCH --recv-type "
          "ELEMENT_CENTROID\n"
          "mpirun --np 1 stk_transfer --input-file transfer.yaml"),
      m_transferParserStatus(TransferParserStatus::PARSE_ERROR)
{
  add_options_to_parser();
}

void TransferMainParser::add_options_to_parser()
{
  stk::CommandLineOption sendMeshFilename{m_optionNames.sendMesh, "s",
      "Mesh to read from, i.e., send mesh. (required, unless specified in the input file)"};
  stk::CommandLineOption recvMeshFilename{m_optionNames.recvMesh, "r",
      "Mesh to transfer data to, i.e., receive mesh. (required, unless specified in the input file)"};
  stk::CommandLineOption outputMeshFilename{m_optionNames.outputMesh, "o",
          "Mesh to write transferred data to (optional). (default: transferred_[recvMesh])"};
  stk::CommandLineOption transferType{m_optionNames.transferType, "x",
          "Set type of transfer. Options: INTERP, COPY, PATCH  (default: INTERP)"};
  stk::CommandLineOption fieldList{m_optionNames.fieldList, "f",
          "Comma-separated list of send/receive field names to transfer. To transfer between fields of different names, use a colon separator (send_field_1:recv_field_1). To transfer a field of the same name, specify the field name once (field_2). Note: When only one field name is specified, a field with that name will be created on the receive mesh if it does not exist already. Syntax example: field_1, field_2:recv_field_2. If this option is not specified, all fields will be transferred. "};
  stk::CommandLineOption partList{m_optionNames.partList, "p",
          "Colon-separated list of send/receive part names to transfer. Syntax example: send_part1,send_part2:recv_part1,recv_part2. If this option is not specified, all element block parts will be transferred (ex. :part1 will use all send mesh parts and receive mesh part1). "};
  stk::CommandLineOption ExtrapolateOption{m_optionNames.ExtrapolateOption, "e",
          "Behavior for receive mesh points that are not inside any receive mesh element. "
          "Options are: IGNORE, EXTRAPOLATE, TRUNCATE, PROJECT, and ABORT. "
          "If this option is not specified, it will default to IGNORE."};

#ifdef STK_BUILT_FOR_SIERRA
  stk::CommandLineOption useMasterElements{m_optionNames.UseMasterElements, "u",
          "Use specified MasterElements implementation, INTREPID2 or SIERRA. (default: "+default_master_element_name+")"};
#else
  stk::CommandLineOption useMasterElements{m_optionNames.UseMasterElements, "u",
          "Use specified MasterElements implementation. (default: "+default_master_element_name+") INTREPID2 is the only option unless STK is built inside Sierra."};
#endif

  stk::CommandLineOption recvType{m_optionNames.recvType, "t",
          "Set receive type used for transfer. Options: NODE, [ELEMENT,FACE,EDGE]_[CENTROID,GAUSS_POINT] (default: NODE)"};
  stk::CommandLineOption timeSteps{m_optionNames.timeSteps, "",
          "Time steps to be transferred. Options: ALL, FIRST, LAST, or colon-separated begin:end:stride. Note that begin, end and stride are 1-based integer step numbers, not time-values. (default: ALL)"};
  stk::CommandLineOption mappingType2dTo3d{m_optionNames.mappingType2dTo3d, "",
          "Mapping type to use for transferring from a 2D mesh to a 3D mesh. Valid values are ZPLANE or EXTRUDE. EXTRUDE maps the 2D mesh to all z-values of the 3D mesh. ZPLANE maps the 2D mesh to the Z=0 plane. (default: ZPLANE)"};

  stk::CommandLineOption twoDto3DAxisymmetric{m_optionNames.twoDto3DAxisymmetric, "2d3daxi",
          "Enables a 2D to 3D axisymmetric coordinate transform, and allows the user to specify parameters for the transform."
          " The receive mesh coordinates are transformed "
          "by converting to cylindrical coordinates and then discarding the theta coordinate "
          "so that only (r, z) remain.  The parameters are either 0, 3, or 6 comma separated values: "
          "axis_x,axis_y,axis_z,trans_x,trans_y,trans_z.  These are:\n"
          "  * (axis_x, axis_y,axis_z): vector along the axis of the cylinder, default (0, 0, 1)\n"
          "  * (trans_x,trans_y,trans_z): translation to apply before other transforms, default (0, 0,0)\n"
          "Examples:\n\n"
          "  --" + m_optionNames.twoDto3DAxisymmetric + " 0,1,0"
          "\nwill rotate the coordinates such that the y axis is used for the axis of the cylinder\n\n"
          "  --" + m_optionNames.twoDto3DAxisymmetric + " 0,1,0,4,5,6"
          "\nwill add (4,5,6) to the coordinates of each point before applying the rotation\n\n"
          "The default parameters are 0,0,1,0,0,0 (z axis and no translation)"};

  stk::CommandLineOption threeDto3DAxisymmetric{m_optionNames.threeDto3DAxisymmetric, "3d3daxi",
          "Allows the user to specify parameters for the 3D to 3D axisymmetric coordinate transform."
          " The send mesh should be a wedge and the receive mesh should be a shape obtained by revolving a wedge"
          " about the tip of the wedge.  The parameters for the transform are 8 comma separated values: "
          " theta_min,theta_max,axis_x,axis_y,axis_z,trans_x,trans_y,trans_z.  These are:\n"
          "  * theta_min, theta_max: the range of theta values for the send mesh wedge (in degrees).  Both values mandatory\n"
          "  * (axis_x, axis_y,axis_z): a vector along the tip of the wedge.  Defines the axis used for the conversion"
          " to cylindrical coordinates.  Default (0, 0, 1)\n"
          "  * (trans_x, trans_y, trans_z): an optional translation to apply before the other parts of the transform, default (0, 0,0)"};

  stk::CommandLineOption coordTransformZExpr{m_optionNames.coordTransformZExpr, "",
          "Expression to use for the z-coordinate for a 2D mesh to a 3D mesh transfer. Valid values are numerical constants or expressions containing x and y, e.g. z=x+y (default: z=0)"};
  stk::CommandLineOption inputFile{m_optionNames.inputFile, "i",
      "YAML input file specifying transfer options. The top-level keys are the long names of the "
      "command-line options (e.g. 'send-mesh: send.exo'). Options given on the command line take "
      "precedence over options in the input file. List-type options use YAML sequences and maps:\n"
      "  field-list: [field_1, {send: send_field_2, recv: recv_field_2}]\n"
      "  part-list: {send: [block_1, block_2], recv: [surface_1]}\n"
      "  time-steps: LAST    (or ALL, FIRST, \"1:5:2\", or {begin: 1, end: 5, stride: 2})\n"
      "  2d-to-3d-axisymmetric: {axis: [0, 0, 1], translation: [0, 0, 0]}    (use {} for defaults)\n"
      "  3d-to-3d-axisymmetric: {theta-min: 30, theta-max: 60, axis: [0, 0, 1], translation: [0, 0, 0]}" +
          std::string(
              transfer_input_file_supported()
                  ? ""
                  : "\nNOTE: this build of stk_transfer does not support input files (built without yaml-cpp).")};

  m_cmdLineParser.add_optional<std::string>(inputFile);
  m_cmdLineParser.add_optional<std::string>(sendMeshFilename);
  m_cmdLineParser.add_optional<std::string>(recvMeshFilename);
  m_cmdLineParser.add_optional<std::string>(outputMeshFilename);
  m_cmdLineParser.add_optional<std::string>(transferType);
  m_cmdLineParser.add_optional<std::string>(fieldList);
  m_cmdLineParser.add_optional<std::string>(partList);
  m_cmdLineParser.add_optional<std::string>(ExtrapolateOption);
  m_cmdLineParser.add_optional<std::string>(useMasterElements);
  m_cmdLineParser.add_optional<std::string>(recvType);
  m_cmdLineParser.add_optional<std::string>(timeSteps);
  m_cmdLineParser.add_optional<std::string>(mappingType2dTo3d);
  m_cmdLineParser.add_optional_implicit<std::string>(twoDto3DAxisymmetric, "0,0,1,0,0,0");
  m_cmdLineParser.add_optional<std::string>(threeDto3DAxisymmetric);
  m_cmdLineParser.add_optional<std::string>(coordTransformZExpr);

  m_cmdLineParser.disallow_unrecognized();
}

namespace
{

// STK_ThrowRequireMsg exceptions contain the failed expression, the source location and a
// stack trace, followed by "Error: <message>". For input errors, users only need the message.
std::string user_error_message(const std::string& what)
{
  std::string message = what;
  const std::size_t locationPos = what.find("Error occurred at: ");
  if (locationPos != std::string::npos) {
    const std::size_t messagePos = what.find("\nError: ", locationPos);
    if (messagePos != std::string::npos) {
      message = what.substr(messagePos + 1);
    }
  }
  const std::size_t lastNonSpace = message.find_last_not_of(" \t\n\r");
  return (lastNonSpace == std::string::npos) ? "" : message.substr(0, lastNonSpace + 1);
}

}  // namespace

void TransferMainParser::parse_command_line_options(int argc, const char** argv)
{
  m_commandLineParsed = true;
  stk::CommandLineParser::ParseState parseState = m_cmdLineParser.parse(argc, argv);

  if (parseState == stk::CommandLineParser::ParseComplete) {
    m_transferParserStatus = TransferParserStatus::SUCCESS;

    try {
      read_input_file_if_specified();
      require_mesh_specified(m_optionNames.sendMesh, m_inputFileOptions.sendMesh);
      require_mesh_specified(m_optionNames.recvMesh, m_inputFileOptions.recvMesh);
    } catch (std::exception& e) {
      m_transferParserStatus = TransferParserStatus::PARSE_ERROR;
      if (stk::parallel_machine_rank(m_comm) == 0) {
        std::cerr << user_error_message(e.what()) << std::endl;
      }
      stk::outputP0() << m_quickExample << std::endl;
      stk::outputP0() << "Use 'stk_transfer --help' for more information." << std::endl;
    }
  }
  else {
    switch(parseState) {
      case stk::CommandLineParser::ParseError:
        m_transferParserStatus = TransferParserStatus::PARSE_ERROR;
        stk::outputP0() << m_quickExample << std::endl;
        stk::outputP0() << "Use 'stk_transfer --help' for more information." << std::endl;
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

void TransferMainParser::set_input_file_contents(const std::string& yamlText)
{
  STK_ThrowRequireMsg(!m_commandLineParsed,
      "TransferMainParser::set_input_file_contents must be called before parse_command_line_options.");
  m_inputFileOptions = parse_transfer_input_file_yaml(yamlText);
}

bool TransferMainParser::input_file_specified() const
{
  return m_cmdLineParser.is_option_parsed(m_optionNames.inputFile);
}

std::string TransferMainParser::get_input_file_name() const
{
  return m_cmdLineParser.get_option_value<std::string>(m_optionNames.inputFile);
}

void TransferMainParser::read_input_file_if_specified()
{
  if (input_file_specified()) {
    const std::string fileName = get_input_file_name();
    STK_ThrowRequireMsg(transfer_input_file_supported(),
        "stk_transfer was built without the yaml-cpp library, so --"
            << m_optionNames.inputFile
            << " is not supported. Specify the transfer options on the command line instead.");
    m_inputFileOptions = parse_transfer_input_file_yaml(stk::read_and_broadcast_file(fileName, m_comm));
  }
}

void TransferMainParser::require_mesh_specified(
    const std::string& optionName, const std::optional<std::string>& fileValue) const
{
  const std::string inputFileHint = transfer_input_file_supported()
                                        ? " or as '" + optionName + ":' in the --" + m_optionNames.inputFile + " file"
                                        : "";
  STK_ThrowRequireMsg(get_string_option(optionName, fileValue).has_value(),
      "required option '--" << optionName << "' not found. Specify it on the command line" << inputFileHint << ".");
}

std::optional<std::string> TransferMainParser::get_string_option(
    const std::string& optionName, const std::optional<std::string>& fileValue) const
{
  if (m_cmdLineParser.is_option_parsed(optionName)) {
    return m_cmdLineParser.get_option_value<std::string>(optionName);
  }
  return fileValue;
}

TransferMainSettings TransferMainParser::generate_transfer_settings()
{
  TransferMainSettings settings;

  if (m_transferParserStatus == TransferParserStatus::SUCCESS) {
    set_num_input_processors(settings);
    set_num_output_processors(settings);
    set_sendMesh_filename(settings);
    set_recvMesh_filename(settings);
    set_outputMesh_filename(settings);
    set_transfer_type(settings);
    set_transfer_fields(settings);
    set_transfer_parts(settings);
    set_extrapolate_option(settings);
    set_master_elements_name(settings);
    set_recv_type(settings);
    set_time_steps_spec(settings);
    set_2d_to_3d_mapping_type(settings);
    set_2d_to_3d_axisymmetric_transfer_params(settings);
    set_3d_to_3d_axisymmetric_transfer_params(settings);
    set_coord_transf_z_expr(settings);
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

void TransferMainParser::set_sendMesh_filename(TransferMainSettings& settings)
{
  settings.set_sendMesh_filename(get_string_option(m_optionNames.sendMesh, m_inputFileOptions.sendMesh).value_or(""));
}

void TransferMainParser::set_recvMesh_filename(TransferMainSettings& settings)
{
  settings.set_recvMesh_filename(get_string_option(m_optionNames.recvMesh, m_inputFileOptions.recvMesh).value_or(""));
}

void TransferMainParser::set_outputMesh_filename(TransferMainSettings& settings)
{
  const std::string recvMesh = get_string_option(m_optionNames.recvMesh, m_inputFileOptions.recvMesh).value_or("");
  const std::string outputMesh =
      get_string_option(m_optionNames.outputMesh, m_inputFileOptions.outputMesh).value_or("transferred_" + recvMesh);

  settings.set_outputMesh_filename(outputMesh);
}

void TransferMainParser::set_transfer_type(TransferMainSettings& settings)
{
  const std::string transferType = sierra::make_upper(
      get_string_option(m_optionNames.transferType, m_inputFileOptions.transferType).value_or("INTERP"));

  STK_ThrowRequireMsg(settings.set_transfer_type(transferType),
                      "Invalid option " + transferType + " provided for transfer type.");
}

void TransferMainParser::set_2d_to_3d_mapping_type(TransferMainSettings& settings)
{
  const std::string mappingType = sierra::make_upper(
      get_string_option(m_optionNames.mappingType2dTo3d, m_inputFileOptions.mappingType2dTo3d).value_or("ZPLANE"));

  STK_ThrowRequireMsg(settings.set_2d_to_3d_mapping_type(mappingType),
                      "Invalid option " + mappingType + " provided for 2d-to-3d-mapping-type.");
}

void TransferMainParser::set_coord_transf_z_expr(TransferMainSettings& settings)
{
  const std::string expression =
      get_string_option(m_optionNames.coordTransformZExpr, m_inputFileOptions.coordTransformZExpr).value_or("z=0");

  settings.set_coord_transf_z_expr(expression);
}

void TransferMainParser::set_transfer_fields(TransferMainSettings& settings)
{
  if (m_cmdLineParser.is_option_parsed(m_optionNames.fieldList))
  {
    const std::string fieldList = m_cmdLineParser.get_option_value<std::string>(m_optionNames.fieldList);
    std::vector<std::string> fieldNamePairs = stk::split_csv_string(fieldList);
    for (const std::string & fieldNames : fieldNamePairs) {
      auto fieldNamePair = stk::split_string(fieldNames, ':');
      STK_ThrowRequireMsg(fieldNamePair.size() <= 2 && fieldNamePair[0] != "",
                          "Please check for syntax errors in the provided field names list: " + fieldList);
      if (fieldNamePair.size() == 1) {
        settings.set_transfer_field({fieldNamePair[0], fieldNamePair[0]});
      }
      else {
        settings.set_transfer_field({fieldNamePair[0], fieldNamePair[1]});
      }
    }
  } else if (m_inputFileOptions.fields) {
    for (const std::pair<std::string, std::string>& fieldNamePair : *m_inputFileOptions.fields) {
      settings.set_transfer_field(fieldNamePair);
    }
  }
}

void TransferMainParser::set_transfer_parts(TransferMainSettings& settings)
{
  if (m_cmdLineParser.is_option_parsed(m_optionNames.partList))
  {
    const std::string partList = m_cmdLineParser.get_option_value<std::string>(m_optionNames.partList);

    auto partNames = stk::split_string(partList, ':');
    STK_ThrowRequireMsg(partNames.size() <= 2,
                        "Please check for syntax errors in the provided parts names list: " + partList);

    std::vector<std::string> sendPartNames;
    std::vector<std::string> recvPartNames;

    if (partNames.size() == 1) {
      sendPartNames = stk::split_csv_string(partNames[0]);
    }

    if (partNames.size() == 2) {
      sendPartNames = stk::split_csv_string(partNames[0]);
      recvPartNames = stk::split_csv_string(partNames[1]);
    }

    settings.set_transfer_send_parts(sendPartNames);
    settings.set_transfer_recv_parts(recvPartNames);
  } else if (m_inputFileOptions.parts) {
    settings.set_transfer_send_parts(m_inputFileOptions.parts->send);
    settings.set_transfer_recv_parts(m_inputFileOptions.parts->recv);
  }
}

void TransferMainParser::set_extrapolate_option(TransferMainSettings& settings)
{
  const std::string policy = sierra::make_upper(
      get_string_option(m_optionNames.ExtrapolateOption, m_inputFileOptions.extrapolateOption).value_or("IGNORE"));

  STK_ThrowRequireMsg(settings.set_extrapolate_option(policy),
                      "Invalid option " + policy + " provided for extrapolate option.");
}

void TransferMainParser::set_master_elements_name(TransferMainSettings& settings)
{
  const std::string masterElemName =
      sierra::make_upper(get_string_option(m_optionNames.UseMasterElements, m_inputFileOptions.masterElements)
              .value_or(default_master_element_name));

  settings.set_master_elements_name(masterElemName);
}

void TransferMainParser::set_recv_type(TransferMainSettings& settings)
{
  const std::string recvType =
      sierra::make_upper(get_string_option(m_optionNames.recvType, m_inputFileOptions.recvType).value_or("NODE"));

  STK_ThrowRequireMsg(settings.set_recv_type(recvType),
                      "Invalid option " + recvType + " provided for receive type.");
}

void TransferMainParser::set_time_steps_spec(TransferMainSettings& settings)
{
  const std::string timeSteps =
      sierra::make_upper(get_string_option(m_optionNames.timeSteps, m_inputFileOptions.timeSteps).value_or("ALL"));
  settings.set_time_steps_spec(timeSteps);
}

void TransferMainParser::set_2d_to_3d_axisymmetric_transfer_params(TransferMainSettings& settings)
{
  if (m_cmdLineParser.is_option_parsed(m_optionNames.twoDto3DAxisymmetric))
  {
    std::string params = m_cmdLineParser.get_option_value<std::string>(m_optionNames.twoDto3DAxisymmetric);
    std::vector<std::string> paramsVec = stk::split_csv_string(params);

    TwoDTo3DAxisymmetricParams params_struct;
    STK_ThrowRequireMsg(paramsVec.size() == 0 || paramsVec.size() == 3 || paramsVec.size() == 6,
                        "--" + m_optionNames.twoDto3DAxisymmetric + " takes either 0, 3, or 6 parameters, "
                       + std::to_string(paramsVec.size()) + " provided");
    if (paramsVec.size() >= 3)
    {
      params_struct.axis[0] = stk::stod(paramsVec[0]);
      params_struct.axis[1] = stk::stod(paramsVec[1]);
      params_struct.axis[2] = stk::stod(paramsVec[2]);
    }

    if (paramsVec.size() == 6)
    {
      params_struct.trans[0] = stk::stod(paramsVec[3]);
      params_struct.trans[1] = stk::stod(paramsVec[4]);
      params_struct.trans[2] = stk::stod(paramsVec[5]);
    }

    settings.set_2d_to_3d_axisymmetric_transfer_params(params_struct);
  } else if (m_inputFileOptions.twoDTo3DAxisymmetric) {
    settings.set_2d_to_3d_axisymmetric_transfer_params(*m_inputFileOptions.twoDTo3DAxisymmetric);
  }
}

void TransferMainParser::set_3d_to_3d_axisymmetric_transfer_params(TransferMainSettings& settings)
{
  if (m_cmdLineParser.is_option_parsed(m_optionNames.threeDto3DAxisymmetric))
  {
    std::string params = m_cmdLineParser.get_option_value<std::string>(m_optionNames.threeDto3DAxisymmetric);
    std::vector<std::string> paramsVec = stk::split_csv_string(params);

    ThreeDTo3DAxisymmetricParams params_struct;
    STK_ThrowRequireMsg(paramsVec.size() == 2 || paramsVec.size() == 5 || paramsVec.size() == 8,
                        "--" + m_optionNames.threeDto3DAxisymmetric + " takes either 2, 5, or 8 parameters, "
                       + std::to_string(paramsVec.size()) + " provided");

    double pi = 2*std::atan2(1, 0);
    params_struct.theta_min = (pi/180.0)*stk::stod(paramsVec[0]);
    params_struct.theta_max = (pi/180.0)*stk::stod(paramsVec[1]);
    if (paramsVec.size() >= 5)
    {
      params_struct.axis[0] = stk::stod(paramsVec[2]);
      params_struct.axis[1] = stk::stod(paramsVec[3]);
      params_struct.axis[2] = stk::stod(paramsVec[4]);
    }

    if (paramsVec.size() == 8)
    {
      params_struct.trans[0] = stk::stod(paramsVec[5]);
      params_struct.trans[1] = stk::stod(paramsVec[6]);
      params_struct.trans[2] = stk::stod(paramsVec[7]);
    }

    settings.set_3d_to_3d_axisymmetric_transfer_params(params_struct);
  } else if (m_inputFileOptions.threeDTo3DAxisymmetric) {
    settings.set_3d_to_3d_axisymmetric_transfer_params(*m_inputFileOptions.threeDTo3DAxisymmetric);
  }
}

} // namespace transfer_util
} // namespace stk

