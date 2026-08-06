/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <utility>
#include "TransferMainParser.hpp"
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <stk_transfer_util/LogMessage.hpp>
#include <stk_util/Version.hpp>
#include <stk_util/diag/StringUtil.hpp>
#include "stk_transfer_util/TransferMainSettings.hpp"

namespace stk {
namespace transfer_util {

TransferMainParser::TransferMainParser(MPI_Comm comm)
  : m_comm(comm),
    m_cmdLineParser(comm),
    m_quickExample("Basic example: mpirun --np 1 stk_transfer -s send.exo -r recv.exo -o xfer_result.exo"),
    m_longExamples("Examples:\nmpirun --np 1 stk_transfer -s send.exo -r recv.exo -o xfer_result.exo --xfer-type INTERP --recv-type FACE_CENTROID --part-list block_1:surface_1\n"
                   "mpirun --np 1 stk_transfer -s send.exo -r recv.exo -o xfer_result.exo --xfer-type PATCH --recv-type ELEMENT_CENTROID")
{
  add_options_to_parser();
}

void TransferMainParser::add_options_to_parser()
{

  stk::CommandLineOption sendMeshFilename{m_optionNames.sendMesh, "s",
          "Mesh to read from, i.e., send mesh."};
  stk::CommandLineOption recvMeshFilename{m_optionNames.recvMesh, "r",
          "Mesh to transfer data to, i.e., receive mesh."};
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
  m_cmdLineParser.add_required<std::string>(sendMeshFilename);
  m_cmdLineParser.add_required<std::string>(recvMeshFilename);
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

void TransferMainParser::parse_command_line_options(int argc, const char** argv)
{
  stk::CommandLineParser::ParseState parseState = m_cmdLineParser.parse(argc, argv);

  if (parseState == stk::CommandLineParser::ParseComplete) {
    m_transferParserStatus = TransferParserStatus::SUCCESS;
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
  const std::string sendMesh = m_cmdLineParser.get_option_value<std::string>(m_optionNames.sendMesh);
  settings.set_sendMesh_filename(sendMesh);
}

void TransferMainParser::set_recvMesh_filename(TransferMainSettings& settings)
{
  const std::string recvMesh = m_cmdLineParser.get_option_value<std::string>(m_optionNames.recvMesh);
  settings.set_recvMesh_filename(recvMesh);
}

void TransferMainParser::set_outputMesh_filename(TransferMainSettings& settings)
{
  const std::string recvMesh = m_cmdLineParser.get_option_value<std::string>(m_optionNames.recvMesh);
  std::string outputMesh = "transferred_" + recvMesh;
  if (m_cmdLineParser.is_option_parsed(m_optionNames.outputMesh)) {
    outputMesh = m_cmdLineParser.get_option_value<std::string>(m_optionNames.outputMesh);
  }

  settings.set_outputMesh_filename(outputMesh);
}

void TransferMainParser::set_transfer_type(TransferMainSettings& settings)
{
  const std::string transferType = m_cmdLineParser.is_option_parsed(m_optionNames.transferType) ?
    sierra::make_upper(m_cmdLineParser.get_option_value<std::string>(m_optionNames.transferType)) : "INTERP";

  STK_ThrowRequireMsg(settings.set_transfer_type(transferType),
                      "Invalid option " + transferType + " provided for transfer type.");
}

void TransferMainParser::set_2d_to_3d_mapping_type(TransferMainSettings& settings)
{
  const std::string mappingType = m_cmdLineParser.is_option_parsed(m_optionNames.mappingType2dTo3d) ?
    sierra::make_upper(m_cmdLineParser.get_option_value<std::string>(m_optionNames.mappingType2dTo3d)) : "ZPLANE";

  STK_ThrowRequireMsg(settings.set_2d_to_3d_mapping_type(mappingType),
                      "Invalid option " + mappingType + " provided for 2d-to-3d-mapping-type.");
}

void TransferMainParser::set_coord_transf_z_expr(TransferMainSettings& settings)
{
  const std::string expression = m_cmdLineParser.is_option_parsed(m_optionNames.coordTransformZExpr) ?
    m_cmdLineParser.get_option_value<std::string>(m_optionNames.coordTransformZExpr) : "z=0";

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
  }
}

void TransferMainParser::set_extrapolate_option(TransferMainSettings& settings)
{
  const std::string policy = m_cmdLineParser.is_option_parsed(m_optionNames.ExtrapolateOption) ?
    sierra::make_upper(m_cmdLineParser.get_option_value<std::string>(m_optionNames.ExtrapolateOption)) : "IGNORE";

  STK_ThrowRequireMsg(settings.set_extrapolate_option(policy),
                      "Invalid option " + policy + " provided for extrapolate option.");
}

void TransferMainParser::set_master_elements_name(TransferMainSettings& settings)
{
  const std::string masterElemName = m_cmdLineParser.is_option_parsed(m_optionNames.UseMasterElements) ?
    sierra::make_upper(m_cmdLineParser.get_option_value<std::string>(m_optionNames.UseMasterElements)) : sierra::make_upper(default_master_element_name);

  settings.set_master_elements_name(masterElemName);
}

void TransferMainParser::set_recv_type(TransferMainSettings& settings)
{
  const std::string recvType = m_cmdLineParser.is_option_parsed(m_optionNames.recvType) ?
    sierra::make_upper(m_cmdLineParser.get_option_value<std::string>(m_optionNames.recvType)) : "NODE";

  STK_ThrowRequireMsg(settings.set_recv_type(recvType),
                      "Invalid option " + recvType + " provided for receive type.");
}

void TransferMainParser::set_time_steps_spec(TransferMainSettings& settings)
{
  const std::string timeSteps = m_cmdLineParser.is_option_parsed(m_optionNames.timeSteps) ?
    sierra::make_upper(m_cmdLineParser.get_option_value<std::string>(m_optionNames.timeSteps)) : "ALL";
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
  }
}

  void set_3d_to_3d_axisymmetric_transfer_params(TransferMainSettings& settings);

} // namespace transfer_util
} // namespace stk

