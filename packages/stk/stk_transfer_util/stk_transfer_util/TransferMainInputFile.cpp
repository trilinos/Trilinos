/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_util/stk_config.h>

#include <stdexcept>
#include <stk_transfer_util/TransferMainInputFile.hpp>
#include <stk_util/util/ReportHandler.hpp>

#ifdef STK_HAS_YAML_CPP
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <map>
#include <set>

#include "yaml-cpp/yaml.h"
#endif

namespace stk {
namespace transfer_util {

#ifdef STK_HAS_YAML_CPP

namespace {

// Prefix for all input-file error messages. (STK_ThrowRequireMsg prepends "Error: ".)
const char* const inputFileError = "stk_transfer input file: ";

std::string join_keys(const std::vector<std::string>& keys)
{
  std::string joined;
  for (const std::string& key : keys) {
    joined += (joined.empty() ? "" : ", ") + key;
  }
  return joined;
}

std::string node_type_name(const YAML::Node& node)
{
  if (node.IsNull()) { return "null (no value)"; }
  if (node.IsScalar()) { return "a scalar"; }
  if (node.IsSequence()) { return "a sequence"; }
  if (node.IsMap()) { return "a map"; }
  return "undefined";
}

std::string get_key_string(const YAML::Node& keyNode, const std::string& context)
{
  STK_ThrowRequireMsg(keyNode.IsScalar(), inputFileError << "keys in '" << context << "' must be plain strings.");
  return keyNode.as<std::string>();
}

// Calls func(key, value) for every entry of the map 'node'. Errors on keys that are
// not in allowedKeys and on keys that appear more than once.
void for_each_map_entry(const YAML::Node& node, const std::string& context,
                        const std::vector<std::string>& allowedKeys,
                        const std::function<void(const std::string&, const YAML::Node&)>& func)
{
  std::set<std::string> seenKeys;
  for (auto& entry : node) {
    const std::string key = get_key_string(entry.first, context);

    const bool isAllowed = std::find(allowedKeys.begin(), allowedKeys.end(), key) != allowedKeys.end();
    STK_ThrowRequireMsg(isAllowed, inputFileError << "unrecognized key '" << key << "' in " << context
                                                  << ". Valid keys are: " << join_keys(allowedKeys) << ".");

    const bool isFirstOccurrence = seenKeys.insert(key).second;
    STK_ThrowRequireMsg(
        isFirstOccurrence, inputFileError << "key '" << key << "' is specified more than once in " << context << ".");

    func(key, entry.second);
  }
}

std::string get_scalar_string(const YAML::Node& node, const std::string& name)
{
  STK_ThrowRequireMsg(node.IsScalar(),
      inputFileError << "'" << name << "' must be a single value, but is " << node_type_name(node) << ".");
  const std::string value = node.as<std::string>();
  STK_ThrowRequireMsg(!value.empty(), inputFileError << "'" << name << "' must not be empty.");
  return value;
}

double get_scalar_double(const YAML::Node& node, const std::string& name)
{
  STK_ThrowRequireMsg(
      node.IsScalar(), inputFileError << "'" << name << "' must be a number, but is " << node_type_name(node) << ".");
  double returnValue = 0.0;
  try {
    returnValue = node.as<double>();
  }
  catch (const YAML::Exception&) {
    STK_ThrowErrorMsg("'" + name + "' must be a number, but is '" + node.as<std::string>() + "'.");
  }
  return returnValue;
}

int get_scalar_int(const YAML::Node& node, const std::string& name)
{
  STK_ThrowRequireMsg(
      node.IsScalar(), inputFileError << "'" << name << "' must be an integer, but is " << node_type_name(node) << ".");
  int returnValue = -1;
  try {
    returnValue = node.as<int>();
  }
  catch (const YAML::Exception&) {
    STK_ThrowErrorMsg("'" + name + "' must be an integer, but is '" + node.as<std::string>() + "'.");
  }
  return returnValue;
}

std::array<double, 3> get_vector3(const YAML::Node& node, const std::string& name)
{
  STK_ThrowRequireMsg(node.IsSequence() && node.size() == 3,
      inputFileError << "'" << name << "' must be a sequence of 3 numbers, e.g. [0, 0, 1].");
  std::array<double, 3> values;
  for (unsigned i = 0; i < 3; ++i) {
    values[i] = get_scalar_double(node[i], name + "[" + std::to_string(i) + "]");
  }
  return values;
}

std::vector<std::string> get_string_list(const YAML::Node& node, const std::string& name)
{
  std::vector<std::string> values;
  if (node.IsNull()) {
    return values;
  }
  if (node.IsScalar()) {
    values.push_back(get_scalar_string(node, name));
    return values;
  }
  STK_ThrowRequireMsg(node.IsSequence(), inputFileError << "'" << name
                                                        << "' must be a name or a sequence of names, but is "
                                                        << node_type_name(node) << ".");
  for (std::size_t i = 0; i < node.size(); ++i) {
    values.push_back(get_scalar_string(node[i], name + "[" + std::to_string(i) + "]"));
  }
  return values;
}

//The following 'parse_*' functions parse the YAML syntax that
//this input-file capability supports.
//
// Input file syntax (all keys are top-level and optional):
//
//   send-mesh: send.exo
//   recv-mesh: recv.exo
//   output-mesh: output.exo
//   xfer-type: INTERP
//   extrapolate-option: IGNORE
//   use-master-elements: SIERRA
//   recv-type: NODE
//   2d-to-3d-mapping-type: ZPLANE
//   2d-to-3d-z-coord: z=x+y
//   time-steps: LAST                               # or ALL, FIRST, "1:5:2", or
//   time-steps: {begin: 1, end: 5, stride: 2}      # stride is optional (default 1)
//   field-list:
//     - field_1                                    # same name on send and recv mesh
//     - {send: send_field_2, recv: recv_field_2}
//   part-list:
//     send: [block_1, block_2]                     # omitted or empty => all parts
//     recv: [surface_1]
//   2d-to-3d-axisymmetric:                         # use {} for default parameters
//     axis: [0, 0, 1]
//     translation: [0, 0, 0]
//   3d-to-3d-axisymmetric:
//     theta-min: 30                                # degrees, required
//     theta-max: 60                                # degrees, required
//     axis: [0, 0, 1]
//     translation: [0, 0, 0]
//

std::pair<std::string, std::string> parse_field_entry(const YAML::Node& entry, const std::string& entryName)
{
  if (entry.IsScalar()) {
    const std::string fieldName = get_scalar_string(entry, entryName);
    STK_ThrowRequireMsg(fieldName.find(':') == std::string::npos,
        inputFileError << "'" << entryName << "' (" << fieldName
                       << ") contains ':'. To transfer between "
                          "fields with different names, use the map form {send: send_field, recv: recv_field}.");
    STK_ThrowRequireMsg(fieldName.find(',') == std::string::npos,
        inputFileError << "'" << entryName << "' (" << fieldName
                       << ") contains ','. Specify multiple "
                          "fields as a YAML sequence, e.g. [field_1, field_2].");
    return {fieldName, fieldName};
  }

  if (!entry.IsMap()) {
    STK_ThrowErrorMsg("'" + entryName + "' must be a field name or a {send: name, recv: name} map, but is "
                    + node_type_name(entry) + ".");
  }

  std::optional<std::string> sendField;
  std::optional<std::string> recvField;
  for_each_map_entry(entry, "'" + entryName + "'", {"send", "recv"},
    [&](const std::string& key, const YAML::Node& value) {
      if (key == "send") { sendField = get_scalar_string(value, entryName + ".send"); }
      else               { recvField = get_scalar_string(value, entryName + ".recv"); }
    });
  STK_ThrowRequireMsg(sendField.has_value(), inputFileError << "'" << entryName << "' requires a 'send' field name.");
  return {*sendField, recvField.value_or(*sendField)};
}

// A sequence of entries, each either a field name or a {send, recv} map. A single
// entry may also be given without the enclosing sequence.
std::vector<std::pair<std::string, std::string>> parse_field_list(const YAML::Node& node, const std::string& name)
{
  if (node.IsScalar() || node.IsMap()) {
    return {parse_field_entry(node, name)};
  }

  STK_ThrowRequireMsg(node.IsSequence(), inputFileError << "'" << name
                                                        << "' must be a sequence of field names and/or "
                                                           "{send: name, recv: name} maps, but is "
                                                        << node_type_name(node) << ".");

  std::vector<std::pair<std::string, std::string>> fields;
  for (std::size_t i = 0; i < node.size(); ++i) {
    fields.push_back(parse_field_entry(node[i], name + "[" + std::to_string(i) + "]"));
  }
  return fields;
}

TransferInputFilePartLists parse_part_list(const YAML::Node& node, const std::string& name)
{
  STK_ThrowRequireMsg(node.IsMap(), inputFileError
                                        << "'" << name
                                        << "' must be a map with optional 'send' and 'recv' "
                                           "lists of part names, e.g. {send: [block_1], recv: [block_2]}, but is "
                                        << node_type_name(node) << ".");

  TransferInputFilePartLists parts;
  for_each_map_entry(node, "'" + name + "'", {"send", "recv"},
    [&](const std::string& key, const YAML::Node& value) {
      if (key == "send") { parts.send = get_string_list(value, name + ".send"); }
      else               { parts.recv = get_string_list(value, name + ".recv"); }
    });
  return parts;
}

std::string parse_time_steps(const YAML::Node& node, const std::string& name)
{
  if (node.IsScalar()) {
    return get_scalar_string(node, name);
  }

  STK_ThrowRequireMsg(node.IsMap(), inputFileError << "'" << name
                                                   << "' must be ALL, FIRST, LAST, \"begin:end:stride\", "
                                                      "or a map {begin: 1, end: 5, stride: 1}, but is "
                                                   << node_type_name(node) << ".");

  std::optional<int> begin;
  std::optional<int> end;
  int stride = 1;
  for_each_map_entry(node, "'" + name + "'", {"begin", "end", "stride"},
    [&](const std::string& key, const YAML::Node& value) {
      const int intValue = get_scalar_int(value, name + "." + key);
      if (key == "begin")    { begin = intValue; }
      else if (key == "end") { end = intValue; }
      else                   { stride = intValue; }
    });

  STK_ThrowRequireMsg(begin.has_value() && end.has_value(),
      inputFileError << "'" << name << "' requires both 'begin' and 'end' step numbers.");
  STK_ThrowRequireMsg(*begin > 0 && *end > 0 && stride > 0 && *begin <= *end,
      inputFileError << "'" << name
                     << "' requires positive step numbers with begin <= end and a "
                        "positive stride.");

  return std::to_string(*begin) + ":" + std::to_string(*end) + ":" + std::to_string(stride);
}

TwoDTo3DAxisymmetricParams parse_2d_to_3d_axisymmetric(const YAML::Node& node, const std::string& name)
{
  STK_ThrowRequireMsg(node.IsMap(), inputFileError
                                        << "'" << name
                                        << "' must be a map with optional 'axis' and "
                                           "'translation' entries (use {} for the default parameters), but is "
                                        << node_type_name(node) << ".");

  TwoDTo3DAxisymmetricParams params;
  for_each_map_entry(node, "'" + name + "'", {"axis", "translation"},
    [&](const std::string& key, const YAML::Node& value) {
      if (key == "axis") { params.axis = get_vector3(value, name + ".axis"); }
      else               { params.trans = get_vector3(value, name + ".translation"); }
    });
  return params;
}

ThreeDTo3DAxisymmetricParams parse_3d_to_3d_axisymmetric(const YAML::Node& node, const std::string& name)
{
  STK_ThrowRequireMsg(
      node.IsMap(), inputFileError << "'" << name
                                   << "' must be a map with required 'theta-min' and "
                                      "'theta-max' entries (in degrees) and optional 'axis' and 'translation' entries, "
                                      "but is "
                                   << node_type_name(node) << ".");

  const double pi = 2*std::atan2(1, 0);
  ThreeDTo3DAxisymmetricParams params;
  bool haveThetaMin = false;
  bool haveThetaMax = false;
  for_each_map_entry(node, "'" + name + "'", {"theta-min", "theta-max", "axis", "translation"},
    [&](const std::string& key, const YAML::Node& value) {
      if (key == "theta-min") {
        params.theta_min = (pi/180.0)*get_scalar_double(value, name + ".theta-min");
        haveThetaMin = true;
      }
      else if (key == "theta-max") {
        params.theta_max = (pi/180.0)*get_scalar_double(value, name + ".theta-max");
        haveThetaMax = true;
      }
      else if (key == "axis") {
        params.axis = get_vector3(value, name + ".axis");
      }
      else {
        params.trans = get_vector3(value, name + ".translation");
      }
    });

  STK_ThrowRequireMsg(haveThetaMin && haveThetaMax,
      inputFileError << "'" << name << "' requires both 'theta-min' and 'theta-max' (in degrees).");
  return params;
}

}

bool transfer_input_file_supported()
{
  return true;
}

TransferInputFileOptions parse_transfer_input_file_yaml(const std::string& yamlText)
{
  YAML::Node root;
  try {
    root = YAML::Load(yamlText);
  }
  catch (const YAML::Exception& e) {
    STK_ThrowErrorMsg(std::string("invalid YAML syntax: ") + e.what());
  }

  STK_ThrowRequireMsg(root.IsMap(), inputFileError << "expected a map of 'option-name: value' entries at the top "
                                                      "level, but the file contains "
                                                   << node_type_name(root) << ".");

  const OptionNames names;
  TransferInputFileOptions options;

  using StringMember = std::optional<std::string> TransferInputFileOptions::*;
  const std::map<std::string, StringMember> scalarOptions = {
    {names.sendMesh,            &TransferInputFileOptions::sendMesh},
    {names.recvMesh,            &TransferInputFileOptions::recvMesh},
    {names.outputMesh,          &TransferInputFileOptions::outputMesh},
    {names.transferType,        &TransferInputFileOptions::transferType},
    {names.ExtrapolateOption,   &TransferInputFileOptions::extrapolateOption},
    {names.UseMasterElements,   &TransferInputFileOptions::masterElements},
    {names.recvType,            &TransferInputFileOptions::recvType},
    {names.mappingType2dTo3d,   &TransferInputFileOptions::mappingType2dTo3d},
    {names.coordTransformZExpr, &TransferInputFileOptions::coordTransformZExpr},
  };

  const std::vector<std::string> allowedKeys = {
    names.sendMesh, names.recvMesh, names.outputMesh, names.transferType, names.fieldList, names.partList,
    names.ExtrapolateOption, names.UseMasterElements, names.recvType, names.timeSteps, names.mappingType2dTo3d,
    names.twoDto3DAxisymmetric, names.threeDto3DAxisymmetric, names.coordTransformZExpr
  };

  try {
    for_each_map_entry(root, "the top level", allowedKeys,
      [&](const std::string& key, const YAML::Node& value) {
        auto scalarIt = scalarOptions.find(key);
        if (scalarIt != scalarOptions.end()) {
          options.*(scalarIt->second) = get_scalar_string(value, key);
        }
        else if (key == names.timeSteps) {
          options.timeSteps = parse_time_steps(value, key);
        }
        else if (key == names.fieldList) {
          options.fields = parse_field_list(value, key);
        }
        else if (key == names.partList) {
          options.parts = parse_part_list(value, key);
        }
        else if (key == names.twoDto3DAxisymmetric) {
          options.twoDTo3DAxisymmetric = parse_2d_to_3d_axisymmetric(value, key);
        }
        else if (key == names.threeDto3DAxisymmetric) {
          options.threeDTo3DAxisymmetric = parse_3d_to_3d_axisymmetric(value, key);
        }
      });
  }
  catch (const YAML::Exception& e) {
    STK_ThrowErrorMsg(e.what());
  }

  return options;
}

#else // STK_HAS_YAML_CPP

bool transfer_input_file_supported()
{
  return false;
}

TransferInputFileOptions parse_transfer_input_file_yaml(const std::string& /*yamlText*/)
{
  throw std::runtime_error("Error: stk_transfer was built without the yaml-cpp library, so --input-file is not supported. "
                           "Specify the transfer options on the command line instead.");
}

#endif // STK_HAS_YAML_CPP

} // namespace transfer_util
} // namespace stk
