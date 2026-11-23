/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_TRANSFER_UTIL_TRANSFERMAININPUTFILE_HPP
#define STK_TRANSFER_UTIL_TRANSFERMAININPUTFILE_HPP

#include <optional>
#include <string>
#include <utility>
#include <vector>
#include <stk_transfer_util/TransferMainOptionNames.hpp>
#include <stk_transfer_util/TransferMainSettings.hpp>

namespace stk {
namespace transfer_util {

// Options read from a stk_transfer YAML input file.
// Values are stored as given (no case conversion or validation of enumerated values);
// TransferMainParser applies the same processing to them as to the corresponding
// command-line values.
// See the .cpp file for a comment describing the yaml syntax that
// this input-file capability supports.
//
struct TransferInputFilePartLists
{
  std::vector<std::string> send;
  std::vector<std::string> recv;
};

// One entry of the grouped form of 'field-list'. Each group is transferred with its own
// transfer type, receive type, etc.; options not given in the group default to the
// corresponding top-level options.
struct TransferInputFileFieldGroup
{
  std::vector<std::pair<std::string, std::string>> fields;
  std::optional<std::string> transferType;
  std::optional<std::string> recvType;
  std::optional<std::string> extrapolateOption;
  std::optional<TransferInputFilePartLists> parts;
};

struct TransferInputFileOptions
{
  std::optional<std::string> sendMesh;
  std::optional<std::string> recvMesh;
  std::optional<std::string> outputMesh;
  std::optional<std::string> transferType;
  std::optional<std::string> extrapolateOption;
  std::optional<std::string> masterElements;
  std::optional<std::string> recvType;
  std::optional<std::string> timeSteps;
  std::optional<std::string> mappingType2dTo3d;
  std::optional<std::string> coordTransformZExpr;
  // At most one of 'fields' and 'fieldGroups' is set, depending on which form of
  // 'field-list' the file uses (a list of field names, or a list of field groups).
  std::optional<std::vector<std::pair<std::string, std::string>>> fields;
  std::optional<std::vector<TransferInputFileFieldGroup>> fieldGroups;
  std::optional<TransferInputFilePartLists> parts;
  std::optional<TwoDTo3DAxisymmetricParams> twoDTo3DAxisymmetric;
  std::optional<ThreeDTo3DAxisymmetricParams> threeDTo3DAxisymmetric; // theta values in radians
};

// True if stk_transfer_util was built with yaml-cpp, i.e. input files can be parsed.
bool transfer_input_file_supported();

// Parse the contents (not the file name) of a stk_transfer YAML input file.
// Throws an exception derived from std::exception, with a descriptive message, on any
// syntax error, or if stk_transfer_util was built without yaml-cpp.
TransferInputFileOptions parse_transfer_input_file_yaml(const std::string& yamlText);

} // namespace transfer_util
} // namespace stk

#endif // STK_TRANSFER_UTIL_TRANSFERMAININPUTFILE_HPP
