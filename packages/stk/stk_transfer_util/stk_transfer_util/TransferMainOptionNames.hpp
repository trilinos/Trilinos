/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_TRANSFER_UTIL_TRANSFERMAINOPTIONNAMES_HPP
#define STK_TRANSFER_UTIL_TRANSFERMAINOPTIONNAMES_HPP

#include <string>

namespace stk {
namespace transfer_util {

// Long names of the stk_transfer options. These are used both as command-line
// option names (--send-mesh) and as top-level keys in a YAML input file (send-mesh:).
struct OptionNames
{
  const std::string sendMesh = "send-mesh";
  const std::string recvMesh = "recv-mesh";
  const std::string outputMesh = "output-mesh";
  const std::string fieldList = "field-list";
  const std::string partList = "part-list";
  const std::string ExtrapolateOption = "extrapolate-option";
  const std::string UseMasterElements = "use-master-elements";
  const std::string recvType = "recv-type";
  const std::string transferType = "xfer-type";
  const std::string timeSteps = "time-steps";
  const std::string mappingType2dTo3d = "2d-to-3d-mapping-type";
  const std::string twoDto3DAxisymmetric = "2d-to-3d-axisymmetric";
  const std::string threeDto3DAxisymmetric = "3d-to-3d-axisymmetric";
  const std::string coordTransformZExpr = "2d-to-3d-z-coord";
  const std::string inputFile = "input-file";
};

} // namespace transfer_util
} // namespace stk

#endif // STK_TRANSFER_UTIL_TRANSFERMAINOPTIONNAMES_HPP
