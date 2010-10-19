#ifndef Ioi_Types_H
#define Ioi_Types_H
#include <vector>
#include <string>
#include <stk_util/util/string_case_compare.hpp>
#include <Ioi_VariableNamePair.h>

namespace Ioi {

  /** Ouput mesh type -- default is 'REFINED'.  'UNREFINED'
means that variables will be output on genesis-mesh
(the mesh read from the input file).  'SURFACE' means
that only the surface of the model will be output.
  */
  enum OutputMesh{REFINED=0,UNREFINED=1,BLOCK_SURFACE=2,EXPOSED_SURFACE=4};
  
  typedef std::vector<VariableNamePair> VarContainer;
  
  enum InputOutputType {RESULTS=1, RESTART=2, HEARTBEAT=4, HISTORY=8, MESH_MODEL=16, RESTART_IN=32};
  enum RestartNeed {NONE=1, REQUIRED=2, OPTIONAL=3};

  enum FieldRoleType {UNKNOWN = -1,
		MODEL,                 // One-time read of non-transient field (i.e., shell-thickness)
		TRANSIENT,              // One-time read of transient data (temperature @ time=3.2)
		TIME_INTERPOLATION,      // Interpolate data at multiple times
		COPY_TRANSFER_INTERPOLATION,  // Interpolate data at multiple times from another model
		FIXED_TIME_INTERPOLATION, // One-time interpolation at specified time
		CLOSEST_TIME,              // Multiple read, no interpolation, pick closest time.
		THREE_STATE};

  enum FieldType { INVALID_TYPE, NODE, EDGE, FACE, ELEMENT, GLOBAL, NODESET };

  inline FieldType get_field_type(const std::string &type) {
    switch(type.c_str()[0]) {
    case 'g':
    case 'G':
return GLOBAL;
    case 'n':
    case 'N':
if (stk::equal_case(type.c_str(), "nodeset") == 0) {
  return NODESET;
} else {
  return NODE;
}
    case 'e':
    case 'E':
switch(type.c_str()[1]) {
case 'd':
case 'D':
  return EDGE;
  break;
case 'l':
case 'L':
  return ELEMENT;
}
break;
    case 'f':
    case 'F':
return FACE;
    default:
return INVALID_TYPE;
    }
    return INVALID_TYPE;  // Should not get here.  Quiet the compiler.
  }
}
#endif
