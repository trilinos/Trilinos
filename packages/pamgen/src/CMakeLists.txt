

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")

SET(PE_HEADERS "")
SET(PE_SOURCES "")

# Set headers (for IDEs and installing)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

# NOTE: Above, we add the XXX_config.h file so tht MS Visual Studio
# will know to include this file in the dependency analysis

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/mesh_spec_lt)
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/extra_functions)
INCLUDE_DIRECTORIES(${PACKAGE_SOURCE_DIR}/rtcompiler)

#
# Core Files
#

APPEND_SET(HEADERS
  create_inline_mesh.h
  )

# Set sources

APPEND_SET(SOURCES
  brick_inline_mesh_desc.C
  calc_decomp_cuts.C
  create_inline_mesh.C
  element_density_function.C
  geometry_transform.C
  inline_mesh_desc.C
  inline_mesh_driver_lt.C
  legacy_inline_mesh_desc.C
  parse_routines.C
  radial_inline_mesh_desc.C
  radial_trisection_inline_mesh_desc.C
  uns_inline_decomp.C
  )

#
# Parser Files
#

APPEND_SET(SOURCES
  ../parser/InputBlock.C
  ../parser/lexer.C
  ../parser/parse_table.C
  ../parser/parser.C
  ../parser/token.C
  ../parser/token_stream.C
  )

#
# MS_LT Files
#

APPEND_SET(HEADERS
  ../mesh_spec_lt/pamgen_im_exodusII.h
  ../mesh_spec_lt/pamgen_im_exodusII_l.h
  ../mesh_spec_lt/pamgen_im_exodusII_ext.h
  ../mesh_spec_lt/pamgen_im_ne_nemesisI.h
  ../mesh_spec_lt/pamgen_im_ne_nemesisI_l.h
  ../mesh_spec_lt/pamgen_mesh_specification.h
  ../mesh_spec_lt/pamgen_fudges.h
  ../mesh_spec_lt/pamgen_element_dictionary.h
  ../mesh_spec_lt/pamgen_code_types.h
  )

APPEND_SET(SOURCES
  ../mesh_spec_lt/im_ex_c_interface.C
  ../mesh_spec_lt/im_ex_c_interface_l.C
  ../mesh_spec_lt/im_ne_c_interface.C
  ../mesh_spec_lt/im_ne_c_interface_l.C
  ../mesh_spec_lt/pamgen_mesh_specification.C
  ../mesh_spec_lt/pamgen_fudges.C
  ../mesh_spec_lt/pamgen_element_dictionary.C
  )

#
# RTCOMPILER Files
#

APPEND_SET(HEADERS
  ../rtcompiler/RTC_ArrayIndexRTC.hh
  ../rtcompiler/RTC_ArrayNumberRTC.hh
  ../rtcompiler/RTC_ArrayVarRTC.hh
  ../rtcompiler/RTC_Bessel_I_RTC.hh
  ../rtcompiler/RTC_BlockRTC.hh
  ../rtcompiler/RTC_ConditionalBlockRTC.hh
  ../rtcompiler/RTC_ExecutableRTC.hh
  ../rtcompiler/RTC_ForBlockRTC.hh
  ../rtcompiler/RTC_FunctionRTC.hh
  ../rtcompiler/RTC_IfElseifElseBlockRTC.hh
  ../rtcompiler/RTC_LineRTC.hh
  ../rtcompiler/RTC_NormalBlockRTC.hh
  ../rtcompiler/RTC_ObjectRTC.hh
  ../rtcompiler/RTC_OperatorRTC.hh
  ../rtcompiler/RTC_RegistrarRTC.hh
  ../rtcompiler/RTC_ScalarNumberRTC.hh
  ../rtcompiler/RTC_ScalarVarRTC.hh
  ../rtcompiler/RTC_TokenizerRTC.hh
  ../rtcompiler/RTC_ValueRTC.hh
  ../rtcompiler/RTC_VariableRTC.hh
  ../rtcompiler/RTC_WhileBlockRTC.hh
  ../rtcompiler/RTC_commonRTC.hh
  )

APPEND_SET(SOURCES
  ../rtcompiler/Bessel_I.C
  ../rtcompiler/Block.C
  ../rtcompiler/ConditionalBlock.C
  ../rtcompiler/ForBlock.C
  ../rtcompiler/Function.C
  ../rtcompiler/IfElseifElseBlock.C
  ../rtcompiler/Line.C
  ../rtcompiler/NormalBlock.C
  ../rtcompiler/Operator.C
  ../rtcompiler/Registrar.C
  ../rtcompiler/Tokenizer.C
  ../rtcompiler/WhileBlock.C
  ../rtcompiler/common.C
  )

#
# Extra Functions
#
APPEND_SET(PE_SOURCES
  ../extra_functions/pamgen_extras.cpp
  )


APPEND_SET(PE_HEADERS
  ../extra_functions/pamgen_extras.h
  )



#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  pamgen
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

TRIBITS_ADD_LIBRARY(
  pamgen_extras
  HEADERS ${PE_HEADERS}
  SOURCES ${PE_SOURCES}
  )

#
# D) Export the dependancy variables of this package for other packages to use
#

#PACKAGE_EXPORT_DEPENDENCY_VARIABLES()

if (WIN32)
  target_compile_definitions (pamgen PRIVATE _USE_MATH_DEFINES)
endif ()

