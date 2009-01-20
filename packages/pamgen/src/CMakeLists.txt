
INCLUDE(TrilinosPackageLibraryHelpers)

#
# A) Package-specific configuration options
#

TRILINOS_PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")

# Set headers (for IDEs and installing)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

# NOTE: Above, we add the XXX_config.h file so tht MS Visual Studio
# will know to include this file in the dependency analysis

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# Core Files
#

APPEND_SET(HEADERS
  Random.h
  StrLoopLimits.h
  Vector.h
  bc_specification.h
  brick_inline_mesh_desc.h
  calc_decomp_cuts.h
  element_density_function.h
  geometry_transform.h
  inline_geometries.h
  inline_mes_driver_lt.h
  inline_mesh_desc.h
  inline_mesh_driver_lt.h
  legacy_inline_mesh_desc.h
  parse_routines.h
  radial_inline_mesh_desc.h
  radial_trisection_inline_mesh_desc.h
  topology_enum.h
  uns_inline_decomp.h
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

APPEND_SET(HEADERS
  ../parser/InputBlock.h
  ../parser/keyword.h
  ../parser/parse_table.h
  ../parser/parser.h
  ../parser/token.h
  ../parser/token_enum.h
  ../parser/token_stream.h
  ../parser/token_value.h
  )

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
  ../mesh_spec_lt/im_exodusII.h
  ../mesh_spec_lt/im_exodusII_ext.h
  ../mesh_spec_lt/im_ne_nemesisI.h
  ../mesh_spec_lt/mesh_specification.h
  )

APPEND_SET(SOURCES
  ../mesh_spec_lt/im_ex_c_interface.C
  ../mesh_spec_lt/im_ne_c_interface.C
  ../mesh_spec_lt/mesh_specification.C
  )

#
# RTCOMPILER Files
#

APPEND_SET(HEADERS
  ../rtcompiler/ArrayIndexRTC.hh
  ../rtcompiler/ArrayNumberRTC.hh
  ../rtcompiler/ArrayVarRTC.hh
  ../rtcompiler/Bessel_I.hh
  ../rtcompiler/BlockRTC.hh
  ../rtcompiler/ConditionalBlockRTC.hh
  ../rtcompiler/ExecutableRTC.hh
  ../rtcompiler/ForBlockRTC.hh
  ../rtcompiler/FunctionRTC.hh
  ../rtcompiler/IfElseifElseBlockRTC.hh
  ../rtcompiler/LineRTC.hh
  ../rtcompiler/NormalBlockRTC.hh
  ../rtcompiler/ObjectRTC.hh
  ../rtcompiler/OperatorRTC.hh
  ../rtcompiler/RegistrarRTC.hh
  ../rtcompiler/ScalarNumberRTC.hh
  ../rtcompiler/ScalarVarRTC.hh
  ../rtcompiler/TokenizerRTC.hh
  ../rtcompiler/ValueRTC.hh
  ../rtcompiler/VariableRTC.hh
  ../rtcompiler/WhileBlockRTC.hh
  ../rtcompiler/commonRTC.hh
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
# ASRC Files
#

APPEND_SET(HEADERS
  ../asrc/code_types.h
  ../asrc/element_dictionary.h
  ../asrc/fudges.h
  )

APPEND_SET(SOURCES
  ../asrc/element_dictionary.C
  ../asrc/fudges.C
  )

#
# C) Define the targets for package's library(s)
#

TRILINOS_PACKAGE_ADD_LIBRARY(
  pamgen
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

#
# D) Export the dependancy variables of this package for other packages to use
#

TRILINOS_PACKAGE_EXPORT_DEPENDENCY_VARIABLES()
