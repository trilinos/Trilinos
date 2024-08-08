include(CMakeFindDependencyMacro)
find_dependency(Intrepid2 REQUIRED)
find_dependency(MPI REQUIRED)
find_dependency(SEACAS REQUIRED)
find_dependency(sierra_common REQUIRED)
find_dependency(stk REQUIRED)
find_dependency(yaml-cpp REQUIRED)

include(${CMAKE_CURRENT_LIST_DIR}/createParserTarget.cmake)
create_parser_target(TARGET krino_commands_xmldb SOURCES krino_commands.xmldb)

include("${CMAKE_CURRENT_LIST_DIR}/krinoTargets.cmake")

