#
# B) Define the header and source files (and include directories)
#

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS
  "${DIR}/Teuchos_FiniteAutomaton.hpp"
  "${DIR}/Teuchos_Grammar.hpp"
  "${DIR}/Teuchos_Graph.hpp"
  "${DIR}/Teuchos_Language.hpp"
  "${DIR}/Teuchos_Parser.hpp"
  "${DIR}/Teuchos_Reader.hpp"
  "${DIR}/Teuchos_ReaderTables.hpp"
  "${DIR}/Teuchos_regex.hpp"
  "${DIR}/Teuchos_string.hpp"
  "${DIR}/Teuchos_TableDecl.hpp"
  "${DIR}/Teuchos_XML.hpp"
  "${DIR}/Teuchos_YAML.hpp"
  "${DIR}/Teuchos_MathExpr.hpp"
  )

SET(SOURCES
  "${DIR}/Teuchos_make_lalr1_parser.cpp"
  "${DIR}/Teuchos_chartab.cpp"
  "${DIR}/Teuchos_FiniteAutomaton.cpp"
  "${DIR}/Teuchos_Grammar.cpp"
  "${DIR}/Teuchos_Graph.cpp"
  "${DIR}/Teuchos_Language.cpp"
  "${DIR}/Teuchos_Parser.cpp"
  "${DIR}/Teuchos_Reader.cpp"
  "${DIR}/Teuchos_regex.cpp"
  "${DIR}/Teuchos_XML.cpp"
  "${DIR}/Teuchos_YAML.cpp"
  "${DIR}/Teuchos_MathExpr.cpp"
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  teuchosparser
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEFINES -DTEUCHOSPARSER_LIB_EXPORTS_MODE
  )

