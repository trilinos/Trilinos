
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  XMLParser_test  
  SOURCES XMLParser_test.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  XMLParser_ordered  
  SOURCES XMLParser_ordered.cpp
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

FILE(GLOB  BAD_XML_FILES RELATIVE  ${CMAKE_CURRENT_SOURCE_DIR} "b_*.xml")
FILE(GLOB GOOD_XML_FILES RELATIVE  ${CMAKE_CURRENT_SOURCE_DIR} "g_*.xml")

TRIBITS_ADD_EXECUTABLE(
  XMLParser_ParseFileForError
  SOURCES XMLParser_ParseFileForError.cpp
          ${TEUCHOS_STD_UNIT_TEST_MAIN}
)

TRIBITS_ADD_EXECUTABLE(
  XMLParser_ParseWithoutFailing
  SOURCES XMLParser_ParseWithoutFailing
          ${TEUCHOS_STD_UNIT_TEST_MAIN}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(BadXmlTestFiles
  SOURCE_FILES ${BAD_XML_FILES} ${GOOD_XML_FILES}
  EXEDEPS XMLParser_ParseFileForError
)

FOREACH (B ${BAD_XML_FILES})
  TRIBITS_ADD_TEST(
    XMLParser_ParseFileForError
    NAME XMLParse_BadParse_${B}
    ARGS "--filename=${B}"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
ENDFOREACH()

FOREACH (G ${GOOD_XML_FILES})
  TRIBITS_ADD_TEST(
    XMLParser_ParseWithoutFailing
    NAME XMLParse_GoodParse_${G}
    ARGS "--filename=${G}"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
ENDFOREACH()
