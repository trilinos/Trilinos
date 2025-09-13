if (${CMAKE_PROJECT_NAME} STREQUAL "Seacas")
  INCLUDE(InstallSymLink)
endif()

SET( ACCESSDIR ${CMAKE_INSTALL_PREFIX} )
SET( SEPARATOR "_")
INCLUDE(GNUInstallDirs)

ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_SEACASNemspread)
ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_SEACASNemslice)
IF (${PROJECT_NAME}_ENABLE_SEACASNemspread AND ${PROJECT_NAME}_ENABLE_SEACASNemslice)
  IF (${PROJECT_NAME}_ENABLE_Zoltan)
     SET( DECOMP_DEFAULT rib )
  ELSE()
     SET( DECOMP_DEFAULT inertial )
  ENDIF()

  CONFIGURE_FILE(
        ${CMAKE_CURRENT_SOURCE_DIR}/decomp.in
        ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/decomp
        @ONLY
        )
  FILE(COPY ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/decomp
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
       FILE_PERMISSIONS OWNER_READ OWNER_WRITE
         OWNER_EXECUTE GROUP_READ GROUP_EXECUTE)
  INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/decomp DESTINATION bin)
ENDIF()

ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_SEACASBlot)
IF (${PROJECT_NAME}_ENABLE_SEACASBlot)
  CONFIGURE_FILE(
        ${CMAKE_CURRENT_SOURCE_DIR}/blot.in
        ${CMAKE_CURRENT_BINARY_DIR}/blot
        @ONLY
        )
  INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/blot DESTINATION bin)
ENDIF()

ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_SEACASFastq)
IF (${PROJECT_NAME}_ENABLE_SEACASFastq)
  CONFIGURE_FILE(
        ${CMAKE_CURRENT_SOURCE_DIR}/fastq.in
        ${CMAKE_CURRENT_BINARY_DIR}/fastq
        @ONLY
        )
  INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/fastq DESTINATION bin)
ENDIF()

ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_SEACASEpu)
IF (${PROJECT_NAME}_ENABLE_SEACASEpu)
  CONFIGURE_FILE(
        ${CMAKE_CURRENT_SOURCE_DIR}/epup.in
        ${CMAKE_CURRENT_BINARY_DIR}/epup
        @ONLY
        )
  INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/epup DESTINATION bin)
ENDIF()

ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_SEACASConjoin)
IF (${PROJECT_NAME}_ENABLE_SEACASConjoin)
  CONFIGURE_FILE(
        ${CMAKE_CURRENT_SOURCE_DIR}/pconjoin.in
        ${CMAKE_CURRENT_BINARY_DIR}/pconjoin
        @ONLY
        )
  INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/pconjoin DESTINATION bin)
ENDIF()

TRIBITS_INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
)

SET(HEADERS "")

IF (NOT CMAKE_HOST_SYSTEM_NAME STREQUAL "Windows")
   APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/getopt.seacas.c)
   ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_SEACASBlot)
   ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_SEACASFastq)
   ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_SEACASEpu)
   IF (${PROJECT_NAME}_ENABLE_SEACASBlot OR
       ${PROJECT_NAME}_ENABLE_SEACASFastq OR
       ${PROJECT_NAME}_ENABLE_SEACASEpu)
      ADD_EXECUTABLE(getopt.seacas ${SOURCES} )
      install_executable(getopt.seacas)
   endif()
ENDIF()

ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_SEACASExodus)
ASSERT_DEFINED(BUILD_SHARED_LIBS)
IF (${PROJECT_NAME}_ENABLE_SEACASExodus)
   find_package(Python COMPONENTS Interpreter)
   IF (${Python_FOUND})
      message(STATUS "Found python version ${Python_VERSION}")
      IF (${Python_VERSION} VERSION_LESS "3.0")
         SET( EXODUSPY "exodus2.py" )
         SET( EXOMERGE "exomerge2.py")
         MESSAGE(STATUS "A Python-2 version of exodus.py and exomerge.py will be installed.")
      ELSE()
         SET( EXODUSPY "exodus3.py" )
         SET( EXOMERGE "exomerge3.py")
         # IF (BUILD_SHARED_LIBS OR SEACASExodus_ENABLE_SHARED)
         #    IF (NOT ${PROJECT_NAME}_HOSTNAME STREQUAL "APPVEYOR-VM")
         #       ADD_TEST(NAME test_exodus3.py COMMAND ${CMAKE_INSTALL_PREFIX}/lib/tests/test_exodus3.py)
         #       ADD_TEST(NAME test_exomerge3.py COMMAND ${CMAKE_INSTALL_PREFIX}/lib/tests/exomerge_unit_test.py)
         #    ENDIF()
         # ENDIF()
         MESSAGE(STATUS "A Python-3 version of exodus.py and exomerge.py will be installed.")
      ENDIF()

      IF (BUILD_SHARED_LIBS OR SEACASExodus_ENABLE_SHARED)
         IF (NOT "${CMAKE_INSTALL_PREFIX}" STREQUAL "")
            string(FIND ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR} ${CMAKE_INSTALL_PREFIX} POSITION)
            cmake_path(IS_ABSOLUTE ${PROJECT_NAME}_INSTALL_INCLUDE_DIR IS_PATH_ABSOLUTE)

            IF (${POSITION} EQUAL 0)
	       # The include install path is a under the install path.  
               file(RELATIVE_PATH SEACAS_INCLUDEDIR ${CMAKE_INSTALL_PREFIX} ${${PROJECT_NAME}_INSTALL_INCLUDE_DIR} )
            ELSEIF(IS_PATH_ABSOLUTE)
	       # The include install path is different than the install path
               SET( SEACAS_INCLUDEDIR "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}" )
            ELSE()
	       # If neither, just set to include.
               SET( SEACAS_INCLUDEDIR "include" )
            ENDIF()

            string(FIND ${${PROJECT_NAME}_INSTALL_LIB_DIR} ${CMAKE_INSTALL_PREFIX} POSITION)
            cmake_path(IS_ABSOLUTE ${PROJECT_NAME}_INSTALL_LIB_DIR IS_PATH_ABSOLUTE)
            IF (${POSITION} EQUAL 0)
	       # The libary path is a under the install path.  
               file(RELATIVE_PATH SEACAS_LIBDIR ${CMAKE_INSTALL_PREFIX} ${${PROJECT_NAME}_INSTALL_LIB_DIR} )
            ELSEIF(IS_PATH_ABSOLUTE)
	       # The libarary install path is different than the install path
               SET( SEACAS_LIBDIR "${${PROJECT_NAME}_INSTALL_LIB_DIR}" )
            ELSE()
               SET( SEACAS_LIBDIR "lib" )
            ENDIF()
         ENDIF()

         FILE(GLOB_RECURSE EXODUSIIH "${CMAKE_SOURCE_DIR}/packages/seacas/*exodusII.h")
         FILE(STRINGS "${EXODUSIIH}" EXODUS_VERSION REGEX "#define EXODUS_VERSION .*")

         CONFIGURE_FILE(
                   ${CMAKE_CURRENT_SOURCE_DIR}/exodus2.in.py
                   ${CMAKE_CURRENT_BINARY_DIR}/exodus2.py
                   @ONLY
                   )
         CONFIGURE_FILE(
                   ${CMAKE_CURRENT_SOURCE_DIR}/exodus3.in.py
                   ${CMAKE_CURRENT_BINARY_DIR}/exodus3.py
                   @ONLY
                   )
         CONFIGURE_FILE(
                   ${CMAKE_CURRENT_SOURCE_DIR}/tests/exomerge_unit_test.py
                   ${CMAKE_CURRENT_BINARY_DIR}/tests/exomerge_unit_test.py
                   @ONLY
                   )
        CONFIGURE_FILE(
                   ${CMAKE_CURRENT_SOURCE_DIR}/exomerge3.py
                   ${CMAKE_CURRENT_BINARY_DIR}/exomerge3.py
                   @ONLY
                   )
        CONFIGURE_FILE(
                   ${CMAKE_CURRENT_SOURCE_DIR}/exomerge2.py
                   ${CMAKE_CURRENT_BINARY_DIR}/exomerge2.py
                   @ONLY
                   )
         CONFIGURE_FILE(
                   ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_exodus3.py
                   ${CMAKE_CURRENT_BINARY_DIR}/tests/test_exodus3.py
                   @ONLY
                   )
         CONFIGURE_FILE(
                   ${CMAKE_CURRENT_SOURCE_DIR}/tests/test-assembly.exo
                   ${CMAKE_CURRENT_BINARY_DIR}/tests/test-assembly.exo
                   COPYONLY
                   )

         INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/exodus2.py DESTINATION lib)
         INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/exodus3.py DESTINATION lib)
         INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/tests/test_exodus3.py DESTINATION lib/tests/)
         INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/tests/exomerge_unit_test.py DESTINATION lib/tests/)
         INSTALL(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test-assembly.exo DESTINATION lib/tests/)
         INSTALL(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/tests/exomerge_unit_test.e DESTINATION lib/tests/)
         INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/exomerge2.py DESTINATION lib)
         INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/exomerge3.py DESTINATION lib)

         if (${CMAKE_PROJECT_NAME} STREQUAL "Seacas" AND NOT WIN32)
            InstallSymLink(${EXODUSPY} ${CMAKE_INSTALL_PREFIX}/lib/exodus.py)
            InstallSymLink(${EXOMERGE} ${CMAKE_INSTALL_PREFIX}/lib/exomerge.py)
         endif()
      ELSE()
         MESSAGE(STATUS "Exodus.py and exomerge.py NOT installed since shared exodus library is not available.")
      ENDIF()
   else()
      MESSAGE(STATUS "Exodus.py and exomerge.py NOT installed since Python could not be found.")
   endif()
ELSE()
   MESSAGE(STATUS "Exodus.py and exomerge.py NOT installed since Exodus is not enabled.")
endif()
