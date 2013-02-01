# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#
#
# Copyright (2001) Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in this data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
#
# ************************************************************************
# @HEADER

INCLUDE(TribitsGeneralMacros)

#
#  This function will take a list and turn it into a space separated string
#  adding the prefix to the front of every entry. 
#
FUNCTION(TRIBITS_LIST_TO_STRING LIST PREFIX OUTPUT_STRING)
  SET(LIST_STRING "")
  
  FOREACH(ITEM ${LIST})
    SET(LIST_STRING "${LIST_STRING} ${PREFIX}${ITEM}")
  ENDFOREACH()

  SET(${OUTPUT_STRING} ${LIST_STRING} PARENT_SCOPE)
ENDFUNCTION()

#
#  This function will take a list of libraries and turn it into a space
#  separated string. In this case though the prefix is not always added
#  to the front of each entry as libraries can be specified either as a
#  name of a library to find or the absolute path to the library file
#  with any decorations the system uses. When an absolute path is given
#  the entry is used verbatim.
#
FUNCTION(TRIBITS_LIBRARY_LIST_TO_STRING LIST PREFIX OUTPUT_STRING)
  SET(LIST_STRING "")
  
  FOREACH(ITEM ${LIST})
    IF(EXISTS ${ITEM})
      SET(LIST_STRING "${LIST_STRING} ${ITEM}")
    ELSE()
      SET(LIST_STRING "${LIST_STRING} ${PREFIX}${ITEM}")
    ENDIF()
  ENDFOREACH()

  SET(${OUTPUT_STRING} ${LIST_STRING} PARENT_SCOPE)
ENDFUNCTION()

#
# CMAKE_CURRENT_LIST_DIR is not defined in CMake versions < 2.8.3, but the
# Trilinos writes paths that use the value of that variable to this file.
# Make sure it is available at *find_package* time. Note that all variable
# references in the code snippet are escaped. This is to keep them from
# being evaluated until they are actually in the install tree. This is
# done to handle movable install trees.
#
# This function defines the variable 
# DEFINE_CMAKE_CURRENT_LIST_DIR_CODE_CODE_SNIPPET in the caller's scope
# as a string that can be referenced from CONFIGURE_FILE input files
# to ensure that the CMAKE_CURRENT_LIST_DIR will be defined on the installation
# target machine, even if it has an older version of cmake.
#
FUNCTION(TRIBITS_SET_DEFINE_CMAKE_CURRENT_LIST_DIR_CODE_SNIPPET)
  SET(DEFINE_CMAKE_CURRENT_LIST_DIR_CODE_SNIPPET "
IF (NOT DEFINED CMAKE_CURRENT_LIST_DIR)
  GET_FILENAME_COMPONENT(_THIS_SCRIPT_PATH \${CMAKE_CURRENT_LIST_FILE} PATH)
  SET(CMAKE_CURRENT_LIST_DIR \${_THIS_SCRIPT_PATH})
ENDIF()
"
  PARENT_SCOPE )
ENDFUNCTION()

FUNCTION(TRIBITS_WRITE_PACKAGE_CONFIG_FILE PACKAGE_NAME)
  
  TRIBITS_SET_DEFINE_CMAKE_CURRENT_LIST_DIR_CODE_SNIPPET()

  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("For package ${PACKAGE_NAME} creating ${PACKAGE_NAME}Config.cmake")
  ENDIF()

  #Create the full package set only for those packages which are enabled
  SET(FULL_PACKAGE_SET)
  SET(FULL_LIBRARY_SET)
  FOREACH(TRIBITS_PACKAGE ${${PACKAGE_NAME}_FULL_EXPORT_DEP_PACKAGES})
    LIST(FIND ${PACKAGE_NAME}_LIB_OPTIONAL_DEP_PACKAGES ${TRIBITS_PACKAGE} PACKAGE_IS_OPTIONAL_DEP)
    IF(PACKAGE_IS_OPTIONAL_DEP GREATER -1)
      IF(${PACKAGE_NAME}_ENABLE_${TRIBITS_PACKAGE})
        LIST(APPEND FULL_PACKAGE_SET ${TRIBITS_PACKAGE})
        LIST(APPEND FULL_LIBRARY_SET ${${TRIBITS_PACKAGE}_LIBRARIES})
      ENDIF()
    ELSE()
      LIST(APPEND FULL_PACKAGE_SET ${TRIBITS_PACKAGE})
      LIST(APPEND FULL_LIBRARY_SET ${${TRIBITS_PACKAGE}_LIBRARIES})
    ENDIF()
  ENDFOREACH()

  SET(FULL_TPL_SET "")
  FOREACH(TRIBITS_PACKAGE ${FULL_PACKAGE_SET})
    LIST(APPEND FULL_TPL_SET ${${TRIBITS_PACKAGE}_LIB_REQUIRED_DEP_TPLS})

    SET(OPTIONAL_TPLS ${${TRIBITS_PACKAGE}_LIB_OPTIONAL_DEP_TPLS})

    FOREACH(TPL ${OPTIONAL_TPLS})
      IF(${TRIBITS_PACKAGE}_ENABLE_${TPL})
        LIST(APPEND FULL_TPL_SET ${TPL})
      ENDIF()
    ENDFOREACH()
  ENDFOREACH()

  #We will use the complete list of supported tpls for the project
  #to help us create a properly ordered list of tpls.
  IF(FULL_TPL_SET)
    # Reversing the tpl list so that the list of tpls will be produced in
    # order of most dependent to least dependent.
    SET(TPL_LIST ${${PROJECT_NAME}_TPLS})
    LIST(REVERSE TPL_LIST)
    
    SET(ORDERED_FULL_TPL_SET "")
    FOREACH(TPL ${TPL_LIST})
      LIST(FIND FULL_TPL_SET ${TPL} FOUND)
      
      IF(FOUND GREATER -1)
        LIST(APPEND ORDERED_FULL_TPL_SET ${TPL})
      ENDIF()
    ENDFOREACH()
  ENDIF()
  
  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("ORDERED_FULL_TPL_SET after all tpls = ${ORDERED_FULL_TPL_SET}")
  ENDIF()
  
  SET(${PACKAGE_NAME}_TPL_LIBRARIES "")
  SET(${PACKAGE_NAME}_TPL_INCLUDE_DIRS "")
  SET(${PACKAGE_NAME}_TPL_LIBRARY_DIRS "")
  FOREACH(TPL ${ORDERED_FULL_TPL_SET})
    LIST(APPEND ${PACKAGE_NAME}_TPL_LIBRARIES ${TPL_${TPL}_LIBRARIES})
    LIST(APPEND ${PACKAGE_NAME}_TPL_INCLUDE_DIRS ${TPL_${TPL}_INCLUDE_DIRS})
    LIST(APPEND ${PACKAGE_NAME}_TPL_LIBRARY_DIRS ${TPL_${TPL}_LIBRARY_DIRS})
  ENDFOREACH()
  
  # Generate a note discouraging editing of the <package>Config.cmake file
  SET(DISCOURAGE_EDITING "Do not edit: This file was generated automatically by CMake.")
  
  ######
  # Create a configure file for the build tree. Creating this file in the base dir
  # of the package since it is possible that the cmake returned path where the file
  # was found would be useful for a package, and having to dig through the hiding that
  # is done wouldn't be nice.
  ######

  SET(LIBRARY_DIRS ${${PACKAGE_NAME}_LIBRARY_DIRS})
  SET(INCLUDE_DIRS ${${PACKAGE_NAME}_INCLUDE_DIRS})

  # Custom code in configuration file.
  SET(PACKAGE_CONFIG_CODE "")

  # Import build tree targets into applications.
  IF(FULL_LIBRARY_SET)
    SET(PACKAGE_CONFIG_CODE "${PACKAGE_CONFIG_CODE}
# Import ${PROJECT_NAME} targets
IF(NOT ${PROJECT_NAME}_TARGETS_IMPORTED)
  SET(${PROJECT_NAME}_TARGETS_IMPORTED 1)
  INCLUDE(\"${PROJECT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake\")
ENDIF()
")
  ENDIF()

  # Write the specification of the rpath if necessary. This is only needed if we're building shared libraries. 
  IF(BUILD_SHARED_LIBS)
    STRING(REPLACE ";" ":" SHARED_LIB_RPATH_COMMAND "${LIBRARY_DIRS}")
    SET(SHARED_LIB_RPATH_COMMAND ${CMAKE_SHARED_LIBRARY_RUNTIME_CXX_FLAG}${SHARED_LIB_RPATH_COMMAND})
  ENDIF()

  CONFIGURE_FILE(
    ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsPackageConfigTemplate.cmake.in 
    ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}Config.cmake
    )

  IF(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES)
    ######
    # Create a Makefile.export.<package_name> for the build tree. This is the equivalent
    # of the cmake version only slightly changed so that it can be directly imported into
    # a Makefile.
    ######

    TRIBITS_LIST_TO_STRING("${FULL_LIBRARY_SET}" ${CMAKE_LINK_LIBRARY_FLAG} MAKEFILE_FULL_LIBRARY_SET)
    TRIBITS_LIST_TO_STRING("${LIBRARY_DIRS}" ${CMAKE_LIBRARY_PATH_FLAG} MAKEFILE_LIBRARY_DIRS)
    TRIBITS_LIST_TO_STRING("${INCLUDE_DIRS}" "-I" MAKEFILE_INCLUDE_DIRS)
    TRIBITS_LIST_TO_STRING("${${PACKAGE_NAME}_TPL_INCLUDE_DIRS}" "-I" MAKEFILE_${PACKAGE_NAME}_TPL_INCLUDE_DIRS)
    TRIBITS_LIST_TO_STRING("${${PACKAGE_NAME}_TPL_LIBRARY_DIRS}" ${CMAKE_LIBRARY_PATH_FLAG} MAKEFILE_${PACKAGE_NAME}_TPL_LIBRARY_DIRS)
    #the TPL library names have to be treated differently
    TRIBITS_LIBRARY_LIST_TO_STRING("${${PACKAGE_NAME}_TPL_LIBRARIES}" ${CMAKE_LINK_LIBRARY_FLAG} MAKEFILE_${PACKAGE_NAME}_TPL_LIBRARIES)

    TRIBITS_LIBRARY_LIST_TO_STRING("${${TPL_MPI_LIBRARIES}}" ${CMAKE_LINK_LIBRARY_FLAG} "MAKEFILE_TPL_MPI_LIBRARIES")
    TRIBITS_LIST_TO_STRING("${${TPL_MPI_LIBRARY_DIRS}}" ${CMAKE_LIBRARY_PATH_FLAG} "MAKEFILE_TPL_MPI_LIBRARY_DIRS")
    TRIBITS_LIST_TO_STRING("${${TPL_MPI_INCLUDE_DIRS}}" "-I" "MAKEFILE_TPL_MPI_INCLUDE_DIRS")
    
    TRIBITS_LIST_TO_STRING("${FULL_PACKAGE_SET}" "" MAKEFILE_FULL_PACKAGE_SET)
    TRIBITS_LIST_TO_STRING("${ORDERED_FULL_TPL_SET}" "" MAKEFILE_ORDERED_FULL_TPL_SET)

    #create an upper case name of the package so that we can make deprecated versions of them to help people
    #transistioning from the autotools version diagnose any missed variables.
    STRING(TOUPPER ${PACKAGE_NAME} PACKAGE_NAME_UPPER)

    CONFIGURE_FILE(
      ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsPackageConfigTemplate.export.in 
      ${CMAKE_CURRENT_BINARY_DIR}/Makefile.export.${PACKAGE_NAME}
      )
  ENDIF()

  ######
  # Create a configure file for the install tree and set the install target for it. This
  # file isn't generally useful inside the build tree so it is being "hidden" in the 
  # CMakeFiles directory. It will be placed in the base install directory for ${PROJECT_NAME}
  # when installed.
  ######

  # Set the include and library directories relative to the location
  # at which the ${PROJECT_NAME}Config.cmake file is going to be
  # installed. Note the variable reference below is escaped so it
  # won't be replaced until a client project attempts to locate
  # directories using the installed config file. This is to deal with
  # installers that allow relocation of the install tree at *install*
  # time.
  SET(LIBRARY_DIRS "\${CMAKE_CURRENT_LIST_DIR}/../../../${${PROJECT_NAME}_INSTALL_LIB_DIR}")
  SET(INCLUDE_DIRS "\${CMAKE_CURRENT_LIST_DIR}/../../../${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}")

  # Custom code in configuration file.
  SET(PACKAGE_CONFIG_CODE "")

  # Import install tree targets into applications.
  GET_PROPERTY(HAS_INSTALL_TARGETS GLOBAL PROPERTY ${PACKAGE_NAME}_HAS_INSTALL_TARGETS)
  IF(HAS_INSTALL_TARGETS)
    SET(PACKAGE_CONFIG_CODE "${PACKAGE_CONFIG_CODE}
# Import ${PROJECT_NAME} targets
IF(NOT ${PROJECT_NAME}_TARGETS_IMPORTED)
  SET(${PROJECT_NAME}_TARGETS_IMPORTED 1)
  INCLUDE(\"\${CMAKE_CURRENT_LIST_DIR}/../${PROJECT_NAME}/${PROJECT_NAME}Targets.cmake\")
ENDIF()
")
  ENDIF()

  # Write the specification of the rpath if necessary. This is only needed if we're building shared libraries. 
  IF(BUILD_SHARED_LIBS)
    SET(SHARED_LIB_RPATH_COMMAND ${CMAKE_SHARED_LIBRARY_RUNTIME_CXX_FLAG}${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR})
  ENDIF()

  CONFIGURE_FILE(
    ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsPackageConfigTemplate.cmake.in 
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${PACKAGE_NAME}Config_install.cmake
    )

  INSTALL(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${PACKAGE_NAME}Config_install.cmake
    DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake/${PACKAGE_NAME}"
    RENAME ${PACKAGE_NAME}Config.cmake
  )

  IF(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES)
    ######
    # Create a Makefile.export.<package_name> for the install tree. This is the equivalent
    # of the cmake version only slightly changed so that it can be directly imported into
    # a Makefile.
    ######

    # Generated Make imports must use CMAKE_INSTALL_PREFIX, rather
    # than the more platform friendly method of locating the libraries
    # and includes using the config file path above. The underlying
    # assumption here is that a generator that uses
    # CMAKE_INSTALL_PREFIX is being used.
    SET(LIBRARY_DIRS ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR})
    SET(INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

    TRIBITS_LIST_TO_STRING("${LIBRARY_DIRS}" ${CMAKE_LIBRARY_PATH_FLAG} MAKEFILE_LIBRARY_DIRS)
    TRIBITS_LIST_TO_STRING("${INCLUDE_DIRS}" "-I" MAKEFILE_INCLUDE_DIRS)

    CONFIGURE_FILE(
      ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsPackageConfigTemplate.export.in 
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Makefile.export.${PACKAGE_NAME}_install
      )

    INSTALL(
      FILES ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Makefile.export.${PACKAGE_NAME}_install
      DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}"
      RENAME Makefile.export.${PACKAGE_NAME}
    )
  ENDIF()
ENDFUNCTION()

FUNCTION(TRIBITS_WRITE_CONFIG_FILE)

  TRIBITS_SET_DEFINE_CMAKE_CURRENT_LIST_DIR_CODE_SNIPPET()

  # Reversing the package list so that libraries will be produced in order of
  # most dependent to least dependent.
  SET(PACKAGE_LIST ${${PROJECT_NAME}_SE_PACKAGES})
  LIST(REVERSE PACKAGE_LIST)
  

  # Loop over all packages to determine which were enabled. Then build a list
  # of all their libraries/includes in the proper order for linking
  SET(FULL_PACKAGE_SET "")
  SET(FULL_LIBRARY_SET "")
  SET(FULL_INCLUDE_DIRS_SET "")
  SET(FULL_LIBRARY_DIRS_SET "")
  FOREACH(TRIBITS_PACKAGE ${PACKAGE_LIST})
    IF(${PROJECT_NAME}_ENABLE_${TRIBITS_PACKAGE})
      LIST(APPEND FULL_PACKAGE_SET ${TRIBITS_PACKAGE})
      LIST(APPEND FULL_LIBRARY_SET ${${TRIBITS_PACKAGE}_LIBRARIES})
      LIST(APPEND FULL_INCLUDE_DIRS_SET ${${TRIBITS_PACKAGE}_INCLUDE_DIRS})
      LIST(APPEND FULL_LIBRARY_DIRS_SET ${${TRIBITS_PACKAGE}_LIBRARY_DIRS})
    ENDIF()
  ENDFOREACH()
  
  SET(${PROJECT_NAME}_CONFIG_LIBRARIES ${FULL_LIBRARY_SET})
  
  # Reversing the tpl list so that the list of tpls will be produced in
  # order of most dependent to least dependent.
  SET(TPL_LIST ${${PROJECT_NAME}_TPLS})
  LIST(REVERSE TPL_LIST)
  
  # Loop over all TPLs to determine which were enabled. Then build a list
  # of all their libraries/includes in the proper order for linking
  SET(FULL_TPL_SET "")
  SET(FULL_TPL_LIBRARY_SET "")
  SET(FULL_TPL_INCLUDE_DIRS_SET "")
  SET(FULL_TPL_LIBRARY_DIRS_SET "")
  FOREACH(TPL ${TPL_LIST})
    IF(TPL_ENABLE_${TPL})
      LIST(APPEND FULL_TPL_SET ${TPL})
      LIST(APPEND FULL_TPL_LIBRARY_SET ${TPL_${TPL}_LIBRARIES})
      LIST(APPEND FULL_TPL_INCLUDE_DIRS_SET ${TPL_${TPL}_INCLUDE_DIRS})
      LIST(APPEND FULL_TPL_LIBRARY_DIRS_SET ${TPL_${TPL}_LIBRARY_DIRS})
    ENDIF()
  ENDFOREACH()
  
  # it is possible that tpls are in the same directory, to keep from 
  # having a very long include path or library path we will strip out
  # any duplicates. This shouldn't affect which include or library is
  # found since the first instance of any path will be the one that is
  # kept.
  LIST(REMOVE_DUPLICATES FULL_TPL_INCLUDE_DIRS_SET)
  LIST(REMOVE_DUPLICATES FULL_TPL_LIBRARY_DIRS_SET)
  
  SET(${PROJECT_NAME}_CONFIG_TPL_INCLUDE_DIRS ${FULL_TPL_INCLUDE_DIRS_SET})
  SET(${PROJECT_NAME}_CONFIG_TPL_LIBRARY_DIRS ${FULL_TPL_LIBRARY_DIRS_SET})
  SET(${PROJECT_NAME}_CONFIG_TPL_LIBRARIES ${FULL_TPL_LIBRARY_SET})
  
  #
  # Configure two files for finding ${PROJECT_NAME}. One for the build tree and one for installing
  #

  # Generate a note discouraging editing of the <package>Config.cmake file
  SET(DISCOURAGE_EDITING "Do not edit: This file was generated automatically by CMake.")
  
  #Config file for setting variables and finding include/library paths from the build directory
  SET(${PROJECT_NAME}_CONFIG_INCLUDE_DIRS ${FULL_INCLUDE_DIRS_SET})
  SET(${PROJECT_NAME}_CONFIG_LIBRARY_DIRS ${FULL_LIBRARY_DIRS_SET})

  # Write the specification of the rpath if necessary. This is only needed if we're building shared libraries. 
  IF(BUILD_SHARED_LIBS)
    STRING(REPLACE ";" ":" SHARED_LIB_RPATH_COMMAND "${${PROJECT_NAME}_CONFIG_LIBRARY_DIRS}")
    SET(SHARED_LIB_RPATH_COMMAND ${CMAKE_SHARED_LIBRARY_RUNTIME_CXX_FLAG}${SHARED_LIB_RPATH_COMMAND})
  ENDIF()

  # Custom code in configuration file.
  SET(PROJECT_CONFIG_CODE "")

  # Export targets from the build tree.
  IF(FULL_LIBRARY_SET)
    LIST(SORT FULL_LIBRARY_SET)
    LIST(REMOVE_DUPLICATES FULL_LIBRARY_SET)
    EXPORT(TARGETS ${FULL_LIBRARY_SET} FILE "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake")
    # Import the targets in applications.
    SET(PROJECT_CONFIG_CODE "${PROJECT_CONFIG_CODE}
# Import ${PROJECT_NAME} targets
IF(NOT ${PROJECT_NAME}_TARGETS_IMPORTED)
  SET(${PROJECT_NAME}_TARGETS_IMPORTED 1)
  INCLUDE(\"${PROJECT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake\")
ENDIF()
")
  ENDIF()

  # Appending the logic to include each package's config file.
  SET(LOAD_CODE "# Load configurations from enabled packages\n")
  FOREACH(TRIBITS_PACKAGE ${FULL_PACKAGE_SET})
    SET(LOAD_CODE "${LOAD_CODE}include(\"${${TRIBITS_PACKAGE}_BINARY_DIR}/${TRIBITS_PACKAGE}Config.cmake\")\n")
  ENDFOREACH()
  SET(PROJECT_CONFIG_CODE "${PROJECT_CONFIG_CODE}\n${LOAD_CODE}")

  CONFIGURE_FILE( 
    ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsProjectConfigTemplate.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake )

  IF(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES)
    ######
    # Create a Makefile.export.<project_name> for the build tree. This is the equivalent
    # of the cmake version only slightly changed so that it can be directly imported into
    # a Makefile.
    ######

    TRIBITS_LIST_TO_STRING("${${PROJECT_NAME}_CONFIG_LIBRARIES}" ${CMAKE_LINK_LIBRARY_FLAG} MAKEFILE_${PROJECT_NAME}_CONFIG_LIBRARIES)
    TRIBITS_LIST_TO_STRING("${${PROJECT_NAME}_CONFIG_LIBRARY_DIRS}" ${CMAKE_LIBRARY_PATH_FLAG} MAKEFILE_${PROJECT_NAME}_CONFIG_LIBRARY_DIRS)
    TRIBITS_LIST_TO_STRING("${${PROJECT_NAME}_CONFIG_INCLUDE_DIRS}" "-I" MAKEFILE_${PROJECT_NAME}_CONFIG_INCLUDE_DIRS)
    TRIBITS_LIST_TO_STRING("${${PROJECT_NAME}_CONFIG_TPL_INCLUDE_DIRS}" "-I" MAKEFILE_${PROJECT_NAME}_CONFIG_TPL_INCLUDE_DIRS)
    TRIBITS_LIST_TO_STRING("${${PROJECT_NAME}_CONFIG_TPL_LIBRARY_DIRS}" ${CMAKE_LIBRARY_PATH_FLAG} MAKEFILE_${PROJECT_NAME}_CONFIG_TPL_LIBRARY_DIRS)
    #the TPL library names have to be treated differently
    TRIBITS_LIBRARY_LIST_TO_STRING("${${PROJECT_NAME}_CONFIG_TPL_LIBRARIES}" ${CMAKE_LINK_LIBRARY_FLAG} MAKEFILE_${PROJECT_NAME}_CONFIG_TPL_LIBRARIES)

    TRIBITS_LIBRARY_LIST_TO_STRING("${${TPL_MPI_LIBRARIES}}" ${CMAKE_LINK_LIBRARY_FLAG} "MAKEFILE_TPL_MPI_LIBRARIES")
    TRIBITS_LIST_TO_STRING("${${TPL_MPI_LIBRARY_DIRS}}" ${CMAKE_LIBRARY_PATH_FLAG} "MAKEFILE_TPL_MPI_LIBRARY_DIRS")
    TRIBITS_LIST_TO_STRING("${${TPL_MPI_INCLUDE_DIRS}}" "-I" "MAKEFILE_TPL_MPI_INCLUDE_DIRS")
    
    TRIBITS_LIST_TO_STRING("${FULL_PACKAGE_SET}" "" MAKEFILE_FULL_PACKAGE_SET)
    TRIBITS_LIST_TO_STRING("${FULL_TPL_SET}" "" MAKEFILE_FULL_TPL_SET)

    CONFIGURE_FILE( 
      ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsProjectConfigTemplate.export.in
      ${CMAKE_CURRENT_BINARY_DIR}/Makefile.export.${PROJECT_NAME})
  ENDIF()  

  ######
  # Create a configure file for the install tree and set the install target for it. This
  # file isn't generally useful inside the build tree. It will be placed in the base
  # install directory for ${PROJECT_NAME} when installed.
  ######
  
  # Set the include and library directories relative to the location
  # at which the ${PROJECT_NAME}Config.cmake file is going to be
  # installed. Note the variable reference below is escaped so it
  # won't be replaced until a client project attempts to locate
  # directories using the installed config file. This is to deal with
  # installers that allow relocation of the install tree at *install*
  # time.
  SET(${PROJECT_NAME}_CONFIG_INCLUDE_DIRS "\${CMAKE_CURRENT_LIST_DIR}/../../../${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}")
  SET(${PROJECT_NAME}_CONFIG_LIBRARY_DIRS "\${CMAKE_CURRENT_LIST_DIR}/../../../${${PROJECT_NAME}_INSTALL_LIB_DIR}")

  # Write the specification of the rpath if necessary. This is only needed if we're building shared libraries. 
  IF(BUILD_SHARED_LIBS)
    SET(SHARED_LIB_RPATH_COMMAND ${CMAKE_SHARED_LIBRARY_RUNTIME_CXX_FLAG}${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR})
  ENDIF()

  # Custom code in configuration file.
  SET(PROJECT_CONFIG_CODE "")

  # Export targets from the install tree.
  GET_PROPERTY(HAS_INSTALL_TARGETS GLOBAL PROPERTY ${PROJECT_NAME}_HAS_INSTALL_TARGETS)
  IF(HAS_INSTALL_TARGETS)
    INSTALL(
      EXPORT ${PROJECT_NAME}
      DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake/${PROJECT_NAME}"
      FILE ${PROJECT_NAME}Targets.cmake
      )
    # Import the targets in applications.
    SET(PROJECT_CONFIG_CODE "${PROJECT_CONFIG_CODE}
# Import ${PROJECT_NAME} targets
IF(NOT ${PROJECT_NAME}_TARGETS_IMPORTED)
  SET(${PROJECT_NAME}_TARGETS_IMPORTED 1)
  INCLUDE(\"\${CMAKE_CURRENT_LIST_DIR}/${PROJECT_NAME}Targets.cmake\")
ENDIF()
")
  ENDIF()

  # Appending the logic to include each package's config file.
  SET(LOAD_CODE "# Load configurations from enabled packages\n")
  FOREACH(TRIBITS_PACKAGE ${FULL_PACKAGE_SET})
    SET(LOAD_CODE "${LOAD_CODE}include(\"\${CMAKE_CURRENT_LIST_DIR}/../${TRIBITS_PACKAGE}/${TRIBITS_PACKAGE}Config.cmake\")\n")
  ENDFOREACH()
  SET(PROJECT_CONFIG_CODE "${PROJECT_CONFIG_CODE}\n${LOAD_CODE}")

  CONFIGURE_FILE( 
    ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsProjectConfigTemplate.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config_install.cmake )

  INSTALL(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config_install.cmake
    DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake/${PROJECT_NAME}"
    RENAME ${PROJECT_NAME}Config.cmake
  )
  
  IF(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES)
    ######
    # Create a Makefile.export.<project_name> for the install tree. This is the equivalent
    # of the cmake version only slightly changed so that it can be directly imported into
    # a Makefile.
    ######

    # Generated Make imports must use CMAKE_INSTALL_PREFIX, rather
    # than the more platform friendly method of locating the libraries
    # and includes using the config file path above. The underlying
    # assumption here is that a generator that uses
    # CMAKE_INSTALL_PREFIX is being used.
    SET(${PROJECT_NAME}_CONFIG_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})
    SET(${PROJECT_NAME}_CONFIG_LIBRARY_DIRS ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_LIB_DIR})

    TRIBITS_LIST_TO_STRING("${${PROJECT_NAME}_CONFIG_LIBRARY_DIRS}" ${CMAKE_LIBRARY_PATH_FLAG} MAKEFILE_${PROJECT_NAME}_CONFIG_LIBRARY_DIRS)
    TRIBITS_LIST_TO_STRING("${${PROJECT_NAME}_CONFIG_INCLUDE_DIRS}" "-I" MAKEFILE_${PROJECT_NAME}_CONFIG_INCLUDE_DIRS)

    CONFIGURE_FILE( 
      ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsProjectConfigTemplate.export.in
      ${CMAKE_CURRENT_BINARY_DIR}/Makefile.export.${PROJECT_NAME}_install )

    INSTALL(
      FILES ${CMAKE_CURRENT_BINARY_DIR}/Makefile.export.${PROJECT_NAME}_install
      DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}"
      RENAME Makefile.export.${PROJECT_NAME}
    )
  ENDIF()
  
  #
  # Configure the version file for ${PROJECT_NAME}
  #
  
  CONFIGURE_FILE(
    ${${PROJECT_NAME}_TRIBITS_DIR}/${TRIBITS_CMAKE_INSTALLATION_FILES_DIR}/TribitsConfigVersionTemplate.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake @ONLY)

  INSTALL(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake/${PROJECT_NAME}"
  )

ENDFUNCTION()

