# - Try to find PETSc
# Once done this will define
#
#  PETSC_FOUND             - system has PETSc
#  TPL_PETSC_INCLUDE_DIRS  - the PETSc include directories
#  TPL_PETSC_LIBRARIES     - Link these to use PETSc
#  PETSC_COMPILER          - Compiler used by PETSc, helpful to find a compatible MPI
#  PETSC_DEFINITIONS       - Compiler switches for using PETSc
#  PETSC_MPIEXEC           - Executable for running MPI programs
#  PETSC_VERSION           - Version string (MAJOR.MINOR.SUBMINOR)
#
#  Usage:
#  find_package(PETSc COMPONENTS CXX)  - required if build --with-clanguage=C++ --with-c-support=0
#  find_package(PETSc COMPONENTS C)    - standard behavior of checking build using a C compiler
#  find_package(PETSc)                 - same as above
#
# Setting these changes the behavior of the search
#  PETSC_DIR - directory in which PETSc resides
#  PETSC_ARCH - build architecture
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

#
# Copyright:
#
#   Constantine Khroulev <ckhroulev@alaska.edu>
#   Jed Brown <jed@59A2.org>
#   johnfettig <john.fettig@gmail.com>
#
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
# 
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# 
# * Redistributions in binary form must reproduce the above copyright notice, this
#   list of conditions and the following disclaimer in the documentation and/or
#   other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#


# Allow cmake cache varaibles to set env vars.

ASSERT_DEFINED(${PROJECT_NAME}_VERBOSE_CONFIGURE)

IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
  PRINT_VAR(PETSC_DIR)
  PRINT_VAR(ENV{PETSC_DIR})
  PRINT_VAR(PETSC_ARCH)
  PRINT_VAR(ENV{PETSC_ARCH})
ENDIF()

IF (PETSC_DIR AND NOT ENV{PETSC_DIR})
  SET(ENV{PETSC_DIR} PETSC_DIR)
  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("-- " "Setting ENV{PETSC_DIR} from PETSC_DIR=${PETSC_DIR}") 
  ENDIF()
ENDIF()

IF (PETSC_ARCH AND NOT ENV{PETSC_ARCH})
  SET(ENV{PETSC_ARCH} PETSC_ARCH)
  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("-- " "Setting ENV{PETSC_ARCH} from PETSC_ARCH=${PETSC_ARCH}") 
  ENDIF()
ENDIF()


set(PETSC_VALID_COMPONENTS
  C
  CXX)
  
if(NOT PETSc_FIND_COMPONENTS)
  set(PETSC_LANGUAGE_BINDINGS "C")
else()
  # Right now, this is designed for compatability with the --with-clanguage option, so
  # only allow one item in the components list.
  list(LENGTH ${PETSc_FIND_COMPONENTS} components_length)
  if(${components_length} GREATER 1)
    message(FATAL_ERROR "Only one component for PETSc is allowed to be specified")
  endif()
  # This is a stub for allowing multiple components should that time ever come. Perhaps
  # to also test Fortran bindings?
  foreach(component ${PETSc_FIND_COMPONENTS})
    list(FIND PETSC_VALID_COMPONENTS ${component} component_location)
    if(${component_location} EQUAL -1)
      message(FATAL_ERROR "\"${component}\" is not a valid PETSc component.")
    else()
      list(APPEND PETSC_LANGUAGE_BINDINGS ${component})
    endif()
  endforeach()
endif()

function (petsc_get_version)
  if (EXISTS "${PETSC_DIR}/include/petscversion.h")
    file (STRINGS "${PETSC_DIR}/include/petscversion.h" vstrings REGEX "#define PETSC_VERSION_(RELEASE|MAJOR|MINOR|SUBMINOR|PATCH) ")
    foreach (line ${vstrings})
      string (REGEX REPLACE " +" ";" fields ${line}) # break line into three fields (the first is always "#define")
      list (GET fields 1 var)
      list (GET fields 2 val)
      set (${var} ${val} PARENT_SCOPE)
      set (${var} ${val})         # Also in local scope so we have access below
    endforeach ()
    if (PETSC_VERSION_RELEASE)
      set (PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}p${PETSC_VERSION_PATCH}" PARENT_SCOPE)
    else ()
      # make dev version compare higher than any patch level of a released version
      set (PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}.99" PARENT_SCOPE)
    endif ()
  else ()
    message (SEND_ERROR "PETSC_DIR can not be used, ${PETSC_DIR}/include/petscversion.h does not exist")
  endif ()
endfunction ()

find_path (PETSC_DIR include/petsc.h
  HINTS ENV PETSC_DIR
  PATHS
  /usr/lib/petscdir/3.1 /usr/lib/petscdir/3.0.0 /usr/lib/petscdir/2.3.3 /usr/lib/petscdir/2.3.2 # Debian
  $ENV{HOME}/petsc
  DOC "PETSc Directory")

find_program (MAKE_EXECUTABLE NAMES make gmake)

if (PETSC_DIR AND NOT PETSC_ARCH)
  set (_petsc_arches
    $ENV{PETSC_ARCH}                   # If set, use environment variable first
    linux-gnu463-omp16-release linux-gnu463-mpi16-release # sunspear arches
    linux-gnu463-omp16-debug linux-gnu463-omp16-debug
    linux-intel121-mpi16-release
    linux-intel121-mpi16-debug
    linux-gnu-c-debug linux-gnu-c-opt  # Debian defaults
    x86_64-unknown-linux-gnu i386-unknown-linux-gnu)
  set (petscconf "NOTFOUND" CACHE FILEPATH "Cleared" FORCE)
  foreach (arch ${_petsc_arches})
    if (NOT PETSC_ARCH)
      find_path (petscconf petscconf.h
	HINTS ${PETSC_DIR}
	PATH_SUFFIXES ${arch}/include bmake/${arch}
	NO_DEFAULT_PATH)
      if (petscconf)
	set (PETSC_ARCH "${arch}" CACHE STRING "PETSc build architecture")
      endif (petscconf)
    endif (NOT PETSC_ARCH)
  endforeach (arch)
  set (petscconf "NOTFOUND" CACHE INTERNAL "Scratch variable" FORCE)
endif (PETSC_DIR AND NOT PETSC_ARCH)

set (petsc_slaves LIBRARIES_SYS LIBRARIES_VEC LIBRARIES_MAT LIBRARIES_DM LIBRARIES_KSP LIBRARIES_SNES LIBRARIES_TS
  INCLUDE_DIR INCLUDE_CONF)
include (${${PROJECT_NAME}_TRIBITS_DIR}/tpls/FindPackageMultipass.cmake)
find_package_multipass (PETSc petsc_config_current
  STATES DIR ARCH
  DEPENDENTS INCLUDES LIBRARIES COMPILER MPIEXEC ${petsc_slaves})

# Determine whether the PETSc layout is old-style (through 2.3.3) or
# new-style (>= 3.0.0)
IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
  MESSAGE("-- " "Determing PETSC layout for different versions for PETSC_DIR=${PETSC_DIR} ...") 
ENDIF()
if (EXISTS "${PETSC_DIR}/include/petscconf.h")   # > 2.3.3 (sunspear)
  set (petsc_conf_rules "${PETSC_DIR}/conf/rules")
  set (petsc_conf_variables "${PETSC_DIR}/conf/variables")
  set (petsc_conf_petscvariables "${PETSC_DIR}/conf/petscvariables")
elseif (EXISTS "${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h")   # > 2.3.3
  set (petsc_conf_rules "${PETSC_DIR}/conf/rules")
  set (petsc_conf_variables "${PETSC_DIR}/conf/variables")
  set (petsc_conf_petscvariables "${PETSC_DIR}/${PETSC_ARCH}/conf/petscvariables")
elseif (EXISTS "${PETSC_DIR}/bmake/${PETSC_ARCH}/petscconf.h") # <= 2.3.3
  set (petsc_conf_rules "${PETSC_DIR}/bmake/common/rules")
  set (petsc_conf_variables "${PETSC_DIR}/bmake/common/variables")
elseif (PETSC_DIR)
  message (SEND_ERROR "The pair PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} do not specify a valid PETSc installation")
endif ()

if (petsc_conf_rules AND petsc_conf_variables)
  petsc_get_version()

  # Put variables into environment since they are needed to get
  # configuration (petscvariables) in the PETSc makefile
  set (ENV{PETSC_DIR} "${PETSC_DIR}")
  set (ENV{PETSC_ARCH} "${PETSC_ARCH}")

  # A temporary makefile to probe the PETSc configuration
  set (petsc_config_makefile "${PROJECT_BINARY_DIR}/Makefile.petsc")
  file (WRITE "${petsc_config_makefile}"
"## This file was autogenerated by FindPETSc.cmake
PETSC_DIR  = $ENV{PETSC_DIR}
PETSC_ARCH = $ENV{PETSC_ARCH}
include ${petsc_conf_rules}
include ${petsc_conf_variables}
include ${petsc_conf_petscvariables}
show :
	-@echo -n \${\${VARIABLE}}
")



  macro (PETSC_GET_VARIABLE name var)
    set (${var} "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)
    execute_process (COMMAND ${MAKE_EXECUTABLE} --no-print-directory -f ${petsc_config_makefile} show VARIABLE=${name}
      OUTPUT_VARIABLE ${var}
      RESULT_VARIABLE petsc_return)
  endmacro (PETSC_GET_VARIABLE)
  petsc_get_variable (PETSC_LIB_DIR            petsc_lib_dir)
  petsc_get_variable (PETSC_EXTERNAL_LIB_BASIC petsc_libs_external)
  petsc_get_variable (PETSC_CCPPFLAGS          petsc_cpp_line)
  petsc_get_variable (PETSC_FC_INCLUDES        petsc_include)
  petsc_get_variable (PCC                      petsc_cc)
  petsc_get_variable (MPIEXEC                  petsc_mpiexec)
  petsc_get_variable (PETSC_LIB                petsc_lib)
  
  # We are done with the temporary Makefile, calling PETSC_GET_VARIABLE after this point is invalid!
  file (REMOVE ${petsc_config_makefile})

  include (${${PROJECT_NAME}_TRIBITS_DIR}/tpls/ResolveCompilerPaths.cmake)
  # Extract include paths and libraries from compile command line
  resolve_includes (petsc_includes_all "${petsc_cpp_line}")

  macro (PETSC_FIND_LIBRARY suffix name)
    set (PETSC_LIBRARY_${suffix} "NOTFOUND" CACHE INTERNAL "Cleared" FORCE) # Clear any stale value, if we got here, we need to find it again
    find_library (PETSC_LIBRARY_${suffix} NAMES ${name} HINTS ${petsc_lib_dir} NO_DEFAULT_PATH)
    set (PETSC_LIBRARIES_${suffix} "${PETSC_LIBRARY_${suffix}}")
    mark_as_advanced (PETSC_LIBRARY_${suffix})
  endmacro (PETSC_FIND_LIBRARY suffix name)

  # Look for petscvec first, if it doesn't exist, we must be using single-library
  petsc_find_library (VEC petscvec)
  if (PETSC_LIBRARY_VEC)
    petsc_find_library (SYS  "petscsys;petsc") # libpetscsys is called libpetsc prior to 3.1 (when single-library was introduced)
    petsc_find_library (MAT  petscmat)
    petsc_find_library (DM   petscdm)
    petsc_find_library (KSP  petscksp)
    petsc_find_library (SNES petscsnes)
    petsc_find_library (TS   petscts)
    macro (PETSC_JOIN libs deps)
      list (APPEND PETSC_LIBRARIES_${libs} ${PETSC_LIBRARIES_${deps}})
    endmacro (PETSC_JOIN libs deps)
    petsc_join (VEC  SYS)
    petsc_join (MAT  VEC)
    petsc_join (DM   MAT)
    petsc_join (KSP  DM)
    petsc_join (SNES KSP)
    petsc_join (TS   SNES)
    petsc_join (ALL  TS)
  else ()
    set (PETSC_LIBRARY_VEC "NOTFOUND" CACHE INTERNAL "Cleared" FORCE) # There is no libpetscvec
    petsc_find_library (SINGLE petsc)
    foreach (pkg SYS VEC MAT DM KSP SNES TS ALL)
      set (PETSC_LIBRARIES_${pkg} "${PETSC_LIBRARY_SINGLE}")
    endforeach ()
  endif ()
  if (PETSC_LIBRARY_TS)
    message (STATUS "Recognized PETSc install with separate libraries for each package")
  else ()
    message (STATUS "Recognized PETSc install with single library for all packages")
  endif ()

  find_path (PETSC_INCLUDE_DIR petscts.h HINTS "${PETSC_DIR}" PATH_SUFFIXES include NO_DEFAULT_PATH)
  find_path (PETSC_INCLUDE_CONF petscconf.h HINTS "${PETSC_DIR}" PATH_SUFFIXES "${PETSC_ARCH}/include" "bmake/${PETSC_ARCH}" NO_DEFAULT_PATH)
  mark_as_advanced (PETSC_INCLUDE_DIR PETSC_INCLUDE_CONF)
  
    # We do an out-of-source build so __FILE__ will be an absolute path, hence __INSDIR__ is superfluous
  if (${PETSC_VERSION} VERSION_LESS 3.1)
    set (PETSC_DEFINITIONS "-D__SDIR__=\"\"" CACHE STRING "PETSc definitions" FORCE)
  else ()
    set (PETSC_DEFINITIONS "-D__INSDIR__=" CACHE STRING "PETSc definitions" FORCE)
  endif ()
  # Sometimes this can be used to assist FindMPI.cmake
  set (PETSC_MPIEXEC ${petsc_mpiexec} CACHE FILEPATH "Executable for running PETSc MPI programs" FORCE)
  string (REGEX REPLACE " -I" ";" PETSC_INCLUDE_DIRS "${petsc_include}")
  string (REGEX REPLACE "-I" "" PETSC_INCLUDE_DIRS "${PETSC_INCLUDE_DIRS}")
  set (TPL_PETSC_INCLUDE_DIRS "${PETSC_INCLUDE_DIRS}" CACHE STRING "PETSc include path" FORCE)
  set (TPL_PETSC_LIBRARIES ${PETSC_LIBRARIES_ALL} CACHE STRING "PETSc libraries" FORCE)
  set (PETSC_COMPILER ${petsc_cc} CACHE FILEPATH "PETSc compiler" FORCE)
  set (TPL_PETSC_LIBRARY_DIRS "${PETSC_DIR}/${PETSC_ARCH}/lib" CACHE STRING "PETSc library paths" FORCE)
  MESSAGE(STATUS "  TPL_PETSC_INCLUDE_DIRS='${TPL_PETSC_INCLUDE_DIRS}'")
  MESSAGE(STATUS "  TPL_PETSC_LIBRARIES='${TPL_PETSC_LIBRARIES}'")
  # Note that we have forced values for all these choices.  If you
  # change these, you are telling the system to trust you that they
  # work.  It is likely that you will end up with a broken build.
  mark_as_advanced (PETSC_INCLUDES PETSC_LIBRARIES PETSC_COMPILER PETSC_DEFINITIONS PETSC_MPIEXEC)

  #check if cmake version > 2.8.5
  SET(major_ver_req 2)
  SET(minor_ver_req 8)
  SET(patch_ver_req 5)
  IF(${CMAKE_MAJOR_VERSION} GREATER major_ver_req OR ${CMAKE_MAJOR_VERSION} EQUAL major_ver_req)
    IF(${CMAKE_MINOR_VERSION} GREATER minor_ver_req OR ${CMAKE_MINOR_VERSION} EQUAL minor_ver_req)
      IF(${CMAKE_PATCH_VERSION} GREATER patch_ver_req OR ${CMAKE_PATCH_VERSION} EQUAL patch_ver_req)
        IF (NOT "${petsc_lib}" STREQUAL "")
          # Look for certain flags that may need to be added/moved to the end of the linking process:
          
          STRING(FIND ${petsc_lib} "-lX11" position)
          IF (position GREATER 0)
            SET(${PROJECT_NAME}_EXTRA_LINK_FLAGS ${${PROJECT_NAME}_EXTRA_LINK_FLAGS} "-lX11")
          ENDIF ()
          
          STRING(FIND ${petsc_lib} "lparmetis" position)
          IF (position GREATER 0)
            SET(${PROJECT_NAME}_EXTRA_LINK_FLAGS ${${PROJECT_NAME}_EXTRA_LINK_FLAGS} "-lparmetis")
          ENDIF () 
          
          STRING(FIND ${petsc_lib} "lmetis" position) 
          IF (position GREATER 0)
            SET(${PROJECT_NAME}_EXTRA_LINK_FLAGS ${${PROJECT_NAME}_EXTRA_LINK_FLAGS} "-lmetis")
          ENDIF ()
          
          STRING(FIND ${petsc_lib} "lflapack" position)
          IF (position GREATER 0)
            SET(${PROJECT_NAME}_EXTRA_LINK_FLAGS ${${PROJECT_NAME}_EXTRA_LINK_FLAGS} "-lflapack")
          ENDIF () 
          
          STRING(FIND ${petsc_lib} "lfblas" position)
          IF (position GREATER 0)
            SET(${PROJECT_NAME}_EXTRA_LINK_FLAGS ${${PROJECT_NAME}_EXTRA_LINK_FLAGS} "-lfblas")
          ENDIF () 
          
          STRING(FIND ${petsc_lib} "lmkl_intel_lp64" position)
          IF (position GREATER 0)
            SET(${PROJECT_NAME}_EXTRA_LINK_FLAGS ${${PROJECT_NAME}_EXTRA_LINK_FLAGS} "-lmkl_intel_lp64")
          ENDIF ()
          
          STRING(FIND ${petsc_lib} "lmkl_intel_thread" position)
          IF (position GREATER 0)
            SET(${PROJECT_NAME}_EXTRA_LINK_FLAGS ${${PROJECT_NAME}_EXTRA_LINK_FLAGS} "-lmkl_intel_thread")
          ENDIF ()    
          
          STRING(FIND ${petsc_lib} "lmkl_core" position)
          IF (position GREATER 0)
            SET(${PROJECT_NAME}_EXTRA_LINK_FLAGS ${${PROJECT_NAME}_EXTRA_LINK_FLAGS} "-lmkl_core")
          ENDIF ()    
          
          STRING(FIND ${petsc_lib} "liomp5" position)
          IF (position GREATER 0)
            SET(${PROJECT_NAME}_EXTRA_LINK_FLAGS ${${PROJECT_NAME}_EXTRA_LINK_FLAGS} "-liomp5")
          ENDIF ()    

        ENDIF ()  
      ELSE()
        MESSAGE(STATUS "STRING FIND is not available with this cmake version (2.8.5 required).  There may be issues with link flags.")
      ENDIF()
    ELSE()
      MESSAGE(STATUS "STRING FIND is not available with this cmake version (2.8.5 required).  There may be issues with link flags.")
    ENDIF()
  ELSE()
    MESSAGE(STATUS "STRING FIND is not available with this cmake version (2.8.5 required).  There may be issues with link flags.")
  ENDIF()

endif ()
