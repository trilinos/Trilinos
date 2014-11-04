# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
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
# ************************************************************************
# @HEADER



INCLUDE(ParseVariableArguments)

SET(MOCK_PROJECT_NAME DummyProj)


#####################################################################
#
# Unit test TribitsProcessEnabledTpl.cmake and
# TribitsTplFindIncludeDirsAndLibraries.cmake
#
#####################################################################


FUNCTION(CREATE_PROCESS_ENABLED_TPLS_TEST_CASE  TEST_NAME)

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "TPL_NAME;TPL_FINDMOD;INCLUDE_DIRS;LIBRARY_DIRS;PASS_REGULAR_EXPRESSION_ALL;EXTRA_ARGS"
     #options
     "TPL_FIND_SHARED_LIBS;VERBOSE_DEBUG"
     ${ARGN}
     )

  SET(EXTENDED_TEST_NAME ProcessEnabledTpl_${TEST_NAME})
  SET(FULL_TEST_NAME ${PACKAGE_NAME}_${EXTENDED_TEST_NAME})

  SET(TEST_OUTPUT_REL_DIR ${FULL_TEST_NAME})
  SET(TEST_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/${TEST_OUTPUT_REL_DIR})

  TRIBITS_ADD_ADVANCED_TEST( ${EXTENDED_TEST_NAME}
    OVERALL_WORKING_DIRECTORY ${TEST_OUTPUT_REL_DIR}
    OVERALL_NUM_MPI_PROCS 1
    TEST_0 CMND "${CMAKE_COMMAND}"
      ARGS
        -DPROJECT_NAME=${MOCK_PROJECT_NAME}
        -D${MOCK_PROJECT_NAME}_TRIBITS_DIR=${${PROJECT_NAME}_TRIBITS_DIR}
        -DTPL_NAME=${PARSE_TPL_NAME}
        -D${PARSE_TPL_NAME}_FINDMOD=${PARSE_TPL_FINDMOD}
        -D${PARSE_TPL_NAME}_INCLUDE_DIRS=${PARSE_INCLUDE_DIRS}
        -D${PARSE_TPL_NAME}_LIBRARY_DIRS=${PARSE_LIBRARY_DIRS}
        -DTPL_FIND_SHARED_LIBS=${PARSE_TPL_FIND_SHARED_LIBS}
        -DTRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_VERBOSE=${PARSE_VERBOSE_DEBUG}
        ${PARSE_EXTRA_ARGS}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/TribitsProcessEnabledTplTester.cmake
      PASS_REGULAR_EXPRESSION_ALL ${PARSE_PASS_REGULAR_EXPRESSION_ALL}
    XHOSTTYPE Windows
    )

ENDFUNCTION()

#
# Test with HeadersOnlyTpl
#


SET(HeaderOnlyTpl_BASE_DIR ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/tpls/HeaderOnlyTpl)
SET(HeaderOnlyTpl_FINDMOD
  ${${PROJECT_NAME}_TRIBITS_DIR}/doc/examples/TribitsExampleProject/cmake/tpls/FindTPLHeaderOnlyTpl.cmake)


CREATE_PROCESS_ENABLED_TPLS_TEST_CASE(
  HeaderOnlyTpl_FindAll
  TPL_NAME  HeaderOnlyTpl
  TPL_FINDMOD  ${HeaderOnlyTpl_FINDMOD}
  INCLUDE_DIRS  ${HeaderOnlyTpl_BASE_DIR}
  PASS_REGULAR_EXPRESSION_ALL
    "Processing enabled TPL: HeaderOnlyTpl .enabled explicitly, disable with -DTPL_ENABLE_HeaderOnlyTpl=OFF."
    "-- Must find at least one header in each of the header sets .HeaderOnlyTpl_stuff.hpp."
    "-- Searching for headers in HeaderOnlyTpl_INCLUDE_DIRS='.+/HeaderOnlyTpl'"
    "-- Searching for a header file in the set .HeaderOnlyTpl_stuff.hpp.:"
    "--   Searching for header 'HeaderOnlyTpl_stuff.hpp' ..."
    "--     Found header '.+/HeaderOnlyTpl/HeaderOnlyTpl_stuff.hpp'"
    "-- Found TPL 'HeaderOnlyTpl' include dirs '.+/HeaderOnlyTpl'"
    "-- TPL_HeaderOnlyTpl_INCLUDE_DIRS='.+/HeaderOnlyTpl'"
    "Exported TPL_ENABLE_HeaderOnlyTpl='ON'"
    "Exported TPL_HeaderOnlyTpl_NOT_FOUND=''"
    "Exported TPL_HeaderOnlyTpl_INCLUDE_DIRS='.+/doc/examples/tpls/HeaderOnlyTpl'"
    "Exported TPL_HeaderOnlyTpl_LIBRARIES=''"
  )


CREATE_PROCESS_ENABLED_TPLS_TEST_CASE(
  HeaderOnlyTpl_FindAll_TentativeEnable
  TPL_NAME  HeaderOnlyTpl
  TPL_FINDMOD  ${HeaderOnlyTpl_FINDMOD}
  INCLUDE_DIRS  ${HeaderOnlyTpl_BASE_DIR}
  EXTRA_ARGS  -DTPL_TENTATIVE_ENABLE_HeaderOnlyTpl=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Processing enabled TPL: HeaderOnlyTpl .enabled explicitly, disable with -DTPL_ENABLE_HeaderOnlyTpl=OFF."
    "-- Attempting to tentatively enable TPL 'HeaderOnlyTpl' [.][.][.]"
    "-- TPL_HeaderOnlyTpl_INCLUDE_DIRS='.+/doc/examples/tpls/HeaderOnlyTpl'"
    "-- Attempt to tentatively enable TPL 'HeaderOnlyTpl' passed!"
    "Exported TPL_ENABLE_HeaderOnlyTpl='ON'"
    "Exported TPL_HeaderOnlyTpl_NOT_FOUND=''"
    "Exported TPL_HeaderOnlyTpl_INCLUDE_DIRS='.+/doc/examples/tpls/HeaderOnlyTpl'"
    "Exported TPL_HeaderOnlyTpl_LIBRARIES=''"
  )

CREATE_PROCESS_ENABLED_TPLS_TEST_CASE(
  HeaderOnlyTpl_TplIncludeDirsOverride
  TPL_NAME  HeaderOnlyTpl
  TPL_FINDMOD  ${HeaderOnlyTpl_FINDMOD}
  EXTR_ARGS  -DTPL_HeaderOnlyTpl_INCLUDE_DIRS=/some/dir/somewhere
  PASS_REGULAR_EXPRESSION_ALL
    "Processing enabled TPL: HeaderOnlyTpl .enabled explicitly, disable with -DTPL_ENABLE_HeaderOnlyTpl=OFF."
    "-- TPL_HeaderOnlyTpl_INCLUDE_DIRS='/some/dir/somewhere'"
    "Exported TPL_HeaderOnlyTpl_NOT_FOUND=''"
    "Exported TPL_HeaderOnlyTpl_INCLUDE_DIRS='/some/dir/somewhere'"
    "Exported TPL_HeaderOnlyTpl_LIBRARIES=''"
  )


CREATE_PROCESS_ENABLED_TPLS_TEST_CASE(
  HeaderOnlyTpl_BadPath_ImplicitEnable
  TPL_NAME  HeaderOnlyTpl
  TPL_FINDMOD  ${HeaderOnlyTpl_FINDMOD}
  INCLUDE_DIRS  "/path_does_not_exist"
  EXTR_ARGS  -DTPL_HeaderOnlyTpl_ENABLING_PKG=SimpleCxx
  PASS_REGULAR_EXPRESSION_ALL
    "Processing enabled TPL: HeaderOnlyTpl .enabled by SimpleCxx, disable with -DTPL_ENABLE_HeaderOnlyTpl=OFF."
    "-- Must find at least one header in each of the header sets .HeaderOnlyTpl_stuff.hpp."
    "-- Searching for headers in HeaderOnlyTpl_INCLUDE_DIRS='/path_does_not_exist'"
    "-- Searching for a header file in the set .HeaderOnlyTpl_stuff.hpp.:"
    "--   Searching for header 'HeaderOnlyTpl_stuff.hpp' ..."
    "-- ERROR: Could not find a header file in the set .HeaderOnlyTpl_stuff.hpp."
    "-- ERROR: Could not find the include directories for TPL 'HeaderOnlyTpl'!"
    "-- TIP: If the TPL 'HeaderOnlyTpl' is on your system then you can set:"
    "     -DHeaderOnlyTpl_INCLUDE_DIRS='<dir0>[;]<dir1>[;]...'"
    "   to point to directories where these header files may be found."
    "   Or, just set:"
    "     -DTPL_HeaderOnlyTpl_INCLUDE_DIRS='<dir0>[;]<dir1>[;]...'"
    "   to point to the include directories which will bypass any search for"
    "   header files and these include directories will be used without"
    "   question in the build.  .But this will result in a build-time error"
    "   obviously if the necessary header files are not found in these"
    "   include directories.."
    "-- ERROR: Failed finding all of the parts of TPL 'HeaderOnlyTpl' .see above., Aborting!"
    "TIP: One way to get past the configure failure for the"
    "TPL 'HeaderOnlyTpl' is to simply disable it with:"
    "  -DTPL_ENABLE_HeaderOnlyTpl=OFF"
    "which will disable it and will recursively disable all of the"
    "downstream packages that have required dependencies on it, including"
    "the package 'SimpleCxx' which triggered its enable."
    "When you reconfigure, just grep the cmake stdout for 'HeaderOnlyTpl'"
    "and then follow the disables that occur as a result to see what impact"
    "this TPL disable has on the configuration of DummyProj."
    "CMake Error at .+/TribitsProcessEnabledTpl[.]cmake:[0-9]+ .MESSAGE.:"
    "ERROR: TPL_HeaderOnlyTpl_NOT_FOUND=TRUE, aborting!"
  )


CREATE_PROCESS_ENABLED_TPLS_TEST_CASE(
  HeaderOnlyTpl_BadPath_ExplicitEnable
  TPL_NAME  HeaderOnlyTpl
  TPL_FINDMOD  ${HeaderOnlyTpl_FINDMOD}
  INCLUDE_DIRS  "/path_does_not_exist"
  PASS_REGULAR_EXPRESSION_ALL
    "Processing enabled TPL: HeaderOnlyTpl .enabled explicitly, disable with -DTPL_ENABLE_HeaderOnlyTpl=OFF."
    "TIP: Even though the TPL 'HeaderOnlyTpl' was explicitly enabled in input,"
    "it can be disabled with:"
    "  -DTPL_ENABLE_HeaderOnlyTpl=OFF"
    "which will disable it and will recursively disable all of the"
    "downstream packages that have required dependencies on it."
    "When you reconfigure, just grep the cmake stdout for 'HeaderOnlyTpl'"
    "and then follow the disables that occur as a result to see what impact" 
    "this TPL disable has on the configuration of ${MOCK_PROJECT_NAME}."
    "-- ERROR: Failed finding all of the parts of TPL 'HeaderOnlyTpl' .see above., Aborting!"
    "CMake Error at .+/TribitsProcessEnabledTpl[.]cmake:[0-9]+ .MESSAGE.:"
    "ERROR: TPL_HeaderOnlyTpl_NOT_FOUND=TRUE, aborting!"
  )


CREATE_PROCESS_ENABLED_TPLS_TEST_CASE(
  HeaderOnlyTpl_BadPath_TentativeEnable
  TPL_NAME  HeaderOnlyTpl
  TPL_FINDMOD  ${HeaderOnlyTpl_FINDMOD}
  INCLUDE_DIRS  "/path_does_not_exist"
  EXTRA_ARGS  -DTPL_TENTATIVE_ENABLE_HeaderOnlyTpl=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Processing enabled TPL: HeaderOnlyTpl .enabled explicitly, disable with -DTPL_ENABLE_HeaderOnlyTpl=OFF."
    "-- Attempting to tentatively enable TPL 'HeaderOnlyTpl' ..."
    "-- Must find at least one header in each of the header sets .HeaderOnlyTpl_stuff.hpp."
    "-- Searching for headers in HeaderOnlyTpl_INCLUDE_DIRS='/path_does_not_exist'"
    "-- Searching for a header file in the set .HeaderOnlyTpl_stuff.hpp.:"
    "--   Searching for header 'HeaderOnlyTpl_stuff.hpp' ..."
    "-- ERROR: Could not find a header file in the set .HeaderOnlyTpl_stuff.hpp."
    "-- ERROR: Could not find the include directories for TPL 'HeaderOnlyTpl'!"
    "-- TIP: If the TPL 'HeaderOnlyTpl' is on your system then you can set:"
    "-- ERROR: Failed finding all of the parts of TPL 'HeaderOnlyTpl' .see above., Aborting!"
    "-- TPL_HeaderOnlyTpl_INCLUDE_DIRS=''"
    "-- Attempt to tentatively enable TPL 'HeaderOnlyTpl' failed!  Setting TPL_ENABLE_HeaderOnlyTpl=OFF"
    "Exported TPL_ENABLE_HeaderOnlyTpl='OFF'"
    "Exported TPL_HeaderOnlyTpl_NOT_FOUND='TRUE'"
    "Exported TPL_HeaderOnlyTpl_LIBRARIES=''"
    "Exported TPL_HeaderOnlyTpl_INCLUDE_DIRS=''"
  )


#
# Test with HeadersAndLibsTpl
#


SET(HeadersAndLibsTpl_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tpls/HeadersAndLibsTpl)
SET(HeadersAndLibsTpl_FINDMOD ${HeadersAndLibsTpl_BASE_DIR}/FindTPLHeadersAndLibsTpl.cmake)


CREATE_PROCESS_ENABLED_TPLS_TEST_CASE(
  HeadersAndLibsTpl_Shared_FindAll
  TPL_NAME  HeadersAndLibsTpl
  TPL_FINDMOD  ${HeadersAndLibsTpl_FINDMOD}
  INCLUDE_DIRS  ${HeadersAndLibsTpl_BASE_DIR}/include
  LIBRARY_DIRS  ${HeadersAndLibsTpl_BASE_DIR}/lib
  TPL_FIND_SHARED_LIBS
  PASS_REGULAR_EXPRESSION_ALL
    "Processing enabled TPL: HeadersAndLibsTpl .enabled explicitly, disable with -DTPL_ENABLE_HeadersAndLibsTpl=OFF." 
    "-- Must find at least one lib in each of the lib sets .haltpl1[;]haltpl2."
    "-- Searching for libs in HeadersAndLibsTpl_LIBRARY_DIRS='.+/HeadersAndLibsTpl/lib'"
    "-- Searching for a lib in the set .haltpl1.:"
    "--   Searching for lib 'haltpl1' [.][.][.]"
    "--     Found lib '.+/HeadersAndLibsTpl/lib/libhaltpl1.so'"
    "-- Searching for a lib in the set .haltpl2.:"
    "--   Searching for lib 'haltpl2' [.][.][.]"
    "--     Found lib '.+/HeadersAndLibsTpl/lib/libhaltpl2.so'"
    "-- TPL_HeadersAndLibsTpl_LIBRARIES='.+/HeadersAndLibsTpl/lib/libhaltpl1.so[;].+/HeadersAndLibsTpl/lib/libhaltpl2.so'"
    "-- Must find at least one header in each of the header sets .HeadersAndLibsTpl_header1.hpp[;]HeadersAndLibsTpl_header2.hpp."
    "-- Searching for headers in HeadersAndLibsTpl_INCLUDE_DIRS='.+/HeadersAndLibsTpl/include'"
    "-- Searching for a header file in the set .HeadersAndLibsTpl_header1.hpp.:"
    "--   Searching for header 'HeadersAndLibsTpl_header1.hpp' [.][.][.]"
    "--     Found header '.+/HeadersAndLibsTpl/include/HeadersAndLibsTpl_header1.hpp'"
    "-- Searching for a header file in the set .HeadersAndLibsTpl_header2.hpp.:"
    "--   Searching for header 'HeadersAndLibsTpl_header2.hpp' [.][.][.]"
    "--     Found header '.+/HeadersAndLibsTpl/include/HeadersAndLibsTpl_header2.hpp'"
    "-- Found TPL 'HeadersAndLibsTpl' include dirs '.+/HeadersAndLibsTpl/include'"
    "-- TPL_HeadersAndLibsTpl_INCLUDE_DIRS='.+/HeadersAndLibsTpl/include'"
    "Exported TPL_HeadersAndLibsTpl_NOT_FOUND=''"
    "Exported TPL_HeadersAndLibsTpl_LIBRARIES='.+/HeadersAndLibsTpl/lib/libhaltpl1.so[;].+/HeadersAndLibsTpl/lib/libhaltpl2.so'"
    "Exported TPL_HeadersAndLibsTpl_INCLUDE_DIRS='.+/HeadersAndLibsTpl/include'"
  )


CREATE_PROCESS_ENABLED_TPLS_TEST_CASE(
  HeadersAndLibsTpl_Shared_FindAll_sets1
  TPL_NAME  HeadersAndLibsTpl
  TPL_FINDMOD  ${HeadersAndLibsTpl_FINDMOD}
  INCLUDE_DIRS  ${HeadersAndLibsTpl_BASE_DIR}/include
  LIBRARY_DIRS  ${HeadersAndLibsTpl_BASE_DIR}/lib
  TPL_FIND_SHARED_LIBS
  EXTRA_ARGS -DFIND_ONE_IN_SET_OF_HEADERS_1=ON -DFIND_ONE_IN_SET_OF_LIBS_1=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Processing enabled TPL: HeadersAndLibsTpl .enabled explicitly, disable with -DTPL_ENABLE_HeadersAndLibsTpl=OFF."
    "-- Must find at least one lib in each of the lib sets .missinglib1 haltpl1[;]missinglib2 haltpl2."
    "-- Searching for libs in HeadersAndLibsTpl_LIBRARY_DIRS='.+/HeadersAndLibsTpl/lib'"
    "-- Searching for a lib in the set .missinglib1 haltpl1.:"
    "--   Searching for lib 'missinglib1' [.][.][.]"
    "--   Searching for lib 'haltpl1' [.][.][.]"
    "--     Found lib '.+/HeadersAndLibsTpl/lib/libhaltpl1.so'"
    "-- Searching for a lib in the set .missinglib2 haltpl2.:"
    "--   Searching for lib 'missinglib2' [.][.][.]"
    "--   Searching for lib 'haltpl2' [.][.][.]"
    "--     Found lib '.+/HeadersAndLibsTpl/lib/libhaltpl2.so'"
    "-- TPL_HeadersAndLibsTpl_LIBRARIES='.+/HeadersAndLibsTpl/lib/libhaltpl1.so[;].+/HeadersAndLibsTpl/lib/libhaltpl2.so'"
    "-- Must find at least one header in each of the header sets .MissingHeader1.hpp HeadersAndLibsTpl_header1.hpp[;]MissingHeader2.hpp HeadersAndLibsTpl_header2.hpp."
    "-- Searching for headers in HeadersAndLibsTpl_INCLUDE_DIRS='.+/HeadersAndLibsTpl/include'"
    "-- Searching for a header file in the set .MissingHeader1.hpp HeadersAndLibsTpl_header1.hpp.:"
    "--   Searching for header 'MissingHeader1.hpp' [.][.][.]"
    "--   Searching for header 'HeadersAndLibsTpl_header1.hpp' [.][.][.]"
    "--     Found header '.+/HeadersAndLibsTpl/include/HeadersAndLibsTpl_header1.hpp'"
    "-- Searching for a header file in the set .MissingHeader2.hpp HeadersAndLibsTpl_header2.hpp.:"
    "--   Searching for header 'MissingHeader2.hpp' [.][.][.]"
    "--   Searching for header 'HeadersAndLibsTpl_header2.hpp' [.][.][.]"
    "--     Found header '.+/HeadersAndLibsTpl/include/HeadersAndLibsTpl_header2.hpp'"
    "-- Found TPL 'HeadersAndLibsTpl' include dirs '.+/HeadersAndLibsTpl/include'"
    "-- TPL_HeadersAndLibsTpl_INCLUDE_DIRS='.+/HeadersAndLibsTpl/include'"
    "Exported TPL_HeadersAndLibsTpl_NOT_FOUND=''"
    "Exported TPL_HeadersAndLibsTpl_LIBRARIES='.+/HeadersAndLibsTpl/lib/libhaltpl1.so[;].+/HeadersAndLibsTpl/lib/libhaltpl2.so'"
    "Exported TPL_HeadersAndLibsTpl_INCLUDE_DIRS='.+/HeadersAndLibsTpl/include'"
  )


CREATE_PROCESS_ENABLED_TPLS_TEST_CASE(
  HeadersAndLibsTpl_Static_FindAll
  TPL_NAME  HeadersAndLibsTpl
  TPL_FINDMOD  ${HeadersAndLibsTpl_FINDMOD}
  INCLUDE_DIRS  ${HeadersAndLibsTpl_BASE_DIR}/include
  LIBRARY_DIRS  ${HeadersAndLibsTpl_BASE_DIR}/lib
  PASS_REGULAR_EXPRESSION_ALL
    "Processing enabled TPL: HeadersAndLibsTpl .enabled explicitly, disable with -DTPL_ENABLE_HeadersAndLibsTpl=OFF."
    "-- Must find at least one lib in each of the lib sets .haltpl1[;]haltpl2."
    "-- Searching for libs in HeadersAndLibsTpl_LIBRARY_DIRS='.+/HeadersAndLibsTpl/lib'"
    "-- Searching for a lib in the set .haltpl1.:"
    "--   Searching for lib 'haltpl1' [.][.][.]"
    "--     Found lib '.+/HeadersAndLibsTpl/lib/libhaltpl1.a'"
    "-- Searching for a lib in the set .haltpl2.:"
    "--   Searching for lib 'haltpl2' [.][.][.]"
    "--     Found lib '.+/HeadersAndLibsTpl/lib/libhaltpl2.a'"
    "-- TPL_HeadersAndLibsTpl_LIBRARIES='.+/HeadersAndLibsTpl/lib/libhaltpl1.a[;].+/HeadersAndLibsTpl/lib/libhaltpl2.a"
    "Exported TPL_HeadersAndLibsTpl_LIBRARIES='.+/HeadersAndLibsTpl/lib/libhaltpl1.a[;].+/HeadersAndLibsTpl/lib/libhaltpl2.a'"
    "Exported TPL_HeadersAndLibsTpl_INCLUDE_DIRS='.+/HeadersAndLibsTpl/include'"
  )


CREATE_PROCESS_ENABLED_TPLS_TEST_CASE(
  HeadersAndLibsTpl_Static_OverrideLibNames
  TPL_NAME  HeadersAndLibsTpl
  TPL_FINDMOD  ${HeadersAndLibsTpl_FINDMOD}
  INCLUDE_DIRS  ${HeadersAndLibsTpl_BASE_DIR}/include
  LIBRARY_DIRS  ${HeadersAndLibsTpl_BASE_DIR}/lib
  EXTRA_ARGS  -DHeadersAndLibsTpl_LIBRARY_NAMES=haltpl2
  PASS_REGULAR_EXPRESSION_ALL
    "Exported TPL_HeadersAndLibsTpl_LIBRARIES='.+/HeadersAndLibsTpl/lib/libhaltpl2.a'"
    "Exported TPL_HeadersAndLibsTpl_INCLUDE_DIRS='.+/HeadersAndLibsTpl/include'"
  )


CREATE_PROCESS_ENABLED_TPLS_TEST_CASE(
  HeadersAndLibsTpl_Static_SkipLibsOverride
  TPL_NAME  HeadersAndLibsTpl
  TPL_FINDMOD  ${HeadersAndLibsTpl_FINDMOD}
  INCLUDE_DIRS  ${HeadersAndLibsTpl_BASE_DIR}/include
  EXTRA_ARGS  -DHeadersAndLibsTpl_LIBRARY_NAMES=
  PASS_REGULAR_EXPRESSION_ALL
    "Exported TPL_HeadersAndLibsTpl_LIBRARIES=''"
    "Exported TPL_HeadersAndLibsTpl_INCLUDE_DIRS='.+/HeadersAndLibsTpl/include'"
  )


CREATE_PROCESS_ENABLED_TPLS_TEST_CASE(
  HeadersAndLibsTpl_OverrideAll
  TPL_NAME  HeadersAndLibsTpl
  TPL_FINDMOD  ${HeadersAndLibsTpl_FINDMOD}
  EXTRA_ARGS  -DTPL_HeadersAndLibsTpl_LIBRARIES=/path/to/lib/lib1.so
    -DTPL_HeadersAndLibsTpl_INCLUDE_DIRS=/path/to/include
  PASS_REGULAR_EXPRESSION_ALL
    "-- TPL_HeadersAndLibsTpl_LIBRARIES='/path/to/lib/lib1.so'"
    "-- TPL_HeadersAndLibsTpl_INCLUDE_DIRS='/path/to/include'"
    "Exported TPL_HeadersAndLibsTpl_LIBRARIES='/path/to/lib/lib1.so'"
    "Exported TPL_HeadersAndLibsTpl_INCLUDE_DIRS='/path/to/include'"
  )


CREATE_PROCESS_ENABLED_TPLS_TEST_CASE(
  HeadersAndLibsTpl_Static_BadLibPath
  TPL_NAME  HeadersAndLibsTpl
  TPL_FINDMOD  ${HeadersAndLibsTpl_FINDMOD}
  INCLUDE_DIRS  ${HeadersAndLibsTpl_BASE_DIR}/include
  LIBRARY_DIRS  /path/does/not/exist
  EXTR_ARGS  -DTPL_HeadersAndLibsTpl_ENABLING_PKG=SomePackage
  PASS_REGULAR_EXPRESSION_ALL
    "Processing enabled TPL: HeadersAndLibsTpl .enabled by SomePackage, disable with -DTPL_ENABLE_HeadersAndLibsTpl=OFF."
    "-- Must find at least one lib in each of the lib sets .haltpl1[;]haltpl2."
    "-- Searching for libs in HeadersAndLibsTpl_LIBRARY_DIRS='/path/does/not/exist'"
    "-- Searching for a lib in the set .haltpl1.:"
    "--   Searching for lib 'haltpl1' ..."
    "-- ERROR: Did not find a lib in the lib set .haltpl1. for the TPL 'HeadersAndLibsTpl'!"
    "-- Searching for a lib in the set .haltpl2.:"
    "--   Searching for lib 'haltpl2' ..."
    "-- ERROR: Did not find a lib in the lib set .haltpl2. for the TPL 'HeadersAndLibsTpl'!"
    "-- ERROR: Could not find the libraries for the TPL 'HeadersAndLibsTpl'!"
    "-- TIP: If the TPL 'HeadersAndLibsTpl' is on your system then you can set:"
    "     -DHeadersAndLibsTpl_LIBRARY_DIRS='<dir0>[;]<dir1>[;]...'"
    "   to point to the directories where these libraries may be found."
    "   Or, just set:"
    "     -DTPL_HeadersAndLibsTpl_LIBRARIES='<path-to-libs0>[;]<path-to-libs1>[;]...'"
    "   to point to the full paths for the libraries which will"
    "   bypass any search for libraries and these libraries will be used without"
    "   question in the build.  .But this will result in a build-time error"
    "   obviously if all all of the necssary symbols are not found.."
    "-- ERROR: Failed finding all of the parts of TPL 'HeadersAndLibsTpl' .see above., Aborting!"
    "TIP: One way to get past the configure failure for the"
    "TPL 'HeadersAndLibsTpl' is to simply disable it with:"
    "  -DTPL_ENABLE_HeadersAndLibsTpl=OFF"
    "which will disable it and will recursively disable all of the"
    "downstream packages that have required dependencies on it, including"
    "the package 'SomePackage' which triggered its enable."
    "When you reconfigure, just grep the cmake stdout for 'HeadersAndLibsTpl'"
    "and then follow the disables that occur as a result to see what impact"
    "this TPL disable has on the configuration of DummyProj."
   "CMake Error at .+/TribitsProcessEnabledTpl.cmake:[0-9]+ .MESSAGE.:"
    "ERROR: TPL_HeadersAndLibsTpl_NOT_FOUND=TRUE, aborting!"
  )


CREATE_PROCESS_ENABLED_TPLS_TEST_CASE(
  HeadersAndLibsTpl_Static_BadLibPath_NotMustFindAllLibs
  TPL_NAME  HeadersAndLibsTpl
  TPL_FINDMOD  ${HeadersAndLibsTpl_FINDMOD}
  INCLUDE_DIRS  ${HeadersAndLibsTpl_BASE_DIR}/include
  LIBRARY_DIRS  /path/does/not/exist
  EXTR_ARGS  -DTPL_HeadersAndLibsTpl_ENABLING_PKG=SomePackage
    -DNOT_MUST_FIND_ALL_LIBS=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Processing enabled TPL: HeadersAndLibsTpl .enabled by SomePackage, disable with -DTPL_ENABLE_HeadersAndLibsTpl=OFF."
    "-- Searching for libs in HeadersAndLibsTpl_LIBRARY_DIRS='/path/does/not/exist'"
    "-- Searching for a lib in the set .haltpl1.:"
    "-- ERROR: Did not find a lib in the lib set .haltpl1. for the TPL 'HeadersAndLibsTpl'!"
    "-- ERROR: Could not find the libraries for the TPL 'HeadersAndLibsTpl'!"
    "-- TIP: If the TPL 'HeadersAndLibsTpl' is on your system then you can set:"
    "-- ERROR: Failed finding all of the parts of TPL 'HeadersAndLibsTpl' .see above., Aborting!"
    "TIP: One way to get past the configure failure for the"
    "CMake Error at .+/TribitsProcessEnabledTpl.cmake:[0-9]+ .MESSAGE.:"
    "ERROR: TPL_HeadersAndLibsTpl_NOT_FOUND=TRUE, aborting!"
  )


CREATE_PROCESS_ENABLED_TPLS_TEST_CASE(
  HeadersAndLibsTpl_Static_BadLibName
  TPL_NAME  HeadersAndLibsTpl
  TPL_FINDMOD  ${HeadersAndLibsTpl_FINDMOD}
  INCLUDE_DIRS  ${HeadersAndLibsTpl_BASE_DIR}/include
  LIBRARY_DIRS  ${HeadersAndLibsTpl_BASE_DIR}/lib
  EXTRA_ARGS  -DHeadersAndLibsTpl_LIBRARY_NAMES="badlibname"
    -DTPL_HeadersAndLibsTpl_ENABLING_PKG=SomePackage
  PASS_REGULAR_EXPRESSION_ALL
    "Processing enabled TPL: HeadersAndLibsTpl .enabled by SomePackage, disable with -DTPL_ENABLE_HeadersAndLibsTpl=OFF."
    "-- Must find at least one lib in each of the lib sets .badlibname."
    "-- Searching for libs in HeadersAndLibsTpl_LIBRARY_DIRS='.+/HeadersAndLibsTpl/lib'"
    "-- Searching for a lib in the set .badlibname.:"
    "-- ERROR: Did not find a lib in the lib set .badlibname. for the TPL 'HeadersAndLibsTpl'!"
    "-- ERROR: Could not find the libraries for the TPL 'HeadersAndLibsTpl'!"
    "-- TIP: If the TPL 'HeadersAndLibsTpl' is on your system then you can set:"
    "-- ERROR: Failed finding all of the parts of TPL 'HeadersAndLibsTpl' .see above., Aborting!"
    "TIP: One way to get past the configure failure for the"
    "CMake Error at .+/TribitsProcessEnabledTpl.cmake:[0-9]+ .MESSAGE.:"
    "ERROR: TPL_HeadersAndLibsTpl_NOT_FOUND=TRUE, aborting!"
  )


CREATE_PROCESS_ENABLED_TPLS_TEST_CASE(
  HeadersAndLibsTpl_Shared_BadIncludeDir_NotMustFindAllHeaders
  TPL_NAME  HeadersAndLibsTpl
  TPL_FINDMOD  ${HeadersAndLibsTpl_FINDMOD}
  INCLUDE_DIRS  /include/path/does/not/exist
  LIBRARY_DIRS  ${HeadersAndLibsTpl_BASE_DIR}/lib
  TPL_FIND_SHARED_LIBS
  EXTR_ARGS  -DNOT_MUST_FIND_ALL_HEADERS=ON
  PASS_REGULAR_EXPRESSION_ALL
    "Processing enabled TPL: HeadersAndLibsTpl .enabled explicitly, disable with -DTPL_ENABLE_HeadersAndLibsTpl=OFF."
    "-- TPL_HeadersAndLibsTpl_LIBRARIES='.+/HeadersAndLibsTpl/lib/libhaltpl1.so[;].+/HeadersAndLibsTpl/lib/libhaltpl2.so'"
    "-- Searching for headers in HeadersAndLibsTpl_INCLUDE_DIRS='/include/path/does/not/exist'"
    "-- ERROR: Could not find a header file in the set .HeadersAndLibsTpl_header1.hpp."
    "-- ERROR: Could not find a header file in the set .HeadersAndLibsTpl_header2.hpp."
    "-- ERROR: Could not find the include directories for TPL 'HeadersAndLibsTpl'!"
    "-- TIP: If the TPL 'HeadersAndLibsTpl' is on your system then you can set:"
    "-- ERROR: Failed finding all of the parts of TPL 'HeadersAndLibsTpl' .see above., Aborting!"
    "TIP: Even though the TPL 'HeadersAndLibsTpl' was explicitly enabled in input,"
    "CMake Error at .+/TribitsProcessEnabledTpl.cmake:[0-9]+ .MESSAGE.:"
    "ERROR: TPL_HeadersAndLibsTpl_NOT_FOUND=TRUE, aborting!"
  )

# ToDo: Add test case for only finding one of the headers in REQUIRED_HEADERS ...
# ToDo: Add test case for only finding one of the libs in REQUIRED_LIBRARIES ...
