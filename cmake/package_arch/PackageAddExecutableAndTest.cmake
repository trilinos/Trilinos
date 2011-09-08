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


INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddTest)


MACRO(FWD_PARSE_ARG  VAR_TO_SET_OUT  ARGNAME)
  IF (PARSE_${ARGNAME})
    SET(${VAR_TO_SET_OUT} ${${VAR_TO_SET_OUT}} ${ARGNAME} ${PARSE_${ARGNAME}}) 
  ENDIF()
ENDMACRO()


MACRO(FWD_PARSE_OPT  VAR_TO_SET_OUT  OPTNAME)
  IF (PARSE_${OPTNAME})
    SET(${VAR_TO_SET_OUT} ${${VAR_TO_SET_OUT}} ${OPTNAME}) 
  ENDIF()
ENDMACRO()


FUNCTION(PACKAGE_ADD_EXECUTABLE_WRAPPER)
  IF (PACKAGE_ADD_EXECUTABLE_AND_TEST_TEST_MODE)
    SET(PACKAGE_ADD_EXECUTABLE_CAPTURE_ARGS ${ARGN} CACHE INTERNAL "")
  ELSE()
    PACKAGE_ADD_EXECUTABLE(${ARGN})
  ENDIF()
ENDFUNCTION()


FUNCTION(PACKAGE_ADD_TEST_WRAPPER)
  IF (PACKAGE_ADD_EXECUTABLE_AND_TEST_TEST_MODE)
    SET(PACKAGE_ADD_TEST_CAPTURE_ARGS ${ARGN} CACHE INTERNAL "")
  ELSE()
    PACKAGE_ADD_TEST(${ARGN})
  ENDIF()
ENDFUNCTION()

#
# Add an executable and a test at the same time.
#
# PACKAGE_ADD_EXECUTABLE_AND_TEST(
#   <execName>
#   SOURCES <src1> <src2> ...
#   [NAME <testName> | NAME_POSTFIX <testNamePostfix>]
#   [CATEGORIES <category1>  <category2> ...]
#   [HOST <host1> <host2> ...]
#   [XHOST <host1> <host2> ...]
#   [HOSTTYPE <hosttype1> <hosttype2> ...]
#   [XHOSTTYPE <hosttype1> <hosttype2> ...]
#   [NOEXEPREFIX ]
#   [DIRECTORY <dir> ]
#   [DEPLIBS <lib1> <lib2> ... ]
#   [COMM [serial] [mpi] ]
#   [NUM_MPI_PROCS <numProcs>]
#   [LINKER_LANGUAGE [C|CXX|Fortran] ]
#   [ADD_DIR_TO_NAME ]
#   [DEFINES <-DSOMEDEFINE>]
#   [KEYWORDS <keyword1> <keyword2> ...]
#   [STANDARD_PASS_OUTPUT
#     | PASS_REGULAR_EXPRESSION "<regex1>;<regex2>;..."]
#   [FAIL_REGULAR_EXPRESSION "<regex1>;<regex2>;..."]
#   [WILL_FAIL]
#   [INSTALLABLE]
#   )
#
# This function takes a fairly common set of arguments to
# PACKAGE_ADD_EXECUTABLE(...) and PACAKGE_ADD_TEST(...) but not the full set
# passed to PACAKGE_ADD_TEST(...).  See the documentation for
# PACKAGE_ADD_EXECUTABLE(...) and PACKAGE_ADD_TEST(...) to see which arguments
# are accpeted by which functions.
#

# Arguments that are specific to this function and not contained in
# PACKAGE_ADD_EXECUTABLE(...) or PACAKGE_ADD_TEST(...):
#
#   XHOST_TEST <host1> <host2> ...
#
#     When specified, this disables just running the tests for the named hosts
#     <host1>, <host2> etc. but still builds the executables for the test.
#
#   XHOSTTYPE_TEST <hosttype1> <hosttype2> ...
#
#     When specified, this disables just running the tests for the named host
#     types <hosttype1>, <hosttype2> etc. but still builds the executables for
#     the test.
#

FUNCTION(PACKAGE_ADD_EXECUTABLE_AND_TEST EXE_NAME)
   
  #
  # A) Parse the input arguments
  #

  PARSE_ARGUMENTS(
     #prefix
     PARSE
     #lists
     "SOURCES;DEPLIBS;NAME;NAME_POSTFIX;NUM_MPI_PROCS;DIRECTORY;KEYWORDS;COMM;ARGS;NAME;PASS_REGULAR_EXPRESSION;CATEGORIES;HOST;XHOST;XHOST_TEST;HOSTTYPE;XHOSTTYPE;XHOSTTYPE_TEST;FAIL_REGULAR_EXPRESSION;LINKER_LANGUAGE;DEFINES"
     #options
     "STANDARD_PASS_OUTPUT;WILL_FAIL;ADD_DIR_TO_NAME;INSTALLABLE"
     ${ARGN}
     )

  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("")
    MESSAGE("PACKAGE_ADD_EXECUTABLE_AND_TEST: ${EXE_NAME} ${ARGN}")
  ENDIF()

  #
  # B) Arguments common to both
  #

  SET(COMMON_CALL_ARGS "")
  FWD_PARSE_ARG(COMMON_CALL_ARGS COMM)
  FWD_PARSE_ARG(COMMON_CALL_ARGS CATEGORIES)
  FWD_PARSE_ARG(COMMON_CALL_ARGS HOST)
  FWD_PARSE_ARG(COMMON_CALL_ARGS XHOST)
  FWD_PARSE_ARG(COMMON_CALL_ARGS HOSTTYPE)
  FWD_PARSE_ARG(COMMON_CALL_ARGS XHOSTTYPE)

  #
  # C) PackageAddExecutable(...)
  #

  SET(CALL_ARGS "")
  FWD_PARSE_ARG(CALL_ARGS SOURCES)
  FWD_PARSE_ARG(CALL_ARGS DEPLIBS)
  FWD_PARSE_ARG(CALL_ARGS DIRECTORY)
  FWD_PARSE_OPT(CALL_ARGS ADD_DIR_TO_NAME)
  FWD_PARSE_ARG(CALL_ARGS LINKER_LANGUAGE)
  FWD_PARSE_ARG(CALL_ARGS DEFINES)
  FWD_PARSE_OPT(CALL_ARGS INSTALLABLE)

  PACKAGE_ADD_EXECUTABLE_WRAPPER(${EXE_NAME} ${COMMON_CALL_ARGS} ${CALL_ARGS})

  #
  # D) PackageAddTest(...)
  #

  SET(CALL_ARGS "")
  FWD_PARSE_ARG(CALL_ARGS NAME)
  FWD_PARSE_ARG(CALL_ARGS NAME_POSTFIX)
  FWD_PARSE_ARG(CALL_ARGS DIRECTORY)
  FWD_PARSE_ARG(CALL_ARGS KEYWORDS)
  FWD_PARSE_ARG(CALL_ARGS NUM_MPI_PROCS)
  FWD_PARSE_ARG(CALL_ARGS ARGS)
  FWD_PARSE_ARG(CALL_ARGS PASS_REGULAR_EXPRESSION)
  FWD_PARSE_ARG(CALL_ARGS FAIL_REGULAR_EXPRESSION)
  FWD_PARSE_OPT(CALL_ARGS STANDARD_PASS_OUTPUT)
  FWD_PARSE_OPT(CALL_ARGS WILL_FAIL)
  FWD_PARSE_OPT(CALL_ARGS ADD_DIR_TO_NAME)
  IF (PARSE_XHOST_TEST)
    SET(CALL_ARGS ${CALL_ARGS} XHOST ${PARSE_XHOST_TEST})
  ENDIF()
  IF (PARSE_XHOSTTYPE_TEST)
    SET(CALL_ARGS ${CALL_ARGS} XHOSTTYPE ${PARSE_XHOSTTYPE_TEST})
  ENDIF()

  PACKAGE_ADD_TEST_WRAPPER(${EXE_NAME} ${COMMON_CALL_ARGS} ${CALL_ARGS})

ENDFUNCTION()
