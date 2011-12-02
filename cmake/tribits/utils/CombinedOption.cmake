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

INCLUDE(ParseVariableArguments)
INCLUDE(MultilineSet)
INCLUDE(ConcatStrings)

FUNCTION(COMBINED_OPTION  COMBINED_OPTION_NAME)

  PARSE_ARGUMENTS(
    #prefix
    PARSE
    #lists
    "DEP_OPTIONS_NAMES;DOCSTR"
    #options
    ""
    ${ARGN}
    )

  # ToDo: Assert the the right input was passed in!

  SET(DEFAULT_VAL ON)
  FOREACH( DEP_OPTION_NAME ${PARSE_DEP_OPTIONS_NAMES})
    IF (NOT ${DEP_OPTION_NAME})
      SET(DEFAULT_VAL OFF)
    ENDIF()
  ENDFOREACH()

  CONCAT_STRINGS(DOCSTR ${PARSE_DOCSTR})

  OPTION(${COMBINED_OPTION_NAME} "${DOCSTR}"
    ${DEFAULT_VAL} )

  #PRINT_VAR(${COMBINED_OPTION_NAME})

  SET(ALL_ENABLED TRUE)
  IF (${COMBINED_OPTION_NAME})
    FOREACH( DEP_OPTION_NAME ${PARSE_DEP_OPTIONS_NAMES})
      #PRINT_VAR(${DEP_OPTION_NAME})
      IF (NOT ${DEP_OPTION_NAME})
        SET(ALL_ENABLED FALSE)
        BREAK()
      ENDIF()
    ENDFOREACH()
  ENDIF()

  IF (NOT ALL_ENABLED)

    SET(OPTION_NAMES "")
    SET(OPTION_NAMES_AND_VALUES "")
    FOREACH( DEP_OPTION_NAME ${PARSE_DEP_OPTIONS_NAMES})
      IF (NOT OPTION_NAMES)
        SET(OPTION_NAMES "${DEP_OPTION_NAME}")
      ELSE()
        SET(OPTION_NAMES "${OPTION_NAMES}, ${DEP_OPTION_NAME}")
      ENDIF()
      SET(OPTION_NAMES_AND_VALUES
        "${OPTION_NAMES_AND_VALUES}  ${DEP_OPTION_NAME}='${${DEP_OPTION_NAME}}'\n")
    ENDFOREACH()

    MESSAGE(FATAL_ERROR 
      "Error: you can not enable the option ${COMBINED_OPTION_NAME} unless"
      " you also enable the options ${OPTION_NAMES}.  The current option"
      "values are:\n${OPTION_NAMES_AND_VALUES}" )

  ENDIF()
  
ENDFUNCTION()
