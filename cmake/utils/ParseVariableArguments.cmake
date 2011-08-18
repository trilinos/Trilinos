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


# Set PARSE_ARGUMENTS_DUMP_OUTPUT_ENABLED to TRUE to see output from parsing.

FUNCTION(PARSE_ARGUMENTS_DUMP_OUTPUT)
  IF (PARSE_ARGUMENTS_DUMP_OUTPUT_ENABLED)
    MESSAGE(${ARGN})
  ENDIF()
ENDFUNCTION()

#
# Parse a set of input arguments into different lists
#

MACRO(PARSE_ARGUMENTS prefix arg_names option_names)

  PARSE_ARGUMENTS_DUMP_OUTPUT("PARSE_ARGUMENTS: prefix='${prefix}'")
  PARSE_ARGUMENTS_DUMP_OUTPUT("PARSE_ARGUMENTS: arg_names='${arg_names}'")
  PARSE_ARGUMENTS_DUMP_OUTPUT("PARSE_ARGUMENTS: option_names='${option_names}'")
  PARSE_ARGUMENTS_DUMP_OUTPUT("PARSE_ARGUMENTS: ARGN='${ARGN}'")

  SET(DEFAULT_ARGS)
  FOREACH(arg_name ${arg_names})    
    SET(${prefix}_${arg_name})
  ENDFOREACH(arg_name)
  FOREACH(option ${option_names})
    SET(${prefix}_${option} FALSE)
  ENDFOREACH(option)

  SET(current_arg_name DEFAULT_ARGS)
  SET(current_arg_list)

  FOREACH(arg ${ARGN})            
    SET(larg_names ${arg_names})    
    LIST(FIND larg_names "${arg}" is_arg_name)                   
    IF (is_arg_name GREATER -1)
      SET(${prefix}_${current_arg_name} "${current_arg_list}")
      PARSE_ARGUMENTS_DUMP_OUTPUT("PARSE_ARGUMENTS:"
        " ${prefix}_${current_arg_name} = '${${prefix}_${current_arg_name}}'" )
      SET(current_arg_name "${arg}")
      SET(current_arg_list)
    ELSE()
      SET(loption_names "${option_names}")
      LIST(FIND loption_names "${arg}" is_option)            
      IF (is_option GREATER -1)
        SET(${prefix}_${arg} TRUE)
        PARSE_ARGUMENTS_DUMP_OUTPUT( "PARSE_ARGUMENTS:"
          " ${prefix}_${arg} = '${${prefix}_${arg}}'" )
      ELSE()
        LIST(APPEND current_arg_list "${arg}")
      ENDIF()
    ENDIF()
  ENDFOREACH()

  SET(${prefix}_${current_arg_name} "${current_arg_list}")
  PARSE_ARGUMENTS_DUMP_OUTPUT( "PARSE_ARGUMENTS:"
    " ${prefix}_${current_arg_name} = '${${prefix}_${current_arg_name}}'" )

ENDMACRO(PARSE_ARGUMENTS)
