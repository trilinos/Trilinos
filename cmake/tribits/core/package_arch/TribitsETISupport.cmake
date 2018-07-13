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

include(AppendSet)
include(AppendGlobalSet)
include(Split)
include(Join)
include(PrintVar)

# this macro expands produces a list of all combinations of a variable number of lists:
# called like TRIBITS_ETI_TYPE_EXPANSION(setvar "T=t1|t2|t3" "S=s1|s2" "V=v1")
# returns setvar="T=T1 S=s1 V=v1;T=T1 S=s2 V=v1;
#                 T=T2 S=s1 V=v1;T=T2 S=s2 V=v1;
#                 T=T3 S=s1 V=v1;T=T3 S=s2 V=v1"
FUNCTION(TRIBITS_ETI_TYPE_EXPANSION outvar first_list)
  # extract the field name and the list of types
  IF (NOT "${first_list}" MATCHES "([^=]*)=(.*)")
    SET(${outvar} "TRIBITS_ETI_BAD_ARGUMENTS" PARENT_SCOPE)
    RETURN()
  ELSE()
    SET(label "${CMAKE_MATCH_1}")
    SET(tlist "${CMAKE_MATCH_2}")
  ENDIF()
  SET(eti_result "") 
  SPLIT("${tlist}" "\\|" tlist)
  IF ("${${outvar}}" STREQUAL "")
    SET(accumulate OFF)
  ELSE()
    SET(accumulate ON)
  ENDIF()
  LIST(LENGTH tlist tlist_len)
  IF (${ARGC} GREATER 2)
    TRIBITS_ETI_TYPE_EXPANSION(sub ${ARGN})
    FOREACH(t ${tlist})
      STRING(STRIP "${t}" t)
      SET(t "{${t}}")
      FOREACH(s ${sub})
        LIST(APPEND eti_result "${label}=${t} ${s}")
      ENDFOREACH()
    ENDFOREACH()
  ELSE()
    FOREACH(t ${tlist})
      STRING(STRIP "${t}" t)
      SET(t "{${t}}")
      LIST(APPEND eti_result "${label}=${t}")
    ENDFOREACH()
  ENDIF()
  IF (accumulate)
    SET(${outvar} "${${outvar}};${eti_result}" PARENT_SCOPE)
  ELSE()
    SET(${outvar} "${eti_result}"              PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()

# Explode an ETI set into a variable number of named component fields
FUNCTION(TRIBITS_ETI_EXPLODE fields inst outvar)
  foreach (field ${fields})
    IF(    "${inst}" MATCHES "( |^)${field}={([^{}$=]*)}( |$)")
      LIST(APPEND eti_result "${CMAKE_MATCH_2}")
    ELSEIF("${inst}" MATCHES "( |^)${field}=([^{} $=]*)( |$)")
      LIST(APPEND eti_result "${CMAKE_MATCH_2}")
    ELSEIF(NOT "${inst}" MATCHES "( |^)${field}=")
      LIST(APPEND eti_result "TYPE-MISSING")
    ELSE()
      SET(eti_result "TRIBITS_ETI_BAD_PARSE")
      BREAK()
    ENDIF()
  endforeach()
  JOIN(eti_result "|" FALSE ${eti_result})
  SET(${outvar} "${eti_result}" PARENT_SCOPE)
ENDFUNCTION()

# effectively, a tupled regex, wrapped in a for loop
# given a list of processed excludes (a list of pipe-separated, bracket-delimited tuples of type/regexes),
# determine whether processed_inst (a pipe-separated, bracket tuple of types) is matched
# if the instantiation matches one of the exclusions, result is set true
FUNCTION(TRIBITS_ETI_CHECK_EXCLUSION processed_excludes processed_inst excluded)
  SPLIT("${processed_inst}" "\\|" processed_inst)
  LIST(LENGTH processed_inst numfields)
  MATH(EXPR NFm1 "${numfields}-1")
  SET(${excluded} OFF PARENT_SCOPE)
  # check to see whether this is excluded or not
  list(LENGTH processed_excludes numexcl)
  FOREACH(excl ${processed_excludes})
    SPLIT(${excl} "\\|" excl)
    LIST(LENGTH excl excp_len)
    IF(${excp_len} EQUAL ${numfields})
      SET(lcltest ON)
      FOREACH(i RANGE 0 ${NFm1})
        LIST(GET processed_inst ${i} f)
        LIST(GET excl           ${i} e)
        IF (NOT "${f}" STREQUAL "TYPE-MISSING" AND NOT "${f}" MATCHES "${e}")
          SET(lcltest OFF)
          BREAK()
        ENDIF()
      ENDFOREACH()
      IF(lcltest)
        SET(${excluded} ON PARENT_SCOPE)
        IF (${PROJECT}_VERBOSE_CONFIGURE)
          MESSAGE(STATUS "-- Instantiation excluded by ${excl}")
        ENDIF()
        RETURN()
      ENDIF()
    ENDIF()
  ENDFOREACH()
ENDFUNCTION()

# given a list of field names and list of macrofields, translate the
# field names in macrofields into indices corresponding to the position in
# the list of field names
FUNCTION(TRIBITS_ETI_INDEX_MACRO_FIELDS etifields macrofields indexvar)
  FOREACH(mf ${macrofields})
    STRING(STRIP "${mf}" mf)
    LIST(FIND etifields ${mf} ind)
    IF(${ind} EQUAL -1)
      MESSAGE(FATAL_ERROR "Macro variable ${mf} not found in list of fields: ${etifields}")
    ENDIF()
    LIST(APPEND inds ${ind})
  ENDFOREACH()
  SET(${indexvar} ${inds} PARENT_SCOPE)
ENDFUNCTION()

# given a macro name and a list of tuples, generate a macro string
FUNCTION(TRIBITS_ETI_BUILD_MACRO_STRING macroname tuplelist outvar)
  SET(str "#define ${macroname}(INSTMACRO)")
  FOREACH(tuple ${tuplelist})
    STRING(REPLACE "|" " , " tuple "${tuple}")
    SET(str "${str}\\\n\tINSTMACRO( ${tuple} )")
  ENDFOREACH()
  SET(str "${str}\n")
  SET(${outvar} ${str} PARENT_SCOPE)
ENDFUNCTION()

# utility for adding eti support to another package's library set
FUNCTION(TRIBITS_ADD_ETI_INSTANTIATIONS PACKAGE)
  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE(STATUS "Adding instantiations to ${PACKAGE} library set...")
  ENDIF()
  APPEND_GLOBAL_SET(${PACKAGE}_ETI_LIBRARYSET ${ARGN})
ENDFUNCTION()

# utility for mangling type names to make them safe for the C preprocessor
FUNCTION(TRIBITS_ETI_MANGLE_SYMBOL
         mangledvar
         input)
  SET(newvar "${input}")
  STRING(REPLACE "::" "_" newvar "${newvar}")
  STRING(REPLACE ","  "_" newvar "${newvar}")
  SET(num_brackets 0)
  WHILE("${newvar}" MATCHES "(.*)<([^<>]*)>(.*)")
    SET(newvar "${CMAKE_MATCH_1}${num_brackets}${CMAKE_MATCH_2}${num_brackets}${CMAKE_MATCH_3}")
    MATH(EXPR num_brackets "${num_brackets}+1")
  ENDWHILE()
  STRING(REPLACE " "  "" newvar "${newvar}")
  SET(${mangledvar} "${newvar}" PARENT_SCOPE)
ENDFUNCTION()

# utility for mangling type names to make them safe for the C preprocessor
# upon mangling, it automatically inserts an ifdef into the mangling macro
FUNCTION(TRIBITS_ETI_MANGLE_SYMBOL_AUGMENT_MACRO
         typedefvar
         symbolvar
         manglistvar)
  SET(oldsymbol "${${symbolvar}}")
  IF(NOT "${oldsymbol}" MATCHES "::|[<>, ]")
    RETURN()
  ENDIF()
  TRIBITS_ETI_MANGLE_SYMBOL(newsymbol "${oldsymbol}")
  SET(${symbolvar} "${newsymbol}" PARENT_SCOPE)
  LIST(FIND ${manglistvar} ${newsymbol} already_mangled)
  IF(${already_mangled} EQUAL -1)
    LIST(APPEND ${manglistvar} ${newsymbol})
    SET(${manglistvar} ${${manglistvar}} PARENT_SCOPE)
    LIST(APPEND ${typedefvar}  "typedef ${oldsymbol} ${newsymbol}")
    SET(${typedefvar}  ${${typedefvar}}  PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()

# generate the macros
FUNCTION(TRIBITS_ETI_GENERATE_MACROS etifields etisetvar etiexcludelist manglinglistvar typedeflistvar)
  SET(manglinglist  "${${manglinglistvar}}")
  SET(typedeflist   "${${typedeflistvar}}")
  SPLIT("${etifields}" "\\|" etifields)
  IF(${PROJECT}_VERBOSE_CONFIGURE)
    MESSAGE(STATUS "ETI fields: ${etifields}")
    MESSAGE(STATUS "ETI set: ${etisetvar}")
    MESSAGE(STATUS "ETI excludes: ${etiexcludelist}")
  ENDIF()
  # we make lists of tuples first, because we want to make sure they don't have duplicates
  # this algorithm is O(N^2) in the number of instantiations in etisetvar
  MATH(EXPR num_macros "(${ARGC}-5)/2")
  IF(${num_macros} EQUAL 0)
    RETURN()
  ENDIF()
  # process macro fields into a list of indices
  FOREACH(m RANGE 1 ${num_macros})
    MATH(EXPR m2   "(${m}-1)*2")
    MATH(EXPR m2p1 "${m2}+1")
    LIST(GET ARGN ${m2}   macroarg)
    LIST(GET ARGN ${m2p1} macrovar${m})
    STRING(REGEX REPLACE "^(.*)\\(.*"     "\\1" macroname${m} ${macroarg})
    STRING(REGEX REPLACE "^.*\\((.*)\\)$" "\\1" macrofields${m} ${macroarg})
    SPLIT("${macrofields${m}}" "," macrofields)
    TRIBITS_ETI_INDEX_MACRO_FIELDS("${etifields}" "${macrofields}" macroindex${m})
    IF(${PROJECT}_VERBOSE_CONFIGURE)
      MESSAGE(STATUS "Parsed macro ${macroname${m}}(${macrofields${m}}) into variable ${macrovar${m}} (index list ${macroindex${m}})")
    ENDIF()
  ENDFOREACH()
  # process the exclusions once
  FOREACH(excl ${etiexcludelist})
    TRIBITS_ETI_EXPLODE("${etifields}" "${excl}" e_excl)
    IF("${e_excl}" STREQUAL "TRIBITS_ETI_BAD_PARSE")
      MESSAGE(FATAL_ERROR "TRIBITS_GENERATE_ETI_MACROS: exclusion did not parse: ${excl}")
    ENDIF()
    LIST(APPEND processed_excludes "${e_excl}")
  ENDFOREACH()
  LIST(LENGTH etifields numfields)
  MATH(EXPR NFm1 "${numfields}-1")
  FOREACH(inst ${etisetvar})
    IF (${PROJECT}_VERBOSE_CONFIGURE)
      MESSAGE(STATUS "Processing instantiation: ${inst}") # comment
    ENDIF()
    TRIBITS_ETI_EXPLODE("${etifields}" "${inst}" tmp)
    IF("${tmp}" STREQUAL "TRIBITS_ETI_BAD_PARSE")
      MESSAGE(FATAL_ERROR "TRIBITS_GENERATE_ETI_MACROS: instantiation did not parse: ${inst}")
    ELSE()
      SET(inst "${tmp}")
    ENDIF()
    # check whether it is on the exclude list
    TRIBITS_ETI_CHECK_EXCLUSION("${processed_excludes}" "${inst}" excluded)
    IF(NOT excluded)
      SPLIT("${inst}" "\\|" inst)
      # append tuple to list
      FOREACH(m RANGE 1 ${num_macros})
        SET(tuple "")
        FOREACH(ind ${macroindex${m}})
          LIST(GET inst ${ind} t)
          IF("${t}" STREQUAL "TYPE-MISSING")
            SET(tuple "SKIP-TUPLE")
            BREAK()
          ENDIF()
          # mangle the types in the instantiation
          TRIBITS_ETI_MANGLE_SYMBOL_AUGMENT_MACRO(typedeflist t manglinglist)
          LIST(APPEND tuple ${t})
        ENDFOREACH()
        IF(NOT "${tuple}" STREQUAL "SKIP-TUPLE")
          JOIN(tuple "|" FALSE "${tuple}")
          LIST(APPEND macrotuples${m} ${tuple})
        ENDIF()
      ENDFOREACH()
    ENDIF()
  ENDFOREACH()
  # remove duplicates from lists
  FOREACH(m RANGE 1 ${num_macros})
    IF(DEFINED macrotuples${m})
      LIST(REMOVE_DUPLICATES macrotuples${m})
    ENDIF()
  ENDFOREACH()
  # build the macro strings
  FOREACH(m RANGE 1 ${num_macros})
    TRIBITS_ETI_BUILD_MACRO_STRING("${macroname${m}}" "${macrotuples${m}}" mac)
    SET(${macrovar${m}} "${mac}" PARENT_SCOPE)
  ENDFOREACH()
  # build the typedef string
  SET(${manglinglistvar} ${manglinglist} PARENT_SCOPE)
  SET(${typedeflistvar}  ${typedeflist}  PARENT_SCOPE)
ENDFUNCTION()

# generate the typedef macro
FUNCTION(TRIBITS_ETI_GENERATE_TYPEDEF_MACRO outputvar macroname typedeflist)
  SET(mac "#define ${macroname}() ")
  FOREACH(td ${typedeflist})
    SET(mac "${mac} \\\n\t${td};")
  ENDFOREACH()
  SET(${outputvar} "${mac}" PARENT_SCOPE)
ENDFUNCTION()
