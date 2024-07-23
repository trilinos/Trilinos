# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(AppendSet)
include(AppendGlobalSet)
include(Split)
include(Join)
include(PrintVar)

# this macro expands produces a list of all combinations of a variable number of lists:
# called like tribits_eti_type_expansion(setvar "T=t1|t2|t3" "S=s1|s2" "V=v1")
# returns setvar="T=T1 S=s1 V=v1;T=T1 S=s2 V=v1;
#                 T=T2 S=s1 V=v1;T=T2 S=s2 V=v1;
#                 T=T3 S=s1 V=v1;T=T3 S=s2 V=v1"
function(tribits_eti_type_expansion outvar first_list)
  # extract the field name and the list of types
  if (NOT "${first_list}" MATCHES "([^=]*)=(.*)")
    set(${outvar} "TRIBITS_ETI_BAD_ARGUMENTS" PARENT_SCOPE)
    return()
  else()
    set(label "${CMAKE_MATCH_1}")
    set(tlist "${CMAKE_MATCH_2}")
  endif()
  set(eti_result "") 
  split("${tlist}" "\\|" tlist)
  if ("${${outvar}}" STREQUAL "")
    set(accumulate OFF)
  else()
    set(accumulate ON)
  endif()
  list(LENGTH tlist tlist_len)
  if (${ARGC} GREATER 2)
    tribits_eti_type_expansion(sub ${ARGN})
    foreach(t ${tlist})
      string(STRIP "${t}" t)
      set(t "{${t}}")
      foreach(s ${sub})
        list(APPEND eti_result "${label}=${t} ${s}")
      endforeach()
    endforeach()
  else()
    foreach(t ${tlist})
      string(STRIP "${t}" t)
      set(t "{${t}}")
      list(APPEND eti_result "${label}=${t}")
    endforeach()
  endif()
  if (accumulate)
    set(${outvar} "${${outvar}};${eti_result}" PARENT_SCOPE)
  else()
    set(${outvar} "${eti_result}"              PARENT_SCOPE)
  endif()
endfunction()

# Explode an ETI set into a variable number of named component fields
function(tribits_eti_explode fields inst outvar)
  foreach (field ${fields})
    if(    "${inst}" MATCHES "( |^)${field}={([^{}$=]*)}( |$)")
      list(APPEND eti_result "${CMAKE_MATCH_2}")
    elseif("${inst}" MATCHES "( |^)${field}=([^{} $=]*)( |$)")
      list(APPEND eti_result "${CMAKE_MATCH_2}")
    elseif(NOT "${inst}" MATCHES "( |^)${field}=")
      list(APPEND eti_result "TYPE-MISSING")
    else()
      set(eti_result "TRIBITS_ETI_BAD_PARSE")
      break()
    endif()
  endforeach()
  join(eti_result "|" FALSE ${eti_result})
  set(${outvar} "${eti_result}" PARENT_SCOPE)
endfunction()

# effectively, a tupled regex, wrapped in a for loop
# given a list of processed excludes (a list of pipe-separated, bracket-delimited tuples of type/regexes),
# determine whether processed_inst (a pipe-separated, bracket tuple of types) is matched
# if the instantiation matches one of the exclusions, result is set true
function(tribits_eti_check_exclusion processed_excludes processed_inst excluded)
  split("${processed_inst}" "\\|" processed_inst)
  list(LENGTH processed_inst numfields)
  math(EXPR NFm1 "${numfields}-1")
  set(${excluded} OFF PARENT_SCOPE)
  # check to see whether this is excluded or not
  list(LENGTH processed_excludes numexcl)
  foreach(excl ${processed_excludes})
    split(${excl} "\\|" excl)
    list(LENGTH excl excp_len)
    if(${excp_len} EQUAL ${numfields})
      set(lcltest ON)
      foreach(i RANGE 0 ${NFm1})
        list(GET processed_inst ${i} f)
        list(GET excl           ${i} e)
        if (NOT "${f}" STREQUAL "TYPE-MISSING" AND NOT "${f}" MATCHES "${e}")
          set(lcltest OFF)
          break()
        endif()
      endforeach()
      if(lcltest)
        set(${excluded} ON PARENT_SCOPE)
        if (${PROJECT}_VERBOSE_CONFIGURE)
          message(STATUS "-- Instantiation excluded by ${excl}")
        endif()
        return()
      endif()
    endif()
  endforeach()
endfunction()

# given a list of field names and list of macrofields, translate the
# field names in macrofields into indices corresponding to the position in
# the list of field names
function(tribits_eti_index_macro_fields etifields macrofields indexvar)
  foreach(mf ${macrofields})
    string(STRIP "${mf}" mf)
    list(FIND etifields ${mf} ind)
    if(${ind} EQUAL -1)
      message(FATAL_ERROR "Macro variable ${mf} not found in list of fields: ${etifields}")
    endif()
    list(APPEND inds ${ind})
  endforeach()
  set(${indexvar} ${inds} PARENT_SCOPE)
endfunction()

# given a macro name and a list of tuples, generate a macro string
function(tribits_eti_build_macro_string macroname tuplelist outvar)
  set(str "#define ${macroname}(INSTMACRO)")
  foreach(tuple ${tuplelist})
    string(REPLACE "|" " , " tuple "${tuple}")
    set(str "${str}\\\n\tINSTMACRO( ${tuple} )")
  endforeach()
  set(str "${str}\n")
  set(${outvar} ${str} PARENT_SCOPE)
endfunction()

# utility for adding eti support to another package's library set
function(tribits_add_eti_instantiations PACKAGE)
  if(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    message(STATUS "Adding instantiations to ${PACKAGE} library set...")
  endif()
  append_global_set(${PACKAGE}_ETI_LIBRARYSET ${ARGN})
endfunction()

# utility for mangling type names to make them safe for the C preprocessor
function(tribits_eti_mangle_symbol
         mangledvar
         input)
  set(newvar "${input}")
  string(REPLACE "::" "_" newvar "${newvar}")
  string(REPLACE ","  "_" newvar "${newvar}")
  set(num_brackets 0)
  while("${newvar}" MATCHES "(.*)<([^<>]*)>(.*)")
    set(newvar "${CMAKE_MATCH_1}${num_brackets}${CMAKE_MATCH_2}${num_brackets}${CMAKE_MATCH_3}")
    math(EXPR num_brackets "${num_brackets}+1")
  endwhile()
  string(REPLACE " "  "" newvar "${newvar}")
  set(${mangledvar} "${newvar}" PARENT_SCOPE)
endfunction()

# utility for mangling type names to make them safe for the C preprocessor
# upon mangling, it automatically inserts an ifdef into the mangling macro
function(tribits_eti_mangle_symbol_augment_macro
         typedefvar
         symbolvar
         manglistvar)
  set(oldsymbol "${${symbolvar}}")
  if(NOT "${oldsymbol}" MATCHES "::|[<>, ]")
    return()
  endif()
  tribits_eti_mangle_symbol(newsymbol "${oldsymbol}")
  set(${symbolvar} "${newsymbol}" PARENT_SCOPE)
  list(FIND ${manglistvar} ${newsymbol} already_mangled)
  if(${already_mangled} EQUAL -1)
    list(APPEND ${manglistvar} ${newsymbol})
    set(${manglistvar} ${${manglistvar}} PARENT_SCOPE)
    list(APPEND ${typedefvar}  "typedef ${oldsymbol} ${newsymbol}")
    set(${typedefvar}  ${${typedefvar}}  PARENT_SCOPE)
  endif()
endfunction()

# generate the macros
function(tribits_eti_generate_macros etifields etisetvar etiexcludelist manglinglistvar typedeflistvar)
  set(manglinglist  "${${manglinglistvar}}")
  set(typedeflist   "${${typedeflistvar}}")
  split("${etifields}" "\\|" etifields)
  if(${PROJECT}_VERBOSE_CONFIGURE)
    message(STATUS "ETI fields: ${etifields}")
    message(STATUS "ETI set: ${etisetvar}")
    message(STATUS "ETI excludes: ${etiexcludelist}")
  endif()
  # we make lists of tuples first, because we want to make sure they don't have duplicates
  # this algorithm is o(N^2) in the number of instantiations in etisetvar
  math(EXPR num_macros "(${ARGC}-5)/2")
  if(${num_macros} EQUAL 0)
    return()
  endif()
  # process macro fields into a list of indices
  foreach(m RANGE 1 ${num_macros})
    math(EXPR m2   "(${m}-1)*2")
    math(EXPR m2p1 "${m2}+1")
    list(GET ARGN ${m2}   macroarg)
    list(GET ARGN ${m2p1} macrovar${m})
    string(REGEX REPLACE "^(.*)\\(.*"     "\\1" macroname${m} ${macroarg})
    string(REGEX REPLACE "^.*\\((.*)\\)$" "\\1" macrofields${m} ${macroarg})
    split("${macrofields${m}}" "," macrofields)
    tribits_eti_index_macro_fields("${etifields}" "${macrofields}" macroindex${m})
    if(${PROJECT}_VERBOSE_CONFIGURE)
      message(STATUS "Parsed macro ${macroname${m}}(${macrofields${m}}) into variable ${macrovar${m}} (index list ${macroindex${m}})")
    endif()
  endforeach()
  # process the exclusions once
  foreach(excl ${etiexcludelist})
    tribits_eti_explode("${etifields}" "${excl}" e_excl)
    if("${e_excl}" STREQUAL "TRIBITS_ETI_BAD_PARSE")
      message(FATAL_ERROR "TRIBITS_GENERATE_ETI_MACROS: exclusion did not parse: ${excl}")
    endif()
    list(APPEND processed_excludes "${e_excl}")
  endforeach()
  list(LENGTH etifields numfields)
  math(EXPR NFm1 "${numfields}-1")
  foreach(inst ${etisetvar})
    if (${PROJECT}_VERBOSE_CONFIGURE)
      message(STATUS "Processing instantiation: ${inst}") # comment
    endif()
    tribits_eti_explode("${etifields}" "${inst}" tmp)
    if("${tmp}" STREQUAL "TRIBITS_ETI_BAD_PARSE")
      message(FATAL_ERROR "TRIBITS_GENERATE_ETI_MACROS: instantiation did not parse: ${inst}")
    else()
      set(inst "${tmp}")
    endif()
    # check whether it is on the exclude list
    tribits_eti_check_exclusion("${processed_excludes}" "${inst}" excluded)
    if(NOT excluded)
      split("${inst}" "\\|" inst)
      # append tuple to list
      foreach(m RANGE 1 ${num_macros})
        set(tuple "")
        foreach(ind ${macroindex${m}})
          list(GET inst ${ind} t)
          if("${t}" STREQUAL "TYPE-MISSING")
            set(tuple "SKIP-TUPLE")
            break()
          endif()
          # mangle the types in the instantiation
          tribits_eti_mangle_symbol_augment_macro(typedeflist t manglinglist)
          list(APPEND tuple ${t})
        endforeach()
        if(NOT "${tuple}" STREQUAL "SKIP-TUPLE")
          join(tuple "|" FALSE "${tuple}")
          list(APPEND macrotuples${m} ${tuple})
        endif()
      endforeach()
    endif()
  endforeach()
  # remove duplicates from lists
  foreach(m RANGE 1 ${num_macros})
    if(DEFINED macrotuples${m})
      list(REMOVE_DUPLICATES macrotuples${m})
    endif()
  endforeach()
  # build the macro strings
  foreach(m RANGE 1 ${num_macros})
    tribits_eti_build_macro_string("${macroname${m}}" "${macrotuples${m}}" mac)
    set(${macrovar${m}} "${mac}" PARENT_SCOPE)
  endforeach()
  # build the typedef string
  set(${manglinglistvar} ${manglinglist} PARENT_SCOPE)
  set(${typedeflistvar}  ${typedeflist}  PARENT_SCOPE)
endfunction()

# generate the typedef macro
function(tribits_eti_generate_typedef_macro outputvar macroname typedeflist)
  set(mac "#define ${macroname}() ")
  foreach(td ${typedeflist})
    set(mac "${mac} \\\n\t${td};")
  endforeach()
  set(${outputvar} "${mac}" PARENT_SCOPE)
endfunction()
