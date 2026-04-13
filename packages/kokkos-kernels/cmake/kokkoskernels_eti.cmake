#
# @FUNCTION: KOKKOSKERNELS_ETI_MAKE_LIST
#
# Create combinatorial sets of all enable ETI options.
# Consider a template T<A,B> where A is an index type and B is a floating type.
# If we have two lists INDEX=INT;UINT64_T and FLOAT=FLOAT;DOUBLE,
# we can invoke the function to generate ETI for all combinations as
# KOKKOSKERNELS_ETI_MAKE_LIST(ETI_FOR_T TYPE_LISTS INDEX FLOAT)
# Upon returning from the function, the variable ETI_FOR_T
# will be a list containing four entries:
# ${ETI_FOR_T}=T_INT_FLOAT;T_INT_DOUBLE;T_UINT64_T_FLOAT;T_UINT64_T_DOUBLE;
# Additionally, each of entries in the list is itself a variable name
# containing the C++ ETI type list, e.g.
# ${T_INT_FLOAT}=int,float
#
# Usage::
#
#   KOKKOSKERNELS_ETI_MAKE_LIST(
#     <ETI_LIST_NAME>
#     [TYPE_LISTS list1 [list2 ...]]
#   )
#   ``<ETI_LIST_NAME>``
#
#   The name of the list output variable that will contain all generated ETI combinations
#
#   ``[TYPE_LISTS list1 [[list2...]]``
#
#   The names of the lists containing ETI types. For a template T<A,B>,
#   then A will take every value in list1 and B will take every value in list2.
#   The types listed here should be the CMake names like DOUBLE and EXECSPACE_SERIAL
function(kokkoskernels_eti_make_list ETI_LIST_NAME)
  cmake_parse_arguments(ETI "" "" "TYPE_LISTS" ${ARGN})
  list(LENGTH ETI_TYPE_LISTS ETI_LIST_LENGTH)
  math(EXPR RANGE_VARIABLE "${ETI_LIST_LENGTH} - 1")
  foreach(IDX RANGE ${RANGE_VARIABLE})
    list(GET ETI_TYPE_LISTS ${IDX} LIST_NAME)
    set(LIST${IDX}_NAME ${LIST_NAME})
  endforeach()
  foreach(TYPE0 ${${LIST0_NAME}})
    if(KOKKOSKERNELS_INST_${TYPE0})
      set(NAME0 ${ETI_LIST_NAME}_${TYPE0})
      set(LIST0 ${TYPE0})
      if(ETI_LIST_LENGTH GREATER 1)
        foreach(TYPE1 ${${LIST1_NAME}})
          if(KOKKOSKERNELS_INST_${TYPE1})
            set(NAME1 ${NAME0}_${TYPE1})
            set(LIST1 ${LIST0})
            list(APPEND LIST1 ${TYPE1})
            if(ETI_LIST_LENGTH GREATER 2)
              foreach(TYPE2 ${${LIST2_NAME}})
                if(KOKKOSKERNELS_INST_${TYPE2})
                  set(NAME2 ${NAME1}_${TYPE2})
                  set(LIST2 ${LIST1})
                  list(APPEND LIST2 ${TYPE2})
                  if(ETI_LIST_LENGTH GREATER 3)
                    foreach(TYPE3 ${${LIST3_NAME}})
                      if(KOKKOSKERNELS_INST_${TYPE3})
                        set(NAME3 ${NAME2}_${TYPE3})
                        set(LIST3 ${LIST2})
                        list(APPEND LIST3 ${TYPE3})
                        if(ETI_LIST_LENGTH GREATER 4)
                          foreach(TYPE4 ${${LIST4_NAME}})
                            if(KOKKOSKERNELS_INST_${TYPE4})
                              set(NAME4 ${NAME3}_${TYPE4})
                              set(LIST4 ${LIST3})
                              list(APPEND LIST4 ${TYPE4})
                              if(ETI_LIST_LENGTH GREATER 5)
                                foreach(TYPE4 ${${LIST4_NAME}})
                                  if(KOKKOSKERNELS_INST_${TYPE5})
                                    set(NAME5 ${NAME4}_${TYPE5})
                                    set(LIST5 ${LIST4})
                                    list(APPEND LIST5 ${TYPE5})
                                    if(ETI_LIST_LENGTH GREATER 6)
                                      message(FATAL_ERROR "Do not support ETI with more than 6 types")
                                    else()
                                      #end of the eti list
                                      list(APPEND ${ETI_LIST_NAME} ${NAME5})
                                      set(${NAME5} ${LIST5} PARENT_SCOPE)
                                    endif()
                                  endif()
                                endforeach()
                              else()
                                #end of the eti list
                                list(APPEND ${ETI_LIST_NAME} ${NAME4})
                                set(${NAME4} ${LIST4} PARENT_SCOPE)
                              endif()
                            endif()
                          endforeach()
                        else()
                          #end of the eti list
                          list(APPEND ${ETI_LIST_NAME} ${NAME3})
                          set(${NAME3} ${LIST3} PARENT_SCOPE)
                        endif()
                      endif()
                    endforeach()
                  else()
                    #end of the eti list
                    list(APPEND ${ETI_LIST_NAME} ${NAME2})
                    set(${NAME2} ${LIST2} PARENT_SCOPE)
                  endif()
                endif()
              endforeach()
            else()
              #end of the eti list
              list(APPEND ${ETI_LIST_NAME} ${NAME1})
              set(${NAME1} ${LIST1} PARENT_SCOPE)
            endif()
          endif()
        endforeach()
      else()
        #end of the eti list
        list(APPEND ${ETI_LIST_NAME} ${NAME0})
        set(${NAME0} ${LIST0} PARENT_SCOPE)
      endif()
    endif()
  endforeach()
  set(${ETI_LIST_NAME} ${${ETI_LIST_NAME}} PARENT_SCOPE)
endfunction()

macro(kokkoskernels_generate_eti FUNCTION_NAME SUBFOLDER)
  cmake_parse_arguments(ETI "" "HEADER_LIST;SOURCE_LIST" "TYPE_LISTS;COMPONENTS" ${ARGN})

  string(TOUPPER "${FUNCTION_NAME}" UPPER_NAME)
  set(ETI_AVAIL_MACRO "KOKKOS${UPPER_NAME}_ETI_SPEC_AVAIL")
  set(ETI_DECL_MACRO "KOKKOS${UPPER_NAME}_ETI_SPEC_DECL")
  set(ETI_INST_MACRO "KOKKOS${UPPER_NAME}_ETI_SPEC_INST")

  # if this is tied to particular components
  # see whether those components are enabled
  kokkoskernels_is_enabled(COMPONENTS ${ETI_COMPONENTS} OUTPUT_VARIABLE ETI_COMP_IS_ENABLED)

  if(ETI_COMP_IS_ENABLED)
    message(STATUS "Creating ETI files for ${FUNCTION_NAME}")
    kokkoskernels_eti_make_list(${FUNCTION_NAME}_eti TYPE_LISTS ${ETI_TYPE_LISTS})
    foreach(ETI ${${FUNCTION_NAME}_eti})
      set(MACRO_STRING "(")
      foreach(TYPE_NAME ${${ETI}})
        string(APPEND MACRO_STRING "${${TYPE_NAME}_CPP_TYPE},")
      endforeach()
      string(APPEND MACRO_STRING ")")
      string(REPLACE ",)" ")" MACRO_STRING ${MACRO_STRING})
      #Make a single header file for all instances
      list(APPEND ${UPPER_NAME}_ETI_AVAIL_LIST "${ETI_AVAIL_MACRO}${MACRO_STRING}")
      list(APPEND ${UPPER_NAME}_ETI_DECL_LIST "${ETI_DECL_MACRO}${MACRO_STRING}")
      #Make a different source file for each instance
      set(INST_SOURCE "${ETI_COMPONENTS}/eti/generated_specializations_cpp/${SUBFOLDER}/${ETI}.cpp")
      set(INST_TEMPLATE "${ETI_COMPONENTS}/eti/generated_specializations_cpp/${SUBFOLDER}/Kokkos${FUNCTION_NAME}_eti_spec_inst.cpp.in")
      set(${UPPER_NAME}_ETI_INST_BLOCK "${ETI_INST_MACRO}${MACRO_STRING}")
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${INST_TEMPLATE} ${CMAKE_CURRENT_BINARY_DIR}/${INST_SOURCE})
      list(APPEND ${ETI_SOURCE_LIST} ${CMAKE_CURRENT_BINARY_DIR}/${INST_SOURCE})
    endforeach()
  else()
    message(STATUS "Skipping ETI files for ${FUNCTION_NAME} because not all components are enabled")
  endif()

  set(AVAIL_HEADER "${ETI_COMPONENTS}/eti/generated_specializations_hpp/Kokkos${FUNCTION_NAME}_eti_spec_avail.hpp")
  set(AVAIL_TEMPLATE "${AVAIL_HEADER}.in")

  string(REPLACE ";" "\n" ${UPPER_NAME}_ETI_INST_BLOCK  "${${UPPER_NAME}_ETI_INST_LIST}")
  string(REPLACE ";" "\n" ${UPPER_NAME}_ETI_AVAIL_BLOCK "${${UPPER_NAME}_ETI_AVAIL_LIST}")

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${AVAIL_TEMPLATE} ${CMAKE_CURRENT_BINARY_DIR}/${AVAIL_HEADER})

  list(APPEND ${ETI_HEADER_LIST} ${CMAKE_CURRENT_BINARY_DIR}/${AVAIL_HEADER})

  set(DECL_HEADER "${ETI_COMPONENTS}/eti/generated_specializations_hpp/Kokkos${FUNCTION_NAME}_eti_spec_decl.hpp")
  set(DECL_TEMPLATE "${DECL_HEADER}.in")

  string(REPLACE ";" "\n" ${UPPER_NAME}_ETI_DECL_BLOCK "${${UPPER_NAME}_ETI_DECL_LIST}")

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${DECL_TEMPLATE} ${CMAKE_CURRENT_BINARY_DIR}/${DECL_HEADER})

  list(APPEND ${ETI_HEADER_LIST} ${CMAKE_CURRENT_BINARY_DIR}/${DECL_HEADER})
endmacro()
