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
FUNCTION(KOKKOSKERNELS_ETI_MAKE_LIST ETI_LIST_NAME)
  CMAKE_PARSE_ARGUMENTS(ETI
    ""
    ""
    "TYPE_LISTS"
    ${ARGN}
  )
  LIST(LENGTH ETI_TYPE_LISTS ETI_LIST_LENGTH)
  MATH(EXPR RANGE_VARIABLE "${ETI_LIST_LENGTH} - 1")
  FOREACH(IDX RANGE ${RANGE_VARIABLE})
    LIST(GET ETI_TYPE_LISTS ${IDX} LIST_NAME)
    SET(LIST${IDX}_NAME ${LIST_NAME})
  ENDFOREACH()
  FOREACH(TYPE0 ${${LIST0_NAME}})
   IF (KOKKOSKERNELS_INST_${TYPE0})
    SET(NAME0 ${ETI_LIST_NAME}_${TYPE0})
    SET(LIST0 ${TYPE0})
    IF (ETI_LIST_LENGTH GREATER 1)
     FOREACH(TYPE1 ${${LIST1_NAME}})
      IF (KOKKOSKERNELS_INST_${TYPE1})
       SET(NAME1 ${NAME0}_${TYPE1})
       SET(LIST1 ${LIST0}) 
       LIST(APPEND LIST1 ${TYPE1})
       IF (ETI_LIST_LENGTH GREATER 2)
        FOREACH(TYPE2 ${${LIST2_NAME}})
         IF (KOKKOSKERNELS_INST_${TYPE2})
          SET(NAME2 ${NAME1}_${TYPE2})
          SET(LIST2 ${LIST1}) 
          LIST(APPEND LIST2 ${TYPE2})
          IF (ETI_LIST_LENGTH GREATER 3)
            FOREACH(TYPE3 ${${LIST3_NAME}})
             IF (KOKKOSKERNELS_INST_${TYPE3})
              SET(NAME3 ${NAME2}_${TYPE3})
              SET(LIST3 ${LIST2}) 
              LIST(APPEND LIST3 ${TYPE3})
              IF (ETI_LIST_LENGTH GREATER 4)
               FOREACH(TYPE4 ${${LIST4_NAME}})
                IF (KOKKOSKERNELS_INST_${TYPE4})
                 SET(NAME4 ${NAME3}_${TYPE4})
                 SET(LIST4 ${LIST3}) 
                 LIST(APPEND LIST4 ${TYPE4})
                 IF (ETI_LIST_LENGTH GREATER 5)
                  FOREACH(TYPE4 ${${LIST4_NAME}})
                   IF (KOKKOSKERNELS_INST_${TYPE5})
                    SET(NAME5 ${NAME4}_${TYPE5})
                    SET(LIST5 ${LIST4}) 
                    LIST(APPEND LIST5 ${TYPE5})
                    IF (ETI_LIST_LENGTH GREATER 6)
                      MESSAGE(FATAL_ERROR "Do not support ETI with more than 6 types")
                    ELSE()
                      #end of the eti list
                      LIST(APPEND ${ETI_LIST_NAME} ${NAME5})
                      SET(${NAME5} ${LIST5} PARENT_SCOPE)
                    ENDIF()
                   ENDIF()
                  ENDFOREACH()
                 ELSE()
                   #end of the eti list
                   LIST(APPEND ${ETI_LIST_NAME} ${NAME4})
                   SET(${NAME4} ${LIST4} PARENT_SCOPE)
                 ENDIF()
                ENDIF()
               ENDFOREACH()
              ELSE()
                #end of the eti list
                LIST(APPEND ${ETI_LIST_NAME} ${NAME3})
                SET(${NAME3} ${LIST3} PARENT_SCOPE)
              ENDIF()
             ENDIF()
            ENDFOREACH()
          ELSE()
            #end of the eti list
            LIST(APPEND ${ETI_LIST_NAME} ${NAME2})
            SET(${NAME2} ${LIST2} PARENT_SCOPE)
          ENDIF()
         ENDIF()
        ENDFOREACH()
       ELSE()
         #end of the eti list
         LIST(APPEND ${ETI_LIST_NAME} ${NAME1})
         SET(${NAME1} ${LIST1} PARENT_SCOPE)
       ENDIF()
      ENDIF()
     ENDFOREACH()
    ELSE()
     #end of the eti list
     LIST(APPEND ${ETI_LIST_NAME} ${NAME0})
     SET(${NAME0} ${LIST0} PARENT_SCOPE)
    ENDIF()
   ENDIF()
  ENDFOREACH()
  SET(${ETI_LIST_NAME} ${${ETI_LIST_NAME}} PARENT_SCOPE)
ENDFUNCTION(KOKKOSKERNELS_ETI_MAKE_LIST)

MACRO(KOKKOSKERNELS_GENERATE_ETI FUNCTION_NAME SUBFOLDER)
  CMAKE_PARSE_ARGUMENTS(ETI
    ""
    "HEADER_LIST;SOURCE_LIST"
    "TYPE_LISTS;COMPONENTS"
    ${ARGN})

  STRING(TOUPPER "${FUNCTION_NAME}" UPPER_NAME)
  SET(ETI_AVAIL_MACRO "KOKKOS${UPPER_NAME}_ETI_SPEC_AVAIL")
  SET(ETI_DECL_MACRO "KOKKOS${UPPER_NAME}_ETI_SPEC_DECL")
  SET(ETI_INST_MACRO  "KOKKOS${UPPER_NAME}_ETI_SPEC_INST")

  # if this is tied to particular components
  # see whether those components are enabled
  KOKKOSKERNELS_IS_ENABLED(
    COMPONENTS ${ETI_COMPONENTS}
    OUTPUT_VARIABLE ETI_COMP_IS_ENABLED
  )

  IF (ETI_COMP_IS_ENABLED)
    MESSAGE(STATUS "Creating ETI files for ${FUNCTION_NAME}")
    KOKKOSKERNELS_ETI_MAKE_LIST(${FUNCTION_NAME}_eti TYPE_LISTS ${ETI_TYPE_LISTS})
    FOREACH(ETI ${${FUNCTION_NAME}_eti})
      SET(MACRO_STRING "(")
      FOREACH(TYPE_NAME ${${ETI}})
        STRING(APPEND MACRO_STRING "${${TYPE_NAME}_CPP_TYPE},")
      ENDFOREACH()
      STRING(APPEND MACRO_STRING ")")
      STRING(REPLACE ",)" ")" MACRO_STRING ${MACRO_STRING})
      #Make a single header file for all instances
      LIST(APPEND ${UPPER_NAME}_ETI_AVAIL_LIST "${ETI_AVAIL_MACRO}${MACRO_STRING}")
      LIST(APPEND ${UPPER_NAME}_ETI_DECL_LIST "${ETI_DECL_MACRO}${MACRO_STRING}")
      #Make a different source file for each instance
      SET(INST_SOURCE   "${ETI_COMPONENTS}/eti/generated_specializations_cpp/${SUBFOLDER}/${ETI}.cpp")
      SET(INST_TEMPLATE "${ETI_COMPONENTS}/eti/generated_specializations_cpp/${SUBFOLDER}/Kokkos${FUNCTION_NAME}_eti_spec_inst.cpp.in")
      SET(${UPPER_NAME}_ETI_INST_BLOCK "${ETI_INST_MACRO}${MACRO_STRING}")
      CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${INST_TEMPLATE}
          ${CMAKE_CURRENT_BINARY_DIR}/${INST_SOURCE})
      LIST(APPEND ${ETI_SOURCE_LIST} ${CMAKE_CURRENT_BINARY_DIR}/${INST_SOURCE})
    ENDFOREACH()
  ELSE()
    MESSAGE(STATUS "Skipping ETI files for ${FUNCTION_NAME} because not all components are enabled")
  ENDIF()

  SET(AVAIL_HEADER   "${ETI_COMPONENTS}/eti/generated_specializations_hpp/Kokkos${FUNCTION_NAME}_eti_spec_avail.hpp")
  SET(AVAIL_TEMPLATE "${AVAIL_HEADER}.in")

  STRING(REPLACE ";" "\n" ${UPPER_NAME}_ETI_INST_BLOCK  "${${UPPER_NAME}_ETI_INST_LIST}")
  STRING(REPLACE ";" "\n" ${UPPER_NAME}_ETI_AVAIL_BLOCK "${${UPPER_NAME}_ETI_AVAIL_LIST}")

  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${AVAIL_TEMPLATE}
      ${CMAKE_CURRENT_BINARY_DIR}/${AVAIL_HEADER})

  LIST(APPEND ${ETI_HEADER_LIST} ${CMAKE_CURRENT_BINARY_DIR}/${AVAIL_HEADER})

  SET(DECL_HEADER   "${ETI_COMPONENTS}/eti/generated_specializations_hpp/Kokkos${FUNCTION_NAME}_eti_spec_decl.hpp")
  SET(DECL_TEMPLATE "${DECL_HEADER}.in")

  STRING(REPLACE ";" "\n" ${UPPER_NAME}_ETI_DECL_BLOCK "${${UPPER_NAME}_ETI_DECL_LIST}")

  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${DECL_TEMPLATE}
      ${CMAKE_CURRENT_BINARY_DIR}/${DECL_HEADER})

  LIST(APPEND ${ETI_HEADER_LIST} ${CMAKE_CURRENT_BINARY_DIR}/${DECL_HEADER})
ENDMACRO(KOKKOSKERNELS_GENERATE_ETI)
