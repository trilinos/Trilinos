INCLUDE(PrintVar)
INCLUDE(AppendSet)

#
# Function that does an in-place sort of a list of items according to the
# ordering in a master list
#
# NOTE: This function has wost-case N^2 complexity as the number of packages N
# or TPLs increases.  It actually has N * n complexity where N is the total
# number of packages/TPLs and n is the number of passed-in packages/TPLs.
# However, since N is not likely to ever be more than a few hundred, this is
# likely not going to be a big performance problem.  If this does become a
# performance problem, LIST(SORT ...) could be used but would require some
# work to build up the datastructures to make this very efficient.
#

FUNCTION(TRIBITS_SORT_LIST_ACCORDING_TO_MASTER_LIST  MASTER_LIST  LIST_VAR_INOUT)

  #MESSAGE("TRIBITS_SORT_LIST_ACCORDING_TO_MASTER_LIST:")
  #PRINT_VAR(MASTER_LIST)
  #PRINT_VAR(LIST_VAR_INOUT)
  #PRINT_VAR(${LIST_VAR_INOUT})

  SET(SORTED_LIST)

  FOREACH(ITEM ${MASTER_LIST})
    LIST(FIND ${LIST_VAR_INOUT} ${ITEM} ITEM_IDX)
     IF (NOT ITEM_IDX EQUAL -1)
      APPEND_SET(SORTED_LIST ${ITEM})
    ENDIF()
  ENDFOREACH()

  #PRINT_VAR(SORTED_LIST)

  SET(${LIST_VAR_INOUT} ${SORTED_LIST} PARENT_SCOPE)

ENDFUNCTION()
