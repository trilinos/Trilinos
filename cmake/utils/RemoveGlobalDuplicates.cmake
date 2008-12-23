INCLUDE(AssertDefined)
INCLUDE(GlobalSet)

FUNCTION(REMOVE_GLOBAL_DUPLICATES VARNAME)
  ASSERT_DEFINED(${VARNAME})
  IF (${VARNAME})
    SET(TMP ${${VARNAME}})
    LIST(REMOVE_DUPLICATES TMP)
    GLOBAL_SET(${VARNAME} ${TMP})
  ENDIF()
ENDFUNCTION()

# 2008/11/21: rabartl: The above function is necessary in order to
# preserve the "global" natrue of the variable.  If you just call
# LIST(REMOVE_DUPLICATES ...) it will actually create a local variable
# of the same name and shadow the global varible.  It took me something
# like two hours to track down that bug!
