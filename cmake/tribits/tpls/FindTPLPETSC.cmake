SET(FIND_TPL_NAME "FindTPLPETSC.cmake")
SET(NEW_TPL_DIR "common_tpls")

MESSAGE(WARNING "WARNING: The file tpls/${FIND_TPL_NAME}"
  " has been moved to ${NEW_TPL_DIR}/${FIND_TPL_NAME}!"
  "Please use the moved copy as this deprecated copy will be removed soon!")
INCLUDE("${CMAKE_CURRENT_LIST_DIR}/../${NEW_TPL_DIR}/${FIND_TPL_NAME}")
