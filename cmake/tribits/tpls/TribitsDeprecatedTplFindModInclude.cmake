macro(TPL_DEPRECATED_TPL_FIND_MOD_INCLUDE  TPL_NAME   NEW_TPL_DIR)
  set(FIND_TPL_NAME "FindTPL${TPL_NAME}.cmake")
  message(WARNING "WARNING: The file tpls/${FIND_TPL_NAME}"
    " for the TPL='${TPL_NAME}' has been moved to ${NEW_TPL_DIR}/${FIND_TPL_NAME}!"
    "  Please use the moved copy as this deprecated copy will be removed soon!"
    "  Make this change in the file:\n"
    "  ${${TPL_NAME}_TPLS_LIST_FILE}\n")
  include("${CMAKE_CURRENT_LIST_DIR}/../${NEW_TPL_DIR}/${FIND_TPL_NAME}")
endmacro()
