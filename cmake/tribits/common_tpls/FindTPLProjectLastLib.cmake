# This is a hack to handle the ${PROJECT_NAME}_EXTRA_LINK_FLAGS options as a
# TPL that every downstream TPL and package will depend on.

separate_arguments(TPL_${PROJECT_NAME}TribitsLastLib_LIBRARIES  NATIVE_COMMAND
   "${${PROJECT_NAME}_EXTRA_LINK_FLAGS}")

tribits_tpl_find_include_dirs_and_libraries( ${PROJECT_NAME}TribitsLastLib
  REQUIRED_LIBS_NAMES  NEVER_USED
  )
