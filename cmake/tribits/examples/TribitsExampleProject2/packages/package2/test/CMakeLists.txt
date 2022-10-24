if (${PACAKGE_NAME}_ENABLE_Tpl3)
  set(tpl3_deps_string ", Tpl3{Tpl2a{tpl1}, Tpl2b{no deps}}")
else()
  set(tpl3_deps_string "")
endif()


tribits_add_test(package2-prg  NOEXEPREFIX  NOEXESUFFIX
  NAME Prg
  DIRECTORY ${PACKAGE_BINARY_DIR}/src
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Package2 Deps: Package1{tpl1}${tpl3_deps_string}"
  )
