if (${PACKAGE_NAME}_ENABLE_Package2)
  set(package2_deps_str
    "Package2{Package1{tpl1}, Tpl3{Tpl2a{tpl1}, Tpl2b{no deps}}}, ")
else()
  set(package2_deps_str "")
endif()

if (${PACKAGE_NAME}_ENABLE_Tpl4)
  set(tpl4_deps_str
    "Tpl4{Tpl3{Tpl2a{tpl1}, Tpl2b{no deps}}, Tpl2a{tpl1}, Tpl2b{no deps}}, ")
else()
  set(tpl4_deps_str "")
endif()

tribits_add_test(package3-prg  NOEXEPREFIX  NOEXESUFFIX
  NAME Prg
  DIRECTORY ${PACKAGE_BINARY_DIR}/src
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "Package3 Deps: ${package2_deps_str}Package1{tpl1}, ${tpl4_deps_str}Tpl2a{tpl1}, Tpl2b{no deps}"
  )
