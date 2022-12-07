tribits_repository_define_tpls(
  Tpl1      "cmake/tpls/"      PT
  Tpl2      "cmake/tpls/"      PT
  Tpl3      "${CMAKE_CURRENT_LIST_DIR}/cmake/tpls/"      PT
  Tpl4      "cmake/tpls/"      PT
  )

# NOTE: Above we are setting the findmod path to an absolute path just to test
# that case with TPL dependencies (see trilinos/Trilinos#10774).
