tribits_repository_define_tpls(
  Tpl1      "cmake/tpls/"      PT
  Tpl2      "cmake/tpls/"      PT
  Tpl3      "cmake/tpls/"      PT
  Tpl4      "cmake/tpls/"      PT
  )

# Temp hack for setting up TPL dependencies
set(Tpl2_LIB_ENABLED_DEPENDENCIES Tpl1 CACHE STRING "Tpl2 upstream dependencies")
set(Tpl3_LIB_ENABLED_DEPENDENCIES Tpl2 CACHE STRING "Tpl3 upstream dependencies")
set(Tpl4_LIB_ENABLED_DEPENDENCIES "Tpl2;Tpl3" CACHE STRING "Tpl4 upstream dependencies")
