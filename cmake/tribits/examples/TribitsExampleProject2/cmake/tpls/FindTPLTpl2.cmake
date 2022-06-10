tribits_tpl_find_include_dirs_and_libraries( Tpl2
  REQUIRED_HEADERS Tpl2a.hpp  # Only look for one header file to find include dir
  REQUIRED_LIBS_NAMES tpl2b tpl2a
  )
