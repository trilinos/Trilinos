if(CMAKE_PROJECT_NAME STREQUAL "Seacas" OR CMAKE_PROJECT_NAME STREQUAL "SEACAS" )
TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  LIB_REQUIRED_TPLS fmt
)
else()
# Typically for Trilinos since don't have fmt as TPL, but instead have embedded versions.
TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(

)
endif()
