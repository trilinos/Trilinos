# Do a shared library build (BinUtils needs -fPIC to work with this)
SET(TPL_ENABLE_BinUtils OFF CACHE BOOL "")
SET(BUILD_SHARED_LIBS ON CACHE BOOL "")
