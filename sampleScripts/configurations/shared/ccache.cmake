# Compiler caching configuration via ccache
# Caches compilation results to speed up rebuilds (especially effective with few changes)

# Wrap C++ and C compilers with ccache
set(CMAKE_CXX_COMPILER_LAUNCHER "ccache" CACHE STRING "")
set(CMAKE_C_COMPILER_LAUNCHER "ccache" CACHE STRING "")
