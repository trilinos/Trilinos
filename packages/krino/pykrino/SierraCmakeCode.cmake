# Standalone build of the pykrino bindings (e.g. the py-krino Spack package):
# build against an already-installed krino. Nothing here depends on the krino
# source tree -- krino::krino_lib supplies the krino/stk/Trilinos headers and
# (via krino_lib's closure) krino_rol transitively.

cmake_minimum_required(VERSION 3.23)
# C is needed because krino's package config chain resolves MPI's C component.
project(pykrino LANGUAGES C CXX)

find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 REQUIRED)
find_package(krino REQUIRED)

set(_install_dir
    lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages/pykrino)

# Pure-python package files: src/*.py
file(GLOB _py_files ${CMAKE_CURRENT_SOURCE_DIR}/src/*.py)
install(FILES ${_py_files} DESTINATION ${_install_dir})

# pybind11 extension module: src/bindings/*.cpp
file(GLOB _sources ${CMAKE_CURRENT_SOURCE_DIR}/src/bindings/*.cpp)
pybind11_add_module(_pykrino ${_sources})
target_include_directories(_pykrino PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/bindings)
target_compile_features(_pykrino PRIVATE cxx_std_17)
target_link_libraries(_pykrino PRIVATE krino::krino_lib)
set_target_properties(_pykrino PROPERTIES SUFFIX ".so")
install(TARGETS _pykrino DESTINATION ${_install_dir})
