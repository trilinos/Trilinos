cmake_minimum_required(VERSION 3.0)

# Disable Kokkos warning about not supporting C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)

project(find_package_Trilinos C CXX)

find_package(Trilinos REQUIRED)
message("Trilinos_FOUND = '${Trilinos_FOUND}'")
message("Trilinos_SELECTED_PACAKGES_LIST = '${Trilinos_SELECTED_PACAKGES_LIST}'")
message("Trilinos_LIBRARIES = '${Trilinos_LIBRARIES}'")
