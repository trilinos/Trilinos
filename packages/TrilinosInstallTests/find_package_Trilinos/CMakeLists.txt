cmake_minimum_required(VERSION 3.17)

project(find_package_Trilinos NONE)

find_package(Trilinos REQUIRED)
message("Trilinos_FOUND = '${Trilinos_FOUND}'")
message("Trilinos_SELECTED_PACAKGES_LIST = '${Trilinos_SELECTED_PACAKGES_LIST}'")
message("Trilinos_LIBRARIES = '${Trilinos_LIBRARIES}'")
