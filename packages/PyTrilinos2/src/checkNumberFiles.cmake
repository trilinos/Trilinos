# @HEADER
# *****************************************************************************
#          PyTrilinos2: Automatic Python Interfaces to Trilinos Packages
#
# Copyright 2022 NTESS and the PyTrilinos2 contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

# Designed to be run as a script with "cmake -P"
#
# Set variable NUMBER_FILE when running this.

if(NOT NUMBER_FILE)
    message(FATAL_ERROR "NUMBER_FILE must be specified")
endif()

if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/../binder/PyTrilinos2_${NUMBER_FILE}.cpp)
    MATH(EXPR INDEX "${NUMBER_FILE}+1")
    while (EXISTS ${CMAKE_CURRENT_BINARY_DIR}/../binder/PyTrilinos2_${INDEX}.cpp)
        MATH(EXPR INDEX "${INDEX}+1")
    endwhile()
    MATH(EXPR INDEX "${INDEX}-1")
    message(FATAL_ERROR "File PyTrilinos2_${NUMBER_FILE}.cpp exists; please rerun the configuration with PyTrilinos2_BINDER_NUM_FILES at least equal to ${INDEX}.")
endif()
