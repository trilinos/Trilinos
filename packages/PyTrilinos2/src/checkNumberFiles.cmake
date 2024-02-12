# @HEADER
# ***********************************************************************
#
#          PyTrilinos2: Automatic Python Interfaces to Trilinos Packages
#                 Copyright (2022) Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia
# Corporation, the U.S. Government retains certain rights in this
# software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Kim Liegeois (knliege@sandia.gov)
#
# ***********************************************************************
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
