# -*- cmake -*-

# @HEADER
# ***********************************************************************
#
#          PyTrilinos: Python Interfaces to Trilinos Packages
#                 Copyright (2014) Sandia Corporation
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
# Questions? Contact William F. Spotz (wfspotz@sandia.gov)
#
# ***********************************************************************
# @HEADER

SET(ADDITIONAL_CLEAN_FILES "")

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/docstrings.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/docstrings.py)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/doxy2swig.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/doxy2swig.py)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/makedeps.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/makedeps.py)

FOREACH(Package ${PyTrilinos_PACKAGES})

  IF(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${Package}_dox.i)

    IF(PyTrilinos_DOCSTRINGS)

      IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile_${Package}.in)

        CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile_${Package}.in
          ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile_${Package}
          @ONLY
          )
        EXECUTE_PROCESS(COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/docstrings.py
          --doxygen=${DOXYGEN_EXECUTABLE} ${Package}
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
          ERROR_VARIABLE Doxy_CONVERT_ERROR)
        IF(Doxy_CONVERT_ERROR)
          MESSAGE(FATAL_ERROR
            "Error converting ${Package} C++ Doxygen comments to Python docstrings:\n  "
            ${Doxy_CONVERT_ERROR})
        ENDIF(Doxy_CONVERT_ERROR)

      ENDIF()

    ELSE(PyTrilinos_DOCSTRINGS)

      EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E touch
        ${CMAKE_CURRENT_BINARY_DIR}/${Package}_dox.i
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
      MESSAGE(STATUS
        "PyTrilinos.${Package} will not have doxygen-generated docstrings")

    ENDIF(PyTrilinos_DOCSTRINGS)

    IF(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${Package}_dox.i)
      LIST(APPEND ADDITIONAL_CLEAN_FILES ${Package}_dox.i)
    ENDIF()

  ENDIF()

ENDFOREACH(Package)

# Set this directory's clean files
GET_DIRECTORY_PROPERTY(clean_files ADDITIONAL_MAKE_CLEAN_FILES)
LIST(APPEND            clean_files ${ADDITIONAL_CLEAN_FILES})
LIST(REMOVE_DUPLICATES clean_files)
LIST(REMOVE_ITEM       clean_files "")
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${clean_files}")
