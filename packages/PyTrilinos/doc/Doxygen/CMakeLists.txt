# -*- cmake -*-

# @HEADER
# ************************************************************************
#
#                PyTrilinos: Python Interface to Trilinos
#                   Copyright (2009) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
# USA
# Questions? Contact Bill Spotz (wfspotz@sandia.gov)
#
# ************************************************************************
# @HEADER

SET(ADDITIONAL_CLEAN_FILES "")

FOREACH(Package ${PyTrilinos_PACKAGES})

  IF(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${Package}_dox.i)

    IF(PyTrilinos_DOCSTRINGS)

      IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile_${Package}.in)

        CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile_${Package}.in
          ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile_${Package}
          @ONLY
          )
        EXECUTE_PROCESS(COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/docstrings.py
          --doxygen=${DOXYGEN_EXECUTABLE} ${Package}
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

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
