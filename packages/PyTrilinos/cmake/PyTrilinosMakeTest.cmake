# @HEADER
# ************************************************************************
#
#                PyTrilinos: Python Interface to Trilinos
#                   Copyright (2010) Sandia Corporation
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Bill Spotz (wfspotz@sandia.gov)
#
# ************************************************************************
# @HEADER

INCLUDE(TribitsAddTest)

MACRO(PyTrilinos_MAKE_TEST TEST_NAME)

  ADD_CUSTOM_COMMAND(OUTPUT ${TEST_NAME}.py
    COMMAND ${CMAKE_BINARY_DIR}/packages/PyTrilinos/util/copyWithCMakeSubstitutions.py
            ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}.py.in
	    ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.py
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}.py.in)

  ADD_CUSTOM_TARGET(PyTrilinos_${TEST_NAME} ALL
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.py)

  TRIBITS_ADD_TEST(
    ${PYTHON_EXECUTABLE}
    NOEXEPREFIX
    NOEXESUFFIX
    NAME ${TEST_NAME}
    ARGS "${TEST_NAME}.py --testharness"
    STANDARD_PASS_OUTPUT
    ${ARGN}
    )

ENDMACRO(PyTrilinos_MAKE_TEST TEST_NAME)
