# -*- cmake -*-

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

# Create the PyTrilinos_LIST string
SET(PyTrilinos_LIST "")
FOREACH(PACKAGE ${PyTrilinos_PACKAGES})
  SET(PyTrilinos_LIST "${PyTrilinos_LIST}    ${PACKAGE}\n")
ENDFOREACH(PACKAGE ${PyTrilinos_PACKAGES})

# Create the PyTrilinos_ALL string
SET(PyTrilinos_ALL "__all__ = [")
FOREACH(PACKAGE ${PyTrilinos_PACKAGES})
  SET(PyTrilinos_ALL "${PyTrilinos_ALL}'${PACKAGE}',\n          ")
ENDFOREACH(PACKAGE ${PyTrilinos_PACKAGES})
SET(PyTrilinos_ALL "${PyTrilinos_ALL}]")

# Define a trigger target
ADD_CUSTOM_TARGET(Pure_Python_Modules ALL)

# Define the list of local modules
SET(LOCAL_MODULES __init__ PropertyBase typed_tuple typed_list typed_dict)

# Initialize the additional "make clean" files
SET(ADDITIONAL_CLEAN_FILES "")

# Create each local module, its byte-compiled file and its
# installation rules
FOREACH(MODULE ${LOCAL_MODULES})

  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}.py.in
    ${CMAKE_CURRENT_BINARY_DIR}/${MODULE}.py)
  ADD_CUSTOM_COMMAND(TARGET Pure_Python_Modules
    COMMAND ${PYTHON_EXECUTABLE} -c "import py_compile; py_compile.compile('${MODULE}.py')"
    COMMENT "Byte compiling ${MODULE}.py"
    VERBATIM)
  INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/${MODULE}.py
    ${CMAKE_CURRENT_BINARY_DIR}/${MODULE}.pyc
    DESTINATION ${PyTrilinos_INSTALL_DIR})
  LIST(APPEND ADDITIONAL_CLEAN_FILES ${MODULE}.pyc)
  
ENDFOREACH(MODULE ${LOCAL_MODULES})

# Add the additional "make clean" files
GET_DIRECTORY_PROPERTY(clean_files ADDITIONAL_MAKE_CLEAN_FILES)
LIST(APPEND            clean_files ${ADDITIONAL_CLEAN_FILES})
LIST(REMOVE_DUPLICATES clean_files)
LIST(REMOVE_ITEM       clean_files "")
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${clean_files}")

# The Isorropia module employs namespaces, so include the Isorropia
# directory
IF(PyTrilinos_NABLE_Isorropia)
  ADD_SUBDIRECTORY(Isorropia)
ENDIF(PyTrilinos_NABLE_Isorropia)

# The NOX module employs namespaces, so include the NOX directory
IF(PyTrilinos_ENABLE_NOX)
  ADD_SUBDIRECTORY(NOX)
ENDIF(PyTrilinos_ENABLE_NOX)
