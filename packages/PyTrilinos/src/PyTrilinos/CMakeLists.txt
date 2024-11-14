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

# Create the PyTrilinos_LIST string
SET(PyTrilinos_LIST "")
FOREACH(TRIBITS_PACKAGE ${PyTrilinos_PACKAGES})
  SET(PyTrilinos_LIST "${PyTrilinos_LIST}    ${TRIBITS_PACKAGE}\n")
ENDFOREACH(TRIBITS_PACKAGE ${PyTrilinos_PACKAGES})

# Create the PyTrilinos_ALL string
SET(PyTrilinos_ALL "__all__ = [")
FOREACH(TRIBITS_PACKAGE ${PyTrilinos_PACKAGES})
  SET(PyTrilinos_ALL "${PyTrilinos_ALL}'${TRIBITS_PACKAGE}',\n          ")
ENDFOREACH(TRIBITS_PACKAGE ${PyTrilinos_PACKAGES})
SET(PyTrilinos_ALL "${PyTrilinos_ALL}]")

# Define a trigger target
ADD_CUSTOM_TARGET(PyTrilinos_Pure_Python_Modules ALL)

# Define the list of local modules
SET(LOCAL_MODULES __init__ PropertyBase typed_tuple typed_list typed_dict)

# Initialize the additional "make clean" files
SET(ADDITIONAL_CLEAN_FILES "")

# Create each local module, its byte-compiled file and its
# installation rules
FOREACH(MODULE ${LOCAL_MODULES})

  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}.py.in
    ${CMAKE_CURRENT_BINARY_DIR}/${MODULE}.py)
  ADD_CUSTOM_COMMAND(TARGET PyTrilinos_Pure_Python_Modules
    COMMAND ${Python3_EXECUTABLE} -c "import py_compile; py_compile.compile('${MODULE}.py', '${MODULE}.pyc')"
    COMMENT "Byte compiling ${MODULE}.py"
    BYPRODUCTS ${MODULE}.pyc
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

# The Teuchos module employs namespaces, so include the Teuchos directory
ADD_SUBDIRECTORY(Teuchos)

# The Isorropia module employs namespaces, so include the Isorropia
# directory
IF(PyTrilinos_ENABLE_Isorropia)
  ADD_SUBDIRECTORY(Isorropia)
ENDIF(PyTrilinos_ENABLE_Isorropia)

# The NOX module employs namespaces, so include the NOX directory
IF(PyTrilinos_ENABLE_NOX)
  ADD_SUBDIRECTORY(NOX)
ENDIF(PyTrilinos_ENABLE_NOX)

# The LOCA module employs namespaces, so include the LOCA directory
IF(NOX_ENABLE_LOCA)
  ADD_SUBDIRECTORY(LOCA)
ENDIF(NOX_ENABLE_LOCA)
