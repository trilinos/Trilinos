# Create the PyTrilinos_PACKAGES string
SET(${PACKAGE_NAME}_LIST "")
FOREACH(PACKAGE ${${PACKAGE_NAME}_PACKAGES})
  SET(${PACKAGE_NAME}_LIST "${${PACKAGE_NAME}_LIST}    ${PACKAGE}\n")
ENDFOREACH(PACKAGE ${${PACKAGE_NAME}_PACKAGES})

# Create the PyTrilinos_ALL string
SET(${PACKAGE_NAME}_ALL "__all__ = [")
FOREACH(PACKAGE ${${PACKAGE_NAME}_PACKAGES})
  SET(${PACKAGE_NAME}_ALL "${${PACKAGE_NAME}_ALL}'${PACKAGE}',\n          ")
ENDFOREACH(PACKAGE ${${PACKAGE_NAME}_PACKAGES})
SET(${PACKAGE_NAME}_ALL "${${PACKAGE_NAME}_ALL}]")

# Create the __init__.py file
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/__init__.src
  ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)

# Create the PropertyBase.py file
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/PropertyBase.src
  ${CMAKE_CURRENT_BINARY_DIR}/PropertyBase.py)

# Create the typed_tuple.py file
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/typed_tuple.src
  ${CMAKE_CURRENT_BINARY_DIR}/typed_tuple.py)

# Create the typed_list.py file
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/typed_list.src
  ${CMAKE_CURRENT_BINARY_DIR}/typed_list.py)

# Create the typed_dict.py file
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/typed_dict.src
  ${CMAKE_CURRENT_BINARY_DIR}/typed_dict.py)

# The NOX module employs namespaces, so include the NOX directory
IF(${PACKAGE_NAME}_ENABLE_NOX)
  ADD_SUBDIRECTORY(NOX)
ENDIF(${PACKAGE_NAME}_ENABLE_NOX)
