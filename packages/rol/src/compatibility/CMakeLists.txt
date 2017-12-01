
SET(DIR ${CMAKE_CURRENT_SOURCE_DIR})

IF(${PACKAGE_NAME}_ENABLE_STD_SHARED_PTR)
  INCLUDE_DIRECTORIES(${DIR}/std/shared_ptr)
  APPEND_GLOB(HEADERS ${DIR}/std/shared_ptr/ROL_Ptr.hpp)
ELSE()
  INCLUDE_DIRECTORIES(${DIR}/teuchos/rcp)
  APPEND_GLOB(HEADERS ${DIR}/teuchos/rcp/ROL_Ptr.hpp)
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_BOOST_PROPERTY_TREE)
  INCLUDE_DIRECTORIES(${DIR}/boost/property_tree)
  APPEND_GLOB(HEADERS ${DIR}/boost/property_tree/ROL_PropertyList.hpp)
ELSE()
  INCLUDE_DIRECTORIES(${DIR}/teuchos/)
  APPEND_GLOB(HEADERS ${DIR}/teuchos/parameterlist/ROL_SharedPointer.hpp)
ENDIF()
