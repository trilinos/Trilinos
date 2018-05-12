INCLUDE(BuildOptionFunctions)
INCLUDE(Deprecated)

SET_PROPERTY( GLOBAL PROPERTY TEUCHOS_PARAMETERLIST "Teuchos::ParameterList" )
SET_PROPERTY( GLOBAL PROPERTY TEUCHOS_RCP           "Teuchos::RCP"           )
SET_PROPERTY( GLOBAL PROPERTY BOOST_PROPERTY_TREE   "boost::property_tree"   )
SET_PROPERTY( GLOBAL PROPERTY STD_SHARED_PTR        "std::shared_ptr"        )

SET_DEFAULT( ROL_Ptr           "Teuchos::RCP" )
SET_DEFAULT( ROL_ParameterList "Teuchos::ParameterList" )

# Compatibility Directory 
SET( DIR ${${PACKAGE_NAME}_SOURCE_DIR}/src/compatibility )

# Override if specified
IF( ${PACKAGE_NAME}_ENABLE_STD_SHARED_PTR ) 
  SET_PROPERTY( GLOBAL PROPERTY PTR_IMPL "std::shared_ptr"       )
  SET_PROPERTY( GLOBAL PROPERTY PTR_DIR  "${DIR}/std/shared_ptr" )
ELSE()
  SET_PROPERTY( GLOBAL PROPERTY PTR_IMPL "Teuchos::RCP"       )
  SET_PROPERTY( GLOBAL PROPERTY PTR_DIR "${DIR}/teuchos/rcp" )
ENDIF()

IF( ${PACKAGE_NAME}_ENABLE_BOOST_PROPERTY_TREE ) 
  SET_PROPERTY( GLOBAL PROPERTY PARAMETERLIST_IMPL "boost::property_tree"       )
  SET_PROPERTY( GLOBAL PROPERTY PARAMETERLIST_DIR  "${DIR}/boost/property_tree" )
ELSE()
  SET_PROPERTY( GLOBAL PROPERTY PARAMETERLIST_IMPL "Teuchos::ParameterList"       )
  SET_PROPERTY( GLOBAL PROPERTY PARAMETERLIST_DIR  "${DIR}/teuchos/parameterlist" )
ENDIF()



