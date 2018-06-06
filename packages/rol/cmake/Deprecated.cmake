# Deprecated Parameters                   Use Instead
#------------------------------------------------------------------------------------------
# ROL_ENABLE_STD_SHARED_PTR:BOOL          ROL_Ptr:STRING="std::shared_ptr"
# ROL_ENABLE_TEUCHOS_RCP:BOOL             ROL_Ptr:STRING="Teuchos::RCP"
#
# ROL_ENABLE_BOOST_PROPERTY_TREE:BOOL     ROL_ParameterList:STRING="boost::property_tree"
# ROL_ENABLE_TEUCHOS_PARAMETERLIST:BOOL   ROL_ParameterList:STRING="Teuchos::ParameterList"
#

# If the deprecated BOOL parameter is defined, issue warning to use the new STRING
# parameter with the corresponding value. If the parameter is ON, the set the 
# new parameter accordingly
FUNCTION( SET_DEPRECATION_WARNING BOOL_KEY STRING_KEY STRING_VALUE )
  MESSAGE( WARNING "${PACKAGE_NAME}_ENABLE_${BOOL_KEY}:BOOL is deprecated."
  "Use ${PACKAGE_NAME}_${STRING_KEY}:STRING=\"${STRING_VALUE}\" instead. " )
  IF( ${${PACKAGE_NAME}_ENABLE_${BOOL_KEY}} )
    SET( ${PACKAGE_NAME}_${STRING_KEY} ${STRING_VALUE} )
  ENDIF()
ENDFUNCTION()

IF( DEFINED ${PACKAGE_NAME}_ENABLE_STD_SHARED_PTR )
SET_DEPRECATION_WARNING( "STD_SHARED_PTR" "Ptr" "std::shared_ptr" )
ENDIF()

IF( DEFINED ${PACKAGE_NAME}_ENABLE_TEUCHOS_RCP )
SET_DEPRECATION_WARNING( "TEUCHOS_RCP" "Ptr" "Teuchos::RCP" )
ENDIF()

IF( DEFINED ${PACKAGE_NAME}_ENABLE_BOOST_PROPERTY_TREE )
SET_DEPRECATION_WARNING( "BOOST_PROPERTY_TREE" "ParameterList" "boost::property_tree" )
ENDIF()

IF( DEFINED ${PACKAGE_NAME}_ENABLE_TEUCHOS_PARAMETERLIST )
SET_DEPRECATION_WARNING( "TEUCHOS_PARAMETERLIST" "ParameterList" "Teuchos::ParameterList" )
ENDIF()

