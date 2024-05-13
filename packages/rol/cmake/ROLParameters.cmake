if( ROL_ENABLE_PARAMETERLIST_VALIDATION )

	message("Enabling automated ParameterList detection and validation")

	if(NOT DEFINED ${PYTHON_EXECUTABLE})
		find_program(PYTHON_EXECUTABLE NAMES python3 python REQUIRED)
	endif()

	set( ROL_SOURCE_DIR "${PROJECT_SOURCE_DIR}/packages/rol" )
	set( ROL_PARAMETERS_SOURCE_DIR "${ROL_SOURCE_DIR}/rol_parameters" )

  set( ROL_BINARY_DIR "${PROJECT_BINARY_DIR}/packages/rol" )
	set( ROL_PARAMETERS_BINARY_DIR "${ROL_BINARY_DIR}/rol_parameters" )

  set( REQUIREMENTS_FILE "${ROL_PARAMETERS_SOURCE_DIR}/requirements.txt" )
	set( VENV_PATH "${ROL_PARAMETERS_BINARY_DIR}/venv" )

  # Set up Python virtual environment
	add_custom_target( setup_venv 
		                 COMMAND ${CMAKE_COMMAND} -E env ${PYTHON_EXECUTABLE} -m venv ${VENV_PATH}
	                   COMMAND ${CMAKE_COMMAND} -E env ${VENV_PATH}/bin/python -m pip install -r ${REQUIREMENTS_FILE} 
										 COMMENT "Setting up virtual environment and installing required Python packages"
										 WORKING_DIRECTORY ${CMAKE_BINARY_DIR} )

	message( "Python virtual environment path: ${VENV_PATH}" )
	message( STATUS "Run 'make setup_venv` or your equivalent build system command (e.g. ninja setup_venv') to setup the Python virtual environment before building rol_parameters")

	add_custom_target( rol_parameters 
 		                 COMMAND ${CMAKE_COMMAND} -E env PYTHONPYCACHEPREFIX=${CMAKE_BINARY_DIR}/pycache 
	        									 ${VENV_PATH}/bin/python ${ROL_PARAMETERS_SOURCE_DIR}/rol_parameters.py ${ROL_SOURCE_DIR} ${ROL_PARAMETERS_BINARY_DIR}
                     DEPENDS setup_venv
	              		 WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
 										 COMMENT "Running rol_parameters.py using the virtual environment")

   message( STATUS "Run 'make rol_parameters` or your equivalent build system command (e.g. ninja rol_parameters') to build the hierarchical parameter list from the ROL source tree")

endif()

