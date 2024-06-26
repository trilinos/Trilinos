find_package(ROCSOLVER)

if(ROCSOLVER_FOUND)
	tribits_extpkg_create_imported_all_libs_target_and_config_file( ROCSOLVER
		INNER_FIND_PACKAGE_NAME  ROCSOLVER
        IMPORTED_TARGETS_FOR_ALL_LIBS  roc::rocsolver )
endif()

