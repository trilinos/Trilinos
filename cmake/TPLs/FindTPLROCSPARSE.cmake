find_package(ROCSPARSE)

if(ROCSPARSE_FOUND)
	tribits_extpkg_create_imported_all_libs_target_and_config_file( ROCSPARSE
		INNER_FIND_PACKAGE_NAME  ROCSPARSE
        IMPORTED_TARGETS_FOR_ALL_LIBS  roc::rocsparse )
endif()

