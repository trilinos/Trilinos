find_package(magistrate)

if(magistrate_FOUND)
    tribits_extpkg_create_imported_all_libs_target_and_config_file( magistrate
	INNER_FIND_PACKAGE_NAME  magistrate
        IMPORTED_TARGETS_FOR_ALL_LIBS  vt::lib::magistrate )
endif()
