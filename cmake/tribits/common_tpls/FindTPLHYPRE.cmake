find_package(HYPRE)

if(HYPRE_FOUND)
	tribits_extpkg_create_imported_all_libs_target_and_config_file( HYPRE
		INNER_FIND_PACKAGE_NAME  HYPRE
        IMPORTED_TARGETS_FOR_ALL_LIBS  HYPRE::HYPRE )
endif()
