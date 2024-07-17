find_package(ROCBLAS)

if(ROCBLAS_FOUND)
    tribits_extpkg_create_imported_all_libs_target_and_config_file( ROCBLAS
	INNER_FIND_PACKAGE_NAME  ROCBLAS
        IMPORTED_TARGETS_FOR_ALL_LIBS  roc::rocblas )
endif()

