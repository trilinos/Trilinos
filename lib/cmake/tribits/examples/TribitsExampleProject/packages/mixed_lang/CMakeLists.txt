tribits_package( MixedLang  ENABLE_SHADOWING_WARNINGS  CLEANED )

tribits_add_debug_option()

tribits_add_show_deprecated_warnings_option()

add_subdirectory(src)

tribits_add_test_directories(test)

tribits_package_postprocess()
