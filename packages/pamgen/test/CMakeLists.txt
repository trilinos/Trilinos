INCLUDE(PackageAddExecutableAndTest)

PACKAGE_COPY_FILES_TO_BINARY_DIR(pamgenTests
     SOURCE_FILES run_test_script_01
     SOURCE_DIR test_files
     DEST_FILES run_test_script_01.exe
     EXEDEPS libpamgen.a
     )


