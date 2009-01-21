INCLUDE(PackageAddExecutableAndTest)

PACKAGE_COPY_FILES_TO_BINARY_DIR(pamgenTests
     SOURCE_FILES run_test_script_01
     SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files
     DEST_FILES run_test_script_01.exe
     )

PACKAGE_COPY_FILES_TO_BINARY_DIR(test_01_copy
     SOURCE_FILES mesh.txt
                  run_test_script_01.res
     SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files
     )


