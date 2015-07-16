# command that copies files during build
TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_input_files_for_driver
    SOURCE_FILES
        input_07132015.xml
        rcbTest.xml
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)
